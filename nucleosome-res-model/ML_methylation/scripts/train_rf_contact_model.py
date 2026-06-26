#!/usr/bin/env python3
"""
train_rf_contact_model.py
=========================
Train Random Forest models that predict 200 bp Micro-C contact enrichment
from one-dimensional epigenomic tracks.

This public version keeps the analysis logic used for the manuscript but
collects the user-editable settings near the top of the file. It uses
checkpoint files so that completed stages are reused when the script is
restarted.
"""

import os, json, warnings, sys, platform
import numpy as np
import pandas as pd
import h5py
import cooler
import pyBigWig
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pathlib import Path
import joblib
from joblib import dump, load, Parallel, delayed
from scipy.stats import pearsonr
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.model_selection import GroupKFold
from skopt import BayesSearchCV
from skopt.space import Integer, Real, Categorical
from scipy.ndimage import convolve, gaussian_filter1d
import logging

warnings.filterwarnings('ignore')

# ============================================================
# STAGE 0 — Configuration
# ============================================================

# -------- Logging Setup --------
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    handlers=[
        logging.FileHandler("microc_rf_pipeline.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

plt.rcParams.update({
    'figure.dpi': 130,
    'savefig.dpi': 300,
    'font.size': 11,
    'axes.titlesize': 12,
    'axes.labelsize': 11,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'axes.spines.top': False,
    'axes.spines.right': False,
})

# Redirect all standard prints to the logger to capture state & progress
def print(*args, **kwargs):
    logger.info(" ".join(map(str, args)))

# -------- Paths --------
# These paths can be changed by editing this block or by setting environment variables:
#   export ML_CONTACT_DATA_DIR=/path/to/input_data
#   export ML_CONTACT_CHECKPOINT_DIR=checkpoints
#   export ML_CONTACT_PLOT_DIR=plots_regression
DATA_DIR       = Path(os.environ.get("ML_CONTACT_DATA_DIR", "/scratch/sk12344/ml_contact_prediction/data"))
CHECKPOINT_DIR = Path(os.environ.get("ML_CONTACT_CHECKPOINT_DIR", "checkpoints"))
PLOT_DIR       = Path(os.environ.get("ML_CONTACT_PLOT_DIR", "plots_regression"))

MCOOL_FILE = str(DATA_DIR / "GSE130275_mESC_WT_combined_2.6B.mcool")

# Grouped by context for plotting alignment
TRACK_FILES = {
    # Active marks
    "H3K27ac":  str(DATA_DIR / "ENCFF230RNU_h3k27ac.bigWig"),
    "H3K4me1":  str(DATA_DIR / "ENCFF282RLA_h3k4me1.bigWig"),
    "H3K4me3":  str(DATA_DIR / "ENCFF363LCX_h3k4me3.bigWig"),
    "H3K36me3": str(DATA_DIR / "ENCFF300QWW_h3k36me3.bigWig"),
    # Inactive marks
    "H3K9me3":  str(DATA_DIR / "ENCFF549VBE_h3k9me3.bigWig"),
    "H3K27me3": str(DATA_DIR / "ENCFF595SIA_h3k27me3.bigWig"),
    # Other features
    "ATAC":     str(DATA_DIR / "GSE172392_ATAC_seq_Control.bw"),
    "CTCF":     str(DATA_DIR / "GSE96107_ES_CTCF.bw"),
}

MARK_NAMES = list(TRACK_FILES.keys())
N_MARKS    = len(MARK_NAMES)

# Matplotlib colors assigned based on context
MARK_COLORS = {
    "H3K4me3": "forestgreen",
    "H3K27ac": "limegreen",
    "H3K4me1": "mediumseagreen",
    "H3K36me3": "darkgreen",
    "H3K9me3": "firebrick",
    "H3K27me3": "indianred",
    "ATAC": "royalblue",
    "CTCF": "darkorchid"
}

# -------- Resolution & Window --------
REQUESTED_RESOLUTION = 200
WINDOW_SIZE_BP       = 100_000

MCOOL_LEVELS = {
    11: 6400, 12: 3200, 13: 1600, 14: 800,
    15: 400,  16: 200,  17: 100
}

# -------- Task Definition --------
# 'regression' for continuous Log2 O/E prediction, 'classification' for enriched loop detection
ML_TASK = 'regression'
# Target clipping bounds for regression (prevents the model from getting skewed by extremes)
TARGET_MIN, TARGET_MAX = -5.0, 5.0
# Loop threshold for classification (Log2 O/E > 1.5 indicates strong enrichment)
CLASSIFICATION_THRESHOLD = 1.5

# -------- Chromosome split (mm10 autosomes) --------
TRAIN_CHROMS = [f'chr{i}' for i in range(1, 16)]   
VAL_CHROMS   = ['chr16', 'chr17']
TEST_CHROMS  = ['chr18', 'chr19']
ALL_CHROMS   = TRAIN_CHROMS + VAL_CHROMS + TEST_CHROMS

# -------- HPC & Reproducibility --------
RANDOM_SEED    = 42
N_VAL_WINDOWS  = 10   
TRACKING_WINDOW_SPLIT_COUNTS = {'train': 8, 'val': 6, 'test': 6}
TRACKING_MANIFEST_VERSION = 'regression_tracking_v2_20_windows'
N_JOBS         = -1  # Uses all available CPU cores for parallel tasks
RF_TRAIN_N_JOBS = 8

# -------- BayesSearchCV --------
N_BAYES_ITER = 20
CV_FOLDS     = 5

# Memory limits to prevent RAM exhaustion
MAX_PAIRS_PER_CHROM_TRAIN = 50_000
MAX_VALS_PER_STRATUM = 100_000
MAX_PAIRS_PER_BUCKET_PER_CHROM_TRAIN = 10_000
MAX_PAIRS_PER_BUCKET_PER_CHROM_VAL = 5_000
PREDICT_BATCH_SIZE = 500_000
N_STABILITY_RETRAINS = 5
N_PERMUTATION_REPEATS = 5
DISTANCE_BASELINE_MODE = "train_stratum_median"
PIPELINE_VERSION = "residual_rf_v5"
DIAG_BUCKETS = [1, 5, 10, 20, 50, 100, 200]
SKIP_UPSTREAM_PLOTS_IF_FEATURES_READY = False
SMOOTH_KERNEL_SIZE = 3
SMOOTH_GROW_THRESHOLD = 1.5
EXACT_DIAG_CORRECTION_SIGMA = 2.0
PRED_LEVEL_CORRECTION_SIGMA = 1.0
PRED_LEVEL_CORRECTION_BINS = 31
PRIMARY_REGRESSION_TARGET_SPACE = "smoothed_oe"
ENABLE_CONTROL_VARIANTS = False
ACTIVE_REGRESSION_VARIANTS = [
    'epi_plus_distance_regression',
    'epi_plus_distance_pairwise_regression',
    'epi_plus_distance_exact_diag_calibrated_regression',
    'epi_plus_distance_pairwise_exact_diag_calibrated_regression',
]
CONTROL_REGRESSION_VARIANTS = [
    'epi_only_regression',
]

# -------- Checkpoint map --------
CKPT = {
    'res_config':      CHECKPOINT_DIR / 'resolution_config.json',
    'chipseq_scalers': CHECKPOINT_DIR / 'chipseq_scalers.npz',
    'stratum_expected': CHECKPOINT_DIR / 'stratum_expected.npz',  
    'window_index':    CHECKPOINT_DIR / 'window_index.csv',
    'val_windows':     CHECKPOINT_DIR / 'validation_windows.npy',
    'features_h5':     CHECKPOINT_DIR / 'features.h5',
    'sampled_features': CHECKPOINT_DIR / 'sampled_features.h5',
    'residual_baseline': CHECKPOINT_DIR / 'stratum_residual_baseline.npz',
    'model':           CHECKPOINT_DIR / 'model.joblib',
    'best_params':     CHECKPOINT_DIR / 'best_params.joblib',
    'model_meta':      CHECKPOINT_DIR / 'model_meta.json',
    'eval_report':     CHECKPOINT_DIR / 'evaluation_report.csv',
    'sample_report':   CHECKPOINT_DIR / 'sample_report_regression.csv',
    'tracking_windows': CHECKPOINT_DIR / 'tracking_windows_regression.csv',
}

# Create required directories
for d in [CHECKPOINT_DIR, PLOT_DIR]:
    d.mkdir(exist_ok=True)


def write_run_environment_metadata() -> None:
    """Write a small metadata file to help reproduce the run."""
    try:
        import sklearn, scipy
        payload = {
            "python": sys.version,
            "platform": platform.platform(),
            "data_dir": str(DATA_DIR),
            "checkpoint_dir": str(CHECKPOINT_DIR),
            "plot_dir": str(PLOT_DIR),
            "random_seed": RANDOM_SEED,
            "train_chroms": TRAIN_CHROMS,
            "val_chroms": VAL_CHROMS,
            "test_chroms": TEST_CHROMS,
            "packages": {
                "numpy": np.__version__,
                "pandas": pd.__version__,
                "scikit-learn": sklearn.__version__,
                "scipy": scipy.__version__,
                "joblib": joblib.__version__,
            },
            "environment_variables": {
                "ML_CONTACT_DATA_DIR": os.environ.get("ML_CONTACT_DATA_DIR"),
                "ML_CONTACT_CHECKPOINT_DIR": os.environ.get("ML_CONTACT_CHECKPOINT_DIR"),
                "ML_CONTACT_PLOT_DIR": os.environ.get("ML_CONTACT_PLOT_DIR"),
            },
        }
        with open(CHECKPOINT_DIR / "run_environment.json", "w") as fh:
            json.dump(payload, fh, indent=2)
    except Exception as exc:
        logger.warning(f"Could not write run_environment.json: {exc}")


write_run_environment_metadata()


# ============================================================
# STAGE 1 — Resolve mcool resolution
# ============================================================
print("\n" + "="*60)
print("STAGE 1 — Resolving mcool resolution")
print("="*60)

if CKPT['res_config'].exists():
    with open(CKPT['res_config']) as f:
        res_cfg = json.load(f)
    print(f"  Loaded resolution config from checkpoint.")
else:
    candidates = {
        lvl: bs for lvl, bs in MCOOL_LEVELS.items()
        if bs <= REQUESTED_RESOLUTION and REQUESTED_RESOLUTION % bs == 0
    }
    if not candidates:
        raise ValueError("No compatible mcool resolution found.")

    native_level = max(candidates, key=lambda k: candidates[k])            
    native_res   = MCOOL_LEVELS[native_level]
    agg_factor   = REQUESTED_RESOLUTION // native_res
    n_bins       = WINDOW_SIZE_BP // REQUESTED_RESOLUTION

    res_cfg = {
        'requested_resolution': REQUESTED_RESOLUTION,
        'native_level':         native_level,
        'native_res':           native_res,
        'agg_factor':           agg_factor,
        'window_size_bp':       WINDOW_SIZE_BP,
        'n_bins':               n_bins,
        'mcool_uri':            f"{MCOOL_FILE}::/{native_level}",
    }
    with open(CKPT['res_config'], 'w') as f:
        json.dump(res_cfg, f, indent=2)

NATIVE_RES = res_cfg['native_res']
AGG_FACTOR = res_cfg['agg_factor']
N_BINS     = res_cfg['n_bins']
MCOOL_URI  = res_cfg['mcool_uri']
NATIVE_BINS_PER_WINDOW = N_BINS * AGG_FACTOR

# We added gaussian-smoothed versions for context
N_FEATURES   = 5 * N_MARKS + 2
FEATURE_NAMES = (
    [f'{m}_i'    for m in MARK_NAMES] +
    [f'{m}_j'    for m in MARK_NAMES] +
    [f'{m}_prod' for m in MARK_NAMES] +
    [f'{m}_int'  for m in MARK_NAMES] +
    [f'{m}_smooth_prod' for m in MARK_NAMES] + # New smoothed context feature
    ['log1p_dist', 'diag_index']
)
ACTIVE_FEATURE_NAMES = FEATURE_NAMES[:-2]
N_ACTIVE_FEATURES = len(ACTIVE_FEATURE_NAMES)
ALL_ACTIVE_FEATURE_NAMES = FEATURE_NAMES[:]
DIAG_BUCKET_EDGES = np.array(DIAG_BUCKETS + [N_BINS], dtype=np.int32)
N_DIAG_BUCKETS = len(DIAG_BUCKET_EDGES) - 1
DIAG_BUCKET_LABELS = [
    f"{DIAG_BUCKET_EDGES[i]}-{DIAG_BUCKET_EDGES[i + 1] - 1}" if i < N_DIAG_BUCKETS - 1 else f"{DIAG_BUCKET_EDGES[i]}+"
    for i in range(N_DIAG_BUCKETS)
]
MARK_FEATURE_GROUPS = {
    mark: [mi, N_MARKS + mi, 2 * N_MARKS + mi, 3 * N_MARKS + mi, 4 * N_MARKS + mi]
    for mi, mark in enumerate(MARK_NAMES)
}
DISTANCE_FEATURE_GROUP = {'Distance': [len(ACTIVE_FEATURE_NAMES), len(ACTIVE_FEATURE_NAMES) + 1]}
PAIRWISE_REGRESSION_FEATURE_NAMES = (
    [f'{m}_min' for m in MARK_NAMES] +
    [f'{m}_absdiff' for m in MARK_NAMES] +
    [f'{m}_contrast_i' for m in MARK_NAMES] +
    [f'{m}_contrast_j' for m in MARK_NAMES]
)
PAIRWISE_VARIANT_SUFFIX = '_pairwise_regression'
TARGET_SPACE_LABELS = {
    'raw_oe': 'Raw O/E',
    'raw_resid': 'Raw residual',
    'smoothed_oe': 'Smoothed O/E',
    'smoothed_resid': 'Smoothed residual',
}


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def ensure_parent(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def shared_plot_path(*parts: str) -> Path:
    return ensure_parent(PLOT_DIR.joinpath(*parts))


def comparison_plot_path(*parts: str) -> Path:
    return shared_plot_path('comparison', *parts)


def shared_plot_dir(name: str) -> Path:
    return ensure_dir(PLOT_DIR / name)


def make_variant_dirs(variant_name: str):
    variant_ckpt_dir = CHECKPOINT_DIR / variant_name
    variant_plot_dir = PLOT_DIR / variant_name
    variant_window_dir = variant_ckpt_dir / 'window_predictions_regression'
    ensure_dir(variant_ckpt_dir)
    ensure_dir(variant_plot_dir)
    ensure_dir(variant_window_dir)
    return variant_ckpt_dir, variant_plot_dir, variant_window_dir


def variant_ckpt_paths(variant_name: str) -> dict:
    variant_ckpt_dir, variant_plot_dir, variant_window_dir = make_variant_dirs(variant_name)
    return {
        'dir': variant_ckpt_dir,
        'plot_dir': variant_plot_dir,
        'window_dir': variant_window_dir,
        'model': variant_ckpt_dir / 'model.joblib',
        'best_params': variant_ckpt_dir / 'best_params.joblib',
        'model_meta': variant_ckpt_dir / 'model_meta.json',
        'eval_report': variant_ckpt_dir / 'evaluation_report.csv',
        'bias_report': variant_ckpt_dir / 'bucket_bias_report.csv',
        'cv_report': variant_ckpt_dir / 'cv_results.csv',
        'sample_report': variant_ckpt_dir / 'sample_report.csv',
        'stability_report': variant_ckpt_dir / 'stability_report.csv',
        'tracking_outputs': variant_ckpt_dir / 'tracking_outputs_regression.csv',
        'correction_summary': variant_ckpt_dir / 'exact_diag_correction_summary.json',
        'pred_level_correction': variant_ckpt_dir / 'pred_level_calibration_curve.csv',
        'exact_diag_report': variant_ckpt_dir / 'exact_diag_report.csv',
        'variant_summary': variant_ckpt_dir / 'train_val_test_summary.csv',
        'bucket_summary': variant_ckpt_dir / 'bucket_summary.csv',
        'comparison_dir': variant_plot_dir / 'comparison',
        'sampling_plot_dir': variant_plot_dir / 'sampling_diagnostics',
        'prediction_plot_dir': variant_plot_dir / 'prediction_maps',
        'distance_plot_dir': variant_plot_dir / 'distance_diagnostics',
        'feature_plot_dir': variant_plot_dir / 'feature_importance',
        'stability_plot_dir': variant_plot_dir / 'stability',
    }


def ordered_bucket_frame(df: pd.DataFrame, bucket_col: str = 'bucket') -> pd.DataFrame:
    if df.empty or bucket_col not in df.columns:
        return df
    out = df.copy()
    bucket_categories = ['overall'] + DIAG_BUCKET_LABELS
    bucket_order = {label: idx for idx, label in enumerate(bucket_categories)}
    out[bucket_col] = out[bucket_col].astype(str)
    out['_bucket_order'] = out[bucket_col].map(bucket_order).fillna(len(bucket_categories)).astype(int)
    sort_cols = ['_bucket_order', bucket_col]
    if 'split' in out.columns:
        sort_cols = ['split'] + sort_cols
    if 'target_space' in out.columns:
        sort_cols = ['target_space'] + sort_cols
    if 'chrom' in out.columns:
        sort_cols = ['chrom', '_bucket_order', bucket_col]
    return out.sort_values(sort_cols).drop(columns=['_bucket_order']).reset_index(drop=True)


def require_matching_row(df: pd.DataFrame, description: str, **filters) -> pd.Series:
    match = df.copy()
    for col, value in filters.items():
        match = match[match[col] == value]
    if match.empty:
        available = {
            col: sorted(df[col].dropna().astype(str).unique().tolist())
            for col in filters
            if col in df.columns
        }
        raise RuntimeError(f"Missing {description}. filters={filters}. available={available}")
    return match.iloc[0]


def get_feature_variant_columns(feature_variant: str):
    if feature_variant.startswith('epi_plus_distance'):
        return np.arange(len(ALL_ACTIVE_FEATURE_NAMES), dtype=np.int32), list(ALL_ACTIVE_FEATURE_NAMES)
    return np.arange(len(ACTIVE_FEATURE_NAMES), dtype=np.int32), list(ACTIVE_FEATURE_NAMES)


def build_pairwise_regression_features(X: np.ndarray) -> np.ndarray:
    x = np.asarray(X, dtype=np.float32)
    x_i = x[:, :N_MARKS]
    x_j = x[:, N_MARKS:2 * N_MARKS]
    x_int = x[:, 3 * N_MARKS:4 * N_MARKS]
    derived = [
        np.minimum(x_i, x_j),
        np.abs(x_i - x_j),
        x_i - x_int,
        x_j - x_int,
    ]
    return np.concatenate([x] + derived, axis=1)


def feature_names_for_variant(feature_variant: str) -> list[str]:
    _, feature_names = get_feature_variant_columns(feature_variant)
    if feature_variant.endswith(PAIRWISE_VARIANT_SUFFIX):
        return feature_names + list(PAIRWISE_REGRESSION_FEATURE_NAMES)
    return feature_names


def feature_groups_for_variant(feature_variant: str) -> dict:
    if feature_variant.startswith('epi_plus_distance'):
        base_groups = {
            mark: idxs[:] for mark, idxs in MARK_FEATURE_GROUPS.items()
        }
        base_groups['Distance'] = DISTANCE_FEATURE_GROUP['Distance'][:]
    else:
        base_groups = {mark: idxs[:] for mark, idxs in MARK_FEATURE_GROUPS.items()}
    if feature_variant.endswith(PAIRWISE_VARIANT_SUFFIX):
        offset = len(ALL_ACTIVE_FEATURE_NAMES)
        for mi, mark in enumerate(MARK_NAMES):
            base_groups[mark] = base_groups.get(mark, []) + [
                offset + mi,
                offset + N_MARKS + mi,
                offset + 2 * N_MARKS + mi,
                offset + 3 * N_MARKS + mi,
            ]
    return base_groups


def select_feature_matrix(X: np.ndarray, feature_variant: str) -> np.ndarray:
    cols, _ = get_feature_variant_columns(feature_variant)
    out = X[:, cols]
    if feature_variant.endswith(PAIRWISE_VARIANT_SUFFIX):
        out = build_pairwise_regression_features(out)
    return out

clr = cooler.Cooler(MCOOL_URI)

# Helper functions for loading matrices
def aggregate_matrix(mat: np.ndarray, k: int) -> np.ndarray:
    if k == 1: return mat
    n = (mat.shape[0] // k) * k
    mat = mat[:n, :n]
    return mat.reshape(n // k, k, n // k, k).sum(axis=(1, 3))

def load_microc_window(chrom: str, start: int, end: int) -> np.ndarray:
    mat = clr.matrix(balance=True).fetch((chrom, start, end))
    # KEEP NANS: Do not replace NaNs with 0. NaNs indicate unmappable regions which must be masked!
    if AGG_FACTOR > 1: 
        mat = aggregate_matrix(mat, AGG_FACTOR)
    return mat.astype(np.float32)

def load_bigwig_window(bw_path: str, chrom: str, start: int, end: int, n_bins: int) -> np.ndarray:
    with pyBigWig.open(bw_path) as bw:
        vals = bw.stats(chrom, start, end, type="mean", nBins=n_bins)
    out = np.array([0.0 if v is None else v for v in vals], dtype=np.float32)
    return np.nan_to_num(out, nan=0.0)


# ============================================================
# STAGE 2 — Genome-wide normalization stats
# ============================================================
print("\n" + "="*60)
print("STAGE 2 — Genome-wide normalization stats")
print("="*60)

chipseq_p99, stratum_expected = {}, np.array([])
chipseq_done = CKPT['chipseq_scalers'].exists()
expected_done = CKPT['stratum_expected'].exists()

if chipseq_done and expected_done:
    print("  Both stat files found — loading from disk.")
    cs = np.load(CKPT['chipseq_scalers'])
    chipseq_p99 = {m: float(cs[m]) for m in MARK_NAMES}
    stratum_expected = np.load(CKPT['stratum_expected'])['expected']
else:
    # Calculate 99th percentile for ChIP-seq scaling
    if not chipseq_done:
        print("  Calculating ChIP-seq percentiles...")
        mark_vals = {m: [] for m in MARK_NAMES}
        for chrom in TRAIN_CHROMS:
            if chrom not in clr.chromnames: continue
            n_bins_chrom = int(clr.chromsizes[chrom]) // REQUESTED_RESOLUTION
            if n_bins_chrom < 1: continue
            
            for mark, bw_path in TRACK_FILES.items():
                vals = load_bigwig_window(bw_path, chrom, 0, n_bins_chrom * REQUESTED_RESOLUTION, n_bins_chrom)
                mark_vals[mark].append(np.log1p(vals))
        
        for mark in MARK_NAMES:
            if mark_vals[mark]:
                all_v = np.concatenate(mark_vals[mark])
                chipseq_p99[mark] = float(np.percentile(all_v[all_v > 0], 99))
            else:
                chipseq_p99[mark] = 1.0
        np.savez(CKPT['chipseq_scalers'], **{m: np.float32(v) for m, v in chipseq_p99.items()})

    # Calculate distance-dependent expected values (stratum means) using randomized block sampling
    if not expected_done:
        print("  Calculating expected contact frequencies (mean per distance stratum)...")
        stratum_sums = np.zeros(N_BINS, dtype=np.float64)
        stratum_counts = np.zeros(N_BINS, dtype=np.float64)
        rng = np.random.default_rng(RANDOM_SEED)

        for chrom in TRAIN_CHROMS:
            if chrom not in clr.chromnames: continue
            try:
                chrom_len = clr.chromsizes[chrom]
                valid_max_start = chrom_len - WINDOW_SIZE_BP
                if valid_max_start <= 0: continue
                
                # Sample 100 random windows per chromosome for robust mean estimation
                starts = rng.integers(0, valid_max_start, size=100)
                for st in starts:
                    st = int(st // NATIVE_RES) * NATIVE_RES
                    
                    mat_raw = clr.matrix(balance=True).fetch((chrom, st, st + WINDOW_SIZE_BP))
                    if AGG_FACTOR > 1:
                        mat_raw = aggregate_matrix(mat_raw, AGG_FACTOR)
                        
                    for s in range(N_BINS):
                        diag = np.diag(mat_raw, k=s)
                        # np.isfinite drops the NaNs (unmappable bins) from the mean calculation
                        valid_diag = diag[np.isfinite(diag)] 
                        if len(valid_diag) > 0:
                            stratum_sums[s] += np.sum(valid_diag)
                            stratum_counts[s] += len(valid_diag)
            except Exception as e:
                logger.warning(f"Error sampling {chrom} for expected values: {e}")
                
        stratum_expected = np.zeros(N_BINS, dtype=np.float32)
        for s in range(N_BINS):
            val = stratum_sums[s] / stratum_counts[s] if stratum_counts[s] > 0 else 1e-9
            stratum_expected[s] = max(float(val), 1e-9)
        stratum_expected[0] = 1.0  
        np.savez(CKPT['stratum_expected'], expected=stratum_expected)

# Generate Stage 2 normalization tracking plots
print("  Generating Stage 2 normalization plots...")
fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(stratum_expected)
ax.set_title("Micro-C Stratum Means (Expected Contacts by Distance)")
ax.set_xlabel("Genomic Distance (bins)")
ax.set_ylabel("Mean Contact Frequency")
ax.set_yscale('log')
plt.tight_layout()
plt.savefig(shared_plot_path('normalization_diagnostics', 'stratum_expected.png'))
plt.close()

fig, ax = plt.subplots(figsize=(8, 4))
ax.bar(chipseq_p99.keys(), chipseq_p99.values(), color=[MARK_COLORS[m] for m in chipseq_p99.keys()])
ax.set_title("ChIP-seq 99th Percentile Values")
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(shared_plot_path('normalization_diagnostics', 'chipseq_p99.png'))
plt.close()

# Normalization functions
def normalize_chipseq(vals: np.ndarray, mark: str) -> np.ndarray:
    v = np.log1p(vals.astype(np.float32))
    p99 = chipseq_p99[mark]
    # No longer clipping at 1.0; divide by p99 to scale, but allow to exceed 1
    return (v / p99).astype(np.float32) if p99 > 0 else v

def normalize_microc(mat: np.ndarray) -> np.ndarray:
    n, eps = mat.shape[0], 1e-9
    oe = np.full((n, n), np.nan, dtype=np.float32) # Init with NaNs
    for s in range(1, n):
        exp = float(stratum_expected[s]) if s < len(stratum_expected) else eps
        exp = max(exp, eps)
        # Using NaNs in np.diag will safely propagate to log2_oe
        log2_oe = np.log2((np.diag(mat, k=s).astype(np.float32) / exp) + eps)
        idx = np.arange(n - s)
        oe[idx, idx + s] = oe[idx + s, idx] = log2_oe
    
    # Target Clipping: Bound the extreme O/E values to prevent ML confusion
    oe = np.clip(oe, TARGET_MIN, TARGET_MAX)
    return oe

# ============================================================
# STAGE 3 — Valid window extraction
# ============================================================
print("\n" + "="*60)
print("STAGE 3 — Valid window extraction")
print("="*60)

if CKPT['window_index'].exists():
    window_df = pd.read_csv(CKPT['window_index'])
    val_window_indices = np.load(CKPT['val_windows'])
else:
    chrom_to_split = {c: 'train' for c in TRAIN_CHROMS}
    chrom_to_split.update({c: 'val' for c in VAL_CHROMS})
    chrom_to_split.update({c: 'test' for c in TEST_CHROMS})
    records = []

    for chrom in ALL_CHROMS:
        if chrom not in clr.chromnames: continue
        cb = clr.bins().fetch(chrom).reset_index(drop=True)
        is_valid = cb['weight'].notna()
        block_ids = (is_valid != is_valid.shift()).cumsum()

        for _, group in cb.groupby(block_ids):
            if not group['weight'].notna().iloc[0] or len(group) < NATIVE_BINS_PER_WINDOW: continue
            grp = group.reset_index(drop=True)
            for w in range(len(grp) // NATIVE_BINS_PER_WINDOW):
                lo, hi = w * NATIVE_BINS_PER_WINDOW, w * NATIVE_BINS_PER_WINDOW + NATIVE_BINS_PER_WINDOW - 1
                g_start, g_end = int(grp.loc[lo, 'start']), int(grp.loc[hi, 'end'])
                if g_end - g_start == WINDOW_SIZE_BP:
                    records.append({'chrom': chrom, 'start': g_start, 'end': g_end, 'split': chrom_to_split[chrom]})

    window_df = pd.DataFrame(records).reset_index(drop=True)
    window_df.to_csv(CKPT['window_index'], index=False)
    
    rng = np.random.default_rng(RANDOM_SEED)
    picks = []
    for split, n_pick in [('train', 4), ('val', 3), ('test', 3)]:
        pool = window_df[window_df['split'] == split].index.tolist()
        picks.extend(rng.choice(pool, min(n_pick, len(pool)), replace=False).tolist())
    val_window_indices = np.array(picks[:N_VAL_WINDOWS])
    np.save(CKPT['val_windows'], val_window_indices)

val_windows = window_df.iloc[val_window_indices]


def build_tracking_windows(window_df: pd.DataFrame, manifest_path: Path) -> pd.DataFrame:
    expected_counts = dict(TRACKING_WINDOW_SPLIT_COUNTS)
    expected_total = int(sum(expected_counts.values()))
    if manifest_path.exists():
        tracking_df = pd.read_csv(manifest_path)
        tracking_df['track_index'] = tracking_df['track_index'].astype(int)
        split_counts = tracking_df['split'].value_counts().to_dict() if 'split' in tracking_df.columns else {}
        manifest_version = str(tracking_df['manifest_version'].iloc[0]) if 'manifest_version' in tracking_df.columns and len(tracking_df) else ''
        if (
            len(tracking_df) == expected_total and
            all(int(split_counts.get(split, 0)) == count for split, count in expected_counts.items()) and
            manifest_version == TRACKING_MANIFEST_VERSION
        ):
            return tracking_df.sort_values('track_index').reset_index(drop=True)
        logger.info(
            "  Existing tracking manifest is stale/incomplete (%s rows, counts=%s, version=%s); regenerating %s.",
            len(tracking_df),
            split_counts,
            manifest_version,
            manifest_path,
        )

    rng = np.random.default_rng(RANDOM_SEED)
    desired = list(TRACKING_WINDOW_SPLIT_COUNTS.items())
    picked_rows = []
    track_index = 0
    for split, n_pick in desired:
        pool = window_df[window_df['split'] == split].reset_index(drop=True)
        if pool.empty:
            logger.warning("No windows available for tracking split=%s", split)
            continue
        n_take = min(n_pick, len(pool))
        chosen = pool.iloc[rng.choice(len(pool), n_take, replace=False)].reset_index(drop=True)
        for _, row in chosen.iterrows():
            track_id = f"track_{track_index:02d}_{split}_{row.chrom}"
            picked_rows.append({
                'track_index': track_index,
                'track_id': track_id,
                'manifest_version': TRACKING_MANIFEST_VERSION,
                'split': split,
                'chrom': row.chrom,
                'start': int(row.start),
                'end': int(row.end),
            })
            track_index += 1
    tracking_df = pd.DataFrame(picked_rows).sort_values('track_index').reset_index(drop=True)
    tracking_df.to_csv(manifest_path, index=False)
    return tracking_df


tracking_windows_stage = build_tracking_windows(window_df, CKPT['tracking_windows'])
prediction_windows_eval = tracking_windows_stage[tracking_windows_stage['split'].isin(['val', 'test'])].reset_index(drop=True)


# ============================================================
# STAGE 4 — Feature extraction (Optimized & Parallelized)
# ============================================================
print("\n" + "="*60)
print("STAGE 4 — Feature extraction → HDF5")
print("="*60)

_ii, _jj = np.triu_indices(N_BINS) 
N_PAIRS  = len(_ii)
OFFDIAG_MASK = (_jj - _ii) > 0
_ii_offdiag = _ii[OFFDIAG_MASK]
_jj_offdiag = _jj[OFFDIAG_MASK]
N_OFFDIAG_PAIRS = len(_ii_offdiag)
DIAG_MATRIX_FULL = np.abs(np.arange(N_BINS)[:, None] - np.arange(N_BINS)[None, :]).astype(np.int16)


def upper_vector_to_symmetric_matrix(values: np.ndarray, fill_value=np.nan, dtype=np.float32) -> np.ndarray:
    if len(values) != N_OFFDIAG_PAIRS:
        raise ValueError(f"Expected {N_OFFDIAG_PAIRS} off-diagonal values, got {len(values)}")
    mat = np.full((N_BINS, N_BINS), fill_value, dtype=dtype)
    mat[_ii_offdiag, _jj_offdiag] = values
    mat[_jj_offdiag, _ii_offdiag] = values
    return mat


def smooth_oe_matrix(oe_mat: np.ndarray, min_diag: int = DIAG_BUCKET_EDGES[0], kernel_size: int = SMOOTH_KERNEL_SIZE) -> np.ndarray:
    if kernel_size == 3:
        kernel = np.array([[1, 2, 1], [2, 4, 2], [1, 2, 1]], dtype=np.float32)
    elif kernel_size == 5:
        kernel = np.array(
            [
                [1, 2, 3, 2, 1],
                [2, 4, 6, 4, 2],
                [3, 6, 9, 6, 3],
                [2, 4, 6, 4, 2],
                [1, 2, 3, 2, 1],
            ],
            dtype=np.float32,
        )
    else:
        raise ValueError("kernel_size must be 3 or 5")
    valid = np.isfinite(oe_mat) & (DIAG_MATRIX_FULL >= min_diag)
    weighted = np.where(valid, oe_mat, 0.0).astype(np.float32)
    denom = convolve(valid.astype(np.float32), kernel, mode='constant', cval=0.0)
    numer = convolve(weighted, kernel, mode='constant', cval=0.0)
    smoothed = np.full_like(oe_mat, np.nan, dtype=np.float32)
    np.divide(numer, np.maximum(denom, 1e-6), out=smoothed, where=denom > 0)
    smoothed[~valid] = np.nan
    np.fill_diagonal(smoothed, np.nan)
    return smoothed


def derive_window_regression_targets(targets: np.ndarray, diag_idx: np.ndarray) -> dict:
    raw_oe = upper_vector_to_symmetric_matrix(targets.astype(np.float32, copy=False))
    smoothed_oe = smooth_oe_matrix(raw_oe)
    raw_vec = raw_oe[_ii_offdiag, _jj_offdiag].astype(np.float32, copy=False)
    smooth_vec = smoothed_oe[_ii_offdiag, _jj_offdiag].astype(np.float32, copy=False)
    valid = np.isfinite(raw_vec) & np.isfinite(smooth_vec) & (diag_idx > 0)
    return {
        'raw_oe_mat': raw_oe,
        'smoothed_oe_mat': smoothed_oe,
        'raw_vec': raw_vec[valid],
        'smoothed_vec': smooth_vec[valid],
        'diag_vec': diag_idx[valid],
        'valid_vec_mask': valid,
    }

def extract_window_features(chrom: str, start: int, end: int, w_idx: int):
    """
    Worker function to extract features for a single window.
    Filters out any bin-pair containing NaN or Inf to prevent training crashes.
    """
    try:
        raw_mat = load_microc_window(chrom, start, end)
        oe_mat  = normalize_microc(raw_mat)

        tracks = np.zeros((N_MARKS, N_BINS), dtype=np.float32)
        tracks_smooth = np.zeros((N_MARKS, N_BINS), dtype=np.float32)
        
        for mi, (mark, bw_path) in enumerate(TRACK_FILES.items()):
            t = normalize_chipseq(load_bigwig_window(bw_path, chrom, start, end, N_BINS), mark)
            tracks[mi] = t
            # Neighborhood context: Gaussian smoothed track (sigma=2 bins)
            tracks_smooth[mi] = gaussian_filter1d(t, sigma=2.0)

        cumsum_ext = np.zeros((N_MARKS, N_BINS + 1), dtype=np.float32)
        cumsum_ext[:, 1:] = np.cumsum(tracks, axis=1)

        feats = np.empty((N_PAIRS, N_FEATURES), dtype=np.float32)
        feats[:, :N_MARKS]            = tracks[:, _ii].T
        feats[:, N_MARKS:2*N_MARKS]   = tracks[:, _jj].T
        feats[:, 2*N_MARKS:3*N_MARKS] = (tracks[:, _ii] * tracks[:, _jj]).T

        span = (_jj - _ii + 1).astype(np.float32)
        int_sum = cumsum_ext[:, _jj + 1] - cumsum_ext[:, _ii]
        feats[:, 3*N_MARKS:4*N_MARKS] = (int_sum / np.where(span > 0, span, 1.0)).T

        # Neighborhood Context Feature: Interaction of smoothed tracks
        feats[:, 4*N_MARKS:5*N_MARKS] = (tracks_smooth[:, _ii] * tracks_smooth[:, _jj]).T

        diag_idx = (_jj - _ii).astype(np.int32)
        feats[:, 5*N_MARKS]     = np.log1p(diag_idx.astype(np.float32) * REQUESTED_RESOLUTION)
        feats[:, 5*N_MARKS + 1] = diag_idx.astype(np.float32)

        if ML_TASK == 'regression':
            targets = oe_mat[_ii, _jj]
        else: # Classification: Encode target as 1 if it exceeds threshold, 0 otherwise
            targets = (oe_mat[_ii, _jj] > CLASSIFICATION_THRESHOLD).astype(np.int32)

        # Explicitly drop NaNs/Infs that would crash scikit-learn
        valid_mask = np.isfinite(targets) & np.all(np.isfinite(feats), axis=1)
        
        return {
            'w_idx': w_idx, 'start': start, 'end': end,
            'features': feats[valid_mask], 'targets': targets[valid_mask], 'diag_idx': diag_idx[valid_mask],
            'success': True
        }
    except Exception as e:
        return {'w_idx': w_idx, 'success': False, 'error': str(e)}

features_checkpoint_preexisting = CKPT['features_h5'].exists()

if features_checkpoint_preexisting:
    print("  features.h5 already exists — skipping extraction.")
else:
    with h5py.File(CKPT['features_h5'], 'w') as hf:
        hf.attrs['n_features'], hf.attrs['n_bins'] = N_FEATURES, N_BINS
        hf.attrs['feature_names'] = np.array(FEATURE_NAMES, dtype='S')
        hf.attrs['ml_task'] = ML_TASK

        for chrom, chrom_wins in window_df.groupby('chrom'):
            cg = hf.create_group(chrom)
            cg.attrs['split'] = chrom_wins.iloc[0]['split']

            # OPTIMIZATION: Extract windows for this chromosome in parallel
            print(f"    Extracting {len(chrom_wins)} windows for {chrom} in parallel...")
            results = Parallel(n_jobs=N_JOBS)(
                delayed(extract_window_features)(row.chrom, row.start, row.end, i) 
                for i, row in chrom_wins.reset_index(drop=True).iterrows()
            )

            # Write results sequentially to avoid HDF5 concurrent write corruption
            for res in results:
                if res['success']:
                    wg = cg.create_group(f"w{res['w_idx']:05d}")
                    wg.attrs['start'], wg.attrs['end'] = res['start'], res['end']
                    wg.create_dataset('features', data=res['features'], compression='gzip', compression_opts=4)
                    wg.create_dataset('targets', data=res['targets'], compression='gzip', compression_opts=4)
                    wg.create_dataset('diag_idx', data=res['diag_idx'], compression='gzip', compression_opts=4)

# Generate feature alignment tracking plots for validation windows
if features_checkpoint_preexisting and SKIP_UPSTREAM_PLOTS_IF_FEATURES_READY:
    print("  Reusing upstream checkpoints — skipping Stage 1/Stage 4 validation plot regeneration.")
else:
    print("  Generating Stage 1 (unnormalized) & Stage 4 (normalized) feature plots for validation windows...")

def plot_multitrack(row, out_path, normalized=True):
    try:
        raw_mat = load_microc_window(row.chrom, row.start, row.end)
        if normalized:
            oe_mat = normalize_microc(raw_mat)
            smoothed_mat = smooth_oe_matrix(oe_mat)
            matrix_panels = [
                (oe_mat, 'Normalized log2 O/E', 'RdBu_r', TARGET_MIN, TARGET_MAX),
                (smoothed_mat, 'Smoothed log2 O/E', 'RdBu_r', TARGET_MIN, TARGET_MAX),
            ]
            track_ylim_default = 3
        else:
            mat_plot = np.log1p(raw_mat)
            finite_positive = mat_plot[np.isfinite(mat_plot) & (mat_plot > 0)]
            vmax_val = np.percentile(finite_positive, 99) if finite_positive.size else 1.0
            matrix_panels = [(mat_plot, 'Raw Micro-C log1p observed', 'Reds', 0, vmax_val)]
            track_ylim_default = None
            
        n_matrix_rows = len(matrix_panels)
        fig = plt.figure(figsize=(9, 1.35 * N_MARKS + 3.8 * n_matrix_rows))
        gs = gridspec.GridSpec(
            N_MARKS + n_matrix_rows,
            2,
            width_ratios=[24, 1],
            height_ratios=[4] * n_matrix_rows + [1] * N_MARKS,
        )
        shared_ax = None
        title_prefix = f"{row.track_id} | {row.split} | {row.chrom}:{int(row.start):,}-{int(row.end):,}"
        for panel_i, (mat_plot, panel_title, cmap, vmin, vmax) in enumerate(matrix_panels):
            ax_hic = fig.add_subplot(gs[panel_i, 0], sharex=shared_ax)
            cax = fig.add_subplot(gs[panel_i, 1])
            if shared_ax is None:
                shared_ax = ax_hic
            im = ax_hic.imshow(mat_plot, cmap=cmap, vmin=vmin, vmax=vmax, aspect='auto')
            ax_hic.set_title(f"{title_prefix}\n{panel_title}")
            fig.colorbar(im, cax=cax)
            ax_hic.set_xlim(-0.5, N_BINS - 0.5)
            if panel_i < n_matrix_rows - 1:
                ax_hic.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
        
        x_axis = np.arange(N_BINS)
        
        for mi, (mark, bw_path) in enumerate(TRACK_FILES.items()):
            raw_track = load_bigwig_window(bw_path, row.chrom, row.start, row.end, N_BINS)
            if normalized:
                track_plot = normalize_chipseq(raw_track, mark)
                ylim_max = track_ylim_default
            else:
                track_plot = raw_track
                ylim_max = np.percentile(track_plot[track_plot > 0], 99) if np.any(track_plot > 0) else 1
                
            ax = fig.add_subplot(gs[mi + n_matrix_rows, 0], sharex=shared_ax)
            color = MARK_COLORS[mark]
            ax.fill_between(x_axis, track_plot, color=color, alpha=0.8, step='mid')
            
            ax.set_ylim(0, ylim_max)
            
            ax.set_ylabel(mark, rotation=0, labelpad=45, ha='right', va='center', fontsize=9, fontweight='bold')
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_visible(True)
            ax.tick_params(axis='y', left=True, labelleft=True, labelsize=8)
            
            if mi < N_MARKS - 1:
                ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
            else:
                ax.tick_params(axis='x', which='both', bottom=True, labelbottom=True)
                
        plt.tight_layout()
        plt.savefig(out_path)
        plt.close()
    except Exception as e:
        logger.warning(f"Failed to plot validation window {row.chrom}:{row.start}-{row.end} (norm={normalized}): {e}")

if not (features_checkpoint_preexisting and SKIP_UPSTREAM_PLOTS_IF_FEATURES_READY):
    for _, row in tracking_windows_stage.reset_index(drop=True).iterrows():
        # Stage 1: Plot the completely unnormalized regions
        plot_multitrack(row, shared_plot_path('raw_window_tracking', f'raw_{row.track_id}.png'), normalized=False)
        # Stage 4: Plot the normalized regions mapped correctly to the ML Input formats
        plot_multitrack(row, shared_plot_path('feature_track_alignment', f'norm_{row.track_id}.png'), normalized=True)

# ============================================================
# Downstream helpers — streamed residualized training/evaluation
# ============================================================
def decode_attr(value):
    if isinstance(value, bytes):
        return value.decode()
    return value

def iter_split_windows(h5_path: str, split: str):
    with h5py.File(h5_path, 'r') as hf:
        for chrom in sorted(hf.keys()):
            cg = hf[chrom]
            if decode_attr(cg.attrs['split']) != split:
                continue
            for win_key in sorted(cg.keys()):
                wg = cg[win_key]
                yield chrom, win_key, wg['features'][:], wg['targets'][:], wg['diag_idx'][:].astype(np.int32)

def assign_diag_buckets(diag_idx: np.ndarray) -> np.ndarray:
    bucket_idx = np.searchsorted(DIAG_BUCKET_EDGES[1:], diag_idx, side='right').astype(np.int16)
    valid = (diag_idx >= DIAG_BUCKET_EDGES[0]) & (diag_idx < DIAG_BUCKET_EDGES[-1])
    bucket_idx[~valid] = -1
    return bucket_idx

def count_diag_rows(h5_path: str, split: str) -> np.ndarray:
    counts = np.zeros(N_BINS, dtype=np.int64)
    for _, _, _, _, diag_idx in iter_split_windows(h5_path, split):
        valid = (diag_idx > 0) & (diag_idx < N_BINS)
        if np.any(valid):
            counts += np.bincount(diag_idx[valid], minlength=N_BINS)[:N_BINS]
    return counts

def compute_train_residual_baseline(h5_path: str, out_path: Path) -> np.ndarray:
    if out_path.exists():
        print(f"  Loaded residual baseline → {out_path}")
        return np.load(out_path)['baseline'].astype(np.float32)

    print("  Computing train-only residual baseline from features.h5...")
    diag_counts = count_diag_rows(h5_path, 'train')
    sampling_prob = np.ones(N_BINS, dtype=np.float64)
    over_cap = diag_counts > MAX_VALS_PER_STRATUM
    sampling_prob[over_cap] = MAX_VALS_PER_STRATUM / diag_counts[over_cap]

    rng = np.random.default_rng(RANDOM_SEED)
    sampled_values = [[] for _ in range(N_BINS)]
    sampled_counts = np.zeros(N_BINS, dtype=np.int64)

    for _, _, _, targets, diag_idx in iter_split_windows(h5_path, 'train'):
        derived = derive_window_regression_targets(targets, diag_idx)
        y_smooth = derived['smoothed_vec']
        d_smooth = derived['diag_vec']
        valid = (d_smooth > 0) & (d_smooth < N_BINS) & np.isfinite(y_smooth)
        if not np.any(valid):
            continue
        y_chunk = y_smooth[valid].astype(np.float32, copy=False)
        d_chunk = d_smooth[valid]
        for diag_val in np.unique(d_chunk):
            diag_mask = d_chunk == diag_val
            vals = y_chunk[diag_mask]
            prob = sampling_prob[diag_val]
            if prob < 1.0:
                vals = vals[rng.random(len(vals)) < prob]
            if len(vals) == 0:
                continue
            sampled_values[diag_val].append(vals)
            sampled_counts[diag_val] += len(vals)

    baseline = np.zeros(N_BINS, dtype=np.float32)
    for diag_val in range(1, N_BINS):
        if not sampled_values[diag_val]:
            continue
        vals = np.concatenate(sampled_values[diag_val]).astype(np.float32, copy=False)
        if len(vals) > MAX_VALS_PER_STRATUM:
            keep_idx = rng.choice(len(vals), MAX_VALS_PER_STRATUM, replace=False)
            vals = vals[keep_idx]
        baseline[diag_val] = np.float32(np.median(vals))
        sampled_counts[diag_val] = len(vals)

    np.savez(
        out_path,
        baseline=baseline,
        diag_counts=diag_counts,
        sampled_counts=sampled_counts,
        pipeline_version=np.array(PIPELINE_VERSION),
        mode=np.array(DISTANCE_BASELINE_MODE),
    )
    print(f"  Saved residual baseline → {out_path}")
    return baseline

def count_bucket_rows_by_chrom(h5_path: str, split: str) -> dict:
    chrom_counts = {}
    for chrom, _, _, _, diag_idx in iter_split_windows(h5_path, split):
        if chrom not in chrom_counts:
            chrom_counts[chrom] = np.zeros(N_DIAG_BUCKETS, dtype=np.int64)
        bucket_idx = assign_diag_buckets(diag_idx)
        valid = bucket_idx >= 0
        if np.any(valid):
            chrom_counts[chrom] += np.bincount(bucket_idx[valid], minlength=N_DIAG_BUCKETS)
    return chrom_counts

def sampled_h5_is_current(path: Path) -> bool:
    if not path.exists():
        return False
    try:
        with h5py.File(path, 'r') as hf:
            return (
                decode_attr(hf.attrs.get('pipeline_version', '')) == PIPELINE_VERSION and
                int(hf.attrs.get('n_features', -1)) == len(ALL_ACTIVE_FEATURE_NAMES) and
                decode_attr(hf.attrs.get('target_mode', '')) == 'smoothed_oe_residual'
            )
    except Exception:
        return False

def build_sampled_features_h5(source_h5: str, sampled_h5: Path, baseline: np.ndarray):
    if sampled_h5_is_current(sampled_h5):
        print(f"  Reusing balanced sampled dataset → {sampled_h5}")
        return

    print("  Building balanced sampled dataset from features.h5...")
    counts_by_split = {
        'train': count_bucket_rows_by_chrom(source_h5, 'train'),
        'val': count_bucket_rows_by_chrom(source_h5, 'val'),
    }
    split_caps = {
        'train': MAX_PAIRS_PER_BUCKET_PER_CHROM_TRAIN,
        'val': MAX_PAIRS_PER_BUCKET_PER_CHROM_VAL,
    }
    rng = np.random.default_rng(RANDOM_SEED)

    with h5py.File(source_h5, 'r') as src, h5py.File(sampled_h5, 'w') as dst:
        dst.attrs['pipeline_version'] = PIPELINE_VERSION
        dst.attrs['target_mode'] = 'smoothed_oe_residual'
        dst.attrs['n_features'] = len(ALL_ACTIVE_FEATURE_NAMES)
        dst.attrs['feature_names'] = np.array(ALL_ACTIVE_FEATURE_NAMES, dtype='S')
        dst.attrs['diag_bucket_edges'] = DIAG_BUCKET_EDGES
        dst.attrs['baseline_mode'] = DISTANCE_BASELINE_MODE
        dst.attrs['ml_task'] = ML_TASK

        for split, cap_per_bucket in split_caps.items():
            sg = dst.create_group(split)
            for chrom in sorted(src.keys()):
                cg = src[chrom]
                if decode_attr(cg.attrs['split']) != split:
                    continue

                chrom_counts = counts_by_split[split].get(chrom, np.zeros(N_DIAG_BUCKETS, dtype=np.int64))
                sample_prob = np.ones(N_DIAG_BUCKETS, dtype=np.float64)
                over_cap = chrom_counts > cap_per_bucket
                sample_prob[over_cap] = cap_per_bucket / chrom_counts[over_cap]

                bucket_store = [{'X': [], 'y_true_raw': [], 'y_true_smooth': [], 'diag': []} for _ in range(N_DIAG_BUCKETS)]

                for win_key in sorted(cg.keys()):
                    feats = cg[win_key]['features'][:].astype(np.float32, copy=False)
                    y_true_raw = cg[win_key]['targets'][:].astype(np.float32, copy=False)
                    diag_idx = cg[win_key]['diag_idx'][:].astype(np.int32)
                    derived = derive_window_regression_targets(y_true_raw, diag_idx)
                    y_true_raw = derived['raw_vec']
                    y_true_smooth = derived['smoothed_vec']
                    diag_idx = derived['diag_vec']
                    feats = feats[derived['valid_vec_mask']]
                    bucket_idx = assign_diag_buckets(diag_idx)
                    valid = (bucket_idx >= 0) & np.isfinite(y_true_raw) & np.isfinite(y_true_smooth)
                    if not np.any(valid):
                        continue

                    feats = feats[valid]
                    y_true_raw = y_true_raw[valid]
                    y_true_smooth = y_true_smooth[valid]
                    diag_idx = diag_idx[valid]
                    bucket_idx = bucket_idx[valid]

                    for bucket_id in range(N_DIAG_BUCKETS):
                        idx = np.flatnonzero(bucket_idx == bucket_id)
                        if len(idx) == 0:
                            continue
                        prob = sample_prob[bucket_id]
                        if prob < 1.0:
                            idx = idx[rng.random(len(idx)) < prob]
                        if len(idx) == 0:
                            continue
                        bucket_store[bucket_id]['X'].append(feats[idx])
                        bucket_store[bucket_id]['y_true_raw'].append(y_true_raw[idx])
                        bucket_store[bucket_id]['y_true_smooth'].append(y_true_smooth[idx])
                        bucket_store[bucket_id]['diag'].append(diag_idx[idx])

                chrom_parts = []
                for bucket_id in range(N_DIAG_BUCKETS):
                    if not bucket_store[bucket_id]['X']:
                        continue
                    X_bucket = np.concatenate(bucket_store[bucket_id]['X'], axis=0)
                    y_true_raw_bucket = np.concatenate(bucket_store[bucket_id]['y_true_raw'], axis=0)
                    y_true_smooth_bucket = np.concatenate(bucket_store[bucket_id]['y_true_smooth'], axis=0)
                    diag_bucket = np.concatenate(bucket_store[bucket_id]['diag'], axis=0)

                    if len(X_bucket) > cap_per_bucket:
                        keep_idx = rng.choice(len(X_bucket), cap_per_bucket, replace=False)
                        X_bucket = X_bucket[keep_idx]
                        y_true_raw_bucket = y_true_raw_bucket[keep_idx]
                        y_true_smooth_bucket = y_true_smooth_bucket[keep_idx]
                        diag_bucket = diag_bucket[keep_idx]

                    if ML_TASK == 'regression':
                        y_train_bucket = y_true_smooth_bucket - baseline[diag_bucket]
                    else:
                        y_train_bucket = y_true_smooth_bucket

                    chrom_parts.append((
                        X_bucket.astype(np.float32, copy=False),
                        y_train_bucket.astype(np.float32, copy=False),
                        y_true_raw_bucket.astype(np.float32, copy=False),
                        y_true_smooth_bucket.astype(np.float32, copy=False),
                        diag_bucket.astype(np.int32, copy=False),
                        np.full(len(diag_bucket), bucket_id, dtype=np.int16),
                    ))

                if not chrom_parts:
                    continue

                X_chrom = np.concatenate([part[0] for part in chrom_parts], axis=0)
                y_train_chrom = np.concatenate([part[1] for part in chrom_parts], axis=0)
                y_true_raw_chrom = np.concatenate([part[2] for part in chrom_parts], axis=0)
                y_true_smooth_chrom = np.concatenate([part[3] for part in chrom_parts], axis=0)
                diag_chrom = np.concatenate([part[4] for part in chrom_parts], axis=0)
                bucket_chrom = np.concatenate([part[5] for part in chrom_parts], axis=0)
                order = rng.permutation(len(X_chrom))

                chrom_group = sg.create_group(chrom)
                chrom_group.create_dataset('features', data=X_chrom[order], compression='gzip', compression_opts=4)
                chrom_group.create_dataset('targets_train', data=y_train_chrom[order], compression='gzip', compression_opts=4)
                chrom_group.create_dataset('targets_true_raw', data=y_true_raw_chrom[order], compression='gzip', compression_opts=4)
                chrom_group.create_dataset('targets_true_smoothed', data=y_true_smooth_chrom[order], compression='gzip', compression_opts=4)
                chrom_group.create_dataset('diag_idx', data=diag_chrom[order], compression='gzip', compression_opts=4)
                chrom_group.create_dataset('bucket_idx', data=bucket_chrom[order], compression='gzip', compression_opts=4)

    print(f"  Saved balanced sampled dataset → {sampled_h5}")


def summarize_sampled_dataset(sampled_h5: str, out_csv: Path, plot_dir: Path):
    records = []
    with h5py.File(sampled_h5, 'r') as hf:
        for split in sorted(hf.keys()):
            for chrom in sorted(hf[split].keys()):
                cg = hf[split][chrom]
                buckets = cg['bucket_idx'][:].astype(np.int16)
                for bucket_id in range(N_DIAG_BUCKETS):
                    mask = buckets == bucket_id
                    if not np.any(mask):
                        continue
                    records.append({
                        'split': split,
                        'chrom': chrom,
                        'bucket': DIAG_BUCKET_LABELS[bucket_id],
                        'bucket_id': bucket_id,
                        'n_total': int(mask.sum()),
                    })
    sample_df = pd.DataFrame(records)
    sample_df.to_csv(out_csv, index=False)
    if not sample_df.empty:
        agg = ordered_bucket_frame(sample_df.groupby(['split', 'bucket'])[['n_total']].sum().reset_index())
        for split in ['train', 'val']:
            sub = ordered_bucket_frame(agg[agg['split'] == split])
            if sub.empty:
                continue
            fig, ax = plt.subplots(figsize=(10, 4))
            ax.bar(sub['bucket'], sub['n_total'], color='slateblue')
            ax.set_title(f"Sampled Regression Dataset Size ({split})")
            ax.set_xlabel("Distance bucket")
            ax.set_ylabel("Count")
            plt.xticks(rotation=45)
            plt.tight_layout()
            plt.savefig(ensure_parent(plot_dir / 'sampling_diagnostics' / f'sampled_balance_{split}.png'))
            plt.close()
    return sample_df

def load_sampled_split(h5_path: str, split: str, feature_variant: str):
    X_list, y_train_list, y_true_raw_list, y_true_smooth_list, diag_list, bucket_list, group_list = [], [], [], [], [], [], []
    with h5py.File(h5_path, 'r') as hf:
        split_group = hf[split]
        for group_id, chrom in enumerate(sorted(split_group.keys())):
            cg = split_group[chrom]
            X = select_feature_matrix(cg['features'][:], feature_variant)
            y_train = cg['targets_train'][:]
            y_true_raw = cg['targets_true_raw'][:]
            y_true_smooth = cg['targets_true_smoothed'][:]
            diag_idx = cg['diag_idx'][:].astype(np.int32)
            bucket_idx = cg['bucket_idx'][:].astype(np.int16)

            X_list.append(X)
            y_train_list.append(y_train)
            y_true_raw_list.append(y_true_raw)
            y_true_smooth_list.append(y_true_smooth)
            diag_list.append(diag_idx)
            bucket_list.append(bucket_idx)
            group_list.append(np.full(len(X), group_id, dtype=np.int32))

    return (
        np.vstack(X_list),
        np.concatenate(y_train_list),
        np.concatenate(y_true_raw_list),
        np.concatenate(y_true_smooth_list),
        np.concatenate(diag_list),
        np.concatenate(bucket_list),
        np.concatenate(group_list),
    )

def make_model(random_state: int, best_params: dict):
    params = dict(best_params)
    params['random_state'] = random_state
    params['n_jobs'] = RF_TRAIN_N_JOBS
    if ML_TASK == 'regression':
        return RandomForestRegressor(**params)
    return RandomForestClassifier(**params)

def model_meta_is_current(variant_paths: dict, feature_variant: str) -> bool:
    if not variant_paths['model_meta'].exists():
        return False
    try:
        with open(variant_paths['model_meta']) as fh:
            meta = json.load(fh)
        return (
            meta.get('pipeline_version') == PIPELINE_VERSION and
            meta.get('ml_task') == ML_TASK and
            meta.get('baseline_mode') == DISTANCE_BASELINE_MODE and
            meta.get('feature_variant') == feature_variant and
            meta.get('target_mode') == 'smoothed_oe_residual'
        )
    except Exception:
        return False

def write_model_meta(variant_paths: dict, best_params: dict, feature_variant: str):
    feature_names = feature_names_for_variant(feature_variant)
    meta = {
        'pipeline_version': PIPELINE_VERSION,
        'ml_task': ML_TASK,
        'baseline_mode': DISTANCE_BASELINE_MODE,
        'target_mode': 'smoothed_oe_residual',
        'feature_variant': feature_variant,
        'n_active_features': len(feature_names),
        'feature_names': feature_names,
        'diag_bucket_edges': DIAG_BUCKET_EDGES.tolist(),
        'rf_train_n_jobs': RF_TRAIN_N_JOBS,
        'best_params': best_params,
    }
    with open(variant_paths['model_meta'], 'w') as fh:
        json.dump(meta, fh, indent=2)

def init_regression_accumulator():
    return {
        'n': 0,
        'sum_y': 0.0,
        'sum_pred': 0.0,
        'sum_y2': 0.0,
        'sum_pred2': 0.0,
        'sum_yp': 0.0,
        'sum_sqerr': 0.0,
    }

def update_regression_accumulator(acc: dict, y_true: np.ndarray, y_pred: np.ndarray):
    if len(y_true) == 0:
        return
    y_true = y_true.astype(np.float64, copy=False)
    y_pred = y_pred.astype(np.float64, copy=False)
    acc['n'] += len(y_true)
    acc['sum_y'] += y_true.sum()
    acc['sum_pred'] += y_pred.sum()
    acc['sum_y2'] += np.square(y_true).sum()
    acc['sum_pred2'] += np.square(y_pred).sum()
    acc['sum_yp'] += (y_true * y_pred).sum()
    acc['sum_sqerr'] += np.square(y_true - y_pred).sum()

def finalize_regression_accumulator(acc: dict) -> dict:
    n = acc['n']
    if n == 0:
        return {'n_pairs': 0, 'MSE': np.nan, 'Pearson_R': np.nan, 'Stratum_Adjusted_R': np.nan}
    denom_left = n * acc['sum_y2'] - acc['sum_y'] ** 2
    denom_right = n * acc['sum_pred2'] - acc['sum_pred'] ** 2
    denom = np.sqrt(max(denom_left, 0.0) * max(denom_right, 0.0))
    pearson_r = np.nan if denom == 0 else (n * acc['sum_yp'] - acc['sum_y'] * acc['sum_pred']) / denom
    mse = acc['sum_sqerr'] / n
    return {'n_pairs': int(n), 'MSE': float(mse), 'Pearson_R': float(pearson_r)}


def init_stratum_adjusted_accumulator() -> dict:
    return {
        'count': np.zeros(N_BINS, dtype=np.int64),
        'sum_y': np.zeros(N_BINS, dtype=np.float64),
        'sum_pred': np.zeros(N_BINS, dtype=np.float64),
        'sum_y2': np.zeros(N_BINS, dtype=np.float64),
        'sum_pred2': np.zeros(N_BINS, dtype=np.float64),
        'sum_yp': np.zeros(N_BINS, dtype=np.float64),
    }


def update_stratum_adjusted_accumulator(acc: dict, diag_idx: np.ndarray, y_true: np.ndarray, y_pred: np.ndarray):
    if len(diag_idx) == 0:
        return
    d = diag_idx.astype(np.int32, copy=False)
    y_true = y_true.astype(np.float64, copy=False)
    y_pred = y_pred.astype(np.float64, copy=False)
    acc['count'] += np.bincount(d, minlength=N_BINS)[:N_BINS]
    acc['sum_y'] += np.bincount(d, weights=y_true, minlength=N_BINS)[:N_BINS]
    acc['sum_pred'] += np.bincount(d, weights=y_pred, minlength=N_BINS)[:N_BINS]
    acc['sum_y2'] += np.bincount(d, weights=np.square(y_true), minlength=N_BINS)[:N_BINS]
    acc['sum_pred2'] += np.bincount(d, weights=np.square(y_pred), minlength=N_BINS)[:N_BINS]
    acc['sum_yp'] += np.bincount(d, weights=y_true * y_pred, minlength=N_BINS)[:N_BINS]


def finalize_stratum_adjusted_accumulator(acc: dict) -> float:
    total_cov = 0.0
    total_var_y = 0.0
    total_var_pred = 0.0
    for diag_val in range(1, N_BINS):
        n = acc['count'][diag_val]
        if n <= 1:
            continue
        mean_y = acc['sum_y'][diag_val] / n
        mean_pred = acc['sum_pred'][diag_val] / n
        cov = acc['sum_yp'][diag_val] - n * mean_y * mean_pred
        var_y = acc['sum_y2'][diag_val] - n * mean_y * mean_y
        var_pred = acc['sum_pred2'][diag_val] - n * mean_pred * mean_pred
        total_cov += cov
        total_var_y += max(var_y, 0.0)
        total_var_pred += max(var_pred, 0.0)
    denom = np.sqrt(total_var_y * total_var_pred)
    if denom == 0:
        return np.nan
    return float(total_cov / denom)

def predict_regression_in_batches(model, X: np.ndarray) -> np.ndarray:
    preds = np.zeros(len(X), dtype=np.float32)
    for start_idx in range(0, len(X), PREDICT_BATCH_SIZE):
        preds[start_idx:start_idx + PREDICT_BATCH_SIZE] = model.predict(X[start_idx:start_idx + PREDICT_BATCH_SIZE])
    return preds

def build_window_pair_features(row):
    raw_mat = load_microc_window(row.chrom, row.start, row.end)
    oe_mat = normalize_microc(raw_mat)

    tracks = np.zeros((N_MARKS, N_BINS), dtype=np.float32)
    tracks_smooth = np.zeros((N_MARKS, N_BINS), dtype=np.float32)
    for mi, (mark, bw_path) in enumerate(TRACK_FILES.items()):
        t = normalize_chipseq(load_bigwig_window(bw_path, row.chrom, row.start, row.end, N_BINS), mark)
        tracks[mi] = t
        tracks_smooth[mi] = gaussian_filter1d(t, sigma=2.0)

    cumsum_ext = np.zeros((N_MARKS, N_BINS + 1), dtype=np.float32)
    cumsum_ext[:, 1:] = np.cumsum(tracks, axis=1)

    feats = np.empty((len(_ii), N_FEATURES), dtype=np.float32)
    feats[:, :N_MARKS] = tracks[:, _ii].T
    feats[:, N_MARKS:2 * N_MARKS] = tracks[:, _jj].T
    feats[:, 2 * N_MARKS:3 * N_MARKS] = (tracks[:, _ii] * tracks[:, _jj]).T

    span = (_jj - _ii + 1).astype(np.float32)
    int_sum = cumsum_ext[:, _jj + 1] - cumsum_ext[:, _ii]
    feats[:, 3 * N_MARKS:4 * N_MARKS] = (int_sum / np.where(span > 0, span, 1.0)).T
    feats[:, 4 * N_MARKS:5 * N_MARKS] = (tracks_smooth[:, _ii] * tracks_smooth[:, _jj]).T

    diag_idx = (_jj - _ii).astype(np.int32)
    feats[:, 5 * N_MARKS] = np.log1p(diag_idx.astype(np.float32) * REQUESTED_RESOLUTION)
    feats[:, 5 * N_MARKS + 1] = diag_idx.astype(np.float32)
    valid_mask = np.isfinite(oe_mat[_ii, _jj]) & np.all(np.isfinite(feats), axis=1)
    return feats, oe_mat, diag_idx, valid_mask

def smooth_exact_diag_values(values: np.ndarray, counts: np.ndarray, sigma: float = EXACT_DIAG_CORRECTION_SIGMA) -> np.ndarray:
    out = np.zeros_like(values, dtype=np.float32)
    valid = counts > 0
    if not np.any(valid):
        return out
    idx = np.arange(len(values), dtype=np.float32)
    interp_vals = np.interp(idx, idx[valid], values[valid]).astype(np.float32)
    smoothed = gaussian_filter1d(interp_vals, sigma=sigma, mode='nearest').astype(np.float32)
    out[valid] = smoothed[valid]
    if np.any(~valid):
        out[~valid] = smoothed[~valid]
    out[0] = 0.0
    return out


def compute_exact_diag_correction(model, baseline: np.ndarray, feature_variant: str) -> tuple[np.ndarray, pd.DataFrame]:
    sums = np.zeros(N_BINS, dtype=np.float64)
    counts = np.zeros(N_BINS, dtype=np.int64)
    for _, _, feats, targets, diag_idx in iter_split_windows(str(CKPT['features_h5']), 'val'):
        derived = derive_window_regression_targets(targets, diag_idx)
        y_true_smooth = derived['smoothed_vec']
        d = derived['diag_vec']
        valid = (d > 0) & (d < N_BINS)
        if not np.any(valid):
            continue
        X = select_feature_matrix(feats[derived['valid_vec_mask']][valid].astype(np.float32, copy=False), feature_variant)
        y_pred = predict_regression_in_batches(model, X) + baseline[d[valid]]
        errors = y_pred - y_true_smooth[valid]
        diag_valid = d[valid]
        sums += np.bincount(diag_valid, weights=errors.astype(np.float64), minlength=N_BINS)[:N_BINS]
        counts += np.bincount(diag_valid, minlength=N_BINS)[:N_BINS]
    raw = np.zeros(N_BINS, dtype=np.float32)
    valid_counts = counts > 0
    raw[valid_counts] = (sums[valid_counts] / counts[valid_counts]).astype(np.float32)
    smooth = smooth_exact_diag_values(raw, counts)
    rows = []
    for diag_val in range(1, N_BINS):
        if counts[diag_val] == 0:
            continue
        rows.append({
            'diag_idx': int(diag_val),
            'target_space': PRIMARY_REGRESSION_TARGET_SPACE,
            'n_pairs': int(counts[diag_val]),
            'mean_error_smoothed_oe': float(raw[diag_val]),
            'smoothed_correction': float(smooth[diag_val]),
        })
    return smooth, pd.DataFrame(rows)


def smooth_pred_level_values(values: np.ndarray, sigma: float = PRED_LEVEL_CORRECTION_SIGMA) -> np.ndarray:
    values = np.asarray(values, dtype=np.float32)
    if values.size == 0:
        return values
    return gaussian_filter1d(values, sigma=sigma, mode='nearest').astype(np.float32)


def compute_pred_level_correction(
    model,
    baseline: np.ndarray,
    feature_variant: str,
    exact_diag_correction: np.ndarray,
) -> tuple[dict | None, pd.DataFrame]:
    pred_chunks = []
    resid_chunks = []
    for _, _, feats, targets, diag_idx in iter_split_windows(str(CKPT['features_h5']), 'val'):
        derived = derive_window_regression_targets(targets, diag_idx)
        y_true_smooth = derived['smoothed_vec']
        d = derived['diag_vec']
        valid = (d > 0) & (d < N_BINS)
        if not np.any(valid):
            continue
        X = select_feature_matrix(feats[derived['valid_vec_mask']][valid].astype(np.float32, copy=False), feature_variant)
        y_pred = predict_regression_in_batches(model, X) + baseline[d[valid]]
        y_pred = y_pred - correction_values_from_diag(d[valid], exact_diag_correction)
        pred_chunks.append(y_pred.astype(np.float32, copy=False))
        resid_chunks.append((y_pred - y_true_smooth[valid]).astype(np.float32, copy=False))

    empty_df = pd.DataFrame(columns=[
        'target_space', 'pred_bin_left', 'pred_bin_right', 'pred_center',
        'n_pairs', 'mean_error_smoothed_oe', 'smoothed_correction',
    ])
    if not pred_chunks:
        return None, empty_df

    pred_all = np.concatenate(pred_chunks)
    resid_all = np.concatenate(resid_chunks)
    if pred_all.size == 0:
        return None, empty_df

    quantiles = np.linspace(0.0, 1.0, PRED_LEVEL_CORRECTION_BINS + 1, dtype=np.float64)
    edges = np.quantile(pred_all.astype(np.float64), quantiles)
    edges = np.unique(edges)
    if edges.size < 3:
        min_pred = float(np.min(pred_all))
        max_pred = float(np.max(pred_all))
        if np.isclose(min_pred, max_pred):
            return None, empty_df
        edges = np.linspace(min_pred, max_pred, min(5, PRED_LEVEL_CORRECTION_BINS) + 1, dtype=np.float64)

    bin_idx = np.digitize(pred_all.astype(np.float64), edges[1:-1], right=False)
    raw_centers = []
    raw_errors = []
    rows = []
    for idx in range(len(edges) - 1):
        mask = bin_idx == idx
        if not np.any(mask):
            continue
        pred_bin = pred_all[mask]
        resid_bin = resid_all[mask]
        center = float(np.mean(pred_bin))
        raw_error = float(np.mean(resid_bin))
        raw_centers.append(center)
        raw_errors.append(raw_error)
        rows.append({
            'target_space': PRIMARY_REGRESSION_TARGET_SPACE,
            'pred_bin_left': float(edges[idx]),
            'pred_bin_right': float(edges[idx + 1]),
            'pred_center': center,
            'n_pairs': int(mask.sum()),
            'mean_error_smoothed_oe': raw_error,
            'smoothed_correction': np.nan,
        })

    centers = np.asarray(raw_centers, dtype=np.float32)
    errors = np.asarray(raw_errors, dtype=np.float32)
    order = np.argsort(centers)
    centers = centers[order]
    errors = errors[order]
    smoothed = smooth_pred_level_values(errors)

    rows_sorted = []
    for idx, row_idx in enumerate(order):
        row = rows[row_idx]
        row['smoothed_correction'] = float(smoothed[idx])
        rows_sorted.append(row)

    return {
        'pred_center': centers,
        'smoothed_correction': smoothed,
    }, pd.DataFrame(rows_sorted)


def correction_values_from_diag(diag_idx: np.ndarray, correction) -> np.ndarray:
    if correction is None:
        return np.zeros(len(diag_idx), dtype=np.float32)
    if isinstance(correction, np.ndarray):
        clipped = np.clip(diag_idx.astype(np.int32), 0, len(correction) - 1)
        return correction[clipped].astype(np.float32, copy=False)
    bucket_idx = assign_diag_buckets(diag_idx)
    return np.array([correction.get(int(b), 0.0) for b in bucket_idx], dtype=np.float32)


def correction_values_from_pred_level(pred_values: np.ndarray, correction: dict | None) -> np.ndarray:
    if correction is None:
        return np.zeros(len(pred_values), dtype=np.float32)
    centers = np.asarray(correction.get('pred_center', []), dtype=np.float32)
    smoothed = np.asarray(correction.get('smoothed_correction', []), dtype=np.float32)
    if centers.size == 0 or smoothed.size == 0:
        return np.zeros(len(pred_values), dtype=np.float32)
    return np.interp(
        pred_values.astype(np.float64),
        centers.astype(np.float64),
        smoothed.astype(np.float64),
        left=float(smoothed[0]),
        right=float(smoothed[-1]),
    ).astype(np.float32)


def apply_prediction_correction(pred_values: np.ndarray, diag_idx: np.ndarray, correction) -> np.ndarray:
    out = pred_values.astype(np.float32, copy=False)
    if correction is None:
        return out
    if isinstance(correction, np.ndarray):
        return out - correction_values_from_diag(diag_idx, correction)
    exact_diag = correction.get('exact_diag') if isinstance(correction, dict) else None
    pred_level = correction.get('pred_level') if isinstance(correction, dict) else None
    if exact_diag is not None:
        out = out - correction_values_from_diag(diag_idx, exact_diag)
    if pred_level is not None:
        out = out - correction_values_from_pred_level(out, pred_level)
    return out


def compute_exact_diag_regression_report(
    split: str,
    diag_idx: np.ndarray,
    y_true_raw: np.ndarray,
    y_true_smoothed: np.ndarray,
    y_pred: np.ndarray,
) -> pd.DataFrame:
    rows = []
    for diag_val in np.unique(diag_idx):
        if diag_val <= 0 or diag_val >= N_BINS:
            continue
        mask = diag_idx == diag_val
        if not np.any(mask):
            continue
        rows.append({
            'split': split,
            'diag_idx': int(diag_val),
            'n_pairs': int(mask.sum()),
            'mean_true_raw_oe': float(np.mean(y_true_raw[mask])),
            'mean_true_smoothed_oe': float(np.mean(y_true_smoothed[mask])),
            'mean_pred_oe': float(np.mean(y_pred[mask])),
            'pred_minus_true_raw': float(np.mean(y_pred[mask] - y_true_raw[mask])),
            'pred_minus_true_smoothed': float(np.mean(y_pred[mask] - y_true_smoothed[mask])),
        })
    return pd.DataFrame(rows)


def evaluate_streamed_regression(model, baseline: np.ndarray, feature_variant: str, split: str, correction: dict | None = None):
    print(f"  Streaming {split} evaluation from features.h5 ({feature_variant})...")
    accs = {
        'raw_oe': init_regression_accumulator(),
        'smoothed_oe': init_regression_accumulator(),
        'raw_resid': init_regression_accumulator(),
        'smoothed_resid': init_regression_accumulator(),
    }
    stratum_accs = {name: init_stratum_adjusted_accumulator() for name in accs.keys()}
    bucket_accs = {
        name: [init_regression_accumulator() for _ in range(N_DIAG_BUCKETS)]
        for name in accs.keys()
    }
    bias_acc = {
        bucket_id: {
            'sum_pred': 0.0,
            'sum_true_raw': 0.0,
            'sum_true_smooth': 0.0,
            'sum_diff_raw': 0.0,
            'sum_diff_smooth': 0.0,
            'n': 0,
        }
        for bucket_id in range(N_DIAG_BUCKETS)
    }
    exact_diag_acc = {
        'sum_pred': np.zeros(N_BINS, dtype=np.float64),
        'sum_true_raw': np.zeros(N_BINS, dtype=np.float64),
        'sum_true_smooth': np.zeros(N_BINS, dtype=np.float64),
        'count': np.zeros(N_BINS, dtype=np.int64),
    }

    for _, _, feats, targets, diag_idx in iter_split_windows(str(CKPT['features_h5']), split):
        derived = derive_window_regression_targets(targets, diag_idx)
        y_true_raw = derived['raw_vec']
        y_true_smooth = derived['smoothed_vec']
        d = derived['diag_vec']
        bucket_idx = assign_diag_buckets(d)
        valid = bucket_idx >= 0
        if not np.any(valid):
            continue

        feats_valid = feats[derived['valid_vec_mask']][valid].astype(np.float32, copy=False)
        X = select_feature_matrix(feats_valid, feature_variant)
        y_pred_resid = predict_regression_in_batches(model, X)
        y_pred = y_pred_resid + baseline[d[valid]]
        if correction is not None:
            y_pred = apply_prediction_correction(y_pred, d[valid], correction)
            y_pred_resid = y_pred - baseline[d[valid]]

        y_raw = y_true_raw[valid]
        y_smooth = y_true_smooth[valid]
        d_valid = d[valid]
        b_valid = bucket_idx[valid]
        exact_diag_acc['sum_pred'] += np.bincount(d_valid, weights=y_pred.astype(np.float64), minlength=N_BINS)[:N_BINS]
        exact_diag_acc['sum_true_raw'] += np.bincount(d_valid, weights=y_raw.astype(np.float64), minlength=N_BINS)[:N_BINS]
        exact_diag_acc['sum_true_smooth'] += np.bincount(d_valid, weights=y_smooth.astype(np.float64), minlength=N_BINS)[:N_BINS]
        exact_diag_acc['count'] += np.bincount(d_valid, minlength=N_BINS)[:N_BINS]

        update_regression_accumulator(accs['raw_oe'], y_raw, y_pred)
        update_regression_accumulator(accs['smoothed_oe'], y_smooth, y_pred)
        update_regression_accumulator(accs['raw_resid'], y_raw - baseline[d_valid], y_pred_resid)
        update_regression_accumulator(accs['smoothed_resid'], y_smooth - baseline[d_valid], y_pred_resid)
        update_stratum_adjusted_accumulator(stratum_accs['raw_oe'], d_valid, y_raw, y_pred)
        update_stratum_adjusted_accumulator(stratum_accs['smoothed_oe'], d_valid, y_smooth, y_pred)
        update_stratum_adjusted_accumulator(stratum_accs['raw_resid'], d_valid, y_raw - baseline[d_valid], y_pred_resid)
        update_stratum_adjusted_accumulator(stratum_accs['smoothed_resid'], d_valid, y_smooth - baseline[d_valid], y_pred_resid)

        for bucket_id, label in enumerate(DIAG_BUCKET_LABELS):
            mask = b_valid == bucket_id
            if not np.any(mask):
                continue
            update_regression_accumulator(bucket_accs['raw_oe'][bucket_id], y_raw[mask], y_pred[mask])
            update_regression_accumulator(bucket_accs['smoothed_oe'][bucket_id], y_smooth[mask], y_pred[mask])
            update_regression_accumulator(bucket_accs['raw_resid'][bucket_id], y_raw[mask] - baseline[d_valid[mask]], y_pred_resid[mask])
            update_regression_accumulator(bucket_accs['smoothed_resid'][bucket_id], y_smooth[mask] - baseline[d_valid[mask]], y_pred_resid[mask])
            bias_acc[bucket_id]['sum_pred'] += float(np.sum(y_pred[mask]))
            bias_acc[bucket_id]['sum_true_raw'] += float(np.sum(y_raw[mask]))
            bias_acc[bucket_id]['sum_true_smooth'] += float(np.sum(y_smooth[mask]))
            bias_acc[bucket_id]['sum_diff_raw'] += float(np.sum(y_pred[mask] - y_raw[mask]))
            bias_acc[bucket_id]['sum_diff_smooth'] += float(np.sum(y_pred[mask] - y_smooth[mask]))
            bias_acc[bucket_id]['n'] += int(mask.sum())

    records = []
    for target_space, acc in accs.items():
        overall_metrics = finalize_regression_accumulator(acc)
        overall_metrics['Stratum_Adjusted_R'] = finalize_stratum_adjusted_accumulator(stratum_accs[target_space])
        records.append({'split': split, 'target_space': target_space, 'bucket': 'overall', **overall_metrics})
        for bucket_id, label in enumerate(DIAG_BUCKET_LABELS):
            records.append({
                'split': split,
                'target_space': target_space,
                'bucket': label,
                **finalize_regression_accumulator(bucket_accs[target_space][bucket_id]),
            })
    bias_records = []
    for bucket_id, label in enumerate(DIAG_BUCKET_LABELS):
        n = bias_acc[bucket_id]['n']
        if n == 0:
            continue
        bias_records.append({
            'split': split,
            'bucket': label,
            'mean_pred_oe': bias_acc[bucket_id]['sum_pred'] / n,
            'mean_true_raw_oe': bias_acc[bucket_id]['sum_true_raw'] / n,
            'mean_true_smoothed_oe': bias_acc[bucket_id]['sum_true_smooth'] / n,
            'pred_minus_true_raw': bias_acc[bucket_id]['sum_diff_raw'] / n,
            'pred_minus_true_smoothed': bias_acc[bucket_id]['sum_diff_smooth'] / n,
            'n_pairs': n,
        })
    exact_diag_records = []
    for diag_val in range(1, N_BINS):
        n = exact_diag_acc['count'][diag_val]
        if n == 0:
            continue
        mean_pred = exact_diag_acc['sum_pred'][diag_val] / n
        mean_true_raw = exact_diag_acc['sum_true_raw'][diag_val] / n
        mean_true_smooth = exact_diag_acc['sum_true_smooth'][diag_val] / n
        exact_diag_records.append({
            'split': split,
            'diag_idx': int(diag_val),
            'n_pairs': int(n),
            'mean_pred_oe': mean_pred,
            'mean_true_raw_oe': mean_true_raw,
            'mean_true_smoothed_oe': mean_true_smooth,
            'pred_minus_true_raw': mean_pred - mean_true_raw,
            'pred_minus_true_smoothed': mean_pred - mean_true_smooth,
        })
    return (
        ordered_bucket_frame(pd.DataFrame(records)),
        ordered_bucket_frame(pd.DataFrame(bias_records)),
        pd.DataFrame(exact_diag_records),
    )


def generate_prediction_heatmaps(
    model,
    baseline: np.ndarray,
    feature_variant: str,
    variant_paths: dict,
    correction: dict | None = None,
    model_feature_variant: str | None = None,
):
    print(f"  Generating prediction heatmaps for tracked validation/test windows ({feature_variant})...")
    model_feature_variant = model_feature_variant or feature_variant
    output_records = []
    for _, row in prediction_windows_eval.iterrows():
        try:
            feats, oe_mat, diag_idx, valid_mask = build_window_pair_features(row)
            derived = derive_window_regression_targets(oe_mat[_ii[valid_mask], _jj[valid_mask]], diag_idx[valid_mask])
            raw_oe = derived['raw_oe_mat']
            smooth_oe = derived['smoothed_oe_mat']
            diag_vec = derived['diag_vec']
            feats_valid = feats[valid_mask][derived['valid_vec_mask']].astype(np.float32, copy=False)
            X = select_feature_matrix(feats_valid, model_feature_variant)
            y_pred_resid = predict_regression_in_batches(model, X)
            y_pred_flat = y_pred_resid + baseline[diag_vec]
            correction_applied = correction is not None
            if correction is not None:
                y_pred_flat = apply_prediction_correction(y_pred_flat, diag_vec, correction)
                y_pred_resid = y_pred_flat - baseline[diag_vec]

            mat_pred = upper_vector_to_symmetric_matrix(y_pred_flat.astype(np.float32))
            pred_minus_raw = mat_pred - raw_oe
            pred_minus_smooth = mat_pred - smooth_oe
            true_resid_raw = raw_oe - baseline[DIAG_MATRIX_FULL]
            true_resid_smoothed = smooth_oe - baseline[DIAG_MATRIX_FULL]
            diag_profile_raw = np.full(N_BINS, np.nan, dtype=np.float32)
            diag_profile_smooth = np.full(N_BINS, np.nan, dtype=np.float32)
            diag_profile_pred = np.full(N_BINS, np.nan, dtype=np.float32)
            diag_profile_bias_smoothed = np.full(N_BINS, np.nan, dtype=np.float32)
            for diag_val in range(1, N_BINS):
                diag_mask = DIAG_MATRIX_FULL == diag_val
                raw_vals = raw_oe[diag_mask]
                smooth_vals = smooth_oe[diag_mask]
                pred_vals = mat_pred[diag_mask]
                raw_vals = raw_vals[np.isfinite(raw_vals)]
                smooth_vals = smooth_vals[np.isfinite(smooth_vals)]
                pred_vals = pred_vals[np.isfinite(pred_vals)]
                if raw_vals.size:
                    diag_profile_raw[diag_val] = float(np.mean(raw_vals))
                if smooth_vals.size:
                    diag_profile_smooth[diag_val] = float(np.mean(smooth_vals))
                if pred_vals.size:
                    diag_profile_pred[diag_val] = float(np.mean(pred_vals))
                if smooth_vals.size and pred_vals.size:
                    diag_profile_bias_smoothed[diag_val] = diag_profile_pred[diag_val] - diag_profile_smooth[diag_val]

            stem = row.track_id
            np.savez_compressed(
                ensure_parent(variant_paths['window_dir'] / f"{stem}.npz"),
                true_oe_raw=raw_oe,
                true_oe_smoothed=smooth_oe,
                true_resid_raw=true_resid_raw,
                true_resid_smoothed=true_resid_smoothed,
                pred_oe=mat_pred,
                pred_minus_true_raw=pred_minus_raw,
                pred_minus_true_smoothed=pred_minus_smooth,
                diag_profile_raw_oe=diag_profile_raw,
                diag_profile_smoothed_oe=diag_profile_smooth,
                diag_profile_pred_oe=diag_profile_pred,
                diag_profile_bias_smoothed=diag_profile_bias_smoothed,
                correction_applied=np.array(correction_applied),
                feature_variant=np.array(feature_variant),
                model_feature_variant=np.array(model_feature_variant),
                split=np.array(row.split),
                track_id=np.array(row.track_id),
                chrom=np.array(row.chrom),
                start=np.array(int(row.start), dtype=np.int64),
                end=np.array(int(row.end), dtype=np.int64),
            )

            fig, axs = plt.subplots(1, 4, figsize=(17, 4.6), constrained_layout=True)
            panels = [
                (raw_oe, "Raw O/E", 'RdBu_r', TARGET_MIN, TARGET_MAX),
                (smooth_oe, "Smoothed O/E", 'RdBu_r', TARGET_MIN, TARGET_MAX),
                (mat_pred, f"Predicted O/E ({'calibrated' if correction_applied else 'uncalibrated'})", 'RdBu_r', TARGET_MIN, TARGET_MAX),
                (pred_minus_smooth, "Pred - Smoothed O/E", 'RdBu_r', -2, 2),
            ]
            fig.suptitle(f"{row.track_id} | {row.split} | {row.chrom}:{int(row.start):,}-{int(row.end):,} | {feature_variant}", y=1.03)
            for ax, (mat_plot, title, cmap, vmin, vmax) in zip(np.ravel(axs), panels):
                im = ax.imshow(mat_plot, cmap=cmap, vmin=vmin, vmax=vmax)
                ax.set_title(title)
                ax.set_xticks([])
                ax.set_yticks([])
                fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

            plt.savefig(ensure_parent(variant_paths['prediction_plot_dir'] / f'{stem}.png'))
            plt.close()
            fig, ax = plt.subplots(figsize=(10, 4))
            diag_axis = np.arange(1, N_BINS)
            ax.plot(diag_axis, diag_profile_raw[1:], label='Raw O/E', linewidth=1.5)
            ax.plot(diag_axis, diag_profile_smooth[1:], label='Smoothed O/E', linewidth=1.5)
            ax.plot(diag_axis, diag_profile_pred[1:], label='Predicted O/E', linewidth=1.5)
            ax.plot(diag_axis, diag_profile_bias_smoothed[1:], label='Pred - Smoothed', linewidth=1.0, linestyle='--')
            for edge in DIAG_BUCKET_EDGES[1:-1]:
                ax.axvline(edge, color='lightgray', linestyle=':', linewidth=0.8)
            ax.set_title(f"Diagonal Profiles vs Smoothed O/E: {stem} ({'calibrated' if correction_applied else 'uncalibrated'})")
            ax.set_xlabel("|i-j| bins")
            ax.set_ylabel("Mean O/E")
            ax.legend(frameon=False, ncol=2)
            plt.tight_layout()
            plt.savefig(ensure_parent(variant_paths['prediction_plot_dir'] / f'{stem}_diag_profile.png'))
            plt.close()
            output_records.append({
                'track_index': int(row.track_index),
                'track_id': row.track_id,
                'split': row.split,
                'chrom': row.chrom,
                'start': int(row.start),
                'end': int(row.end),
                'correction_applied': correction_applied,
                'stage6_plot': str(variant_paths['prediction_plot_dir'] / f'{stem}.png'),
                'diag_profile_plot': str(variant_paths['prediction_plot_dir'] / f'{stem}_diag_profile.png'),
                'stage6_npz': str(variant_paths['window_dir'] / f'{stem}.npz'),
            })
        except Exception as e:
            logger.warning(f"Failed to plot prediction heatmap for {row.track_id}: {e}")
    if output_records:
        pd.DataFrame(output_records).sort_values('track_index').to_csv(variant_paths['tracking_outputs'], index=False)

def compute_grouped_permutation_importance(model, X_val, y_true_val, diag_val, baseline, feature_variant: str):
    base_pred_resid = predict_regression_in_batches(model, X_val)
    base_score, _ = pearsonr(y_true_val, base_pred_resid + baseline[diag_val])

    rng = np.random.default_rng(RANDOM_SEED)
    feature_groups = feature_groups_for_variant(feature_variant)
    results = []
    for mark, group_cols in feature_groups.items():
        perm_scores = []
        for _ in range(N_PERMUTATION_REPEATS):
            perm = rng.permutation(len(X_val))
            X_perm = X_val.copy()
            X_perm[:, group_cols] = X_perm[perm][:, group_cols]
            pred_perm = predict_regression_in_batches(model, X_perm) + baseline[diag_val]
            score, _ = pearsonr(y_true_val, pred_perm)
            perm_scores.append(score)
        results.append({
            'mark': mark,
            'importance': base_score - float(np.mean(perm_scores)),
            'std': float(np.std(perm_scores)),
        })
    return pd.DataFrame(results).sort_values('importance', ascending=False)


def save_bias_plot(bias_df: pd.DataFrame, plot_dir: Path, split: str):
    sub = ordered_bucket_frame(bias_df[bias_df['split'] == split])
    if sub.empty:
        return
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(sub['bucket'], sub['pred_minus_true_smoothed'], marker='o', label='Pred - Smoothed O/E')
    ax.axhline(0.0, color='black', linestyle='--', linewidth=1)
    ax.set_title(f"Bias vs Smoothed O/E by Distance Bucket ({split})")
    ax.set_xlabel("Distance bucket")
    ax.set_ylabel("Mean prediction bias")
    ax.legend(frameon=False)
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(ensure_parent(plot_dir / 'distance_diagnostics' / f'bias_by_bucket_{split}.png'))
    plt.close()

def bootstrap_indices(groups: np.ndarray, bucket_idx: np.ndarray, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    sampled = []
    for group_id in np.unique(groups):
        group_mask = groups == group_id
        group_buckets = bucket_idx[group_mask]
        group_indices = np.flatnonzero(group_mask)
        for bucket_id in np.unique(group_buckets):
            bucket_indices = group_indices[group_buckets == bucket_id]
            if len(bucket_indices) == 0:
                continue
            sampled.append(rng.choice(bucket_indices, size=len(bucket_indices), replace=True))
    out = np.concatenate(sampled)
    rng.shuffle(out)
    return out

def run_stability_analysis(best_params: dict, X_train, y_train, train_buckets, train_groups, X_val, y_true_val, diag_val, val_buckets, baseline):
    print("  Running bootstrap retraining stability analysis...")
    all_run_r = []

    for run_idx in range(N_STABILITY_RETRAINS):
        seed = RANDOM_SEED + run_idx
        idx = bootstrap_indices(train_groups, train_buckets, seed)
        stab_model = make_model(seed, best_params)
        stab_model.fit(X_train[idx], y_train[idx])

        y_pred = predict_regression_in_batches(stab_model, X_val) + baseline[diag_val]
        run_r = []
        for bucket_id in range(N_DIAG_BUCKETS):
            mask = val_buckets == bucket_id
            if mask.sum() > 10:
                r, _ = pearsonr(y_true_val[mask], y_pred[mask])
                run_r.append(r)
            else:
                run_r.append(np.nan)
        all_run_r.append(run_r)
        print(f"  Run {run_idx + 1}/{N_STABILITY_RETRAINS}: {[f'{r:.3f}' if np.isfinite(r) else 'nan' for r in run_r]}")

    return np.array(all_run_r, dtype=np.float32)


def save_regression_cv_summary_plot(cv_df: pd.DataFrame, variant_paths: dict, feature_variant: str):
    if cv_df.empty or 'mean_test_score' not in cv_df.columns:
        return
    fig, ax = plt.subplots(figsize=(8, 4))
    x = np.arange(len(cv_df))
    ax.plot(x, -cv_df['mean_test_score'], marker='o', label='CV validation MSE')
    if 'mean_train_score' in cv_df.columns:
        ax.plot(x, -cv_df['mean_train_score'], marker='o', label='CV train MSE')
    ax.set_title(f"CV Train vs Validation ({feature_variant})")
    ax.set_xlabel("BayesSearch iteration")
    ax.set_ylabel("Training objective MSE (smoothed residual)")
    ax.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(ensure_parent(variant_paths['sampling_plot_dir'] / 'cv_train_vs_validation.png'))
    plt.close()


def save_heatmap(df: pd.DataFrame, index_col: str, column_col: str, value_col: str, out_path: Path, title: str, cmap: str = 'magma', center: float | None = None):
    if df.empty:
        return
    pivot = df.pivot(index=index_col, columns=column_col, values=value_col)
    if pivot.empty:
        return
    pivot = pivot.reindex(columns=[c for c in DIAG_BUCKET_LABELS if c in pivot.columns])
    fig, ax = plt.subplots(figsize=(10, max(3, 0.5 * len(pivot.index))))
    values = pivot.values.astype(float)
    finite = np.isfinite(values)
    if not finite.any():
        return
    vmin = np.nanmin(values)
    vmax = np.nanmax(values)
    if center is not None:
        lim = max(abs(vmin - center), abs(vmax - center))
        vmin, vmax = center - lim, center + lim
    im = ax.imshow(values, aspect='auto', cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_title(title)
    ax.set_xticks(np.arange(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns, rotation=45, ha='right')
    ax.set_yticks(np.arange(len(pivot.index)))
    ax.set_yticklabels(pivot.index)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    plt.tight_layout()
    plt.savefig(ensure_parent(out_path))
    plt.close()


def save_regression_comparison_outputs(summary_df: pd.DataFrame, bucket_df: pd.DataFrame, exact_diag_df: pd.DataFrame, variant_paths_map: dict):
    if summary_df.empty:
        return
    summary_df.to_csv(CHECKPOINT_DIR / 'regression_variant_summary.csv', index=False)
    bucket_df.to_csv(CHECKPOINT_DIR / 'regression_bucket_summary.csv', index=False)
    exact_diag_df.to_csv(CHECKPOINT_DIR / 'regression_exact_diag_summary.csv', index=False)

    primary_label = TARGET_SPACE_LABELS.get(PRIMARY_REGRESSION_TARGET_SPACE, PRIMARY_REGRESSION_TARGET_SPACE)
    overall_primary = summary_df[summary_df['target_space'] == PRIMARY_REGRESSION_TARGET_SPACE].copy()
    if not overall_primary.empty:
        for metric, ylabel, fname in [
            ('Pearson_R', f'Pearson r vs {primary_label}', f'regression_{PRIMARY_REGRESSION_TARGET_SPACE}_pearson.png'),
            ('MSE', f'MSE vs {primary_label}', f'regression_{PRIMARY_REGRESSION_TARGET_SPACE}_mse.png'),
            ('Stratum_Adjusted_R', f'Stratum-adjusted r vs {primary_label}', f'regression_{PRIMARY_REGRESSION_TARGET_SPACE}_stratum_adjusted_r.png'),
        ]:
            fig, ax = plt.subplots(figsize=(10, 4))
            for split in ['train', 'val', 'test']:
                sub = overall_primary[overall_primary['split'] == split]
                if sub.empty:
                    continue
                ax.plot(sub['feature_variant'], sub[metric], marker='o', label=split)
            ax.set_title(f"Regression {ylabel} by Variant and Split")
            ax.set_xlabel("Variant")
            ax.set_ylabel(ylabel)
            ax.legend(frameon=False)
            plt.xticks(rotation=45, ha='right')
            plt.tight_layout()
            plt.savefig(comparison_plot_path('train_val_test', fname))
            plt.close()

    bucket_primary = bucket_df[bucket_df['target_space'] == PRIMARY_REGRESSION_TARGET_SPACE].copy()
    if not bucket_primary.empty:
        for split in ['train', 'val', 'test']:
            sub = ordered_bucket_frame(bucket_primary[bucket_primary['split'] == split])
            save_heatmap(
                sub,
                'feature_variant',
                'bucket',
                'Pearson_R',
                comparison_plot_path('cross_variant', f'regression_bucket_{PRIMARY_REGRESSION_TARGET_SPACE}_pearson_{split}.png'),
                f"Regression {primary_label} Pearson by Bucket ({split})",
            )
            save_heatmap(
                sub,
                'feature_variant',
                'bucket',
                'pred_minus_true_smoothed',
                comparison_plot_path('cross_variant', f'regression_bucket_{PRIMARY_REGRESSION_TARGET_SPACE}_bias_{split}.png'),
                f"Regression {primary_label} Bias by Bucket ({split})",
                cmap='RdBu_r',
                center=0.0,
            )

    if not exact_diag_df.empty:
        for split in ['train', 'val', 'test']:
            fig, ax = plt.subplots(figsize=(10, 4))
            plotted = False
            for variant in summary_df['feature_variant'].unique():
                sub = exact_diag_df[(exact_diag_df['split'] == split) & (exact_diag_df['feature_variant'] == variant)]
                if sub.empty:
                    continue
                ax.plot(sub['diag_idx'], sub['pred_minus_true_smoothed'], linewidth=1.2, label=variant)
                plotted = True
            if not plotted:
                plt.close()
                continue
            ax.axhline(0.0, color='black', linestyle='--', linewidth=1)
            ax.set_title(f"Exact-diagonal {primary_label} Bias ({split})")
            ax.set_xlabel("|i-j| bins")
            ax.set_ylabel(f"Pred - {primary_label}")
            ax.legend(frameon=False)
            plt.tight_layout()
            plt.savefig(comparison_plot_path('cross_variant', f'regression_exact_diag_{PRIMARY_REGRESSION_TARGET_SPACE}_bias_{split}.png'))
            plt.close()

        tracked_ids = sorted({
            path.stem for variant_paths in variant_paths_map.values()
            for path in variant_paths['window_dir'].glob('*.npz')
        })
        for track_id in tracked_ids:
            fig, ax = plt.subplots(figsize=(10, 4))
            plotted = False
            for variant, variant_paths in variant_paths_map.items():
                npz_path = variant_paths['window_dir'] / f'{track_id}.npz'
                if not npz_path.exists():
                    continue
                with np.load(npz_path) as arr:
                    diag_profile = arr['diag_profile_pred_oe']
                ax.plot(np.arange(1, len(diag_profile)), diag_profile[1:], linewidth=1.2, label=variant)
                plotted = True
            if not plotted:
                plt.close()
                continue
            ax.set_title(f"Tracked Window Predicted O/E Profiles ({track_id})")
            ax.set_xlabel("|i-j| bins")
            ax.set_ylabel("Predicted O/E")
            ax.legend(frameon=False)
            plt.tight_layout()
            plt.savefig(comparison_plot_path('cross_variant', f'{track_id}_predicted_diag_overlay.png'))
            plt.close()

# ============================================================
# STAGE 5 — Model training
# ============================================================
print("\n" + "="*60)
print("STAGE 5 — Model training")
print("="*60)

baseline = compute_train_residual_baseline(str(CKPT['features_h5']), CKPT['residual_baseline'])
build_sampled_features_h5(str(CKPT['features_h5']), CKPT['sampled_features'], baseline)
param_space = {
    'n_estimators': Integer(150, 300),
    'max_features': Categorical(['sqrt', 0.3, 0.5]),
    'max_depth': Integer(10, 18),
    'min_samples_split': Integer(10, 50),
    'min_samples_leaf': Integer(5, 25),
    'max_samples': Real(0.05, 0.15),
}

variant_results = []
trained_models = {}
trained_payloads = {}
variant_paths_map = {}
all_eval_frames = []
all_bias_frames = []
all_exact_diag_frames = []

trainable_regression_variants = [
    'epi_plus_distance_regression',
    'epi_plus_distance_pairwise_regression',
]
if ENABLE_CONTROL_VARIANTS:
    trainable_regression_variants.extend(CONTROL_REGRESSION_VARIANTS)

for feature_variant in trainable_regression_variants:
    variant_paths = variant_ckpt_paths(feature_variant)
    variant_paths_map[feature_variant] = variant_paths
    sample_df = summarize_sampled_dataset(str(CKPT['sampled_features']), variant_paths['sample_report'], variant_paths['plot_dir'])
    X_train, y_train, y_true_train_raw, y_true_train_smooth, diag_train, bucket_train, groups_train = load_sampled_split(
        str(CKPT['sampled_features']), 'train', feature_variant
    )
    X_val, y_val_train, y_val_true_raw, y_val_true_smooth, diag_val, bucket_val, groups_val = load_sampled_split(
        str(CKPT['sampled_features']), 'val', feature_variant
    )
    print(f"  Loaded sampled train matrix for {feature_variant}: {X_train.shape}")

    model_ready = variant_paths['model'].exists() and variant_paths['best_params'].exists() and model_meta_is_current(variant_paths, feature_variant)
    best_params_ready = variant_paths['best_params'].exists() and model_meta_is_current(variant_paths, feature_variant)

    if model_ready:
        model = load(variant_paths['model'])
        best_params = load(variant_paths['best_params'])
        print(f"  Loaded regression model checkpoint for {feature_variant}.")
    else:
        if best_params_ready:
            best_params = load(variant_paths['best_params'])
            print(f"  Loaded tuned hyperparameters for {feature_variant}.")
        else:
            estimator = RandomForestRegressor(random_state=RANDOM_SEED, n_jobs=RF_TRAIN_N_JOBS)
            cv = GroupKFold(n_splits=CV_FOLDS)
            print(f"  Running BayesSearchCV for {feature_variant}...")
            bayes = BayesSearchCV(
                estimator=estimator,
                search_spaces=param_space,
                n_iter=N_BAYES_ITER,
                cv=cv,
                scoring='neg_mean_squared_error',
                random_state=RANDOM_SEED,
                n_jobs=1,
                verbose=3,
                return_train_score=True,
            )
            bayes.fit(X_train, y_train, groups=groups_train)
            best_params = dict(bayes.best_params_)
            dump(best_params, variant_paths['best_params'])
            write_model_meta(variant_paths, best_params, feature_variant)
            cv_results_df = pd.DataFrame(bayes.cv_results_)
            cv_results_df.to_csv(variant_paths['cv_report'], index=False)
            save_regression_cv_summary_plot(cv_results_df, variant_paths, feature_variant)
            print(f"  Saved tuned hyperparameters → {variant_paths['best_params']}")

        print(f"  Fitting residualized RF on sampled training set ({feature_variant})...")
        model = make_model(RANDOM_SEED, best_params)
        model.fit(X_train, y_train)
        dump(model, variant_paths['model'])
        write_model_meta(variant_paths, best_params, feature_variant)
        print(f"  Saved trained model → {variant_paths['model']}")

    trained_models[feature_variant] = model
    trained_payloads[feature_variant] = {
        'variant_paths': variant_paths,
        'best_params': best_params,
        'X_train': X_train,
        'y_train': y_train,
        'bucket_train': bucket_train,
        'groups_train': groups_train,
        'X_val': X_val,
        'y_val_true_raw': y_val_true_raw,
        'y_val_true_smooth': y_val_true_smooth,
        'diag_val': diag_val,
        'bucket_val': bucket_val,
    }

    eval_train_df, bias_train_df, exact_train_df = evaluate_streamed_regression(model, baseline, feature_variant, 'train')
    eval_val_df, bias_val_df, exact_val_df = evaluate_streamed_regression(model, baseline, feature_variant, 'val')
    eval_test_df, bias_test_df, exact_test_df = evaluate_streamed_regression(model, baseline, feature_variant, 'test')
    eval_df = pd.concat([eval_train_df, eval_val_df, eval_test_df], ignore_index=True)
    bias_df = pd.concat([bias_train_df, bias_val_df, bias_test_df], ignore_index=True)
    exact_diag_df = pd.concat([exact_train_df, exact_val_df, exact_test_df], ignore_index=True)
    eval_df.insert(0, 'feature_variant', feature_variant)
    bias_df.insert(0, 'feature_variant', feature_variant)
    exact_diag_df.insert(0, 'feature_variant', feature_variant)
    eval_df.to_csv(variant_paths['eval_report'], index=False)
    bias_df.to_csv(variant_paths['bias_report'], index=False)
    exact_diag_df.to_csv(variant_paths['exact_diag_report'], index=False)
    save_bias_plot(bias_df, variant_paths['plot_dir'], 'val')
    save_bias_plot(bias_df, variant_paths['plot_dir'], 'test')
    save_bias_plot(bias_df, variant_paths['plot_dir'], 'train')
    generate_prediction_heatmaps(model, baseline, feature_variant, variant_paths)
    eval_df[eval_df['bucket'] == 'overall'].to_csv(variant_paths['variant_summary'], index=False)
    bucket_summary_df = eval_df[(eval_df['target_space'] == PRIMARY_REGRESSION_TARGET_SPACE) & (eval_df['bucket'] != 'overall')].merge(
        bias_df[['feature_variant', 'split', 'bucket', 'pred_minus_true_raw', 'pred_minus_true_smoothed']],
        on=['feature_variant', 'split', 'bucket'],
        how='left',
    )
    bucket_summary_df.to_csv(variant_paths['bucket_summary'], index=False)
    all_eval_frames.append(eval_df)
    all_bias_frames.append(bias_df)
    all_exact_diag_frames.append(exact_diag_df)

    importance_df = compute_grouped_permutation_importance(model, X_val, y_val_true_smooth, diag_val, baseline, feature_variant)
    importance_df.insert(0, 'target_space', PRIMARY_REGRESSION_TARGET_SPACE)
    importance_df.to_csv(variant_paths['dir'] / 'feature_importances.csv', index=False)
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar(importance_df['mark'], importance_df['importance'], color=[MARK_COLORS.get(m, 'dimgray') for m in importance_df['mark']])
    ax.set_ylabel("Permutation Importance (drop in Pearson r vs smoothed O/E)")
    ax.set_title(f"Grouped Permutation Importances vs Smoothed O/E ({feature_variant})")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(ensure_parent(variant_paths['feature_plot_dir'] / 'feature_importances.png'))
    plt.close()

    all_run_r = run_stability_analysis(
        best_params=best_params,
        X_train=X_train,
        y_train=y_train,
        train_buckets=bucket_train,
        train_groups=groups_train,
        X_val=X_val,
        y_true_val=y_val_true_smooth,
        diag_val=diag_val,
        val_buckets=bucket_val,
        baseline=baseline,
    )
    stability_records = []
    for run_idx, run_values in enumerate(all_run_r, start=1):
        for bucket_label, pearson_r in zip(DIAG_BUCKET_LABELS, run_values):
            stability_records.append({
                'run': run_idx,
                'bucket': bucket_label,
                'target_space': PRIMARY_REGRESSION_TARGET_SPACE,
                'Pearson_R': float(pearson_r) if np.isfinite(pearson_r) else np.nan,
            })
    pd.DataFrame(stability_records).to_csv(variant_paths['stability_report'], index=False)
    fig, ax = plt.subplots(figsize=(10, 5))
    bp = ax.boxplot(all_run_r, labels=DIAG_BUCKET_LABELS, patch_artist=True)
    for patch in bp['boxes']:
        patch.set_facecolor('lightsteelblue')
    ax.set_title(f'Stability vs Smoothed O/E: per-bucket Pearson r across {N_STABILITY_RETRAINS} bootstrap retrains', fontsize=12)
    ax.set_ylabel('Pearson r vs smoothed O/E')
    ax.set_xlabel('Diagonal bucket')
    plt.xticks(rotation=45)
    plt.tight_layout()
    stability_plot = ensure_parent(variant_paths['stability_plot_dir'] / 'stability_stratum_pearson_smoothed_oe.png')
    plt.savefig(stability_plot)
    plt.savefig(ensure_parent(variant_paths['stability_plot_dir'] / 'stability_stratum_pearson.png'))
    plt.close()

    train_row = require_matching_row(eval_df, "regression training overall row", split='train', target_space=PRIMARY_REGRESSION_TARGET_SPACE, bucket='overall')
    val_row = require_matching_row(eval_df, "regression validation overall row", split='val', target_space=PRIMARY_REGRESSION_TARGET_SPACE, bucket='overall')
    test_row = require_matching_row(eval_df, "regression test overall row", split='test', target_space=PRIMARY_REGRESSION_TARGET_SPACE, bucket='overall')
    variant_results.append({
        'feature_variant': feature_variant,
        'primary_target_space': PRIMARY_REGRESSION_TARGET_SPACE,
        f'pearson_train_{PRIMARY_REGRESSION_TARGET_SPACE}': float(train_row['Pearson_R']),
        f'pearson_val_{PRIMARY_REGRESSION_TARGET_SPACE}': float(val_row['Pearson_R']),
        f'pearson_test_{PRIMARY_REGRESSION_TARGET_SPACE}': float(test_row['Pearson_R']),
        f'mse_train_{PRIMARY_REGRESSION_TARGET_SPACE}': float(train_row['MSE']),
        f'mse_val_{PRIMARY_REGRESSION_TARGET_SPACE}': float(val_row['MSE']),
        f'mse_test_{PRIMARY_REGRESSION_TARGET_SPACE}': float(test_row['MSE']),
        f'stratum_adjusted_train_{PRIMARY_REGRESSION_TARGET_SPACE}': float(train_row['Stratum_Adjusted_R']),
        f'stratum_adjusted_val_{PRIMARY_REGRESSION_TARGET_SPACE}': float(val_row['Stratum_Adjusted_R']),
        f'stratum_adjusted_test_{PRIMARY_REGRESSION_TARGET_SPACE}': float(test_row['Stratum_Adjusted_R']),
    })

print("\n" + "="*60)
print("STAGE 6b — Exact-Diagonal Calibrated Evaluation")
print("="*60)


def evaluate_exact_diag_calibrated_variant(base_variant: str, calibrated_variant: str):
    print(f"  Evaluating exact-diagonal calibrated variant: {calibrated_variant} (base={base_variant})")
    base_model = trained_models[base_variant]
    calibrated_paths = variant_ckpt_paths(calibrated_variant)
    variant_paths_map[calibrated_variant] = calibrated_paths
    exact_diag_correction, correction_df = compute_exact_diag_correction(base_model, baseline, base_variant)
    pred_level_correction, pred_level_df = compute_pred_level_correction(base_model, baseline, base_variant, exact_diag_correction)
    combined_correction = {
        'exact_diag': exact_diag_correction,
        'pred_level': pred_level_correction,
    }
    with open(calibrated_paths['correction_summary'], 'w') as fh:
        json.dump(
            {
                'mode': 'exact_diag_plus_pred_level_smoothed_oe',
                'base_variant': base_variant,
                'target_space': PRIMARY_REGRESSION_TARGET_SPACE,
                'exact_diag_sigma': EXACT_DIAG_CORRECTION_SIGMA,
                'pred_level_sigma': PRED_LEVEL_CORRECTION_SIGMA,
                'pred_level_bins': PRED_LEVEL_CORRECTION_BINS,
                'n_diags': int(len(exact_diag_correction)),
                'n_pred_level_points': int(0 if pred_level_correction is None else len(pred_level_correction['pred_center'])),
            },
            fh,
            indent=2,
        )
    summarize_sampled_dataset(str(CKPT['sampled_features']), calibrated_paths['sample_report'], calibrated_paths['plot_dir'])
    eval_train_df, bias_train_df, exact_train_df = evaluate_streamed_regression(base_model, baseline, base_variant, 'train', correction=combined_correction)
    eval_val_df, bias_val_df, exact_val_df = evaluate_streamed_regression(base_model, baseline, base_variant, 'val', correction=combined_correction)
    eval_test_df, bias_test_df, exact_test_df = evaluate_streamed_regression(base_model, baseline, base_variant, 'test', correction=combined_correction)
    eval_df = pd.concat([eval_train_df, eval_val_df, eval_test_df], ignore_index=True)
    bias_df = pd.concat([bias_train_df, bias_val_df, bias_test_df], ignore_index=True)
    exact_diag_df = pd.concat([exact_train_df, exact_val_df, exact_test_df], ignore_index=True)
    eval_df.insert(0, 'feature_variant', calibrated_variant)
    bias_df.insert(0, 'feature_variant', calibrated_variant)
    exact_diag_df.insert(0, 'feature_variant', calibrated_variant)
    eval_df.to_csv(calibrated_paths['eval_report'], index=False)
    bias_df.to_csv(calibrated_paths['bias_report'], index=False)
    exact_diag_df.to_csv(calibrated_paths['exact_diag_report'], index=False)
    correction_df.to_csv(calibrated_paths['dir'] / 'exact_diag_calibration_curve.csv', index=False)
    pred_level_df.to_csv(calibrated_paths['pred_level_correction'], index=False)
    save_bias_plot(bias_df, calibrated_paths['plot_dir'], 'train')
    save_bias_plot(bias_df, calibrated_paths['plot_dir'], 'val')
    save_bias_plot(bias_df, calibrated_paths['plot_dir'], 'test')
    generate_prediction_heatmaps(
        base_model,
        baseline,
        calibrated_variant,
        calibrated_paths,
        correction=combined_correction,
        model_feature_variant=base_variant,
    )
    eval_df[eval_df['bucket'] == 'overall'].to_csv(calibrated_paths['variant_summary'], index=False)
    bucket_summary_df = eval_df[(eval_df['target_space'] == PRIMARY_REGRESSION_TARGET_SPACE) & (eval_df['bucket'] != 'overall')].merge(
        bias_df[['feature_variant', 'split', 'bucket', 'pred_minus_true_raw', 'pred_minus_true_smoothed']],
        on=['feature_variant', 'split', 'bucket'],
        how='left',
    )
    bucket_summary_df.to_csv(calibrated_paths['bucket_summary'], index=False)
    all_eval_frames.append(eval_df)
    all_bias_frames.append(bias_df)
    all_exact_diag_frames.append(exact_diag_df)
    train_row = require_matching_row(eval_df, f"{calibrated_variant} training overall row", split='train', target_space=PRIMARY_REGRESSION_TARGET_SPACE, bucket='overall')
    val_row = require_matching_row(eval_df, f"{calibrated_variant} validation overall row", split='val', target_space=PRIMARY_REGRESSION_TARGET_SPACE, bucket='overall')
    test_row = require_matching_row(eval_df, f"{calibrated_variant} test overall row", split='test', target_space=PRIMARY_REGRESSION_TARGET_SPACE, bucket='overall')
    variant_results.append({
        'feature_variant': calibrated_variant,
        'base_variant': base_variant,
        'primary_target_space': PRIMARY_REGRESSION_TARGET_SPACE,
        f'pearson_train_{PRIMARY_REGRESSION_TARGET_SPACE}': float(train_row['Pearson_R']),
        f'pearson_val_{PRIMARY_REGRESSION_TARGET_SPACE}': float(val_row['Pearson_R']),
        f'pearson_test_{PRIMARY_REGRESSION_TARGET_SPACE}': float(test_row['Pearson_R']),
        f'mse_train_{PRIMARY_REGRESSION_TARGET_SPACE}': float(train_row['MSE']),
        f'mse_val_{PRIMARY_REGRESSION_TARGET_SPACE}': float(val_row['MSE']),
        f'mse_test_{PRIMARY_REGRESSION_TARGET_SPACE}': float(test_row['MSE']),
        f'stratum_adjusted_train_{PRIMARY_REGRESSION_TARGET_SPACE}': float(train_row['Stratum_Adjusted_R']),
        f'stratum_adjusted_val_{PRIMARY_REGRESSION_TARGET_SPACE}': float(val_row['Stratum_Adjusted_R']),
        f'stratum_adjusted_test_{PRIMARY_REGRESSION_TARGET_SPACE}': float(test_row['Stratum_Adjusted_R']),
    })


evaluate_exact_diag_calibrated_variant(
    'epi_plus_distance_regression',
    'epi_plus_distance_exact_diag_calibrated_regression',
)
evaluate_exact_diag_calibrated_variant(
    'epi_plus_distance_pairwise_regression',
    'epi_plus_distance_pairwise_exact_diag_calibrated_regression',
)

pd.DataFrame(variant_results).to_csv(CHECKPOINT_DIR / 'ablation_summary_regression.csv', index=False)
summary_df = pd.concat(all_eval_frames, ignore_index=True)
bias_df = pd.concat(all_bias_frames, ignore_index=True)
exact_diag_df = pd.concat(all_exact_diag_frames, ignore_index=True)
bucket_summary_df = summary_df[(summary_df['target_space'] == PRIMARY_REGRESSION_TARGET_SPACE) & (summary_df['bucket'] != 'overall')].merge(
    bias_df[['feature_variant', 'split', 'bucket', 'pred_minus_true_raw', 'pred_minus_true_smoothed']],
    on=['feature_variant', 'split', 'bucket'],
    how='left',
)
save_regression_comparison_outputs(
    summary_df[summary_df['bucket'] == 'overall'].copy(),
    bucket_summary_df,
    exact_diag_df,
    variant_paths_map,
)

print("\nPIPELINE COMPLETE")

#!/usr/bin/env python3
"""
predict_contact_from_feature_allelic.py
======================================
Stage-wise inference pipeline for predicting 200 bp smoothed O/E contact maps
from allele-specific H3K27me3 bigWigs (Xa and Xi), using a trained Random
Forest model. The script:

- computes shared whole-genome normalization statistics from BOTH allele-
  specific bigWigs, by concatenating log1p-binned values genome-wide the same
  way as the original training code;
- computes a shared H3K27me3 enrichment threshold from the SAME combined Xa+Xi
  whole-genome distribution;
- predicts separate contact maps for Xa and Xi on the same requested region;
- always predicts on a centered 100 kb window, then crops to the requested
  region for downstream intensity / occupancy conversion and fiber assignment;
- gates pair predictions so that both i and j bins must exceed the methylation
  threshold;
- reconstructs an observed-like target intensity from predicted smoothed
  log2(O/E), subtracts a scaled polymer baseline at the probability level, and
  interprets the result as an absolute extra-contact probability;
- assigns extra contacts independently across feasible fibers using Bernoulli
  sampling rather than a fixed global contact budget;
- reads numbered dim files such as dim1.in, dim2.in, ... using natural sorting.

The code keeps the staged design and plotting style close to the original
training pipeline. In the manuscript, the baseline-subtraction scale is denoted
by gamma; in this script the older command-line name --kappa is kept for
backward compatibility.
"""

from __future__ import annotations

import argparse
import glob
import json
import logging
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import joblib
import numpy as np
import pandas as pd

try:
    import pyBigWig
except Exception:
    pyBigWig = None

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle

try:
    import cooler
except Exception:
    cooler = None

from scipy.ndimage import convolve, gaussian_filter1d


# ============================================================
# Constants matching the training pipeline
# ============================================================
REQUESTED_RESOLUTION = 200
WINDOW_SIZE_BP = 100_000
TARGET_MIN, TARGET_MAX = -5.0, 5.0
LOG_COLOR_QUANTILE_LOW = 1.0
LOG_COLOR_QUANTILE_HIGH = 99.5
LOG_FLOOR_INTENSITY = 1e-8
LOG_FLOOR_OCCUPANCY = 1e-7
PREDICT_BATCH_SIZE = 500_000
SMOOTH_KERNEL_SIZE = 3
RANDOM_SEED = 42
LINKER_BP = 9
NUCLEOSOME_BP = 146
FIRST_CORE_OFFSET_BP = 144 // 2  # matches the C code exactly
SCRIPT_VERSION = "allelic_predict_v5_polymer_alpha_cleanup"
STAGE3_INFERENCE_VERSION = "residual_baseline_reconstruction_v1"
DIAG_BUCKET_EDGES_DEFAULT = [1, 5, 10, 20, 50, 100, 200]

MARK_NAMES = [
    "H3K27ac",
    "H3K4me1",
    "H3K4me3",
    "H3K36me3",
    "H3K9me3",
    "H3K27me3",
    "ATAC",
    "CTCF",
]
N_MARKS = len(MARK_NAMES)
MARK_TO_INDEX = {m: i for i, m in enumerate(MARK_NAMES)}

MARK_COLORS = {
    "H3K4me3": "forestgreen",
    "H3K27ac": "limegreen",
    "H3K4me1": "mediumseagreen",
    "H3K36me3": "darkgreen",
    "H3K9me3": "firebrick",
    "H3K27me3": "indianred",
    "ATAC": "royalblue",
    "CTCF": "darkorchid",
}

FEATURE_NAMES = (
    [f"{m}_i" for m in MARK_NAMES]
    + [f"{m}_j" for m in MARK_NAMES]
    + [f"{m}_prod" for m in MARK_NAMES]
    + [f"{m}_int" for m in MARK_NAMES]
    + [f"{m}_smooth_prod" for m in MARK_NAMES]
    + ["log1p_dist", "diag_index"]
)
ACTIVE_FEATURE_NAMES = FEATURE_NAMES[:-2]
ALL_ACTIVE_FEATURE_NAMES = FEATURE_NAMES[:]
PAIRWISE_VARIANT_SUFFIX = "_pairwise_regression"
PAIRWISE_REGRESSION_FEATURE_NAMES = (
    [f"{m}_min" for m in MARK_NAMES]
    + [f"{m}_absdiff" for m in MARK_NAMES]
    + [f"{m}_contrast_i" for m in MARK_NAMES]
    + [f"{m}_contrast_j" for m in MARK_NAMES]
)

plt.rcParams.update({
    "figure.dpi": 130,
    "savefig.dpi": 300,
    "font.size": 11,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
})


# ============================================================
# Dataclasses
# ============================================================
@dataclass
class Region:
    chrom: str
    start: int
    end: int

    @property
    def center(self) -> float:
        return 0.5 * (self.start + self.end)


# ============================================================
# Logging / filesystem helpers
# ============================================================
def setup_logger(log_path: Path) -> logging.Logger:
    logger = logging.getLogger(f"predict_contact_from_feature_allelic::{log_path}")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")

    fh = logging.FileHandler(log_path)
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    sh = logging.StreamHandler()
    sh.setFormatter(formatter)
    logger.addHandler(sh)
    return logger


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def load_json(path: Optional[Path]) -> dict:
    if path is None or not path.exists():
        return {}
    with open(path) as fh:
        return json.load(fh)


def write_json(path: Path, payload: dict) -> None:
    ensure_dir(path.parent)
    with open(path, "w") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True)


def json_safe(value):
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, Region):
        return {"chrom": value.chrom, "start": int(value.start), "end": int(value.end)}
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    return value


def metadata_matches(path: Path, expected: dict) -> bool:
    if not path.exists():
        return False
    observed = load_json(path)
    return observed.get("inputs") == json_safe(expected)


def write_stage_metadata(path: Path, inputs: dict, outputs: Optional[dict] = None) -> None:
    write_json(path, {
        "script_version": SCRIPT_VERSION,
        "inputs": json_safe(inputs),
        "outputs": json_safe(outputs or {}),
    })


def files_exist(paths: Sequence[Path]) -> bool:
    return all(Path(p).exists() for p in paths)


def read_csv_or_empty(path: Path) -> pd.DataFrame:
    try:
        return pd.read_csv(path)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def file_signature(path: Path) -> dict:
    p = Path(path)
    if not p.exists():
        return {"path": str(p), "exists": False}
    st = p.stat()
    return {"path": str(p), "exists": True, "size": int(st.st_size), "mtime_ns": int(st.st_mtime_ns)}


def region_dict(region: Region) -> dict:
    return {"chrom": region.chrom, "start": int(region.start), "end": int(region.end)}


# ============================================================
# Generic helpers
# ============================================================
def infer_feature_variant(model, model_meta: dict) -> str:
    feature_variant = model_meta.get("feature_variant")
    if feature_variant:
        return feature_variant

    n_in = getattr(model, "n_features_in_", None)
    if n_in is None:
        return "epi_plus_distance_regression"

    if n_in == len(ALL_ACTIVE_FEATURE_NAMES):
        return "epi_plus_distance_regression"
    if n_in == len(ACTIVE_FEATURE_NAMES):
        return "epi_only_regression"
    if n_in == len(ALL_ACTIVE_FEATURE_NAMES) + len(PAIRWISE_REGRESSION_FEATURE_NAMES):
        return "epi_plus_distance_pairwise_regression"
    if n_in == len(ACTIVE_FEATURE_NAMES) + len(PAIRWISE_REGRESSION_FEATURE_NAMES):
        return "epi_only_pairwise_regression"
    raise ValueError(f"Could not infer feature variant from model input size: {n_in}")


def get_feature_variant_columns(feature_variant: str):
    if feature_variant.startswith("epi_plus_distance"):
        return np.arange(len(ALL_ACTIVE_FEATURE_NAMES), dtype=np.int32), list(ALL_ACTIVE_FEATURE_NAMES)
    return np.arange(len(ACTIVE_FEATURE_NAMES), dtype=np.int32), list(ACTIVE_FEATURE_NAMES)


def feature_names_for_variant(feature_variant: str) -> List[str]:
    _, feature_names = get_feature_variant_columns(feature_variant)
    if feature_variant.endswith(PAIRWISE_VARIANT_SUFFIX):
        return feature_names + list(PAIRWISE_REGRESSION_FEATURE_NAMES)
    return feature_names


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


def select_feature_matrix(X: np.ndarray, feature_variant: str) -> np.ndarray:
    cols, _ = get_feature_variant_columns(feature_variant)
    out = X[:, cols]
    if feature_variant.endswith(PAIRWISE_VARIANT_SUFFIX):
        out = build_pairwise_regression_features(out)
    return out


def feature_block_and_mark(feature_name: str) -> Tuple[str, str]:
    for mark in MARK_NAMES:
        prefix = f"{mark}_"
        if feature_name.startswith(prefix):
            return feature_name[len(prefix):], mark
    if feature_name in {"log1p_dist", "diag_index"}:
        return "distance", "Distance"
    return "derived", ""


def save_feature_nonzero_summary(X: np.ndarray, feature_names: Sequence[str], out_csv: Path) -> None:
    rows = []
    for col_idx, feature_name in enumerate(feature_names):
        vals = X[:, col_idx]
        finite = vals[np.isfinite(vals)]
        nonzero = np.count_nonzero(np.nan_to_num(vals, nan=0.0))
        block, mark = feature_block_and_mark(feature_name)
        rows.append({
            "feature_index": int(col_idx),
            "feature_name": feature_name,
            "block": block,
            "mark": mark,
            "nonzero_count": int(nonzero),
            "total_count": int(len(vals)),
            "nonzero_fraction": float(nonzero / max(len(vals), 1)),
            "min": float(np.min(finite)) if finite.size else np.nan,
            "mean": float(np.mean(finite)) if finite.size else np.nan,
            "max": float(np.max(finite)) if finite.size else np.nan,
        })
    pd.DataFrame(rows).to_csv(out_csv, index=False)


def predict_in_batches(model, X: np.ndarray, batch_size: int = PREDICT_BATCH_SIZE) -> np.ndarray:
    preds = np.zeros(len(X), dtype=np.float32)
    for start in range(0, len(X), batch_size):
        stop = min(start + batch_size, len(X))
        preds[start:stop] = model.predict(X[start:stop]).astype(np.float32)
    return preds


def natural_sort_key(path_like) -> list:
    s = str(path_like)
    return [int(tok) if tok.isdigit() else tok.lower() for tok in re.split(r"(\d+)", s)]


def parse_mark_bw_args(entries: Optional[Sequence[str]]) -> Dict[str, List[str]]:
    parsed: Dict[str, List[str]] = {}
    if not entries:
        return parsed
    for item in entries:
        if "=" not in item:
            raise ValueError(
                f"Invalid mark-specific BigWig argument '{item}'. Expected MARK=/path/to/file.bw"
            )
        mark, path = item.split("=", 1)
        mark = mark.strip()
        path = path.strip()
        if mark not in MARK_TO_INDEX:
            raise ValueError(
                f"Unsupported mark '{mark}' in mark-specific BigWig input. "
                f"Supported marks: {', '.join(MARK_NAMES)}"
            )
        if not path:
            raise ValueError(f"Empty BigWig path for mark '{mark}'")
        parsed.setdefault(mark, []).append(path)
    return parsed


def ordered_mark_list(mark_map: Dict[str, Sequence[str]]) -> List[str]:
    return [mark for mark in MARK_NAMES if mark in mark_map]


def resolve_allele_mark_bw_map(args) -> Tuple[Dict[str, Dict[str, List[str]]], List[str], str, bool]:
    xa_mark_map = parse_mark_bw_args(getattr(args, "xa_mark_bw", None))
    xi_mark_map = parse_mark_bw_args(getattr(args, "xi_mark_bw", None))
    multi_mark_mode = bool(xa_mark_map or xi_mark_map)

    if multi_mark_mode:
        if not xa_mark_map or not xi_mark_map:
            raise ValueError("Multi-mark mode requires both --xa-mark-bw and --xi-mark-bw")
        if set(xa_mark_map.keys()) != set(xi_mark_map.keys()):
            raise ValueError(
                "Xa and Xi mark-specific inputs must include the same mark set. "
                f"Xa={sorted(xa_mark_map.keys())}, Xi={sorted(xi_mark_map.keys())}"
            )
        included_marks = ordered_mark_list(xa_mark_map)
    else:
        included_marks = [args.feature_slot]
        xa_mark_map = {args.feature_slot: [str(x) for x in args.xa_bw]}
        xi_mark_map = {args.feature_slot: [str(x) for x in args.xi_bw]}

    primary_feature_slot = args.feature_slot if args.feature_slot in included_marks else included_marks[0]
    return {"Xa": xa_mark_map, "Xi": xi_mark_map}, included_marks, primary_feature_slot, multi_mark_mode


def mark_stats_map(shared_stats: dict) -> Dict[str, dict]:
    return shared_stats.get("mark_stats", {})


# ============================================================
# Window / bin geometry
# ============================================================
def resolve_centered_prediction_window(requested: Region, chrom_size: int, resolution: int = REQUESTED_RESOLUTION) -> Region:
    center = int(round(requested.center))
    raw_start = center - WINDOW_SIZE_BP // 2
    raw_end = raw_start + WINDOW_SIZE_BP

    if raw_start < 0:
        raw_start = 0
        raw_end = WINDOW_SIZE_BP
    if raw_end > chrom_size:
        raw_end = chrom_size
        raw_start = max(0, raw_end - WINDOW_SIZE_BP)

    raw_start = (raw_start // resolution) * resolution
    raw_end = raw_start + WINDOW_SIZE_BP
    if raw_end > chrom_size:
        raw_end = (chrom_size // resolution) * resolution
        raw_start = max(0, raw_end - WINDOW_SIZE_BP)
    if (raw_end - raw_start) != WINDOW_SIZE_BP:
        raise ValueError(f"Could not resolve a full 100 kb prediction window on {requested.chrom}")
    return Region(requested.chrom, raw_start, raw_end)


def requested_region_bin_slice(requested: Region, prediction: Region, resolution: int = REQUESTED_RESOLUTION) -> Tuple[int, int]:
    start_offset = max(0, requested.start - prediction.start)
    end_offset = max(0, requested.end - prediction.start)
    b0 = int(math.floor(start_offset / resolution))
    b1 = int(math.ceil(end_offset / resolution))
    max_bins = (prediction.end - prediction.start) // resolution
    b0 = max(0, min(b0, max_bins))
    b1 = max(b0, min(b1, max_bins))
    return b0, b1


# ============================================================
# BigWig / Cooler loading
# ============================================================
def load_bigwig_window_average(bw_paths: Sequence[str], chrom: str, start: int, end: int, n_bins: int) -> np.ndarray:
    if pyBigWig is None:
        raise RuntimeError("pyBigWig is not installed, but BigWig input was requested")
    acc = np.zeros(n_bins, dtype=np.float64)
    for bw_path in bw_paths:
        with pyBigWig.open(str(bw_path)) as bw:
            chrom_sizes = bw.chroms()
            if chrom not in chrom_sizes:
                raise ValueError(
                    f"Chromosome {chrom} is not present in BigWig {bw_path}. "
                    f"Requested interval was {chrom}:{start}-{end}."
                )
            chrom_len = int(chrom_sizes[chrom])
            if start < 0 or end <= start or end > chrom_len:
                raise ValueError(
                    f"Invalid interval bounds for BigWig {bw_path}: "
                    f"requested {chrom}:{start}-{end}, file length={chrom_len}."
                )
            vals = bw.stats(chrom, start, end, type="mean", nBins=n_bins)
        arr = np.array([0.0 if v is None else v for v in vals], dtype=np.float32)
        arr = np.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0)
        acc += arr
    return (acc / max(len(bw_paths), 1)).astype(np.float32)


def collect_shared_bigwig_chrom_sizes(
    allele_mark_bw_map: Dict[str, Dict[str, Sequence[str]]],
    included_marks: Sequence[str],
    logger: Optional[logging.Logger] = None,
) -> Dict[str, int]:
    if pyBigWig is None:
        raise RuntimeError("pyBigWig is not installed, but BigWig input was requested")

    shared_chroms: Optional[set] = None
    min_sizes: Dict[str, int] = {}
    bw_count = 0
    for allele_name, mark_bw_map in allele_mark_bw_map.items():
        for mark in included_marks:
            for bw_path in mark_bw_map[mark]:
                with pyBigWig.open(str(bw_path)) as bw:
                    chrom_map = {str(chrom): int(size) for chrom, size in bw.chroms().items()}
                chrom_set = set(chrom_map.keys())
                if shared_chroms is None:
                    shared_chroms = chrom_set
                    min_sizes = chrom_map.copy()
                else:
                    shared_chroms &= chrom_set
                    for chrom in list(min_sizes.keys()):
                        if chrom in chrom_map:
                            min_sizes[chrom] = min(min_sizes[chrom], chrom_map[chrom])
                        else:
                            min_sizes.pop(chrom, None)
                bw_count += 1

    shared_chroms = shared_chroms or set()
    safe_sizes = {chrom: min_sizes[chrom] for chrom in sorted(shared_chroms) if chrom in min_sizes}
    if not safe_sizes:
        raise ValueError("No common chromosomes were found across the supplied Xa/Xi mark-specific BigWig files.")
    if logger is not None:
        logger.info("Resolved shared chromosome sizes across %d BigWig files (%d common chromosomes)", bw_count, len(safe_sizes))
    return safe_sizes


def compute_combined_whole_genome_feature_stats(
    allele_mark_bw_map: Dict[str, Dict[str, Sequence[str]]],
    included_marks: Sequence[str],
    resolution: int,
    logger: logging.Logger,
    methyl_gate_quantile: float,
    epigenetic_noise_floor_quantile: Optional[float] = None,
    chrom_include_regex: Optional[str] = r"^chr([0-9]+|X|Y|M)$",
) -> Tuple[dict, pd.DataFrame]:
    """Compute shared raw-scale Xa+Xi thresholds before log normalization."""
    logger.info("STAGE 1 — Computing shared whole-genome Xa+Xi statistics")
    if epigenetic_noise_floor_quantile is not None and float(epigenetic_noise_floor_quantile) != float(methyl_gate_quantile):
        logger.warning(
            "--epigenetic-noise-floor-quantile is deprecated and ignored; using methyl_gate_quantile %.2f for both feature denoising and pair gating.",
            methyl_gate_quantile,
        )
    epigenetic_noise_floor_quantile = float(methyl_gate_quantile)
    chrom_pat = re.compile(chrom_include_regex) if chrom_include_regex else None

    chrom_sizes = collect_shared_bigwig_chrom_sizes(allele_mark_bw_map, included_marks, logger)

    chrom_rows = []
    all_vals_by_mark: Dict[str, List[np.ndarray]] = {mark: [] for mark in included_marks}
    mark_summary: Dict[str, dict] = {}
    for mark in included_marks:
        for allele_name, mark_bw_map in allele_mark_bw_map.items():
            bw_paths = mark_bw_map[mark]
            for chrom, chrom_size in chrom_sizes.items():
                if chrom_pat and not chrom_pat.match(chrom):
                    continue
                n_bins = chrom_size // resolution
                if n_bins < 1:
                    continue
                vals = load_bigwig_window_average(bw_paths, chrom, 0, n_bins * resolution, n_bins)
                vals = vals.astype(np.float32, copy=False)
                pos = vals[vals > 0]
                if pos.size:
                    all_vals_by_mark[mark].append(pos)
                    pos_log = np.log1p(pos)
                    chrom_rows.append({
                        "mark": mark,
                        "allele": allele_name,
                        "chrom": chrom,
                        "n_bins": int(n_bins),
                        "n_positive": int(pos.size),
                        "median_raw": float(np.median(pos)),
                        "p99_raw": float(np.percentile(pos, 99)),
                        "median_log1p": float(np.median(pos_log)),
                        "p99_log1p": float(np.percentile(pos_log, 99)),
                    })
                else:
                    chrom_rows.append({
                        "mark": mark,
                        "allele": allele_name,
                        "chrom": chrom,
                        "n_bins": int(n_bins),
                        "n_positive": 0,
                        "median_raw": 0.0,
                        "p99_raw": 0.0,
                        "median_log1p": 0.0,
                        "p99_log1p": 0.0,
                    })

        if not all_vals_by_mark[mark]:
            logger.warning(
                "No positive Xa/Xi feature values found genome-wide for %s; using p99_raw = 1.0 and thresholds = 0.0",
                mark,
            )
            p99_raw = 1.0
            methyl_gate_raw_threshold = 0.0
        else:
            all_concat = np.concatenate(all_vals_by_mark[mark]).astype(np.float32, copy=False)
            p99_raw = float(np.percentile(all_concat, 99)) if all_concat.size else 1.0
            p99_raw = max(p99_raw, 1e-6)
            methyl_gate_raw_threshold = float(np.percentile(all_concat, methyl_gate_quantile)) if all_concat.size else 0.0
        epigenetic_noise_floor_raw_threshold = float(methyl_gate_raw_threshold)
        p99_log1p = max(float(np.log1p(p99_raw)), 1e-6)
        methyl_gate_log1p_threshold = float(np.log1p(methyl_gate_raw_threshold))
        epigenetic_noise_floor_log1p_threshold = float(methyl_gate_log1p_threshold)
        methyl_gate_norm_threshold = float(methyl_gate_log1p_threshold / p99_log1p)
        epigenetic_noise_floor_norm_threshold = float(methyl_gate_norm_threshold)
        mark_summary[mark] = {
            "p99_raw": float(p99_raw),
            "p99_log1p": float(p99_log1p),
            "p99_new": float(p99_log1p),
            "methyl_gate_quantile": float(methyl_gate_quantile),
            "methyl_gate_raw_threshold": float(methyl_gate_raw_threshold),
            "methyl_gate_log1p_threshold": float(methyl_gate_log1p_threshold),
            "methyl_gate_norm_threshold": float(methyl_gate_norm_threshold),
            "epigenetic_noise_floor_quantile": float(epigenetic_noise_floor_quantile),
            "epigenetic_noise_floor_raw_threshold": float(epigenetic_noise_floor_raw_threshold),
            "epigenetic_noise_floor_log1p_threshold": float(epigenetic_noise_floor_log1p_threshold),
            "epigenetic_noise_floor_norm_threshold": float(epigenetic_noise_floor_norm_threshold),
            "threshold_log1p": float(methyl_gate_log1p_threshold),
            "threshold_norm": float(methyl_gate_norm_threshold),
        }

    shared_stats = {
        "included_marks": list(included_marks),
        "feature_marks": list(included_marks),
        "gate_marks": list(included_marks),
        "pair_gate_mode": "any_same_mark_match",
        "shared_feature_and_gate_threshold": True,
        "mark_stats": mark_summary,
    }
    if len(included_marks) == 1:
        primary = mark_summary[included_marks[0]]
        shared_stats.update(primary)
    for mark in included_marks:
        stats = mark_summary[mark]
        logger.info(
            "Computed shared Xa+Xi raw-scale stats for %s: p99_raw = %.6f | shared feature/gate q = %.2f raw %.6f norm %.6f",
            mark,
            stats["p99_raw"],
            methyl_gate_quantile,
            stats["methyl_gate_raw_threshold"],
            stats["methyl_gate_norm_threshold"],
        )
    return shared_stats, pd.DataFrame(chrom_rows)


def normalize_and_denoise_track(raw_track: np.ndarray, p99_log1p: float, epigenetic_noise_floor_raw_threshold: float) -> Tuple[np.ndarray, np.ndarray]:
    raw = np.nan_to_num(raw_track.astype(np.float32, copy=False), nan=0.0, posinf=0.0, neginf=0.0)
    raw_denoised = np.where(raw >= epigenetic_noise_floor_raw_threshold, raw, 0.0).astype(np.float32)
    scale = max(float(p99_log1p), 1e-6)
    norm_denoised = (np.log1p(raw_denoised) / scale).astype(np.float32)
    return raw_denoised, norm_denoised


def load_cooler(uri: Optional[str]):
    if uri is None:
        return None
    if cooler is None:
        raise RuntimeError("cooler is not installed, but --mcool-uri was provided")
    return cooler.Cooler(uri)


def load_microc_window(clr, chrom: str, start: int, end: int) -> np.ndarray:
    mat = clr.matrix(balance=True).fetch((chrom, start, end))
    return mat.astype(np.float32)


# ============================================================
# Matrix helpers matching training logic
# ============================================================
def normalize_microc(mat: np.ndarray, stratum_expected: np.ndarray) -> np.ndarray:
    n, eps = mat.shape[0], 1e-9
    oe = np.full((n, n), np.nan, dtype=np.float32)
    for s in range(1, n):
        exp = float(stratum_expected[s]) if s < len(stratum_expected) else eps
        exp = max(exp, eps)
        log2_oe = np.log2((np.diag(mat, k=s).astype(np.float32) / exp) + eps)
        idx = np.arange(n - s)
        oe[idx, idx + s] = oe[idx + s, idx] = log2_oe
    oe = np.clip(oe, TARGET_MIN, TARGET_MAX)
    return oe


def smooth_oe_matrix(oe_mat: np.ndarray, kernel_size: int = SMOOTH_KERNEL_SIZE) -> np.ndarray:
    n = oe_mat.shape[0]
    diag_matrix = np.abs(np.arange(n)[:, None] - np.arange(n)[None, :]).astype(np.int16)
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
    valid = np.isfinite(oe_mat) & (diag_matrix >= 1)
    weighted = np.where(valid, oe_mat, 0.0).astype(np.float32)
    denom = convolve(valid.astype(np.float32), kernel, mode="constant", cval=0.0)
    numer = convolve(weighted, kernel, mode="constant", cval=0.0)
    smoothed = np.full_like(oe_mat, np.nan, dtype=np.float32)
    np.divide(numer, np.maximum(denom, 1e-6), out=smoothed, where=denom > 0)
    smoothed[~valid] = np.nan
    np.fill_diagonal(smoothed, np.nan)
    return smoothed


def upper_indices(n_bins: int) -> Tuple[np.ndarray, np.ndarray]:
    return np.triu_indices(n_bins)


def vector_to_symmetric_matrix(values: np.ndarray, ii: np.ndarray, jj: np.ndarray, n_bins: int, fill=np.nan) -> np.ndarray:
    mat = np.full((n_bins, n_bins), fill, dtype=np.float32)
    mat[ii, jj] = values
    mat[jj, ii] = values
    np.fill_diagonal(mat, np.nan)
    return mat


# ============================================================
# Feature construction
# ============================================================
def build_feature_tracks(normalized_tracks_by_mark: Dict[str, np.ndarray], n_bins: int) -> Tuple[np.ndarray, np.ndarray]:
    tracks = np.zeros((N_MARKS, n_bins), dtype=np.float32)
    tracks_smooth = np.zeros((N_MARKS, n_bins), dtype=np.float32)
    for mark, normalized_track in normalized_tracks_by_mark.items():
        idx = MARK_TO_INDEX[mark]
        tracks[idx] = normalized_track.astype(np.float32, copy=False)
        tracks_smooth[idx] = gaussian_filter1d(normalized_track.astype(np.float32, copy=False), sigma=2.0)
    return tracks, tracks_smooth


def build_pair_features_for_window(tracks: np.ndarray, tracks_smooth: np.ndarray, resolution: int):
    n_bins = tracks.shape[1]
    ii, jj = upper_indices(n_bins)
    n_pairs = len(ii)

    cumsum_ext = np.zeros((N_MARKS, n_bins + 1), dtype=np.float32)
    cumsum_ext[:, 1:] = np.cumsum(tracks, axis=1)

    feats = np.empty((n_pairs, len(FEATURE_NAMES)), dtype=np.float32)
    feats[:, :N_MARKS] = tracks[:, ii].T
    feats[:, N_MARKS:2 * N_MARKS] = tracks[:, jj].T
    feats[:, 2 * N_MARKS:3 * N_MARKS] = (tracks[:, ii] * tracks[:, jj]).T

    span = (jj - ii + 1).astype(np.float32)
    int_sum = cumsum_ext[:, jj + 1] - cumsum_ext[:, ii]
    feats[:, 3 * N_MARKS:4 * N_MARKS] = (int_sum / np.where(span > 0, span, 1.0)).T
    feats[:, 4 * N_MARKS:5 * N_MARKS] = (tracks_smooth[:, ii] * tracks_smooth[:, jj]).T

    diag_idx = (jj - ii).astype(np.int32)
    feats[:, 5 * N_MARKS] = np.log1p(diag_idx.astype(np.float32) * resolution)
    feats[:, 5 * N_MARKS + 1] = diag_idx.astype(np.float32)

    valid_mask = np.all(np.isfinite(feats), axis=1)
    return feats, diag_idx, valid_mask, ii, jj


# ============================================================
# Calibration / probability / counts / gating
# ============================================================
def resolve_diag_correction_path(path: Optional[Path]) -> Optional[Path]:
    if path is None:
        return None
    if path.exists():
        return path
    fallback = path.with_name("exact_diag_calibration_curve.csv") if path.name == "exact_diag_correction.npy" else None
    if fallback is not None and fallback.exists():
        return fallback
    return path


def resolve_pred_level_correction_path(
    path: Optional[Path],
    diag_correction_path: Optional[Path] = None,
    disabled: bool = False,
) -> Optional[Path]:
    if disabled:
        return None
    if path is not None:
        return path if path.exists() else path
    if diag_correction_path is None:
        return None
    resolved_diag = resolve_diag_correction_path(diag_correction_path)
    if resolved_diag is None:
        return None
    sibling = resolved_diag.with_name("pred_level_calibration_curve.csv")
    return sibling if sibling.exists() else None


def load_diag_correction(path: Optional[Path], logger: Optional[logging.Logger] = None) -> Optional[np.ndarray]:
    requested_path = path
    path = resolve_diag_correction_path(path)
    if path is None:
        if logger is not None:
            logger.info("No exact-diagonal correction provided; using uncalibrated predicted smoothed O/E.")
        return None
    if not path.exists():
        raise FileNotFoundError(
            f"Could not find diagonal correction file: {path}. Accepted formats: "
            ".csv with columns diag_idx and smoothed_correction, .npy, .npz, or .json. "
            "For v5 calibrated models, use exact_diag_calibration_curve.csv from the calibrated variant directory."
        )
    if requested_path is not None and path != requested_path and logger is not None:
        logger.info("Requested correction %s was missing; using v5 CSV fallback %s", requested_path, path)
    suffix = path.suffix.lower()
    if suffix == ".npy":
        arr = np.asarray(np.load(path), dtype=np.float32)
        if logger is not None:
            logger.info("Loaded exact-diagonal correction from %s (%d values)", path, len(arr))
        return arr
    if suffix == ".npz":
        data = np.load(path)
        for key in ["correction", "smoothed_correction", "diag_correction", "arr_0"]:
            if key in data:
                arr = np.asarray(data[key], dtype=np.float32)
                if logger is not None:
                    logger.info("Loaded exact-diagonal correction from %s key=%s (%d values)", path, key, len(arr))
                return arr
        raise KeyError(f"Could not find a correction array in {path}")
    if suffix == ".json":
        payload = load_json(path)
        if isinstance(payload, list):
            arr = np.asarray(payload, dtype=np.float32)
            if logger is not None:
                logger.info("Loaded exact-diagonal correction from %s (%d values)", path, len(arr))
            return arr
        for key in ["correction", "smoothed_correction", "diag_correction"]:
            if key in payload:
                arr = np.asarray(payload[key], dtype=np.float32)
                if logger is not None:
                    logger.info("Loaded exact-diagonal correction from %s key=%s (%d values)", path, key, len(arr))
                return arr
        raise KeyError(f"Could not find a correction array in {path}")
    if suffix == ".csv":
        df = pd.read_csv(path)
        required = {"diag_idx", "smoothed_correction"}
        if not required.issubset(df.columns):
            raise ValueError(
                f"Correction CSV {path} must contain columns {sorted(required)}; found {df.columns.tolist()}"
            )
        diag = df["diag_idx"].astype(int).to_numpy()
        corr = df["smoothed_correction"].astype(np.float32).to_numpy()
        if diag.size == 0:
            raise ValueError(f"Correction CSV {path} contains no rows")
        out = np.zeros(int(np.nanmax(diag)) + 1, dtype=np.float32)
        out[diag] = corr
        out[0] = 0.0
        if logger is not None:
            logger.info("Loaded exact-diagonal correction from %s (%d diagonals, max diag=%d)", path, len(out), int(np.nanmax(diag)))
        return out
    raise ValueError(
        f"Unsupported correction file format: {path}. Accepted formats: .csv, .npy, .npz, .json"
    )


def load_pred_level_correction(
    path: Optional[Path],
    diag_correction_path: Optional[Path] = None,
    disabled: bool = False,
    logger: Optional[logging.Logger] = None,
) -> Optional[dict]:
    if disabled:
        if logger is not None:
            if path is not None:
                logger.info(
                    "Prediction-level correction was explicitly disabled; ignoring --pred-level-correction=%s and using exact-diagonal calibration only.",
                    path,
                )
            else:
                logger.info("Prediction-level correction was explicitly disabled; using exact-diagonal calibration only.")
        return None
    resolved = resolve_pred_level_correction_path(path, diag_correction_path, disabled=False)
    if resolved is None:
        if logger is not None:
            logger.info("No prediction-level correction provided; using exact-diagonal calibration only.")
        return None
    if not resolved.exists():
        raise FileNotFoundError(
            f"Could not find prediction-level correction file: {resolved}. "
            "Expected a CSV with columns pred_center and smoothed_correction."
        )
    df = pd.read_csv(resolved)
    required = {"pred_center", "smoothed_correction"}
    if not required.issubset(df.columns):
        raise ValueError(
            f"Prediction-level correction CSV {resolved} must contain columns {sorted(required)}; "
            f"found {df.columns.tolist()}"
        )
    centers = df["pred_center"].astype(np.float32).to_numpy()
    corr = df["smoothed_correction"].astype(np.float32).to_numpy()
    order = np.argsort(centers)
    centers = centers[order]
    corr = corr[order]
    if logger is not None:
        logger.info(
            "Loaded prediction-level correction from %s (%d support points, pred range %.3f..%.3f)",
            resolved,
            len(centers),
            float(centers[0]),
            float(centers[-1]),
        )
    return {
        "pred_center": centers,
        "smoothed_correction": corr,
        "path": str(resolved),
    }


def apply_diag_correction(pred_vec: np.ndarray, diag_idx: np.ndarray, diag_correction: Optional[np.ndarray]) -> np.ndarray:
    if diag_correction is None:
        return pred_vec.astype(np.float32, copy=False)
    clipped = np.clip(diag_idx.astype(np.int32), 0, len(diag_correction) - 1)
    return pred_vec.astype(np.float32, copy=False) - diag_correction[clipped].astype(np.float32, copy=False)


def apply_pred_level_correction(pred_vec: np.ndarray, pred_level_correction: Optional[dict]) -> np.ndarray:
    pred_vec = pred_vec.astype(np.float32, copy=False)
    if pred_level_correction is None:
        return pred_vec
    centers = np.asarray(pred_level_correction["pred_center"], dtype=np.float32)
    corr = np.asarray(pred_level_correction["smoothed_correction"], dtype=np.float32)
    if centers.size == 0:
        return pred_vec
    interp = np.interp(
        pred_vec.astype(np.float64),
        centers.astype(np.float64),
        corr.astype(np.float64),
        left=float(corr[0]),
        right=float(corr[-1]),
    ).astype(np.float32)
    return pred_vec - interp


def load_residual_baseline(path: Optional[Path], diag_idx: np.ndarray, required: bool, logger: Optional[logging.Logger] = None) -> Optional[np.ndarray]:
    if path is None:
        if required:
            raise ValueError(
                "This model metadata says target_mode=smoothed_oe_residual, so --residual-baseline is required. "
                "For the v5 HPC run, use checkpoints/stratum_residual_baseline.npz."
            )
        return None
    if not path.exists():
        if required:
            raise FileNotFoundError(
                f"Could not find residual baseline file: {path}. "
                "This v5 residual RF needs the train stratum median baseline to reconstruct predicted smoothed O/E."
            )
        if logger is not None:
            logger.warning("Residual baseline path %s does not exist; continuing without baseline because model does not require it", path)
        return None

    data = np.load(path)
    if "baseline" not in data:
        raise KeyError(f"Residual baseline file {path} must contain key 'baseline'")
    baseline = np.asarray(data["baseline"], dtype=np.float32)
    if baseline.ndim != 1:
        raise ValueError(f"Residual baseline in {path} must be 1D; found shape {baseline.shape}")
    max_diag = int(np.nanmax(diag_idx)) if len(diag_idx) else 0
    if len(baseline) <= max_diag:
        raise ValueError(
            f"Residual baseline in {path} has length {len(baseline)}, but predictions require diag_idx up to {max_diag}"
        )
    if logger is not None:
        logger.info("Loaded residual baseline from %s (%d diagonals, max required diag=%d)", path, len(baseline), max_diag)
    return baseline


def add_residual_baseline(pred_resid: np.ndarray, diag_idx: np.ndarray, baseline: Optional[np.ndarray]) -> np.ndarray:
    pred_resid = pred_resid.astype(np.float32, copy=False)
    if baseline is None:
        return pred_resid
    clipped = np.clip(diag_idx.astype(np.int32), 0, len(baseline) - 1)
    return pred_resid + baseline[clipped].astype(np.float32, copy=False)


def baseline_values_from_diag(diag_idx: np.ndarray, baseline: Optional[np.ndarray]) -> np.ndarray:
    if baseline is None:
        return np.zeros(len(diag_idx), dtype=np.float32)
    clipped = np.clip(diag_idx.astype(np.int32), 0, len(baseline) - 1)
    return baseline[clipped].astype(np.float32, copy=False)


def save_stage3_reconstruction_histogram(
    out_path: Path,
    rf_residual: np.ndarray,
    baseline_applied: np.ndarray,
    pred_uncalibrated: np.ndarray,
    pred_calibrated: np.ndarray,
    allele_name: str,
) -> None:
    fig, axes = plt.subplots(1, 4, figsize=(16, 3.8))
    panels = [
        (rf_residual, "RF residual"),
        (baseline_applied, "Baseline[diag]"),
        (pred_uncalibrated, "Residual + baseline"),
        (pred_calibrated, "Calibrated smoothed O/E"),
    ]
    for ax, (vals, title) in zip(axes, panels):
        finite = np.asarray(vals, dtype=np.float32)
        finite = finite[np.isfinite(finite)]
        if finite.size:
            ax.hist(finite, bins=80, color="slateblue", alpha=0.85)
            ax.axvline(float(np.median(finite)), color="firebrick", lw=1.2, ls="--", label="median")
            ax.legend(frameon=False, loc="upper right")
        ax.set_title(f"{allele_name}: {title}")
        ax.set_xlabel("log2 O/E units")
        ax.set_ylabel("Pairs")
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()


def build_pair_gate_from_tracks(
    raw_tracks_by_mark: Dict[str, np.ndarray],
    ii: np.ndarray,
    jj: np.ndarray,
    shared_stats: dict,
    included_marks: Sequence[str],
) -> Tuple[np.ndarray, Dict[str, np.ndarray], Dict[str, np.ndarray], np.ndarray, np.ndarray]:
    mark_summary = mark_stats_map(shared_stats)
    n_bins = len(next(iter(raw_tracks_by_mark.values()))) if raw_tracks_by_mark else 0
    gate_bins_by_mark: Dict[str, np.ndarray] = {}
    pair_gate_by_mark: Dict[str, np.ndarray] = {}
    pair_gate_any = np.zeros(len(ii), dtype=bool)
    for mark in included_marks:
        raw_track = raw_tracks_by_mark[mark]
        threshold = float(mark_summary[mark]["methyl_gate_raw_threshold"])
        gate_bins = raw_track >= threshold
        pair_gate = gate_bins[ii] & gate_bins[jj]
        pair_gate[diag_or_lower(ii, jj)] = False
        gate_bins_by_mark[mark] = gate_bins
        pair_gate_by_mark[mark] = pair_gate
        pair_gate_any |= pair_gate
    combined_gate_bins = (
        np.any(np.stack([gate_bins_by_mark[mark] for mark in included_marks], axis=0), axis=0)
        if included_marks else np.zeros(n_bins, dtype=bool)
    )
    gate_map = vector_to_symmetric_matrix(pair_gate_any.astype(np.float32), ii, jj, n_bins, fill=0.0)
    np.fill_diagonal(gate_map, 0.0)
    return combined_gate_bins, gate_bins_by_mark, pair_gate_by_mark, pair_gate_any, gate_map


def diag_bucket_label(lo: int, hi: int, n_bins: int) -> str:
    return f"{lo}-{hi - 1}" if hi < n_bins else f"{lo}+"


def distance_bucket_summary(values: np.ndarray, ii: np.ndarray, jj: np.ndarray, n_bins: int, value_name: str) -> pd.DataFrame:
    diag = jj - ii
    edges = DIAG_BUCKET_EDGES_DEFAULT + [n_bins]
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        mask = (diag >= lo) & (diag < hi)
        bucket_values = values[mask]
        rows.append({
            "bucket": diag_bucket_label(lo, hi, n_bins),
            "diag_min": int(lo),
            "diag_max_exclusive": int(hi),
            f"{value_name}_sum": float(np.nansum(bucket_values)),
            f"{value_name}_nonzero": int(np.count_nonzero(np.nan_to_num(bucket_values, nan=0.0))),
            "n_pairs": int(mask.sum()),
        })
    return pd.DataFrame(rows)


def summarize_upper_triangle_by_bucket(
    upper_values: np.ndarray,
    diag: np.ndarray,
    n_bins: int,
    value_name: str,
    include_overall: bool = True,
) -> pd.DataFrame:
    upper_values = np.asarray(upper_values, dtype=np.float64)
    diag = np.asarray(diag, dtype=np.int32)
    finite = np.isfinite(upper_values)
    rows = []

    def _row(bucket: str, mask: np.ndarray, lo: int, hi: int):
        vals = upper_values[mask]
        finite_vals = vals[np.isfinite(vals)]
        rows.append({
            "metric": value_name,
            "bucket": bucket,
            "diag_min": int(lo),
            "diag_max_exclusive": int(hi),
            "n_pairs": int(mask.sum()),
            "n_finite": int(finite_vals.size),
            "sum": float(np.nansum(vals)),
            "mean": float(np.nanmean(vals)) if finite_vals.size else np.nan,
            "median": float(np.nanmedian(vals)) if finite_vals.size else np.nan,
            "max": float(np.nanmax(vals)) if finite_vals.size else np.nan,
            "nonzero": int(np.count_nonzero(np.nan_to_num(vals, nan=0.0))),
        })

    if include_overall:
        _row("overall", finite, int(diag[finite].min()) if finite.any() else 0, int(diag[finite].max() + 1) if finite.any() else 0)

    edges = DIAG_BUCKET_EDGES_DEFAULT + [n_bins]
    for lo, hi in zip(edges[:-1], edges[1:]):
        mask = finite & (diag >= lo) & (diag < hi)
        _row(diag_bucket_label(lo, hi, n_bins), mask, lo, hi)
    return pd.DataFrame(rows)


def summarize_by_exact_diag(
    values: np.ndarray,
    diag: np.ndarray,
    allowed: np.ndarray,
    value_name: str,
) -> pd.DataFrame:
    values = np.asarray(values, dtype=np.float64)
    diag = np.asarray(diag, dtype=np.int32)
    allowed = np.asarray(allowed, dtype=bool)
    rows = []
    for d in np.unique(diag[allowed]):
        mask = allowed & (diag == d)
        diag_vals = values[mask]
        finite_vals = diag_vals[np.isfinite(diag_vals)]
        rows.append({
            "metric": value_name,
            "diag_idx": int(d),
            "n_pairs": int(mask.sum()),
            "mean": float(np.mean(finite_vals)) if finite_vals.size else np.nan,
            "median": float(np.median(finite_vals)) if finite_vals.size else np.nan,
            "max": float(np.max(finite_vals)) if finite_vals.size else np.nan,
        })
    return pd.DataFrame(rows)


def statistic_from_name(values: np.ndarray, statistic: str) -> float:
    values = np.asarray(values, dtype=np.float64)
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return float("nan")
    if statistic == "median":
        return float(np.median(finite))
    if statistic == "mean":
        return float(np.mean(finite))
    raise ValueError(f"Unsupported statistic: {statistic}")


def build_feasible_pair_statistics(
    fiber_bin_to_nucs: Optional[Dict[str, Dict[int, List[int]]]],
    n_bins: int,
) -> Tuple[np.ndarray, np.ndarray]:
    feasible_count_map = np.zeros((n_bins, n_bins), dtype=np.int32)
    if not fiber_bin_to_nucs:
        return feasible_count_map > 0, feasible_count_map

    for bin_to_nucs in fiber_bin_to_nucs.values():
        bins = sorted(int(b) for b in bin_to_nucs.keys() if 0 <= int(b) < n_bins)
        for pos, bi in enumerate(bins):
            for bj in bins[pos + 1:]:
                feasible_count_map[bi, bj] += 1
                feasible_count_map[bj, bi] += 1

    np.fill_diagonal(feasible_count_map, 0)
    feasible_mask = feasible_count_map > 0
    np.fill_diagonal(feasible_mask, False)
    return feasible_mask, feasible_count_map


def save_gate_diagnostics(
    table_dir: Path,
    plot_dir: Path,
    raw_tracks_by_mark: Dict[str, np.ndarray],
    normalized_denoised_tracks_by_mark: Dict[str, np.ndarray],
    gate_bins_by_mark: Dict[str, np.ndarray],
    pair_gate_by_mark: Dict[str, np.ndarray],
    combined_gate_bins: np.ndarray,
    pair_gate: np.ndarray,
    ii: np.ndarray,
    jj: np.ndarray,
    prediction_region: Region,
    allele_name: str,
    included_marks: Sequence[str],
):
    n_bins = len(combined_gate_bins)
    bin_idx = np.arange(n_bins, dtype=np.int32)
    starts = prediction_region.start + bin_idx * REQUESTED_RESOLUTION
    gate_df = pd.DataFrame({
        "bin_idx": bin_idx,
        "chrom": prediction_region.chrom,
        "start": starts,
        "end": starts + REQUESTED_RESOLUTION,
        "center": starts + REQUESTED_RESOLUTION // 2,
        "gated_any_mark": combined_gate_bins.astype(bool),
    })
    for mark in included_marks:
        safe = mark.replace(" ", "_")
        gate_df[f"raw_signal_{safe}"] = raw_tracks_by_mark[mark].astype(float)
        gate_df[f"normalized_denoised_signal_{safe}"] = normalized_denoised_tracks_by_mark[mark].astype(float)
        gate_df[f"gated_{safe}"] = gate_bins_by_mark[mark].astype(bool)
    gate_df.to_csv(table_dir / "stage3_gate_bins.csv", index=False)

    pair_frames = [distance_bucket_summary(pair_gate.astype(np.float32), ii, jj, n_bins, "gated_pair").assign(mark="any_same_mark_match")]
    for mark in included_marks:
        pair_frames.append(distance_bucket_summary(pair_gate_by_mark[mark].astype(np.float32), ii, jj, n_bins, "gated_pair").assign(mark=mark))
    pair_summary = pd.concat(pair_frames, ignore_index=True)
    pair_summary.to_csv(table_dir / "stage3_gated_pair_distance_summary.csv", index=False)

    fig, ax = plt.subplots(figsize=(8, 4))
    pivot = pair_summary.pivot(index="bucket", columns="mark", values="gated_pair_nonzero").fillna(0.0)
    x = np.arange(len(pivot.index))
    marks_for_plot = ["any_same_mark_match"] + [mark for mark in included_marks if mark in pivot.columns]
    width = 0.8 / max(len(marks_for_plot), 1)
    for idx, mark in enumerate(marks_for_plot):
        vals = pivot[mark].to_numpy(dtype=float) if mark in pivot.columns else np.zeros(len(x), dtype=float)
        ax.bar(x + (idx - (len(marks_for_plot) - 1) / 2.0) * width, vals, width=width, color=MARK_COLORS.get(mark, "firebrick") if mark != "any_same_mark_match" else "black", alpha=0.85, label=mark)
    ax.set_xticks(x)
    ax.set_xticklabels(pivot.index.tolist(), rotation=45, ha="right")
    ax.set_title(f"{allele_name}: self-similar gated pair counts by distance bucket")
    ax.set_xlabel("Distance bucket (200 bp bins)")
    ax.set_ylabel("Allowed gated pairs")
    ax.legend(frameon=False, fontsize=8, ncol=min(2, len(marks_for_plot)))
    plt.tight_layout()
    plt.savefig(plot_dir / "stage3_gated_pair_distance_summary.png")
    plt.close()


def diag_or_lower(ii: np.ndarray, jj: np.ndarray) -> np.ndarray:
    return jj <= ii


def build_absolute_probability_maps(
    pred_map: np.ndarray,
    stratum_expected: np.ndarray,
    min_diag_bins: int,
    polymer_alpha: float,
    local_diag_min: int,
    local_diag_max: int,
    local_statistic: str,
    p_local: float,
    kappa: float,
):
    # The RF predicts smoothed log2(O/E). We convert that back to an observed-like
    # intensity field T using the same stratum expected curve used during training,
    # then compare it against a scaled polymer baseline S. Both are pushed through
    # the same bounded occupancy map so their difference can be interpreted as an
    # absolute extra-contact probability rather than a globally normalized weight.
    n = pred_map.shape[0]
    ii, jj = np.triu_indices(n, k=max(1, int(min_diag_bins)))
    diag = (jj - ii).astype(np.int32)
    pred_vec = pred_map[ii, jj].astype(np.float64)
    finite = np.isfinite(pred_vec)

    expected_vec = np.full_like(pred_vec, np.nan, dtype=np.float64)
    valid_diag = diag < len(stratum_expected)
    expected_vec[valid_diag] = stratum_expected[diag[valid_diag]].astype(np.float64)
    allowed = finite & valid_diag & np.isfinite(expected_vec) & (expected_vec > 0)

    t_vec = np.zeros_like(pred_vec, dtype=np.float64)
    t_vec[allowed] = expected_vec[allowed] * np.power(2.0, pred_vec[allowed])

    s_raw_vec = np.zeros_like(pred_vec, dtype=np.float64)
    positive_diag = diag >= 1
    s_raw_vec[positive_diag] = np.power(diag[positive_diag].astype(np.float64), -float(polymer_alpha))

    local_lo = max(int(min_diag_bins), int(local_diag_min))
    local_hi = int(local_diag_max)
    local_mask = allowed & (diag >= local_lo) & (diag <= local_hi)
    if not np.any(local_mask):
        local_mask = allowed

    local_t_stat = statistic_from_name(t_vec[local_mask], local_statistic)
    local_s_raw_stat = statistic_from_name(s_raw_vec[local_mask], local_statistic)
    if np.isfinite(local_t_stat) and np.isfinite(local_s_raw_stat) and local_s_raw_stat > 0:
        polymer_scale = float(local_t_stat / local_s_raw_stat)
    else:
        polymer_scale = 0.0

    s_vec = np.zeros_like(pred_vec, dtype=np.float64)
    if polymer_scale > 0:
        s_vec[allowed] = polymer_scale * s_raw_vec[allowed]

    local_s_stat = statistic_from_name(s_vec[local_mask], local_statistic)
    p_local_clipped = min(max(float(p_local), 1e-6), 1.0 - 1e-6)
    if np.isfinite(local_s_stat) and local_s_stat > 0:
        lambda_value = float(-math.log1p(-p_local_clipped) / local_s_stat)
    else:
        lambda_value = 0.0

    q_t_vec = np.zeros_like(pred_vec, dtype=np.float64)
    q_s_vec = np.zeros_like(pred_vec, dtype=np.float64)
    if lambda_value > 0:
        q_t_vec[allowed] = 1.0 - np.exp(-lambda_value * np.clip(t_vec[allowed], 0.0, None))
        q_s_vec[allowed] = 1.0 - np.exp(-lambda_value * np.clip(s_vec[allowed], 0.0, None))

    p_extra_vec = np.zeros_like(pred_vec, dtype=np.float64)
    p_extra_vec[allowed] = np.maximum(0.0, q_t_vec[allowed] - float(kappa) * q_s_vec[allowed])
    subtraction_vec = np.zeros_like(pred_vec, dtype=np.float64)
    subtraction_vec[allowed] = q_t_vec[allowed] - float(kappa) * q_s_vec[allowed]

    allowed_map = vector_to_symmetric_matrix(allowed.astype(np.float32), ii, jj, n, fill=0.0)
    t_map = vector_to_symmetric_matrix(t_vec.astype(np.float32), ii, jj, n, fill=0.0)
    s_raw_map = vector_to_symmetric_matrix(s_raw_vec.astype(np.float32), ii, jj, n, fill=0.0)
    s_map = vector_to_symmetric_matrix(s_vec.astype(np.float32), ii, jj, n, fill=0.0)
    q_t_map = vector_to_symmetric_matrix(q_t_vec.astype(np.float32), ii, jj, n, fill=0.0)
    q_s_map = vector_to_symmetric_matrix(q_s_vec.astype(np.float32), ii, jj, n, fill=0.0)
    p_extra_map = vector_to_symmetric_matrix(p_extra_vec.astype(np.float32), ii, jj, n, fill=0.0)
    subtraction_map = vector_to_symmetric_matrix(subtraction_vec.astype(np.float32), ii, jj, n, fill=0.0)
    diag_map = vector_to_symmetric_matrix(diag.astype(np.float32), ii, jj, n, fill=0.0)
    for mat in [allowed_map, t_map, s_raw_map, s_map, q_t_map, q_s_map, p_extra_map, subtraction_map, diag_map]:
        np.fill_diagonal(mat, 0.0)

    outputs = {
        "ii": ii,
        "jj": jj,
        "diag": diag,
        "allowed": allowed,
        "allowed_map": allowed_map,
        "expected_vec": expected_vec.astype(np.float32),
        "t_vec": t_vec.astype(np.float32),
        "s_raw_vec": s_raw_vec.astype(np.float32),
        "s_vec": s_vec.astype(np.float32),
        "q_t_vec": q_t_vec.astype(np.float32),
        "q_s_vec": q_s_vec.astype(np.float32),
        "p_extra_vec": p_extra_vec.astype(np.float32),
        "subtraction_vec": subtraction_vec.astype(np.float32),
        "diag_map": diag_map,
        "t_map": t_map,
        "s_raw_map": s_raw_map,
        "s_map": s_map,
        "q_t_map": q_t_map,
        "q_s_map": q_s_map,
        "p_extra_map": p_extra_map,
        "subtraction_map": subtraction_map,
        "polymer_scale": polymer_scale,
        "polymer_alpha": float(polymer_alpha),
        "lambda_value": lambda_value,
        "local_diag_min_used": local_lo,
        "local_diag_max_used": local_hi,
        "local_t_stat": local_t_stat,
        "local_s_raw_stat": local_s_raw_stat,
        "local_s_stat": local_s_stat,
        "p_local_clipped": p_local_clipped,
    }
    return outputs


def stage5_bucket_summary(
    maps: Dict[str, np.ndarray],
    diag: np.ndarray,
    allowed: np.ndarray,
    n_bins: int,
) -> pd.DataFrame:
    rows = []
    for metric_name, values in maps.items():
        values = np.asarray(values, dtype=np.float64)
        masked_values = np.full_like(values, np.nan, dtype=np.float64)
        masked_values[allowed] = values[allowed]
        rows.append(summarize_upper_triangle_by_bucket(masked_values, diag, n_bins, metric_name))
    return pd.concat(rows, ignore_index=True)


def p_extra_histogram(p_extra_vec: np.ndarray, allowed: np.ndarray, bins: int = 50) -> pd.DataFrame:
    vals = np.asarray(p_extra_vec, dtype=np.float64)[allowed]
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return pd.DataFrame({"bin_left": [], "bin_right": [], "count": []})
    counts, edges = np.histogram(vals, bins=bins, range=(0.0, max(1e-6, float(vals.max()))))
    return pd.DataFrame({
        "bin_left": edges[:-1],
        "bin_right": edges[1:],
        "count": counts.astype(int),
    })


# ============================================================
# dim.in parsing and contact assignment
# ============================================================
def parse_dim_file(dim_path: Path) -> Tuple[int, List[int]]:
    lines = [line.strip() for line in dim_path.read_text().splitlines() if line.strip()]
    try:
        cores_idx = lines.index("#cores")
        n_cores = int(lines[cores_idx + 1])
        dna_idx = lines.index("#DNA")
    except Exception as exc:
        raise ValueError(f"Could not parse dim.in structure in {dim_path}") from exc

    dna_vals = []
    for line in lines[dna_idx + 1: dna_idx + 1 + n_cores]:
        dna_vals.append(int(line.split()[0]))
    if len(dna_vals) != n_cores:
        raise ValueError(f"Expected {n_cores} DNA entries in {dim_path}, found {len(dna_vals)}")
    return n_cores, dna_vals


def dim_to_core_centers(dim_path: Path, region_start_bp: int) -> np.ndarray:
    n_cores, num_dna = parse_dim_file(dim_path)
    xcore = np.zeros(n_cores, dtype=np.int64)
    xcore[0] = int(region_start_bp + FIRST_CORE_OFFSET_BP)
    for i in range(n_cores - 1):
        xcore[i + 1] = int(xcore[i] + (num_dna[i] + 1) * LINKER_BP + NUCLEOSOME_BP)
    return xcore


def fibers_from_dim_files(dim_paths: Sequence[Path], requested_region: Region, fiber_region_start_bp: int, resolution: int = REQUESTED_RESOLUTION):
    rows = []
    bins_by_fiber: Dict[str, Dict[int, List[int]]] = {}
    centers_by_fiber: Dict[str, np.ndarray] = {}

    for dim_path in dim_paths:
        fiber_id = dim_path.stem
        centers = dim_to_core_centers(dim_path, fiber_region_start_bp)
        centers_by_fiber[fiber_id] = centers
        bins_by_fiber[fiber_id] = {}
        for nuc_id, center in enumerate(centers):
            if requested_region.start <= center < requested_region.end:
                bin_id = int((center - requested_region.start) // resolution)
                rows.append({
                    "fiber_id": fiber_id,
                    "dim_path": str(dim_path),
                    "nuc_id": int(nuc_id),
                    "center_bp": int(center),
                    "bin_id": int(bin_id),
                })
                bins_by_fiber[fiber_id].setdefault(int(bin_id), []).append(int(nuc_id))
    return pd.DataFrame(rows), bins_by_fiber, centers_by_fiber


def build_fiber_assignment_index(
    fiber_bin_to_nucs: Dict[str, Dict[int, List[int]]],
    centers_by_fiber: Dict[str, np.ndarray],
    n_bins: int,
) -> Tuple[List[str], List[Dict[int, List[int]]], List[np.ndarray], Dict[int, List[int]]]:
    fiber_ids = sorted(fiber_bin_to_nucs.keys(), key=natural_sort_key)
    fiber_bins_by_index: List[Dict[int, List[int]]] = []
    fiber_centers_by_index: List[np.ndarray] = []
    bin_to_fiber_indices: Dict[int, List[int]] = {}

    for fiber_idx, fiber_id in enumerate(fiber_ids):
        raw_bin_map = fiber_bin_to_nucs.get(fiber_id, {})
        clean_bin_map: Dict[int, List[int]] = {}
        for bin_id, nuc_ids in raw_bin_map.items():
            bin_i = int(bin_id)
            if 0 <= bin_i < n_bins:
                clean_bin_map[bin_i] = [int(n) for n in nuc_ids]
                bin_to_fiber_indices.setdefault(bin_i, []).append(int(fiber_idx))
        fiber_bins_by_index.append(clean_bin_map)
        fiber_centers_by_index.append(np.asarray(centers_by_fiber[fiber_id], dtype=np.int64))

    return fiber_ids, fiber_bins_by_index, fiber_centers_by_index, bin_to_fiber_indices


def candidate_fiber_indices_for_pair(
    bi: int,
    bj: int,
    bin_to_fiber_indices: Dict[int, List[int]],
) -> List[int]:
    left = bin_to_fiber_indices.get(int(bi), [])
    right = bin_to_fiber_indices.get(int(bj), [])
    if not left or not right:
        return []
    if len(left) > len(right):
        left, right = right, left
    right_set = set(int(idx) for idx in right)
    return [int(idx) for idx in left if int(idx) in right_set]


def apply_feasible_mask(pred_map: np.ndarray, feasible_mask: Optional[np.ndarray]) -> np.ndarray:
    if feasible_mask is None:
        return pred_map
    out = pred_map.copy()
    out[~feasible_mask] = np.nan
    return out


def available_nucs(nucs: Sequence[int], valency: Dict[int, int], vmax: int) -> List[int]:
    return [int(n) for n in nucs if valency.get(int(n), 0) < vmax]


def choose_lowest_valency_nuc(candidates: Sequence[int], valency: Dict[int, int], rng: np.random.Generator) -> int:
    if not candidates:
        raise ValueError("No nucleosome candidates provided")
    loads = np.array([valency.get(int(n), 0) for n in candidates], dtype=np.int32)
    min_load = int(loads.min())
    pool = [int(n) for n, load in zip(candidates, loads) if int(load) == min_load]
    return int(pool[int(rng.integers(0, len(pool)))])


def clipped_logit(p: float, eps: float = 1e-6) -> float:
    p_clip = min(max(float(p), eps), 1.0 - eps)
    return float(np.log(p_clip / (1.0 - p_clip)))


def assign_contacts_capacity_aware(
    p_extra_map: np.ndarray,
    fiber_ids: Sequence[str],
    fiber_bin_to_nucs_by_index: Sequence[Dict[int, List[int]]],
    centers_by_fiber_by_index: Sequence[np.ndarray],
    bin_to_fiber_indices: Dict[int, List[int]],
    vmax: int,
    rng: np.random.Generator,
):
    # P_extra(i, j) is interpreted as the absolute probability that a feasible
    # conformation contains one extra methylation-driven contact for that pair.
    # We therefore sample each feasible fiber independently with Bernoulli(P_extra),
    # then resolve competition within each fiber using stochastic edge priorities
    # so high-probability short-range pairs do not consume capacity globally first.
    n = p_extra_map.shape[0]
    assignments = []
    skipped_rows = []
    sampled_map = np.zeros_like(p_extra_map, dtype=np.int32)
    recon = np.zeros_like(p_extra_map, dtype=np.int32)
    valency: List[Dict[int, int]] = [{} for _ in fiber_ids]
    candidate_edges_by_fiber: Dict[int, List[dict]] = {}
    event_ids_by_fiber: Dict[int, List[int]] = {}
    event_meta: Dict[int, dict] = {}
    pair_stats: Dict[Tuple[int, int], dict] = {}
    fiber_summary: Dict[int, dict] = {
        fiber_idx: {
            "fiber_id": str(fiber_ids[fiber_idx]),
            "sampled_pair_events": 0,
            "candidate_edges": 0,
            "assigned_contacts": 0,
            "conflict_rejected_events": 0,
            "no_valid_edge_events": 0,
        }
        for fiber_idx in range(len(fiber_ids))
    }

    event_id = 0
    for bi, bj in zip(*np.triu_indices(n, k=1)):
        p_extra = float(p_extra_map[bi, bj])
        if p_extra <= 0:
            continue
        candidate_fiber_indices = candidate_fiber_indices_for_pair(bi, bj, bin_to_fiber_indices)
        if not candidate_fiber_indices:
            continue
        sampled_fiber_indices = [
            fiber_idx
            for fiber_idx in candidate_fiber_indices
            if rng.random() < min(max(p_extra, 0.0), 1.0)
        ]
        sampled_count = len(sampled_fiber_indices)
        sampled_map[bi, bj] = sampled_count
        sampled_map[bj, bi] = sampled_count
        if sampled_count == 0:
            continue

        pair_key = (int(bi), int(bj))
        pair_stats[pair_key] = {
            "p_extra": float(p_extra),
            "feasible_fibers": int(len(candidate_fiber_indices)),
            "sampled_contacts": int(sampled_count),
            "assigned_contacts": 0,
            "nucleosome_conflict_after_sampling": 0,
            "no_valid_nucleosome_pair_in_fiber": 0,
        }
        base_score = clipped_logit(p_extra)

        for fiber_idx in sampled_fiber_indices:
            fiber_summary[fiber_idx]["sampled_pair_events"] += 1
            bin_to_nucs = fiber_bin_to_nucs_by_index[fiber_idx]
            left = [int(n) for n in bin_to_nucs.get(int(bi), [])]
            right = [int(n) for n in bin_to_nucs.get(int(bj), [])]
            edge_candidates = []
            event_meta[event_id] = {
                "fiber_idx": int(fiber_idx),
                "pair_key": pair_key,
                "bi": int(bi),
                "bj": int(bj),
                "p_extra": float(p_extra),
                "status": "pending",
                "edge_count": 0,
            }
            if left and right:
                event_score = base_score + float(rng.gumbel())
                for nuc_i in left:
                    for nuc_j in right:
                        if int(nuc_i) == int(nuc_j):
                            continue
                        edge_candidates.append({
                            "event_id": int(event_id),
                            "fiber_idx": int(fiber_idx),
                            "bi": int(bi),
                            "bj": int(bj),
                            "nuc_i": int(nuc_i),
                            "nuc_j": int(nuc_j),
                            "score": float(event_score + 1e-6 * rng.standard_normal()),
                        })
            if not edge_candidates:
                event_meta[event_id]["status"] = "no_valid_edge"
                pair_stats[pair_key]["no_valid_nucleosome_pair_in_fiber"] += 1
                fiber_summary[fiber_idx]["no_valid_edge_events"] += 1
                event_id += 1
                continue

            event_meta[event_id]["edge_count"] = int(len(edge_candidates))
            candidate_edges_by_fiber.setdefault(int(fiber_idx), []).extend(edge_candidates)
            event_ids_by_fiber.setdefault(int(fiber_idx), []).append(int(event_id))
            fiber_summary[fiber_idx]["candidate_edges"] += int(len(edge_candidates))
            event_id += 1

    for fiber_idx, edge_rows in candidate_edges_by_fiber.items():
        edge_rows.sort(key=lambda row: row["score"], reverse=True)
        used_nucs: set[int] = set()
        for edge in edge_rows:
            meta = event_meta[edge["event_id"]]
            if meta["status"] != "pending":
                continue
            nuc_i = int(edge["nuc_i"])
            nuc_j = int(edge["nuc_j"])
            if nuc_i in used_nucs or nuc_j in used_nucs:
                continue
            used_nucs.add(nuc_i)
            used_nucs.add(nuc_j)
            meta["status"] = "assigned"
            fiber_valency = valency[fiber_idx]
            fiber_valency[nuc_i] = fiber_valency.get(nuc_i, 0) + 1
            fiber_valency[nuc_j] = fiber_valency.get(nuc_j, 0) + 1
            bi = int(edge["bi"])
            bj = int(edge["bj"])
            recon[bi, bj] += 1
            recon[bj, bi] += 1
            pair_key = meta["pair_key"]
            pair_stats[pair_key]["assigned_contacts"] += 1
            fiber_summary[fiber_idx]["assigned_contacts"] += 1
            fiber_id = str(fiber_ids[fiber_idx])
            assignments.append({
                "fiber_id": fiber_id,
                "bin_i": bi,
                "bin_j": bj,
                "nuc_i": nuc_i,
                "nuc_j": nuc_j,
                "center_i_bp": int(centers_by_fiber_by_index[fiber_idx][nuc_i]),
                "center_j_bp": int(centers_by_fiber_by_index[fiber_idx][nuc_j]),
            })

        for event_idx in event_ids_by_fiber.get(fiber_idx, []):
            meta = event_meta[event_idx]
            if meta["status"] == "pending":
                meta["status"] = "conflict"
                pair_stats[meta["pair_key"]]["nucleosome_conflict_after_sampling"] += 1
                fiber_summary[fiber_idx]["conflict_rejected_events"] += 1

    for (bi, bj), stats in pair_stats.items():
        for reason in ["nucleosome_conflict_after_sampling", "no_valid_nucleosome_pair_in_fiber"]:
            n_failed = int(stats.get(reason, 0))
            if n_failed <= 0:
                continue
            skipped_rows.append({
                "bin_i": int(bi),
                "bin_j": int(bj),
                "p_extra": float(stats["p_extra"]),
                "feasible_fibers": int(stats["feasible_fibers"]),
                "sampled_contacts": int(stats["sampled_contacts"]),
                "assigned_contacts": int(stats["assigned_contacts"]),
                "failed_contacts": int(n_failed),
                "reason": reason,
            })

    valency_rows = []
    for fiber_idx, vmap in enumerate(valency):
        fiber_id = str(fiber_ids[fiber_idx])
        for nuc_id, v in vmap.items():
            valency_rows.append({"fiber_id": fiber_id, "nuc_id": int(nuc_id), "valency": int(v)})
    fiber_summary_df = pd.DataFrame(fiber_summary.values()).sort_values("fiber_id")
    return (
        pd.DataFrame(assignments),
        sampled_map,
        recon,
        pd.DataFrame(skipped_rows),
        pd.DataFrame(valency_rows),
        fiber_summary_df,
    )


def occupancy_from_counts(count_map: np.ndarray, feasible_count_map: np.ndarray) -> np.ndarray:
    out = np.zeros_like(count_map, dtype=np.float32)
    mask = feasible_count_map > 0
    out[mask] = count_map[mask].astype(np.float32) / feasible_count_map[mask].astype(np.float32)
    np.fill_diagonal(out, 0.0)
    return out


def assignment_efficiency_summary(
    p_extra_map: np.ndarray,
    feasible_count_map: np.ndarray,
    sampled_count_map: np.ndarray,
    recon_map: np.ndarray,
) -> pd.DataFrame:
    n = p_extra_map.shape[0]
    ii, jj = np.triu_indices(n, k=1)
    target = p_extra_map[ii, jj].astype(float)
    feasible = feasible_count_map[ii, jj].astype(float)
    sampled = sampled_count_map[ii, jj].astype(float)
    recon = recon_map[ii, jj].astype(float)
    sampled_occ = np.divide(sampled, feasible, out=np.zeros_like(sampled), where=feasible > 0)
    recon_occ = np.divide(recon, feasible, out=np.zeros_like(recon), where=feasible > 0)
    rows = [{
        "bucket": "overall",
        "mean_target_occupancy": float(np.nanmean(target)),
        "mean_sampled_occupancy": float(np.nanmean(sampled_occ)),
        "mean_reconstructed_occupancy": float(np.nanmean(recon_occ)),
        "occupancy_gap": float(np.nanmean(recon_occ - target)),
        "sampled_contacts": int(np.sum(sampled)),
        "reconstructed_contacts": int(np.sum(recon)),
        "feasible_fibers": int(np.sum(feasible)),
    }]
    diag = (jj - ii).astype(np.int32)
    edges = DIAG_BUCKET_EDGES_DEFAULT + [n]
    for lo, hi in zip(edges[:-1], edges[1:]):
        mask = (diag >= lo) & (diag < hi)
        if not np.any(mask):
            continue
        rows.append({
            "bucket": diag_bucket_label(lo, hi, n),
            "mean_target_occupancy": float(np.nanmean(target[mask])),
            "mean_sampled_occupancy": float(np.nanmean(sampled_occ[mask])),
            "mean_reconstructed_occupancy": float(np.nanmean(recon_occ[mask])),
            "occupancy_gap": float(np.nanmean(recon_occ[mask] - target[mask])),
            "sampled_contacts": int(np.sum(sampled[mask])),
            "reconstructed_contacts": int(np.sum(recon[mask])),
            "feasible_fibers": int(np.sum(feasible[mask])),
        })
    return pd.DataFrame(rows)


def methyl_output_path_for_dim(dim_path: Path, methyl_output_dir: Optional[Path]) -> Path:
    stem = dim_path.stem
    methyl_stem = f"methyl{stem[3:]}" if stem.startswith("dim") else f"methyl_{stem}"
    out_dir = ensure_dir(methyl_output_dir) if methyl_output_dir is not None else dim_path.parent
    return out_dir / f"{methyl_stem}{dim_path.suffix}"


def write_methyl_contact_files(assignments_df: pd.DataFrame, dim_paths: Sequence[Path], methyl_output_dir: Optional[Path]) -> List[str]:
    written = []
    grouped = {fiber_id: sub for fiber_id, sub in assignments_df.groupby("fiber_id")} if not assignments_df.empty else {}
    for dim_path in dim_paths:
        out_path = methyl_output_path_for_dim(dim_path, methyl_output_dir)
        sub = grouped.get(dim_path.stem, pd.DataFrame())
        with open(out_path, "w") as fh:
            if sub.empty:
                fh.write("# no methyl contacts assigned\n")
            else:
                for _, row in sub.iterrows():
                    fh.write(f"{int(row['nuc_i'])} {int(row['nuc_j'])}\n")
        written.append(str(out_path))
    return written


def assign_contacts_randomly(count_map: np.ndarray, fiber_bin_to_nucs: Dict[str, Dict[int, List[int]]], centers_by_fiber: Dict[str, np.ndarray], vmax: int, rng: np.random.Generator, max_attempts_per_contact: int = 50):
    fibers = sorted(fiber_bin_to_nucs.keys(), key=natural_sort_key)
    valency = {fiber_id: {} for fiber_id in fibers}
    assignments = []
    skipped_rows = []
    recon = np.zeros_like(count_map, dtype=np.int32)

    n = count_map.shape[0]
    for bi in range(n):
        for bj in range(bi + 1, n):
            target = int(count_map[bi, bj])
            if target <= 0:
                continue
            assigned = 0
            failed = 0
            for _ in range(target):
                success = False
                for _attempt in range(max_attempts_per_contact):
                    fiber_id = fibers[int(rng.integers(0, len(fibers)))]
                    left = fiber_bin_to_nucs[fiber_id].get(int(bi), [])
                    right = fiber_bin_to_nucs[fiber_id].get(int(bj), [])
                    if not left or not right:
                        continue
                    nuc_i = int(left[int(rng.integers(0, len(left)))])
                    nuc_j = int(right[int(rng.integers(0, len(right)))])
                    if nuc_i == nuc_j:
                        continue
                    v_i = valency[fiber_id].get(nuc_i, 0)
                    v_j = valency[fiber_id].get(nuc_j, 0)
                    if v_i >= vmax or v_j >= vmax:
                        continue
                    valency[fiber_id][nuc_i] = v_i + 1
                    valency[fiber_id][nuc_j] = v_j + 1
                    recon[bi, bj] += 1
                    recon[bj, bi] += 1
                    assignments.append({
                        "fiber_id": fiber_id,
                        "bin_i": int(bi),
                        "bin_j": int(bj),
                        "nuc_i": nuc_i,
                        "nuc_j": nuc_j,
                        "center_i_bp": int(centers_by_fiber[fiber_id][nuc_i]),
                        "center_j_bp": int(centers_by_fiber[fiber_id][nuc_j]),
                    })
                    assigned += 1
                    success = True
                    break
                if not success:
                    failed += 1
            if failed > 0:
                skipped_rows.append({
                    "bin_i": int(bi),
                    "bin_j": int(bj),
                    "target_contacts": int(target),
                    "assigned_contacts": int(assigned),
                    "failed_contacts": int(failed),
                })

    valency_rows = []
    for fiber_id, vmap in valency.items():
        for nuc_id, v in vmap.items():
            valency_rows.append({"fiber_id": fiber_id, "nuc_id": int(nuc_id), "valency": int(v)})
    return pd.DataFrame(assignments), recon, pd.DataFrame(skipped_rows), pd.DataFrame(valency_rows)


# ============================================================
# Plotting helpers
# ============================================================
def save_combined_normalization_plots(out_dir: Path, chrom_stats: pd.DataFrame, shared_stats: dict):
    ensure_dir(out_dir)
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.4))
    included_marks = shared_stats.get("included_marks", [])
    mark_summary = mark_stats_map(shared_stats)
    if not chrom_stats.empty:
        if "mark" in chrom_stats.columns and included_marks:
            x = np.arange(len(included_marks))
            p99_vals = [mark_summary[mark]["p99_raw"] for mark in included_marks]
            gate_vals = [mark_summary[mark]["methyl_gate_raw_threshold"] for mark in included_marks]
            pos_counts = [
                float(chrom_stats.loc[chrom_stats["mark"] == mark, "n_positive"].sum())
                for mark in included_marks
            ]
            axes[0].bar(x, p99_vals, color=[MARK_COLORS.get(mark, "slateblue") for mark in included_marks], alpha=0.85)
            axes[0].scatter(x, gate_vals, color="firebrick", label="feature/gate threshold", zorder=3)
            axes[0].set_xticks(x)
            axes[0].set_xticklabels(included_marks, rotation=30, ha="right")
            axes[0].set_ylabel("Genome-wide raw signal")
            axes[0].set_title("Per-mark p99 raw and gate threshold")
            axes[0].legend(frameon=False)

            axes[1].bar(x, pos_counts, color=[MARK_COLORS.get(mark, "darkorchid") for mark in included_marks], alpha=0.85)
            axes[1].set_xticks(x)
            axes[1].set_xticklabels(included_marks, rotation=30, ha="right")
            axes[1].set_ylabel("Positive 200 bp bins")
            axes[1].set_title("Combined Xa+Xi positive-bin counts by mark")
        else:
            labels = [f"{a}:{c}" for a, c in zip(chrom_stats["allele"], chrom_stats["chrom"])]
            x = np.arange(len(chrom_stats))
            axes[0].bar(x, chrom_stats["p99_raw"].values, color="slateblue")
            axes[0].axhline(
                shared_stats["methyl_gate_raw_threshold"],
                color="firebrick",
                ls="--",
                lw=1.2,
                label=f"feature/gate q{shared_stats['methyl_gate_quantile']:g}",
            )
            axes[0].set_xticks(x)
            axes[0].set_xticklabels(labels, rotation=90)
            axes[0].set_ylabel("Per-chromosome p99 raw signal")
            axes[0].set_title("Combined Xa+Xi raw-scale summaries")
            axes[0].legend(frameon=False)

            axes[1].bar(x, chrom_stats["n_positive"].values.astype(float), color="darkorchid")
            axes[1].set_xticks(x)
            axes[1].set_xticklabels(labels, rotation=90)
            axes[1].set_ylabel("Positive 200 bp bins")
            axes[1].set_title(
                f"p99_raw = {shared_stats['p99_raw']:.4f} | p99_log1p = {shared_stats['p99_log1p']:.4f}\n"
                f"shared feature/gate norm threshold = {shared_stats['methyl_gate_norm_threshold']:.4f}"
            )
    else:
        axes[0].text(0.5, 0.5, "No chromosome stats available", ha="center", va="center")
        axes[1].text(0.5, 0.5, "No thresholds available", ha="center", va="center")
    plt.tight_layout()
    plt.savefig(out_dir / "stage1_combined_whole_genome_normalization.png")
    plt.close()


def save_selected_feature_normalization_plots(
    out_dir: Path,
    raw_tracks_by_mark: Dict[str, np.ndarray],
    norm_tracks_by_mark: Dict[str, np.ndarray],
    gate_bins_by_mark: Dict[str, np.ndarray],
    combined_gate_bins: np.ndarray,
    prediction_region: Region,
    requested_region: Region,
    included_marks: Sequence[str],
    shared_stats: dict,
    allele_name: str,
):
    ensure_dir(out_dir)
    n_bins = len(combined_gate_bins)
    x_bp = np.arange(n_bins) * REQUESTED_RESOLUTION + prediction_region.start
    n_rows = len(included_marks) + 1
    fig, axes = plt.subplots(n_rows, 3, figsize=(15, 2.35 * n_rows), sharex=True, constrained_layout=True)
    axes = np.atleast_2d(axes)
    mark_summary = mark_stats_map(shared_stats)
    for row_idx, mark in enumerate(included_marks):
        raw_track = raw_tracks_by_mark[mark]
        norm_track = norm_tracks_by_mark[mark]
        gate_bins = gate_bins_by_mark[mark]
        stats = mark_summary[mark]
        color = MARK_COLORS.get(mark, "steelblue")
        axes[row_idx, 0].fill_between(x_bp, raw_track, step="mid", color=color, alpha=0.85)
        axes[row_idx, 0].axhline(stats["methyl_gate_raw_threshold"], color="firebrick", ls="--", lw=1.1)
        axes[row_idx, 0].set_ylabel(f"{mark}\nraw")
        axes[row_idx, 0].set_title("Raw track")

        axes[row_idx, 1].fill_between(x_bp, norm_track, step="mid", color=color, alpha=0.85)
        axes[row_idx, 1].axhline(stats["methyl_gate_norm_threshold"], color="firebrick", ls="--", lw=1.1)
        axes[row_idx, 1].set_ylabel("norm")
        axes[row_idx, 1].set_title("q-thresholded normalized")

        axes[row_idx, 2].fill_between(x_bp, gate_bins.astype(float), step="mid", color=color, alpha=0.85)
        axes[row_idx, 2].set_ylabel("gate")
        axes[row_idx, 2].set_yticks([0, 1])
        axes[row_idx, 2].set_title("Per-mark gate bins")

    axes[-1, 0].fill_between(x_bp, combined_gate_bins.astype(float), step="mid", color="black", alpha=0.8)
    axes[-1, 0].set_ylabel("any mark")
    axes[-1, 0].set_yticks([0, 1])
    axes[-1, 0].set_title("Combined bin gate")
    axes[-1, 1].axis("off")
    axes[-1, 2].axis("off")

    for ax in axes.flat:
        if not ax.has_data():
            continue
        ax.axvspan(requested_region.start, requested_region.end, color="gold", alpha=0.18)
        ax.axvline(requested_region.start, color="goldenrod", lw=1)
        ax.axvline(requested_region.end, color="goldenrod", lw=1)
    axes[0, 0].set_title(f"Stage 1 — {allele_name} multi-mark tracks", loc="left")
    for ax in axes[-1, :]:
        if ax.axison:
            ax.set_xlabel(f"{prediction_region.chrom} coordinate (bp)")
    plt.savefig(out_dir / "stage1_prediction_window_track.png")
    plt.close()


def save_feature_diagnostics_plot(out_path: Path, selected_slot: str, included_marks: Sequence[str], tracks: np.ndarray, tracks_smooth: np.ndarray, prod_map: np.ndarray, int_map: np.ndarray, dist_map: np.ndarray, prediction_region: Region, requested_region: Region, allele_name: str):
    n_bins = tracks.shape[1]
    x_bp = np.arange(n_bins) * REQUESTED_RESOLUTION + prediction_region.start
    fig = plt.figure(figsize=(13, 10))
    gs = gridspec.GridSpec(5, 1, height_ratios=[1, 1, 3, 3, 3])

    ax0 = fig.add_subplot(gs[0])
    for mark in included_marks:
        ax0.fill_between(x_bp, tracks[MARK_TO_INDEX[mark]], step="mid", color=MARK_COLORS.get(mark, "steelblue"), alpha=0.35, label=mark)
    ax0.set_ylabel("norm")
    ax0.set_title(f"Stage 2 — {allele_name} training-style input features (map panels shown for {selected_slot})")
    if included_marks:
        ax0.legend(frameon=False, ncol=min(3, len(included_marks)), loc="upper right")

    ax1 = fig.add_subplot(gs[1], sharex=ax0)
    for mark in included_marks:
        ax1.fill_between(x_bp, tracks_smooth[MARK_TO_INDEX[mark]], step="mid", color=MARK_COLORS.get(mark, "steelblue"), alpha=0.35)
    ax1.set_ylabel("smooth")

    for ax in [ax0, ax1]:
        ax.axvspan(requested_region.start, requested_region.end, color="gold", alpha=0.18)
        ax.tick_params(axis="x", labelbottom=False)

    panels = [
        (prod_map, f"{selected_slot}_prod", "viridis"),
        (int_map, f"{selected_slot}_int", "viridis"),
        (dist_map, "diag_index", "magma"),
    ]
    for k, (mat, title, cmap) in enumerate(panels, start=2):
        ax = fig.add_subplot(gs[k])
        im = ax.imshow(mat, cmap=cmap, aspect="auto")
        ax.set_title(title)
        fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01)

    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()


def save_prediction_panels(
    out_path: Path,
    prediction_region: Region,
    requested_region: Region,
    pred_map_uncalibrated_full: np.ndarray,
    pred_map_exactdiag_full: np.ndarray,
    pred_map_fully_calibrated_full: np.ndarray,
    pred_map_fully_calibrated_gated: np.ndarray,
    gate_map: np.ndarray,
    allele_name: str,
    pred_level_correction_applied: bool,
):
    fig, axes = plt.subplots(2, 4, figsize=(22.5, 9.0), constrained_layout=True)
    exact_diag_delta = pred_map_exactdiag_full - pred_map_uncalibrated_full
    pred_level_delta = pred_map_fully_calibrated_full - pred_map_exactdiag_full
    total_delta = pred_map_fully_calibrated_full - pred_map_uncalibrated_full
    pred_level_title = (
        "Prediction-level calibration delta"
        if pred_level_correction_applied
        else "Prediction-level calibration delta (not used)"
    )
    panels = [
        (pred_map_uncalibrated_full, "Uncalibrated predicted smoothed O/E", "RdBu_r", TARGET_MIN, TARGET_MAX),
        (pred_map_exactdiag_full, "Exact-diag calibrated predicted smoothed O/E", "RdBu_r", TARGET_MIN, TARGET_MAX),
        (exact_diag_delta, "Exact-diag calibration delta", "RdBu_r", -1.0, 1.0),
        (pred_map_fully_calibrated_full, "Fully calibrated predicted smoothed O/E", "RdBu_r", TARGET_MIN, TARGET_MAX),
        (pred_level_delta, pred_level_title, "RdBu_r", -1.0, 1.0),
        (total_delta, "Total calibration delta", "RdBu_r", -1.0, 1.0),
        (gate_map, "Self-similar pair gate", "viridis", 0.0, 1.0),
        (pred_map_fully_calibrated_gated, "Gated fully calibrated predicted smoothed O/E", "RdBu_r", TARGET_MIN, TARGET_MAX),
    ]

    for ax, (mat, title, cmap, vmin, vmax) in zip(axes.flat, panels):
        im = ax.imshow(mat, cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
        ax.set_title(f"{allele_name}: {title}")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        b0, b1 = requested_region_bin_slice(requested_region, prediction_region)
        rect = Rectangle((b0 - 0.5, b0 - 0.5), b1 - b0, b1 - b0, fill=False, edgecolor="gold", lw=1.4)
        ax.add_patch(rect)
    for ax in axes.flat[len(panels):]:
        ax.axis("off")

    plt.savefig(out_path)
    plt.close()


def save_truth_comparison_panels(out_path: Path, prediction_region: Region, requested_region: Region, pred_map_calibrated_gated: np.ndarray, true_raw_oe: np.ndarray, true_smooth_oe: np.ndarray, allele_name: str):
    fig, axes = plt.subplots(1, 4, figsize=(17.2, 4.2))
    panels = [
        (true_raw_oe, "Measured raw O/E", "RdBu_r", TARGET_MIN, TARGET_MAX),
        (true_smooth_oe, "Measured smoothed O/E", "RdBu_r", TARGET_MIN, TARGET_MAX),
        (pred_map_calibrated_gated, "Gated predicted smoothed O/E", "RdBu_r", TARGET_MIN, TARGET_MAX),
        (pred_map_calibrated_gated - true_smooth_oe, "Prediction - measured smoothed O/E", "RdBu_r", -2.0, 2.0),
    ]
    for ax, (mat, title, cmap, vmin, vmax) in zip(axes, panels):
        im = ax.imshow(mat, cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
        ax.set_title(f"{allele_name}: {title}")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        b0, b1 = requested_region_bin_slice(requested_region, prediction_region)
        rect = Rectangle((b0 - 0.5, b0 - 0.5), b1 - b0, b1 - b0, fill=False, edgecolor="gold", lw=1.4)
        ax.add_patch(rect)
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()


def save_crop_plot(out_path: Path, full_map: np.ndarray, crop_map: np.ndarray, prediction_region: Region, requested_region: Region, allele_name: str):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    im0 = axes[0].imshow(full_map, cmap="RdBu_r", vmin=TARGET_MIN, vmax=TARGET_MAX, aspect="auto")
    b0, b1 = requested_region_bin_slice(requested_region, prediction_region)
    rect = Rectangle((b0 - 0.5, b0 - 0.5), b1 - b0, b1 - b0, fill=False, edgecolor="gold", lw=1.5)
    axes[0].add_patch(rect)
    axes[0].set_title(f"{allele_name}: full 100 kb predicted smoothed O/E")
    fig.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04)

    im1 = axes[1].imshow(crop_map, cmap="RdBu_r", vmin=TARGET_MIN, vmax=TARGET_MAX, aspect="auto")
    axes[1].set_title(f"{allele_name}: cropped predicted smoothed O/E")
    fig.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04)
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()


def positive_log_norm(
    mat: np.ndarray,
    floor: float,
    lower_q: float = LOG_COLOR_QUANTILE_LOW,
    upper_q: float = LOG_COLOR_QUANTILE_HIGH,
) -> Tuple[Optional[LogNorm], np.ndarray]:
    vals = np.asarray(mat, dtype=np.float32)
    positive = vals[np.isfinite(vals) & (vals > 0)]
    if positive.size == 0:
        return None, vals.astype(np.float32, copy=False)
    vmin = max(float(np.percentile(positive, lower_q)), float(floor))
    vmax = max(float(np.percentile(positive, upper_q)), vmin * 10.0)
    masked = np.where(vals > 0, vals, np.nan).astype(np.float32, copy=False)
    return LogNorm(vmin=vmin, vmax=vmax), masked


def imshow_with_optional_log(
    fig,
    ax,
    mat: np.ndarray,
    title: str,
    cmap: str,
    *,
    use_log: bool = False,
    log_floor: float = LOG_FLOOR_OCCUPANCY,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
):
    title_suffix = ""
    if use_log:
        norm, mat_plot = positive_log_norm(mat, floor=log_floor)
        if norm is not None:
            im = ax.imshow(mat_plot, cmap=cmap, norm=norm, aspect="auto")
            title_suffix = " (log scale)"
        else:
            im = ax.imshow(mat.astype(np.float32), cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
            title_suffix = " (no positive values for log scale)"
    else:
        im = ax.imshow(mat.astype(np.float32), cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
    ax.set_title(f"{title}{title_suffix}")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    return im


def save_probability_count_plots(
    out_dir: Path,
    crop_pred: np.ndarray,
    t_map: np.ndarray,
    s_map: np.ndarray,
    q_t_map: np.ndarray,
    q_s_map: np.ndarray,
    p_extra_map: np.ndarray,
    subtraction_map: np.ndarray,
    feasible_count_map: np.ndarray,
    bucket_df: pd.DataFrame,
    diag_df: pd.DataFrame,
    p_extra_hist_df: pd.DataFrame,
    allele_name: str,
    polymer_alpha: float,
    kappa: float,
    local_diag_min: int,
    local_diag_max: int,
    feasible_mask_applied: bool,
):
    ensure_dir(out_dir)
    fig, axes = plt.subplots(2, 3, figsize=(16.5, 9.0), constrained_layout=True)
    panel_suffix = (
        f"alpha={polymer_alpha:.2f}, kappa={kappa:.2f}, local={local_diag_min}-{local_diag_max}, "
        f"feasible_mask={'yes' if feasible_mask_applied else 'no'}"
    )
    panels = [
        (crop_pred, f"{allele_name}: cropped predicted smoothed O/E\n{panel_suffix}", "RdBu_r", False, LOG_FLOOR_INTENSITY, TARGET_MIN, TARGET_MAX),
        (t_map, f"{allele_name}: target-like intensity T", "viridis", True, LOG_FLOOR_INTENSITY, None, None),
        (s_map, f"{allele_name}: scaled polymer baseline S", "viridis", True, LOG_FLOOR_INTENSITY, None, None),
        (q_t_map, f"{allele_name}: occupancy Q_T", "viridis", True, LOG_FLOOR_OCCUPANCY, None, None),
        (q_s_map, f"{allele_name}: occupancy Q_S", "viridis", True, LOG_FLOOR_OCCUPANCY, None, None),
        (p_extra_map, f"{allele_name}: extra-contact probability P_extra", "viridis", True, LOG_FLOOR_OCCUPANCY, None, None),
    ]
    for ax, (mat, title, cmap, use_log, log_floor, vmin, vmax) in zip(axes.flat, panels):
        imshow_with_optional_log(fig, ax, mat, title, cmap, use_log=use_log, log_floor=log_floor, vmin=vmin, vmax=vmax)
    plt.savefig(out_dir / "stage5_intensity_occupancy_maps.png", bbox_inches="tight")
    plt.close()

    fig, axes = plt.subplots(1, 3, figsize=(16.5, 4.8), constrained_layout=True)
    imshow_with_optional_log(
        fig,
        axes[0],
        feasible_count_map.astype(np.float32),
        f"{allele_name}: feasible fibers per pair",
        "magma",
        use_log=True,
        log_floor=1.0,
    )

    if not p_extra_hist_df.empty:
        centers = 0.5 * (p_extra_hist_df["bin_left"].values + p_extra_hist_df["bin_right"].values)
        widths = p_extra_hist_df["bin_right"].values - p_extra_hist_df["bin_left"].values
        axes[1].bar(centers, p_extra_hist_df["count"].values, width=widths, color="steelblue", align="center")
    axes[1].set_title(f"{allele_name}: P_extra histogram")
    axes[1].set_xlabel("P_extra")
    axes[1].set_ylabel("Pair count")

    mean_df = bucket_df[(bucket_df["bucket"] != "overall") & (bucket_df["metric"].isin(["T", "S", "Q_T", "Q_S", "P_extra"]))]
    for metric_name, sub in mean_df.groupby("metric"):
        axes[2].plot(sub["bucket"], sub["mean"], marker="o", label=metric_name)
    axes[2].set_title(f"{allele_name}: mean target fields by distance bucket")
    axes[2].set_xlabel("Distance bucket")
    axes[2].set_ylabel("Mean value")
    axes[2].tick_params(axis="x", rotation=45)
    axes[2].legend(frameon=False)
    plt.savefig(out_dir / "stage6_bucket_and_feasibility_summaries.png", bbox_inches="tight")
    plt.close()

    fig, axes = plt.subplots(1, 3, figsize=(17.5, 4.8), constrained_layout=True)
    imshow_with_optional_log(
        fig,
        axes[0],
        subtraction_map.astype(np.float32),
        f"{allele_name}: Q_T - kappa * Q_S",
        "RdBu_r",
        use_log=False,
        vmin=-1.0,
        vmax=1.0,
    )

    subtraction_bucket_df = bucket_df[
        (bucket_df["bucket"] != "overall") &
        (bucket_df["metric"] == "Q_T_minus_kappa_Q_S")
    ]
    axes[1].plot(subtraction_bucket_df["bucket"], subtraction_bucket_df["mean"], marker="o", color="firebrick")
    axes[1].axhline(0.0, color="black", lw=1.0, ls="--")
    axes[1].set_title(f"{allele_name}: mean subtraction by bucket")
    axes[1].set_xlabel("Distance bucket")
    axes[1].set_ylabel("Mean Q_T - kappa * Q_S")
    axes[1].tick_params(axis="x", rotation=45)

    subtraction_diag_df = diag_df[diag_df["metric"] == "Q_T_minus_kappa_Q_S"]
    axes[2].plot(subtraction_diag_df["diag_idx"], subtraction_diag_df["mean"], color="darkslateblue", lw=1.4)
    axes[2].axhline(0.0, color="black", lw=1.0, ls="--")
    axes[2].set_title(f"{allele_name}: exact-diagonal subtraction profile")
    axes[2].set_xlabel("Diagonal index")
    axes[2].set_ylabel("Mean Q_T - kappa * Q_S")
    plt.savefig(out_dir / "stage6_subtraction_diagnostics.png", bbox_inches="tight")
    plt.close()


def save_nucleosome_mapping_plots(out_dir: Path, nuc_df: pd.DataFrame, requested_region: Region, allele_name: str):
    ensure_dir(out_dir)
    if nuc_df.empty:
        return
    occ = nuc_df.groupby("bin_id").size().reset_index(name="n_nucleosomes")
    fig, axes = plt.subplots(2, 1, figsize=(11, 6), sharex=False)
    axes[0].bar(occ["bin_id"], occ["n_nucleosomes"], color="slateblue")
    axes[0].set_xlabel("Requested-region 200 bp bin")
    axes[0].set_ylabel("Nucleosomes across fibers")
    axes[0].set_title(f"Stage 7 — {allele_name} nucleosome occupancy by bin")

    sample = nuc_df[nuc_df["fiber_id"].isin(sorted(nuc_df["fiber_id"].unique(), key=natural_sort_key)[:10])].copy()
    fiber_to_y = {fid: i for i, fid in enumerate(sorted(sample["fiber_id"].unique(), key=natural_sort_key))}
    y = sample["fiber_id"].map(fiber_to_y).values
    axes[1].scatter(sample["center_bp"], y, s=10, alpha=0.8)
    axes[1].set_xlabel(f"{requested_region.chrom} coordinate (bp)")
    axes[1].set_ylabel("Fiber index (subset)")
    axes[1].set_title(f"{allele_name} example nucleosome centers in requested region")
    plt.tight_layout()
    plt.savefig(out_dir / "stage7_nucleosome_mapping.png")
    plt.close()


def save_assignment_plots(
    out_dir: Path,
    p_extra_map: np.ndarray,
    sampled_occ_map: np.ndarray,
    recon_occ_map: np.ndarray,
    valency_df: pd.DataFrame,
    assignments_df: pd.DataFrame,
    efficiency_df: pd.DataFrame,
    fiber_summary_df: pd.DataFrame,
    allele_name: str,
):
    ensure_dir(out_dir)
    fig, axes = plt.subplots(1, 4, figsize=(18.5, 4.4), constrained_layout=True)
    mats = [
        (p_extra_map.astype(np.float32), f"{allele_name}: target P_extra", "viridis", True, LOG_FLOOR_OCCUPANCY, None, None),
        (sampled_occ_map.astype(np.float32), f"{allele_name}: sampled occupancy", "viridis", True, LOG_FLOOR_OCCUPANCY, None, None),
        (recon_occ_map.astype(np.float32), f"{allele_name}: reconstructed occupancy", "viridis", True, LOG_FLOOR_OCCUPANCY, None, None),
        ((recon_occ_map - p_extra_map).astype(np.float32), f"{allele_name}: reconstructed - target occupancy", "RdBu_r", False, LOG_FLOOR_OCCUPANCY, -1.0, 1.0),
    ]
    for ax, (mat, title, cmap, use_log, log_floor, vmin, vmax) in zip(axes, mats):
        imshow_with_optional_log(fig, ax, mat, title, cmap, use_log=use_log, log_floor=log_floor, vmin=vmin, vmax=vmax)
    plt.savefig(out_dir / "stage8_reconstruction.png", bbox_inches="tight")
    plt.close()

    fig, axes = plt.subplots(1, 3, figsize=(17.5, 4.6), constrained_layout=True)
    if not valency_df.empty:
        axes[0].hist(valency_df["valency"], bins=np.arange(valency_df["valency"].max() + 2) - 0.5, color="indianred", rwidth=0.85)
    axes[0].set_xlabel("Assigned valency per nucleosome")
    axes[0].set_ylabel("Count")
    axes[0].set_title("Valency distribution")

    if not fiber_summary_df.empty:
        axes[1].hist(fiber_summary_df["assigned_contacts"], bins=20, color="steelblue", alpha=0.85)
    elif not assignments_df.empty:
        per_fiber = assignments_df.groupby("fiber_id").size().reset_index(name="n_contacts")
        axes[1].hist(per_fiber["n_contacts"], bins=20, color="steelblue")
    axes[1].set_xlabel("Assigned methylation contacts per fiber")
    axes[1].set_ylabel("Count")
    axes[1].set_title("Per-fiber assigned contacts")

    if not efficiency_df.empty:
        eff = efficiency_df[efficiency_df["bucket"] != "overall"].copy()
        axes[2].plot(eff["bucket"], eff["mean_target_occupancy"], marker="o", label="Target")
        axes[2].plot(eff["bucket"], eff["mean_sampled_occupancy"], marker="o", label="Sampled")
        axes[2].plot(eff["bucket"], eff["mean_reconstructed_occupancy"], marker="o", label="Reconstructed")
        axes[2].tick_params(axis="x", rotation=45)
        axes[2].legend(frameon=False)
    axes[2].set_xlabel("Distance bucket")
    axes[2].set_ylabel("Mean occupancy")
    axes[2].set_title("Assignment efficiency by bucket")
    plt.savefig(out_dir / "stage8_assignment_summaries.png", bbox_inches="tight")
    plt.close()


# ============================================================
# Core per-allele pipeline
# ============================================================
def run_single_allele_prediction(
    allele_name: str,
    bw_paths_by_mark: Dict[str, Sequence[str]],
    included_marks: Sequence[str],
    primary_feature_slot: str,
    args,
    model,
    model_meta: dict,
    feature_variant: str,
    requested_region: Region,
    prediction_region: Region,
    shared_stats: dict,
    root_out_dir: Path,
    clr=None,
    stratum_expected: Optional[np.ndarray] = None,
):
    allele_out = ensure_dir(root_out_dir / allele_name)
    plot_dir = ensure_dir(allele_out / "plots")
    table_dir = ensure_dir(allele_out / "tables")
    contact_dir = ensure_dir(allele_out / "contacts")
    logger = setup_logger(allele_out / "run.log")
    allele_seed = int(args.seed + (1 if allele_name == "Xi" else 0))
    rng_assign = np.random.default_rng(allele_seed + 10_000)
    resolved_dim_glob = args.xa_dim_glob if allele_name == "Xa" and args.xa_dim_glob else (
        args.xi_dim_glob if allele_name == "Xi" and args.xi_dim_glob else args.dim_glob
    )

    logger.info("Running allele-specific prediction for %s", allele_name)

    # --------------------------------------------------------
    # STAGE 1 — Allele-specific track on centered 100 kb window
    # --------------------------------------------------------
    n_bins_pred = WINDOW_SIZE_BP // REQUESTED_RESOLUTION
    stage1_meta = table_dir / "stage1_metadata.json"
    stage1_inputs = {
        "script_version": SCRIPT_VERSION,
        "allele": allele_name,
        "bw_paths_by_mark": {mark: [str(x) for x in bw_paths_by_mark[mark]] for mark in included_marks},
        "included_marks": list(included_marks),
        "feature_marks": list(included_marks),
        "gate_marks": list(included_marks),
        "pair_gate_mode": "any_same_mark_match",
        "feature_slot": primary_feature_slot,
        "prediction_region": region_dict(prediction_region),
        "requested_region": region_dict(requested_region),
        "normalization_mode": "dataset_specific_shared_feature_gate_threshold",
        "shared_feature_and_gate_threshold": True,
        "mark_stats": {mark: shared_stats["mark_stats"][mark] for mark in included_marks},
        "resolution": REQUESTED_RESOLUTION,
    }
    stage1_files = [
        table_dir / "stage1_mark_tracks_prediction_window.npz",
        plot_dir / "stage1_prediction_window_track.png",
    ]
    if not args.force and metadata_matches(stage1_meta, stage1_inputs) and files_exist(stage1_files):
        logger.info("STAGE 1 — Reusing checkpointed allele-specific prediction-window track")
        with np.load(stage1_files[0]) as data:
            raw_tracks = data["raw_tracks"].astype(np.float32)
            raw_tracks_denoised = data["raw_tracks_denoised"].astype(np.float32)
            norm_tracks = data["norm_tracks"].astype(np.float32)
            gate_bins_matrix = data["gate_bins_matrix"].astype(bool)
            combined_gate_bins = data["combined_gate_bins"].astype(bool)
    else:
        logger.info("STAGE 1 — Loading allele-specific track on centered 100 kb window")
        raw_tracks = np.zeros((N_MARKS, n_bins_pred), dtype=np.float32)
        raw_tracks_denoised = np.zeros((N_MARKS, n_bins_pred), dtype=np.float32)
        norm_tracks = np.zeros((N_MARKS, n_bins_pred), dtype=np.float32)
        gate_bins_matrix = np.zeros((N_MARKS, n_bins_pred), dtype=bool)
        for mark in included_marks:
            raw_track_pred = load_bigwig_window_average(
                bw_paths_by_mark[mark],
                prediction_region.chrom,
                prediction_region.start,
                prediction_region.end,
                n_bins_pred,
            )
            mark_stats = shared_stats["mark_stats"][mark]
            raw_track_denoised_pred, norm_track_pred = normalize_and_denoise_track(
                raw_track_pred,
                mark_stats["p99_log1p"],
                mark_stats["epigenetic_noise_floor_raw_threshold"],
            )
            idx = MARK_TO_INDEX[mark]
            raw_tracks[idx] = raw_track_pred
            raw_tracks_denoised[idx] = raw_track_denoised_pred
            norm_tracks[idx] = norm_track_pred
            gate_bins_matrix[idx] = raw_track_pred >= mark_stats["methyl_gate_raw_threshold"]
        combined_gate_bins = np.any(gate_bins_matrix[[MARK_TO_INDEX[mark] for mark in included_marks]], axis=0)

        np.savez_compressed(
            stage1_files[0],
            raw_tracks=raw_tracks,
            raw_tracks_denoised=raw_tracks_denoised,
            norm_tracks=norm_tracks,
            gate_bins_matrix=gate_bins_matrix.astype(np.int8),
            combined_gate_bins=combined_gate_bins.astype(np.int8),
        )
        raw_tracks_by_mark = {mark: raw_tracks[MARK_TO_INDEX[mark]] for mark in included_marks}
        norm_tracks_by_mark = {mark: norm_tracks[MARK_TO_INDEX[mark]] for mark in included_marks}
        gate_bins_by_mark = {mark: gate_bins_matrix[MARK_TO_INDEX[mark]] for mark in included_marks}
        save_selected_feature_normalization_plots(
            plot_dir,
            raw_tracks_by_mark,
            norm_tracks_by_mark,
            gate_bins_by_mark,
            combined_gate_bins,
            prediction_region,
            requested_region,
            included_marks,
            shared_stats,
            allele_name,
        )
        write_stage_metadata(stage1_meta, stage1_inputs)

    raw_tracks_by_mark = {mark: raw_tracks[MARK_TO_INDEX[mark]] for mark in included_marks}
    raw_tracks_denoised_by_mark = {mark: raw_tracks_denoised[MARK_TO_INDEX[mark]] for mark in included_marks}
    norm_tracks_by_mark = {mark: norm_tracks[MARK_TO_INDEX[mark]] for mark in included_marks}
    gate_bins_by_mark = {mark: gate_bins_matrix[MARK_TO_INDEX[mark]] for mark in included_marks}

    # --------------------------------------------------------
    # STAGE 2 — Build training-compatible pairwise features
    # --------------------------------------------------------
    stage2_meta = table_dir / "stage2_metadata.json"
    stage2_npz = table_dir / "stage2_feature_matrix.npz"
    stage2_summary_csv = table_dir / "stage2_feature_nonzero_summary.csv"
    stage2_inputs = {
        "script_version": SCRIPT_VERSION,
        "allele": allele_name,
        "feature_slot": primary_feature_slot,
        "included_marks": list(included_marks),
        "prediction_region": region_dict(prediction_region),
        "n_bins": int(n_bins_pred),
        "resolution": REQUESTED_RESOLUTION,
        "feature_variant": feature_variant,
    }
    if not args.force and metadata_matches(stage2_meta, stage2_inputs) and files_exist([stage2_npz, stage2_summary_csv, plot_dir / "stage2_feature_diagnostics.png"]):
        logger.info("STAGE 2 — Reusing checkpointed training-compatible pairwise feature matrix")
        with np.load(stage2_npz) as data:
            full_feats = data["full_feats"]
            diag_idx = data["diag_idx"]
            valid_mask = data["valid_mask"].astype(bool)
            ii = data["ii"]
            jj = data["jj"]
    else:
        logger.info("STAGE 2 — Building training-compatible pairwise feature matrix")
        tracks, tracks_smooth = build_feature_tracks(norm_tracks_by_mark, n_bins_pred)
        full_feats, diag_idx, valid_mask, ii, jj = build_pair_features_for_window(tracks, tracks_smooth, REQUESTED_RESOLUTION)
        np.savez_compressed(stage2_npz, full_feats=full_feats, diag_idx=diag_idx, valid_mask=valid_mask, ii=ii, jj=jj)
        write_stage_metadata(stage2_meta, stage2_inputs, {"n_pairs": int(len(diag_idx)), "n_features": int(full_feats.shape[1])})

        selected_idx = MARK_TO_INDEX[primary_feature_slot]
        prod_map = np.full((n_bins_pred, n_bins_pred), np.nan, dtype=np.float32)
        int_map = np.full((n_bins_pred, n_bins_pred), np.nan, dtype=np.float32)
        diag_map = np.full((n_bins_pred, n_bins_pred), np.nan, dtype=np.float32)
        prod_vec = full_feats[:, 2 * N_MARKS + selected_idx]
        int_vec = full_feats[:, 3 * N_MARKS + selected_idx]
        dist_vec = full_feats[:, 5 * N_MARKS + 1]
        prod_map[ii, jj] = prod_vec
        prod_map[jj, ii] = prod_vec
        int_map[ii, jj] = int_vec
        int_map[jj, ii] = int_vec
        diag_map[ii, jj] = dist_vec
        diag_map[jj, ii] = dist_vec
        np.fill_diagonal(prod_map, np.nan)
        np.fill_diagonal(int_map, np.nan)
        np.fill_diagonal(diag_map, np.nan)
        save_feature_diagnostics_plot(plot_dir / "stage2_feature_diagnostics.png", primary_feature_slot, included_marks, tracks, tracks_smooth, prod_map, int_map, diag_map, prediction_region, requested_region, allele_name)

    model_input = select_feature_matrix(full_feats, feature_variant)
    save_feature_nonzero_summary(model_input, feature_names_for_variant(feature_variant), stage2_summary_csv)
    if getattr(model, "n_features_in_", model_input.shape[1]) != model_input.shape[1]:
        raise ValueError(
            f"Model expects {getattr(model, 'n_features_in_', 'unknown')} input features, "
            f"but the constructed matrix has {model_input.shape[1]} features for inferred variant {feature_variant}. "
            f"Base constructed feature count={full_feats.shape[1]}; model_meta={args.model_meta}; model={args.model_joblib}."
        )

    # --------------------------------------------------------
    # STAGE 3 — Predict 100 kb O/E map, then methylation gate
    # --------------------------------------------------------
    stage3_meta = table_dir / "stage3_metadata.json"
    resolved_correction = resolve_diag_correction_path(args.diag_correction)
    resolved_pred_level_correction = resolve_pred_level_correction_path(
        args.pred_level_correction,
        args.diag_correction,
        disabled=bool(args.disable_pred_level_correction),
    )
    target_mode = str(model_meta.get("target_mode", "unknown"))
    requires_residual_baseline = target_mode == "smoothed_oe_residual"
    residual_baseline_used = args.residual_baseline if requires_residual_baseline else None
    stage3_inputs = {
        "script_version": SCRIPT_VERSION,
        "stage3_inference_version": STAGE3_INFERENCE_VERSION,
        "allele": allele_name,
        "feature_variant": feature_variant,
        "model_joblib": str(args.model_joblib),
        "model_meta": str(args.model_meta) if args.model_meta else None,
        "model_target_mode": target_mode,
        "rf_output_space": "smoothed_oe_residual" if requires_residual_baseline else "smoothed_oe",
        "prediction_space": "smoothed_oe",
        "residual_baseline_requested": str(args.residual_baseline) if args.residual_baseline else None,
        "residual_baseline_used": str(residual_baseline_used) if residual_baseline_used else None,
        "diag_correction_requested": str(args.diag_correction) if args.diag_correction else None,
        "diag_correction_used": str(resolved_correction) if resolved_correction else None,
        "pred_level_correction_enabled": not bool(args.disable_pred_level_correction),
        "pred_level_correction_requested": str(args.pred_level_correction) if args.pred_level_correction else None,
        "pred_level_correction_used": str(resolved_pred_level_correction) if resolved_pred_level_correction else None,
        "included_marks": list(included_marks),
        "feature_marks": list(included_marks),
        "gate_marks": list(included_marks),
        "pair_gate_mode": "any_same_mark_match",
        "feature_threshold_source": "allele_dataset_genomewide",
        "methyl_gate_source": "allele_dataset_genomewide",
        "shared_feature_and_gate_threshold": True,
        "mark_stats": {mark: shared_stats["mark_stats"][mark] for mark in included_marks},
        "prediction_region": region_dict(prediction_region),
        "has_truth": bool(clr is not None and stratum_expected is not None),
    }
    stage3_files = [
        table_dir / "stage3_pred_map_raw.npy",
        table_dir / "stage3_pred_map_calibrated.npy",
        table_dir / "stage3_pred_smoothed_oe_uncalibrated_full.npy",
        table_dir / "stage3_pred_smoothed_oe_exactdiag_calibrated_full.npy",
        table_dir / "stage3_pred_smoothed_oe_calibrated_full.npy",
        table_dir / "stage3_pred_smoothed_oe_exactdiag_calibrated_gated.npy",
        table_dir / "stage3_pred_smoothed_oe_calibrated_gated.npy",
        table_dir / "stage3_calibration_delta.npy",
        table_dir / "stage3_exact_diag_calibration_delta.npy",
        table_dir / "stage3_pred_level_calibration_delta.npy",
        table_dir / "stage3_rf_residual.npy",
        table_dir / "stage3_residual_baseline_applied.npy",
        table_dir / "stage3_gate_map.npy",
        table_dir / "stage3_gate_bins.npy",
        table_dir / "stage3_gate_bins_by_mark.npz",
        table_dir / "stage3_gate_bins.csv",
        table_dir / "stage3_gated_pair_distance_summary.csv",
        plot_dir / "stage3_prediction_panels.png",
        plot_dir / "stage3_reconstruction_histograms.png",
        plot_dir / "stage3_gated_pair_distance_summary.png",
    ]
    if stage3_inputs["has_truth"]:
        stage3_files.extend([
            table_dir / "stage3_true_raw_oe.npy",
            table_dir / "stage3_true_smooth_oe.npy",
            plot_dir / "stage3_truth_comparison.png",
        ])

    if not args.force and metadata_matches(stage3_meta, stage3_inputs) and files_exist(stage3_files):
        logger.info("STAGE 3 — Reusing checkpointed 100 kb contact prediction maps")
        pred_map_raw = np.load(table_dir / "stage3_pred_map_raw.npy")
        pred_map_cal = np.load(table_dir / "stage3_pred_map_calibrated.npy")
        pred_map_raw_full = np.load(table_dir / "stage3_pred_smoothed_oe_uncalibrated_full.npy")
        pred_map_exactdiag_full = np.load(table_dir / "stage3_pred_smoothed_oe_exactdiag_calibrated_full.npy")
        pred_map_cal_full = np.load(table_dir / "stage3_pred_smoothed_oe_calibrated_full.npy")
        gate_map = np.load(table_dir / "stage3_gate_map.npy")
        gate_bins = np.load(table_dir / "stage3_gate_bins.npy").astype(bool)
        true_raw_oe = np.load(table_dir / "stage3_true_raw_oe.npy") if (table_dir / "stage3_true_raw_oe.npy").exists() else None
        true_smooth_oe = np.load(table_dir / "stage3_true_smooth_oe.npy") if (table_dir / "stage3_true_smooth_oe.npy").exists() else None
    else:
        logger.info("STAGE 3 — Predicting the 100 kb contact map")
        rf_residual = predict_in_batches(model, model_input)
        residual_baseline = load_residual_baseline(
            args.residual_baseline if requires_residual_baseline else None,
            diag_idx,
            requires_residual_baseline,
            logger,
        )
        baseline_applied = baseline_values_from_diag(diag_idx, residual_baseline)
        pred_vec_raw = add_residual_baseline(rf_residual, diag_idx, residual_baseline)
        diag_correction = load_diag_correction(args.diag_correction, logger)
        pred_vec_exactdiag = apply_diag_correction(pred_vec_raw, diag_idx, diag_correction)
        pred_level_correction = load_pred_level_correction(
            args.pred_level_correction,
            args.diag_correction,
            disabled=bool(args.disable_pred_level_correction),
            logger=logger,
        )
        pred_vec_cal = apply_pred_level_correction(pred_vec_exactdiag, pred_level_correction)
        pred_map_raw_full = vector_to_symmetric_matrix(pred_vec_raw.astype(np.float32), ii, jj, n_bins_pred)
        pred_map_exactdiag_full = vector_to_symmetric_matrix(pred_vec_exactdiag.astype(np.float32), ii, jj, n_bins_pred)
        pred_map_cal_full = vector_to_symmetric_matrix(pred_vec_cal.astype(np.float32), ii, jj, n_bins_pred)
        exact_diag_calibration_delta = pred_map_exactdiag_full - pred_map_raw_full
        pred_level_calibration_delta = pred_map_cal_full - pred_map_exactdiag_full
        calibration_delta = pred_map_cal_full - pred_map_raw_full

        gate_bins, gate_bins_by_mark, pair_gate_by_mark, pair_gate, gate_map = build_pair_gate_from_tracks(
            raw_tracks_by_mark,
            ii,
            jj,
            shared_stats,
            included_marks,
        )
        pred_vec_raw_gated = np.where(pair_gate, pred_vec_raw, np.nan).astype(np.float32)
        pred_vec_exactdiag_gated = np.where(pair_gate, pred_vec_exactdiag, np.nan).astype(np.float32)
        pred_vec_cal_gated = np.where(pair_gate, pred_vec_cal, np.nan).astype(np.float32)
        pred_map_raw = vector_to_symmetric_matrix(pred_vec_raw_gated, ii, jj, n_bins_pred)
        pred_map_exactdiag = vector_to_symmetric_matrix(pred_vec_exactdiag_gated, ii, jj, n_bins_pred)
        pred_map_cal = vector_to_symmetric_matrix(pred_vec_cal_gated, ii, jj, n_bins_pred)

        true_raw_oe = None
        true_smooth_oe = None
        if clr is not None and stratum_expected is not None:
            raw_mat = load_microc_window(clr, prediction_region.chrom, prediction_region.start, prediction_region.end)
            true_raw_oe = normalize_microc(raw_mat, stratum_expected)
            true_smooth_oe = smooth_oe_matrix(true_raw_oe)
            np.save(table_dir / "stage3_true_raw_oe.npy", true_raw_oe)
            np.save(table_dir / "stage3_true_smooth_oe.npy", true_smooth_oe)

        np.save(table_dir / "stage3_pred_map_raw.npy", pred_map_raw)
        np.save(table_dir / "stage3_pred_map_calibrated.npy", pred_map_cal)
        np.save(table_dir / "stage3_pred_smoothed_oe_uncalibrated_full.npy", pred_map_raw_full)
        np.save(table_dir / "stage3_pred_smoothed_oe_exactdiag_calibrated_full.npy", pred_map_exactdiag_full)
        np.save(table_dir / "stage3_pred_smoothed_oe_calibrated_full.npy", pred_map_cal_full)
        np.save(table_dir / "stage3_pred_smoothed_oe_exactdiag_calibrated_gated.npy", pred_map_exactdiag)
        np.save(table_dir / "stage3_pred_smoothed_oe_calibrated_gated.npy", pred_map_cal)
        np.save(table_dir / "stage3_calibration_delta.npy", calibration_delta)
        np.save(table_dir / "stage3_exact_diag_calibration_delta.npy", exact_diag_calibration_delta)
        np.save(table_dir / "stage3_pred_level_calibration_delta.npy", pred_level_calibration_delta)
        np.save(table_dir / "stage3_rf_residual.npy", rf_residual.astype(np.float32))
        np.save(table_dir / "stage3_residual_baseline_applied.npy", baseline_applied.astype(np.float32))
        np.save(table_dir / "stage3_gate_map.npy", gate_map)
        np.save(table_dir / "stage3_gate_bins.npy", gate_bins.astype(np.int8))
        np.savez_compressed(
            table_dir / "stage3_gate_bins_by_mark.npz",
            gate_bins_matrix=np.stack([gate_bins_by_mark[mark].astype(np.int8) for mark in included_marks], axis=0),
        )
        save_gate_diagnostics(
            table_dir,
            plot_dir,
            raw_tracks_by_mark,
            norm_tracks_by_mark,
            gate_bins_by_mark,
            pair_gate_by_mark,
            gate_bins,
            pair_gate,
            ii,
            jj,
            prediction_region,
            allele_name,
            included_marks,
        )
        save_stage3_reconstruction_histogram(
            plot_dir / "stage3_reconstruction_histograms.png",
            rf_residual,
            baseline_applied,
            pred_vec_raw,
            pred_vec_cal,
            allele_name,
        )
        save_prediction_panels(
            plot_dir / "stage3_prediction_panels.png",
            prediction_region,
            requested_region,
            pred_map_raw_full,
            pred_map_exactdiag_full,
            pred_map_cal_full,
            pred_map_cal,
            gate_map,
            allele_name,
            pred_level_correction is not None,
        )
        if true_raw_oe is not None and true_smooth_oe is not None:
            save_truth_comparison_panels(plot_dir / "stage3_truth_comparison.png", prediction_region, requested_region, pred_map_cal, true_raw_oe, true_smooth_oe, allele_name)
        write_stage_metadata(stage3_meta, stage3_inputs, {
            "correction_loaded": str(resolved_correction) if resolved_correction else None,
            "pred_level_correction_enabled": not bool(args.disable_pred_level_correction),
            "pred_level_correction_loaded": str(resolved_pred_level_correction) if resolved_pred_level_correction else None,
            "residual_baseline_loaded": str(args.residual_baseline) if residual_baseline is not None else None,
        })

    # --------------------------------------------------------
    # STAGE 4 — Crop to the requested region for downstream use
    # --------------------------------------------------------
    b0, b1 = requested_region_bin_slice(requested_region, prediction_region)
    stage4_meta = table_dir / "stage4_metadata.json"
    stage4_file = table_dir / "stage4_cropped_pred_map.npy"
    stage4_inputs = {
        "script_version": SCRIPT_VERSION,
        "allele": allele_name,
        "requested_region": region_dict(requested_region),
        "prediction_region": region_dict(prediction_region),
        "bin_slice": {"b0": int(b0), "b1": int(b1)},
        "stage3_inputs": stage3_inputs,
    }
    if not args.force and metadata_matches(stage4_meta, stage4_inputs) and files_exist([stage4_file, plot_dir / "stage4_crop.png"]):
        logger.info("STAGE 4 — Reusing checkpointed cropped prediction map")
        crop_pred = np.load(stage4_file)
    else:
        logger.info("STAGE 4 — Cropping prediction to the requested region")
        crop_pred = pred_map_cal[b0:b1, b0:b1].copy()
        save_crop_plot(plot_dir / "stage4_crop.png", pred_map_cal, crop_pred, prediction_region, requested_region, allele_name)
        np.save(stage4_file, crop_pred)
        write_stage_metadata(stage4_meta, stage4_inputs)

    # --------------------------------------------------------
    # STAGE 4b — Optional dim-file geometry for feasible pair masking
    # --------------------------------------------------------
    dim_paths: List[Path] = []
    dim_signatures: List[dict] = []
    nuc_df = pd.DataFrame()
    fiber_bin_to_nucs = None
    centers_by_fiber = None
    feasible_mask = None
    feasible_count_map = np.zeros_like(crop_pred, dtype=np.int32)

    if resolved_dim_glob:
        dim_paths = sorted([Path(p) for p in glob.glob(resolved_dim_glob)], key=natural_sort_key)
        if not dim_paths:
            logger.warning("No dim files matched resolved dim glob=%s; feasible masking and ensemble assignment will be skipped", resolved_dim_glob)
        elif args.fiber_region_start is None:
            raise ValueError("--fiber-region-start is required when dim files are provided")
        else:
            dim_signatures = [file_signature(p) for p in dim_paths]
            stage7_meta = table_dir / "stage7_metadata.json"
            stage7_table = table_dir / "stage7_nucleosome_table.csv"
            stage7_inputs = {
                "script_version": SCRIPT_VERSION,
                "allele": allele_name,
                "dim_glob": resolved_dim_glob,
                "dim_files": dim_signatures,
                "fiber_region_start": int(args.fiber_region_start),
                "requested_region": region_dict(requested_region),
                "resolution": REQUESTED_RESOLUTION,
            }
            logger.info("STAGE 4b — Loading dim-file geometry for feasible pair masking")
            nuc_df, fiber_bin_to_nucs, centers_by_fiber = fibers_from_dim_files(dim_paths, requested_region, args.fiber_region_start, REQUESTED_RESOLUTION)
            feasible_mask, feasible_count_map = build_feasible_pair_statistics(fiber_bin_to_nucs, crop_pred.shape[0])
            np.save(table_dir / "stage4_feasible_pair_mask.npy", feasible_mask.astype(np.int8))
            if args.force or not (metadata_matches(stage7_meta, stage7_inputs) and stage7_table.exists()):
                nuc_df.to_csv(stage7_table, index=False)
                save_nucleosome_mapping_plots(plot_dir, nuc_df, requested_region, allele_name)
                write_stage_metadata(stage7_meta, stage7_inputs, {
                    "n_nucleosomes": int(len(nuc_df)),
                    "n_feasible_pairs": int(np.count_nonzero(np.triu(feasible_mask, k=1))),
                })
            else:
                logger.info("STAGE 7 — Reusing checkpointed nucleosome mapping table")

    crop_pred_for_sampling = apply_feasible_mask(crop_pred, feasible_mask)
    if feasible_mask is not None:
        np.save(table_dir / "stage5_cropped_pred_feasible_masked.npy", crop_pred_for_sampling)

    # --------------------------------------------------------
    # STAGE 5/6 — Absolute extra-contact probability maps
    # --------------------------------------------------------
    stage56_meta = table_dir / "stage5_6_metadata.json"
    stage56_files = [
        table_dir / "stage5_target_intensity_map.npy",
        table_dir / "stage5_polymer_baseline_raw_map.npy",
        table_dir / "stage5_polymer_baseline_scaled_map.npy",
        table_dir / "stage5_target_occupancy_map.npy",
        table_dir / "stage5_polymer_occupancy_map.npy",
        table_dir / "stage5_extra_probability_map.npy",
        table_dir / "stage5_subtraction_map.npy",
        table_dir / "stage5_allowed_pair_mask.npy",
        table_dir / "stage5_bucket_summary.csv",
        table_dir / "stage5_exact_diag_summary.csv",
        table_dir / "stage5_p_extra_histogram.csv",
        table_dir / "stage6_feasible_fiber_count_map.npy",
        plot_dir / "stage5_intensity_occupancy_maps.png",
        plot_dir / "stage6_bucket_and_feasibility_summaries.png",
        plot_dir / "stage6_subtraction_diagnostics.png",
    ]
    if feasible_mask is not None:
        stage56_files.extend([
            table_dir / "stage4_feasible_pair_mask.npy",
            table_dir / "stage5_cropped_pred_feasible_masked.npy",
        ])
    stage56_inputs = {
        "script_version": SCRIPT_VERSION,
        "allele": allele_name,
        "stage4_inputs": stage4_inputs,
        "feasible_mask_from_dim_files": bool(feasible_mask is not None),
        "dim_files": dim_signatures,
        "min_diag_bins": int(args.min_diag_bins),
        "polymer_alpha": float(args.polymer_alpha),
        "polymer_local_diag_min": int(args.polymer_local_diag_min),
        "polymer_local_diag_max": int(args.polymer_local_diag_max),
        "polymer_local_statistic": str(args.polymer_local_statistic),
        "p_local": float(args.p_local),
        "kappa": float(args.kappa),
        "seed": int(allele_seed),
        "requires_stratum_expected": True,
        "stratum_expected": str(args.stratum_expected) if args.stratum_expected else None,
    }
    if not args.force and metadata_matches(stage56_meta, stage56_inputs) and files_exist(stage56_files):
        logger.info("STAGE 5/6 — Reusing checkpointed absolute intensity / probability maps")
        t_map = np.load(table_dir / "stage5_target_intensity_map.npy")
        s_raw_map = np.load(table_dir / "stage5_polymer_baseline_raw_map.npy")
        s_map = np.load(table_dir / "stage5_polymer_baseline_scaled_map.npy")
        q_t_map = np.load(table_dir / "stage5_target_occupancy_map.npy")
        q_s_map = np.load(table_dir / "stage5_polymer_occupancy_map.npy")
        p_extra_map = np.load(table_dir / "stage5_extra_probability_map.npy")
        subtraction_map = np.load(table_dir / "stage5_subtraction_map.npy")
        allowed_pair_mask = np.load(table_dir / "stage5_allowed_pair_mask.npy").astype(bool)
        bucket_summary_df = pd.read_csv(table_dir / "stage5_bucket_summary.csv")
        exact_diag_summary_df = pd.read_csv(table_dir / "stage5_exact_diag_summary.csv")
        p_extra_hist_df = read_csv_or_empty(table_dir / "stage5_p_extra_histogram.csv")
        feasible_count_map = np.load(table_dir / "stage6_feasible_fiber_count_map.npy")
        stage56_outputs = load_json(stage56_meta).get("outputs", {})
        polymer_scale = float(stage56_outputs.get("polymer_scale", 0.0))
        lambda_value = float(stage56_outputs.get("lambda_value", 0.0))
        total_extra_probability = float(stage56_outputs.get("total_extra_probability", float(np.nansum(p_extra_map))))
    else:
        logger.info("STAGE 5 — Reconstructing target-like intensity and polymer baseline fields")
        if stratum_expected is None:
            raise ValueError(
                "--stratum-expected is required for the absolute-probability inference path because "
                "T_ij = stratum_expected[d_ij] * 2^(predicted smoothed log2 O/E)."
            )
        stage56_maps = build_absolute_probability_maps(
            crop_pred,
            stratum_expected=stratum_expected,
            min_diag_bins=args.min_diag_bins,
            polymer_alpha=args.polymer_alpha,
            local_diag_min=args.polymer_local_diag_min,
            local_diag_max=args.polymer_local_diag_max,
            local_statistic=args.polymer_local_statistic,
            p_local=args.p_local,
            kappa=args.kappa,
        )
        t_map = stage56_maps["t_map"]
        s_raw_map = stage56_maps["s_raw_map"]
        s_map = stage56_maps["s_map"]
        q_t_map = stage56_maps["q_t_map"]
        q_s_map = stage56_maps["q_s_map"]
        p_extra_map = stage56_maps["p_extra_map"]
        subtraction_map = stage56_maps["subtraction_map"]
        allowed_pair_mask = stage56_maps["allowed_map"].astype(bool)
        bucket_summary_df = stage5_bucket_summary(
            {
                "T": stage56_maps["t_vec"],
                "S_raw": stage56_maps["s_raw_vec"],
                "S": stage56_maps["s_vec"],
                "Q_T": stage56_maps["q_t_vec"],
                "Q_S": stage56_maps["q_s_vec"],
                "Q_T_minus_kappa_Q_S": stage56_maps["subtraction_vec"],
                "P_extra": stage56_maps["p_extra_vec"],
                "feasible_fibers": feasible_count_map[stage56_maps["ii"], stage56_maps["jj"]].astype(np.float32),
            },
            stage56_maps["diag"],
            stage56_maps["allowed"],
            crop_pred.shape[0],
        )
        exact_diag_summary_df = pd.concat([
            summarize_by_exact_diag(stage56_maps["q_t_vec"], stage56_maps["diag"], stage56_maps["allowed"], "Q_T"),
            summarize_by_exact_diag(stage56_maps["q_s_vec"], stage56_maps["diag"], stage56_maps["allowed"], "Q_S"),
            summarize_by_exact_diag(stage56_maps["subtraction_vec"], stage56_maps["diag"], stage56_maps["allowed"], "Q_T_minus_kappa_Q_S"),
            summarize_by_exact_diag(stage56_maps["p_extra_vec"], stage56_maps["diag"], stage56_maps["allowed"], "P_extra"),
        ], ignore_index=True)
        p_extra_hist_df = p_extra_histogram(stage56_maps["p_extra_vec"], stage56_maps["allowed"])
        polymer_scale = float(stage56_maps["polymer_scale"])
        lambda_value = float(stage56_maps["lambda_value"])
        total_extra_probability = float(np.nansum(stage56_maps["p_extra_vec"]))

        save_probability_count_plots(
            plot_dir,
            crop_pred,
            t_map,
            s_map,
            q_t_map,
            q_s_map,
            p_extra_map,
            subtraction_map,
            feasible_count_map,
            bucket_summary_df,
            exact_diag_summary_df,
            p_extra_hist_df,
            allele_name,
            args.polymer_alpha,
            args.kappa,
            args.polymer_local_diag_min,
            args.polymer_local_diag_max,
            feasible_mask is not None,
        )
        np.save(table_dir / "stage5_target_intensity_map.npy", t_map)
        np.save(table_dir / "stage5_polymer_baseline_raw_map.npy", s_raw_map)
        np.save(table_dir / "stage5_polymer_baseline_scaled_map.npy", s_map)
        np.save(table_dir / "stage5_target_occupancy_map.npy", q_t_map)
        np.save(table_dir / "stage5_polymer_occupancy_map.npy", q_s_map)
        np.save(table_dir / "stage5_extra_probability_map.npy", p_extra_map)
        np.save(table_dir / "stage5_subtraction_map.npy", subtraction_map)
        np.save(table_dir / "stage5_allowed_pair_mask.npy", allowed_pair_mask.astype(np.int8))
        np.save(table_dir / "stage6_feasible_fiber_count_map.npy", feasible_count_map.astype(np.int32))
        bucket_summary_df.to_csv(table_dir / "stage5_bucket_summary.csv", index=False)
        exact_diag_summary_df.to_csv(table_dir / "stage5_exact_diag_summary.csv", index=False)
        p_extra_hist_df.to_csv(table_dir / "stage5_p_extra_histogram.csv", index=False)
        write_stage_metadata(stage56_meta, stage56_inputs, {
            "polymer_scale": polymer_scale,
            "polymer_alpha": float(stage56_maps["polymer_alpha"]),
            "lambda_value": lambda_value,
            "total_extra_probability": total_extra_probability,
            "local_diag_min_used": int(stage56_maps["local_diag_min_used"]),
            "local_diag_max_used": int(stage56_maps["local_diag_max_used"]),
            "local_t_stat": float(stage56_maps["local_t_stat"]) if np.isfinite(stage56_maps["local_t_stat"]) else None,
            "local_s_raw_stat": float(stage56_maps["local_s_raw_stat"]) if np.isfinite(stage56_maps["local_s_raw_stat"]) else None,
            "local_s_stat": float(stage56_maps["local_s_stat"]) if np.isfinite(stage56_maps["local_s_stat"]) else None,
            "p_local_clipped": float(stage56_maps["p_local_clipped"]),
        })

    # --------------------------------------------------------
    # STAGE 8 — Bernoulli-per-fiber assignment and methyl*.in outputs
    # --------------------------------------------------------
    assignments_df = pd.DataFrame()
    skipped_df = pd.DataFrame()
    valency_df = pd.DataFrame()
    fiber_summary_df = pd.DataFrame()
    sampled_count_map = np.zeros_like(p_extra_map, dtype=np.int32)
    sampled_occ_map = np.zeros_like(p_extra_map, dtype=np.float32)
    recon_map = np.zeros_like(p_extra_map, dtype=np.int32)
    recon_occ_map = np.zeros_like(p_extra_map, dtype=np.float32)
    methyl_files_written: List[str] = []

    if dim_paths and fiber_bin_to_nucs is not None and centers_by_fiber is not None:
        stage8_meta = contact_dir / "stage8_metadata.json"
        stage8_files = [
            contact_dir / "stage8_assigned_contacts.csv",
            contact_dir / "stage8_unassigned_contacts.csv",
            contact_dir / "stage8_valency_table.csv",
            contact_dir / "stage8_fiber_assignment_summary.csv",
            contact_dir / "stage8_assignment_efficiency_by_bucket.csv",
            contact_dir / "stage8_sampled_count_map.npy",
            contact_dir / "stage8_reconstructed_count_map.npy",
            contact_dir / "stage8_sampled_occupancy_map.npy",
            contact_dir / "stage8_reconstructed_occupancy_map.npy",
            plot_dir / "stage8_reconstruction.png",
            plot_dir / "stage8_assignment_summaries.png",
        ]
        stage8_inputs = {
            "script_version": SCRIPT_VERSION,
            "allele": allele_name,
            "dim_files": dim_signatures,
            "methyl_output_dir": str(args.methyl_output_dir) if args.methyl_output_dir else None,
            "stage5_6_inputs": stage56_inputs,
            "vmax": int(args.vmax),
            "seed": int(allele_seed + 10_000),
            "assignment_mode": "per_fiber_stochastic_matching",
            "assignment_score_mode": "gumbel_perturbed_logit_p_extra",
        }
        if not args.force and metadata_matches(stage8_meta, stage8_inputs) and files_exist(stage8_files):
            logger.info("STAGE 8 — Reusing checkpointed per-fiber stochastic matching assignment")
            assignments_df = read_csv_or_empty(contact_dir / "stage8_assigned_contacts.csv")
            skipped_df = read_csv_or_empty(contact_dir / "stage8_unassigned_contacts.csv")
            valency_df = read_csv_or_empty(contact_dir / "stage8_valency_table.csv")
            fiber_summary_df = read_csv_or_empty(contact_dir / "stage8_fiber_assignment_summary.csv")
            sampled_count_map = np.load(contact_dir / "stage8_sampled_count_map.npy")
            recon_map = np.load(contact_dir / "stage8_reconstructed_count_map.npy")
            sampled_occ_map = np.load(contact_dir / "stage8_sampled_occupancy_map.npy")
            recon_occ_map = np.load(contact_dir / "stage8_reconstructed_occupancy_map.npy")
            methyl_files_written = write_methyl_contact_files(assignments_df, dim_paths, args.methyl_output_dir)
            logger.info("STAGE 8 — Wrote %d methyl*.in files", len(methyl_files_written))
        else:
            logger.info("STAGE 8 — Sampling per-fiber candidate contacts and resolving stochastic matching")
            fiber_ids, fiber_bin_to_nucs_by_index, centers_by_fiber_by_index, bin_to_fiber_indices = build_fiber_assignment_index(
                fiber_bin_to_nucs,
                centers_by_fiber,
                p_extra_map.shape[0],
            )
            assignments_df, sampled_count_map, recon_map, skipped_df, valency_df, fiber_summary_df = assign_contacts_capacity_aware(
                p_extra_map,
                fiber_ids,
                fiber_bin_to_nucs_by_index,
                centers_by_fiber_by_index,
                bin_to_fiber_indices,
                args.vmax,
                rng_assign,
            )
            sampled_occ_map = occupancy_from_counts(sampled_count_map, feasible_count_map)
            recon_occ_map = occupancy_from_counts(recon_map, feasible_count_map)
            assignments_df.to_csv(contact_dir / "stage8_assigned_contacts.csv", index=False)
            skipped_df.to_csv(contact_dir / "stage8_unassigned_contacts.csv", index=False)
            valency_df.to_csv(contact_dir / "stage8_valency_table.csv", index=False)
            fiber_summary_df.to_csv(contact_dir / "stage8_fiber_assignment_summary.csv", index=False)
            np.save(contact_dir / "stage8_sampled_count_map.npy", sampled_count_map)
            np.save(contact_dir / "stage8_reconstructed_count_map.npy", recon_map)
            np.save(contact_dir / "stage8_sampled_occupancy_map.npy", sampled_occ_map)
            np.save(contact_dir / "stage8_reconstructed_occupancy_map.npy", recon_occ_map)
            efficiency_df = assignment_efficiency_summary(
                p_extra_map,
                feasible_count_map,
                sampled_count_map,
                recon_map,
            )
            efficiency_df.to_csv(contact_dir / "stage8_assignment_efficiency_by_bucket.csv", index=False)
            methyl_files_written = write_methyl_contact_files(assignments_df, dim_paths, args.methyl_output_dir)
            logger.info("STAGE 8 — Wrote %d methyl*.in files", len(methyl_files_written))
            save_assignment_plots(plot_dir, p_extra_map, sampled_occ_map, recon_occ_map, valency_df, assignments_df, efficiency_df, fiber_summary_df, allele_name)
            write_stage_metadata(stage8_meta, stage8_inputs, {
                "n_assigned_contacts": int(len(assignments_df)),
                "n_unassigned_bin_pairs": int(len(skipped_df)),
                "methyl_files": methyl_files_written,
            })

    metadata = {
        "script_version": SCRIPT_VERSION,
        "stage3_inference_version": STAGE3_INFERENCE_VERSION,
        "allele": allele_name,
        "selected_feature_slot": primary_feature_slot,
        "included_marks": list(included_marks),
        "feature_marks": list(included_marks),
        "gate_marks": list(included_marks),
        "pair_gate_mode": "any_same_mark_match",
        "selected_feature_bw": [str(x) for x in bw_paths_by_mark.get(primary_feature_slot, [])],
        "bw_paths_by_mark": {mark: [str(x) for x in bw_paths_by_mark[mark]] for mark in included_marks},
        "model_joblib": str(args.model_joblib),
        "model_meta": str(args.model_meta) if args.model_meta else None,
        "model_target_mode": str(model_meta.get("target_mode", "unknown")),
        "residual_baseline": str(args.residual_baseline) if args.residual_baseline else None,
        "rf_output_space": "smoothed_oe_residual" if str(model_meta.get("target_mode", "unknown")) == "smoothed_oe_residual" else "smoothed_oe",
        "prediction_space": "smoothed_oe",
        "target_intensity_space": "T_ij = stratum_expected[d_ij] * 2^(predicted smoothed log2 O/E)",
        "polymer_baseline_space": "S_ij = C * d_ij^(-alpha)",
        "extra_probability_space": "P_extra = max(0, Q_T - kappa * Q_S)",
        "diag_correction_requested": str(args.diag_correction) if args.diag_correction else None,
        "diag_correction_used": str(resolve_diag_correction_path(args.diag_correction)) if resolve_diag_correction_path(args.diag_correction) else None,
        "pred_level_correction_enabled": not bool(args.disable_pred_level_correction),
        "pred_level_correction_requested": str(args.pred_level_correction) if args.pred_level_correction else None,
        "pred_level_correction_used": str(resolve_pred_level_correction_path(args.pred_level_correction, args.diag_correction, disabled=bool(args.disable_pred_level_correction))) if resolve_pred_level_correction_path(args.pred_level_correction, args.diag_correction, disabled=bool(args.disable_pred_level_correction)) else None,
        "stratum_expected": str(args.stratum_expected) if args.stratum_expected else None,
        "requested_region": {"chrom": requested_region.chrom, "start": requested_region.start, "end": requested_region.end},
        "prediction_window": {"chrom": prediction_region.chrom, "start": prediction_region.start, "end": prediction_region.end},
        "requested_bin_slice_in_prediction_window": {"b0": int(b0), "b1": int(b1)},
        "resolution_bp": REQUESTED_RESOLUTION,
        "normalization_mode": "dataset_specific_shared_feature_gate_threshold",
        "feature_threshold_source": "allele_dataset_genomewide",
        "methyl_gate_source": "allele_dataset_genomewide",
        "shared_feature_and_gate_threshold": True,
        "mark_stats": shared_stats.get("mark_stats", {}),
        "polymer_alpha": float(args.polymer_alpha),
        "polymer_local_diag_min": int(args.polymer_local_diag_min),
        "polymer_local_diag_max": int(args.polymer_local_diag_max),
        "polymer_local_statistic": str(args.polymer_local_statistic),
        "p_local": float(args.p_local),
        "kappa": float(args.kappa),
        "dim_glob": args.dim_glob,
        "xa_dim_glob": args.xa_dim_glob,
        "xi_dim_glob": args.xi_dim_glob,
        "resolved_dim_glob": resolved_dim_glob,
        "methyl_output_dir": str(args.methyl_output_dir) if args.methyl_output_dir else None,
        "feasible_mask_from_dim_files": bool(feasible_mask is not None),
        "n_feasible_bin_pairs": int(np.count_nonzero(np.triu(feasible_mask, k=1))) if feasible_mask is not None else 0,
        "assignment_seed": int(allele_seed + 10_000),
        "polymer_scale": float(polymer_scale),
        "lambda_value": float(lambda_value),
        "total_extra_probability": float(total_extra_probability),
        "total_sampled_contacts": int(np.sum(sampled_count_map)),
        "vmax": int(args.vmax),
        "assignment_mode": "per_fiber_stochastic_matching",
        "assignment_score_mode": "gumbel_perturbed_logit_p_extra",
        "n_assigned_contacts": int(len(assignments_df)) if not assignments_df.empty else 0,
        "n_unassigned_bin_pairs": int(len(skipped_df)) if not skipped_df.empty else 0,
        "methyl_files_written": methyl_files_written,
    }
    with open(allele_out / "run_metadata.json", "w") as fh:
        json.dump(metadata, fh, indent=2)

    logger.info("Finished allele %s. Outputs written to %s", allele_name, allele_out)


# ============================================================
# Main pipeline
# ============================================================
def main():
    parser = argparse.ArgumentParser(description="Predict Xa and Xi contact maps from allele-specific epigenetic tracks and assign contacts to numbered dim files.")
    parser.add_argument("--model-joblib", required=True, type=Path)
    parser.add_argument("--model-meta", type=Path, default=None)
    parser.add_argument("--residual-baseline", type=Path, default=Path("checkpoints/stratum_residual_baseline.npz"), help="NPZ with key 'baseline' used to reconstruct v5 residual RF outputs into predicted smoothed O/E.")
    parser.add_argument("--diag-correction", type=Path, default=None)
    parser.add_argument("--pred-level-correction", type=Path, default=None, help="Optional second-layer calibration curve CSV. If omitted, the predictor will look for pred_level_calibration_curve.csv next to the exact-diagonal correction.")
    parser.add_argument("--disable-pred-level-correction", action="store_true", help="Disable second-layer prediction-level calibration entirely and do not auto-discover sibling pred_level_calibration_curve.csv files.")
    parser.add_argument(
        "--stratum-expected",
        type=Path,
        default=None,
        help=(
            "NPZ file with key 'expected'. Required for the absolute-probability inference path because "
            "predicted smoothed log2(O/E) is converted back to an observed-like intensity T_ij using the "
            "training stratum expected curve."
        ),
    )
    parser.add_argument("--mcool-uri", type=str, default=None, help="Optional cooler URI for plotting true Raw/Smoothed O/E on the 100 kb prediction window.")
    parser.add_argument("--feature-slot", default="H3K27me3", choices=MARK_NAMES)
    parser.add_argument("--xa-bw", nargs="+", default=None, help="Allele-specific Xa BigWig file(s) for the selected feature in backward-compatible single-mark mode.")
    parser.add_argument("--xi-bw", nargs="+", default=None, help="Allele-specific Xi BigWig file(s) for the selected feature in backward-compatible single-mark mode.")
    parser.add_argument("--xa-mark-bw", action="append", default=None, help="Repeatable Xa mark-specific BigWig input in the form MARK=/path/to/file.bw")
    parser.add_argument("--xi-mark-bw", action="append", default=None, help="Repeatable Xi mark-specific BigWig input in the form MARK=/path/to/file.bw")
    parser.add_argument("--chrom", required=True)
    parser.add_argument("--start", required=True, type=int)
    parser.add_argument("--end", required=True, type=int)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--n-ens", type=int, default=100, help="Deprecated and ignored. The new algorithm does not use a fixed total contact budget.")
    parser.add_argument("--contacts-per-fiber", type=int, default=80, help="Deprecated and ignored. The new algorithm uses Bernoulli sampling per feasible fiber.")
    parser.add_argument("--distance-alpha", type=float, default=1.0, help="Deprecated and ignored. The new algorithm uses T_ij = stratum_expected[d] * 2^(predicted smoothed O/E).")
    parser.add_argument("--beta", type=float, default=1.0, help="Deprecated and ignored. The new algorithm no longer exponentiates predictions into globally normalized pair weights.")
    parser.add_argument("--vmax", type=int, default=1)
    parser.add_argument("--min-diag-bins", type=int, default=1)
    parser.add_argument("--polymer-alpha", type=float, default=1.0, help="Polymer baseline exponent alpha in S_raw = d^(-alpha). Larger alpha weakens long-range subtraction.")
    parser.add_argument("--polymer-local-diag-min", type=int, default=1, help="Minimum exact diagonal used to scale the polymer baseline and calibrate lambda.")
    parser.add_argument("--polymer-local-diag-max", type=int, default=3, help="Maximum exact diagonal used to scale the polymer baseline and calibrate lambda.")
    parser.add_argument("--polymer-local-statistic", choices=["median", "mean"], default="median", help="Statistic used on local diagonals when scaling the polymer baseline and lambda.")
    parser.add_argument("--p-local", type=float, default=0.8, help="Target local baseline occupancy used to calibrate lambda from the scaled polymer baseline.")
    parser.add_argument("--kappa", type=float, default=1.0, help="Baseline subtraction scale used in P_extra = max(0, Q_T - kappa Q_S). This is gamma in the manuscript text; the command-line name is kept for backward compatibility.")
    parser.add_argument("--seed", type=int, default=RANDOM_SEED)
    parser.add_argument(
        "--methyl-gate-quantile",
        "--methyl-threshold-quantile",
        dest="methyl_gate_quantile",
        type=float,
        default=80.0,
        help=(
            "High-signal raw-scale quantile over combined Xa+Xi positive genome-wide bins "
            "used for methyl/contact pair gating. The old --methyl-threshold-quantile "
            "name is kept as an alias; use 75 to reproduce older runs."
        ),
    )
    parser.add_argument(
        "--epigenetic-noise-floor-quantile",
        type=float,
        default=None,
        help="Deprecated and ignored. Feature denoising now always uses --methyl-gate-quantile.",
    )
    parser.add_argument("--dim-glob", type=str, default="dim*.in", help="Glob for numbered dim files, e.g. dim*.in or fibers/dim*.in")
    parser.add_argument("--xa-dim-glob", type=str, default=None, help="Optional Xa-specific dim glob. If omitted, Xa falls back to --dim-glob.")
    parser.add_argument("--xi-dim-glob", type=str, default=None, help="Optional Xi-specific dim glob. If omitted, Xi falls back to --dim-glob.")
    parser.add_argument("--methyl-output-dir", type=Path, default=None, help="Optional directory for methyl*.in outputs. Default writes next to each matched dim*.in file.")
    parser.add_argument("--fiber-region-start", type=int, default=None, help="Genomic start coordinate used for dim file mapping, equivalent to 'start' in the C code.")
    parser.add_argument("--force", action="store_true", help="Recompute all stages even if matching checkpoints already exist.")
    parser.add_argument("--force-shared", action="store_true", help="Recompute only shared Xa+Xi normalization even if its checkpoint is current.")
    args = parser.parse_args()

    root_out_dir = ensure_dir(args.output_dir)
    root_logger = setup_logger(root_out_dir / "run.log")
    if args.epigenetic_noise_floor_quantile is not None:
        root_logger.warning(
            "--epigenetic-noise-floor-quantile is deprecated and ignored; using --methyl-gate-quantile %.2f for both feature denoising and methyl pair gating.",
            args.methyl_gate_quantile,
        )
    if args.n_ens != 100 or args.contacts_per_fiber != 80 or args.distance_alpha != 1.0 or args.beta != 1.0:
        root_logger.warning(
            "--n-ens, --contacts-per-fiber, --distance-alpha, and --beta are deprecated and ignored in the "
            "absolute-probability inference path. Extra methylation contacts now emerge from P_extra and the "
            "available feasible fibers."
        )
    if args.polymer_local_diag_max < args.polymer_local_diag_min:
        raise ValueError("--polymer-local-diag-max must be >= --polymer-local-diag-min")
    if args.polymer_alpha <= 0.0:
        raise ValueError("--polymer-alpha must be > 0")
    if not (0.0 < args.p_local < 1.0):
        raise ValueError("--p-local must be in the open interval (0, 1)")
    if args.stratum_expected is None:
        raise ValueError("--stratum-expected is required for the absolute-probability inference path")

    if not args.xa_mark_bw and not args.xi_mark_bw:
        if not args.xa_bw or not args.xi_bw:
            raise ValueError(
                "Provide either backward-compatible single-mark inputs (--xa-bw/--xi-bw) "
                "or multi-mark inputs (--xa-mark-bw/--xi-mark-bw)."
            )

    allele_mark_bw_map, included_marks, primary_feature_slot, multi_mark_mode = resolve_allele_mark_bw_map(args)

    if not multi_mark_mode and args.feature_slot != "H3K27me3":
        root_logger.warning("This allele-specific mode is intended for H3K27me3. Proceeding with feature slot = %s", args.feature_slot)

    root_logger.info("Loading trained model from %s", args.model_joblib)
    model = joblib.load(args.model_joblib)
    model_meta = load_json(args.model_meta)
    feature_variant = infer_feature_variant(model, model_meta)
    root_logger.info("Using feature variant: %s", feature_variant)
    model_target_mode = str(model_meta.get("target_mode", "unknown"))
    if model_target_mode == "smoothed_oe_residual":
        if not args.residual_baseline.exists():
            raise FileNotFoundError(
                f"Model metadata target_mode={model_target_mode}, so residual baseline is required but was not found: "
                f"{args.residual_baseline}. Use --residual-baseline checkpoints/stratum_residual_baseline.npz."
            )
        root_logger.info("Using residual baseline for smoothed-O/E reconstruction: %s", args.residual_baseline)
    resolved_diag_correction = resolve_diag_correction_path(args.diag_correction)
    if args.diag_correction and resolved_diag_correction != args.diag_correction:
        root_logger.info("Using v5 exact-diagonal correction CSV fallback: %s", resolved_diag_correction)

    if pyBigWig is None:
        raise RuntimeError("pyBigWig is not installed, but BigWig input was requested")
    chrom_sizes = collect_shared_bigwig_chrom_sizes(allele_mark_bw_map, included_marks, root_logger)
    if args.chrom not in chrom_sizes:
        raise ValueError(
            f"Chromosome {args.chrom} is not present in all supplied Xa/Xi BigWig files. "
            "Multi-mark inference requires the requested chromosome to exist across the full input set."
        )

    requested_region = Region(args.chrom, int(args.start), int(args.end))
    prediction_region = resolve_centered_prediction_window(requested_region, int(chrom_sizes[args.chrom]))
    root_logger.info(
        "Requested region: %s:%d-%d | centered 100 kb prediction window: %s:%d-%d",
        requested_region.chrom, requested_region.start, requested_region.end,
        prediction_region.chrom, prediction_region.start, prediction_region.end,
    )

    shared_table_dir = ensure_dir(root_out_dir / "shared_tables")
    shared_plot_dir = ensure_dir(root_out_dir / "shared_plots")
    shared_norm_meta = shared_table_dir / "shared_normalization_metadata.json"
    shared_norm_csv = shared_table_dir / "stage1_combined_whole_genome_feature_stats.csv"
    shared_norm_inputs = {
        "script_version": SCRIPT_VERSION,
        "normalization_mode": "dataset_specific_shared_feature_gate_threshold",
        "shared_feature_and_gate_threshold": True,
        "feature_slot": primary_feature_slot,
        "included_marks": list(included_marks),
        "feature_marks": list(included_marks),
        "gate_marks": list(included_marks),
        "pair_gate_mode": "any_same_mark_match",
        "xa_mark_bw": {mark: [str(x) for x in allele_mark_bw_map["Xa"][mark]] for mark in included_marks},
        "xi_mark_bw": {mark: [str(x) for x in allele_mark_bw_map["Xi"][mark]] for mark in included_marks},
        "resolution": REQUESTED_RESOLUTION,
        "methyl_gate_quantile": float(args.methyl_gate_quantile),
        "epigenetic_noise_floor_quantile": float(args.methyl_gate_quantile),
        "requested_region": region_dict(requested_region),
        "prediction_window": region_dict(prediction_region),
    }
    shared_payload = load_json(shared_norm_meta).get("outputs", {}) if shared_norm_meta.exists() else {}
    if (
        not args.force and
        not args.force_shared and
        metadata_matches(shared_norm_meta, shared_norm_inputs) and
        bool(shared_payload.get("shared_feature_and_gate_threshold")) and
        "mark_stats" in shared_payload and
        files_exist([shared_norm_csv, shared_plot_dir / "stage1_combined_whole_genome_normalization.png"])
    ):
        root_logger.info("STAGE 1 — Reusing checkpointed shared Xa+Xi whole-genome normalization")
        shared_stats = dict(shared_payload)
        chrom_stats = pd.read_csv(shared_norm_csv)
    else:
        shared_stats, chrom_stats = compute_combined_whole_genome_feature_stats(
            allele_mark_bw_map,
            included_marks,
            REQUESTED_RESOLUTION,
            root_logger,
            methyl_gate_quantile=args.methyl_gate_quantile,
            epigenetic_noise_floor_quantile=args.methyl_gate_quantile,
        )
        chrom_stats.to_csv(shared_norm_csv, index=False)
        save_combined_normalization_plots(shared_plot_dir, chrom_stats, shared_stats)
        write_stage_metadata(
            shared_norm_meta,
            shared_norm_inputs,
            shared_stats,
        )

    clr = load_cooler(args.mcool_uri) if args.mcool_uri else None
    stratum_expected = np.load(args.stratum_expected)["expected"].astype(np.float32)

    # Save shared metadata
    shared_meta = {
        "script_version": SCRIPT_VERSION,
        "stage3_inference_version": STAGE3_INFERENCE_VERSION,
        "feature_slot": primary_feature_slot,
        "included_marks": list(included_marks),
        "feature_marks": list(included_marks),
        "gate_marks": list(included_marks),
        "pair_gate_mode": "any_same_mark_match",
        "requested_region": {"chrom": requested_region.chrom, "start": requested_region.start, "end": requested_region.end},
        "prediction_window": {"chrom": prediction_region.chrom, "start": prediction_region.start, "end": prediction_region.end},
        "model_joblib": str(args.model_joblib),
        "model_meta": str(args.model_meta) if args.model_meta else None,
        "model_target_mode": str(model_meta.get("target_mode", "unknown")),
        "residual_baseline": str(args.residual_baseline) if args.residual_baseline else None,
        "rf_output_space": "smoothed_oe_residual" if str(model_meta.get("target_mode", "unknown")) == "smoothed_oe_residual" else "smoothed_oe",
        "prediction_space": "smoothed_oe",
        "target_intensity_space": "T_ij = stratum_expected[d_ij] * 2^(predicted smoothed log2 O/E)",
        "polymer_baseline_space": "S_ij = C * d_ij^(-alpha)",
        "extra_probability_space": "P_extra = max(0, Q_T - kappa * Q_S)",
        "diag_correction_requested": str(args.diag_correction) if args.diag_correction else None,
        "diag_correction_used": str(resolved_diag_correction) if resolved_diag_correction else None,
        "stratum_expected": str(args.stratum_expected),
        "polymer_alpha": float(args.polymer_alpha),
        "polymer_local_diag_min": int(args.polymer_local_diag_min),
        "polymer_local_diag_max": int(args.polymer_local_diag_max),
        "polymer_local_statistic": str(args.polymer_local_statistic),
        "p_local": float(args.p_local),
        "kappa": float(args.kappa),
        "normalization_mode": "dataset_specific_shared_feature_gate_threshold",
        "feature_threshold_source": "allele_dataset_genomewide",
        "methyl_gate_source": "allele_dataset_genomewide",
        "shared_feature_and_gate_threshold": True,
        "mark_stats": shared_stats.get("mark_stats", {}),
        "xa_bw": [str(x) for x in args.xa_bw] if args.xa_bw else [],
        "xi_bw": [str(x) for x in args.xi_bw] if args.xi_bw else [],
        "xa_mark_bw": {mark: [str(x) for x in allele_mark_bw_map["Xa"][mark]] for mark in included_marks},
        "xi_mark_bw": {mark: [str(x) for x in allele_mark_bw_map["Xi"][mark]] for mark in included_marks},
        "dim_glob": args.dim_glob,
        "xa_dim_glob": args.xa_dim_glob,
        "xi_dim_glob": args.xi_dim_glob,
    }
    with open(root_out_dir / "shared_metadata.json", "w") as fh:
        json.dump(shared_meta, fh, indent=2)
    with open(root_out_dir / "shared_run_metadata.json", "w") as fh:
        json.dump(shared_meta, fh, indent=2)

    for allele_name, bw_paths_by_mark in allele_mark_bw_map.items():
        run_single_allele_prediction(
            allele_name=allele_name,
            bw_paths_by_mark=bw_paths_by_mark,
            included_marks=included_marks,
            primary_feature_slot=primary_feature_slot,
            args=args,
            model=model,
            model_meta=model_meta,
            feature_variant=feature_variant,
            requested_region=requested_region,
            prediction_region=prediction_region,
            shared_stats=shared_stats,
            root_out_dir=root_out_dir,
            clr=clr,
            stratum_expected=stratum_expected,
        )

    root_logger.info("Done. Outputs written to %s", root_out_dir)


if __name__ == "__main__":
    main()

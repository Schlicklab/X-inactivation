#!/bin/bash
# Example only. Replace paths with local files.

python3 scripts/predict_xic_h3k27me3_contacts.py   --model-joblib /path/to/model.joblib   --model-meta /path/to/model_meta.json   --residual-baseline /path/to/stratum_residual_baseline.npz   --diag-correction /path/to/exact_diag_calibration_curve.csv   --stratum-expected /path/to/stratum_expected.npz   --feature-slot H3K27me3   --xa-bw /path/to/Xa_H3K27me3.bw   --xi-bw /path/to/Xi_H3K27me3.bw   --chrom chrX   --start 103104661   --end 103194661   --output-dir xic_prediction_output   --dim-glob 'fibers/dim*.in'   --seed 42

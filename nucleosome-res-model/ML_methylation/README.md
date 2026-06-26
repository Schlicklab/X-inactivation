# Xic H3K27me3-associated contact prediction

This repository contains the Python scripts used to train and apply the Random Forest contact-prediction model used in the manuscript. The model predicts 200 bp Micro-C contact enrichment from one-dimensional epigenomic tracks and maps the predicted H3K27me3-associated enrichment to contact restraints for nucleosome-resolution simulations.

The repository is intentionally lightweight. Large public datasets, trained model files, and intermediate checkpoints are not stored in Git. The input datasets are listed in `metadata/data_manifest.tsv`; users should download them from the original repositories and place them in a local data directory.

## Contents

```text
scripts/train_rf_contact_model.py              Train Random Forest contact models from Micro-C and epigenomic tracks
scripts/predict_xic_h3k27me3_contacts.py      Apply a trained model to Xa/Xi H3K27me3 tracks and assign contacts
scripts/slurm_train_example.sh                Example SLURM script for training
scripts/predict_example.sh                    Example prediction command
metadata/data_manifest.tsv                    Dataset list used for training and prediction
metadata/best_params.json                     Final Random Forest hyperparameters
metadata/model_meta_summary.json              Short summary of the final model metadata
metadata/run_environment_example.json         Hardware/software information from the manuscript run
results/*.csv                                 Summary reports used in the manuscript
requirements.txt                              Python package versions from the manuscript run
```

## Data and model files

The raw bigWig and mcool files are not included because they are large public datasets. The trained model and checkpoint files may also be large; for publication, place them in a stable archive such as Zenodo, OSF, Figshare, institutional storage, or a GitHub Release/Git LFS if the files are small enough. The minimal files needed to reproduce allele-specific prediction from a trained model are:

```text
model.joblib
model_meta.json
stratum_residual_baseline.npz
stratum_expected.npz
exact_diag_calibration_curve.csv
pred_level_calibration_curve.csv, if used
```

The manuscript results were generated with the final hyperparameters listed in `metadata/best_params.json`.

## Training

Set the input data directory and run:

```bash
export ML_CONTACT_DATA_DIR=/path/to/input_data
export ML_CONTACT_CHECKPOINT_DIR=checkpoints
export ML_CONTACT_PLOT_DIR=plots_regression
python3 scripts/train_rf_contact_model.py
```

The script uses checkpoints. If a stage has already completed, the saved output is reused when the script is restarted.

## Allele-specific prediction

Use `scripts/predict_example.sh` as a template. The prediction script applies the trained model to Xa and Xi H3K27me3 tracks, computes a shared Xa+Xi normalization and H3K27me3 enrichment threshold, predicts contact enrichment, subtracts a polymer-distance baseline, and samples fiber-specific contact restraints.

## Software environment

The manuscript run used Python 3.13.12 with the package versions listed in `requirements.txt`. The job was run on one HPC node with 32 allocated CPU cores and 1000 GB requested memory; the completed SLURM job used 11:35:38 wall time and 155.84 GB memory.

## Notes

The scripts preserve the analysis logic used for the manuscript. They are provided for reproducibility and review, not as a general-purpose software package.

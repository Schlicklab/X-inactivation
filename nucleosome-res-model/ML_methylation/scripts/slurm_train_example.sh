#!/bin/bash
#SBATCH --job-name=ML_contact_train
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=1000GB
#SBATCH --time=15:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

module purge
module load singularity-ce/4.3.3  # or apptainer, depending on the cluster

export ML_CONTACT_DATA_DIR=/path/to/input_data
export ML_CONTACT_CHECKPOINT_DIR=checkpoints
export ML_CONTACT_PLOT_DIR=plots_regression

singularity exec   --overlay /path/to/overlay.ext3:ro   /path/to/container.sif   /bin/bash -c "source /ext3/env.sh; python3 scripts/train_rf_contact_model.py"

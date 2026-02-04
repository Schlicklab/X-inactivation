#!/bin/bash
#SBATCH --job-name=xi
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --mem=70GB
#SBATCH --time=05:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

RUNDIR=$(pwd)
cd $RUNDIR
# Load any necessary modules if required
# module load intel openmpi

# Run your job using the container
/share/apps/images/run-openmpi-5.0.7-intel.bash mpiexec -np 64 ./chrom_ap1.x dim.in input.1 LH_N0G6C22.in LH_N0G6C22_equil.in >> log_ch


#!/bin/bash
#SBATCH --job-name=snakemake_master
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=long
#SBATCH --time=24:00:00
#SBATCH --output=logs/slurm/snakemake_master_%j.out

# Activate your conda environment if needed
source activate snakemake_env

# Run snakemake
snakemake --profile slurm --use-conda
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --output=logs/%x_%A_%a.out
# One array stage of the simulations. The array index is the TASK argument of the driver.
#
#   sbatch --job-name=sim-pure --array=1-60%16  slurm/run_stage.sbatch 01_pure.R
#   sbatch --job-name=sim-semi --array=1-200%20 slurm/run_stage.sbatch 02_semi_real.R
#   sbatch --job-name=sim-null --array=1-20%8   slurm/run_stage.sbatch 03_null_gwas.R
#   sbatch --job-name=sim-fdr  --array=1-66%16  slurm/run_stage.sbatch 04_false_discovery.R
set -uo pipefail
export FLARE_SIMULATION=${FLARE_SIMULATION:-$PWD}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 JAX_PLATFORMS=cpu
mkdir -p logs
Rscript "$FLARE_SIMULATION/$1" "${SLURM_ARRAY_TASK_ID}"

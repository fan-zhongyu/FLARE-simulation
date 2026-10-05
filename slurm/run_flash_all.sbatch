#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x_%j.out
# FLASH on one family of data sets, one after another in a single job.
#   sbatch --job-name=flash-pure slurm/run_flash_all.sbatch pure 60
#   sbatch --job-name=flash-semi slurm/run_flash_all.sbatch semi 200
#   sbatch --job-name=flash-null slurm/run_flash_all.sbatch null 20
set -uo pipefail
export FLARE_SIMULATION=${FLARE_SIMULATION:-$PWD}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 JAX_PLATFORMS=cpu
mkdir -p logs
for t in $(seq 1 "$2"); do Rscript "$FLARE_SIMULATION/06_flash.R" "$1" "$t" 2>&1 | grep "_DONE\|Error"; done

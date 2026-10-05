#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=08:00:00
#SBATCH --output=logs/%x_%j.out
# GFA on all pure (300) and semi-real (200) data sets, one after another in a single job.
# A fit takes a few seconds, so an array is not needed. Data sets with a result file are skipped.
set -uo pipefail
export FLARE_SIMULATION=${FLARE_SIMULATION:-$PWD}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 JAX_PLATFORMS=cpu
mkdir -p logs
for t in $(seq 1 60);  do Rscript "$FLARE_SIMULATION/05_gfa.R" pure $t 2>&1 | grep "GFA_PURE_DONE\|Error"; done
for t in $(seq 1 200); do Rscript "$FLARE_SIMULATION/05_gfa.R" semi $t 2>&1 | grep "GFA_SEMI_DONE\|Error"; done

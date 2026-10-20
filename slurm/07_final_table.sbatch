#!/bin/bash

#SBATCH --job-name=final_table
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=80G
#SBATCH --time=12:00:00
#SBATCH --qos=medium

set -euo pipefail

source "${ATLASPIPE_ROOT:-${SLURM_SUBMIT_DIR:?}}/slurm/common.sh"
atlaspipe_init_slurm_step "07" "final_table"

echo "Building final contextualized abundance table..."

ml anaconda
conda activate Renv

cd "$PROJECT_ROOT"

Rscript scripts/07_final_table.R

echo "Final contextualized abundance table finished."

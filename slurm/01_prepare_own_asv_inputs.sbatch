#!/bin/bash

#SBATCH --job-name=prepare_own_asv
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --qos=medium

set -euo pipefail

source "${MICROBEATLAS_ROOT:-${SLURM_SUBMIT_DIR:?}}/slurm/common.sh"
microbeatlas_init_slurm_step "01" "prepare_own_asv_inputs"

module load anaconda
conda activate Renv

cd "$PROJECT_ROOT"

Rscript scripts/01_prepare_own_asv_inputs.R "$@"

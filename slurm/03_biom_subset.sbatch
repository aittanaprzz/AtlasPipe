#!/bin/bash

#SBATCH --job-name=biom_subset
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=120G
#SBATCH --time=24:00:00
#SBATCH --qos=medium

set -euo pipefail

source "${ATLASPIPE_ROOT:-${SLURM_SUBMIT_DIR:?}}/slurm/common.sh"
atlaspipe_init_slurm_step "03" "biom_subset"

echo "Starting BIOM subset extraction..."

ml anaconda
conda activate Biom

cd "$PROJECT_ROOT"

bash scripts/03_biom_subset.sh "$@"

echo "BIOM subset extraction finished."

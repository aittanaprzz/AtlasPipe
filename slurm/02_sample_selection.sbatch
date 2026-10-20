#!/bin/bash

#SBATCH --job-name=sample_selection
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --qos=medium

set -euo pipefail

source "${ATLASPIPE_ROOT:-${SLURM_SUBMIT_DIR:?}}/slurm/common.sh"
atlaspipe_init_slurm_step "02" "sample_selection"
cd "$PROJECT_ROOT"

DATA_DIR="${ATLASPIPE_DATA_DIR:-$PROJECT_ROOT/data/intermediate}"
IN_DIR="${MICROBEATLAS_RAW_DATA_DIR:-$PROJECT_ROOT/data/raw}"
BIOM_TABLE="${MICROBEATLAS_BIOM_TABLE:-$IN_DIR/samples-otus.97.metag.minfilter.rod2025companion.biom}"
BIOM_SAMPLE_IDS_FILE="${MICROBEATLAS_BIOM_SAMPLE_IDS:-$DATA_DIR/02_filters/biom_sample_ids.txt}"

python3 scripts/02_write_biom_sample_ids.py \
    --biom "$BIOM_TABLE" \
    --output "$BIOM_SAMPLE_IDS_FILE"

module load anaconda
conda activate Renv

# Example:
#   sbatch slurm/02_sample_selection.sbatch --env-id 30 --exclude-env-id 32 --samples-per-cc 1
args=("$@")
if [ "$#" -eq 0 ]; then
    if [ -n "${MICROBEATLAS_ENV_IDS:-}" ]; then
        args+=(--env-id "$MICROBEATLAS_ENV_IDS")
    fi
    if [ -n "${MICROBEATLAS_EXCLUDE_ENV_IDS:-}" ]; then
        args+=(--exclude-env-id "$MICROBEATLAS_EXCLUDE_ENV_IDS")
    fi
    if [ -n "${MICROBEATLAS_SAMPLES_PER_CC:-}" ]; then
        args+=(--samples-per-cc "$MICROBEATLAS_SAMPLES_PER_CC")
    fi
fi

Rscript scripts/02_sample_selection.R "${args[@]}"

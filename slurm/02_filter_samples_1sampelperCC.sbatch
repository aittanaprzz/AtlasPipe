#!/bin/bash

#SBATCH --job-name=filter_CC_1sample
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --qos=medium

set -euo pipefail

source "${MICROBEATLAS_ROOT:-${SLURM_SUBMIT_DIR:?}}/slurm/common.sh"
microbeatlas_init_slurm_step "02" "filter_samples"
cd "$PROJECT_ROOT"

DATA_DIR="${MICROBEATLAS_DATA_DIR:-$PROJECT_ROOT/data/intermediate}"
IN_DIR="${MICROBEATLAS_RAW_DATA_DIR:-$PROJECT_ROOT/data/raw/microbeatlas}"
BIOM_TABLE="${MICROBEATLAS_BIOM_TABLE:-$IN_DIR/samples-otus.97.mapped.biom}"
BIOM_SAMPLE_IDS_FILE="${MICROBEATLAS_BIOM_SAMPLE_IDS:-$DATA_DIR/02_filters/biom_sample_ids.txt}"

python3 scripts/02_write_biom_sample_ids.py \
    --biom "$BIOM_TABLE" \
    --output "$BIOM_SAMPLE_IDS_FILE"

module load anaconda
conda activate Renv

# Example:
#   sbatch slurm/02_filter_samples_1sampelperCC.sbatch --env-prefix aquatic --samples-per-cc 1
args=("$@")
if [ "$#" -eq 0 ]; then
    if [ -n "${MICROBEATLAS_ENV_PREFIX:-${MICROBEATLAS_ENV_LABEL:-}}" ]; then
        args+=(--env-prefix "${MICROBEATLAS_ENV_PREFIX:-$MICROBEATLAS_ENV_LABEL}")
    fi
    if [ -n "${MICROBEATLAS_SAMPLES_PER_CC:-}" ]; then
        args+=(--samples-per-cc "$MICROBEATLAS_SAMPLES_PER_CC")
    fi
fi

Rscript scripts/02_filter_samples_1sampelperCC.R "${args[@]}"

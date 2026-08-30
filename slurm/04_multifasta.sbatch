#!/bin/bash

#SBATCH --job-name=multifasta_otus97
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=80G
#SBATCH --time=12:00:00
#SBATCH --qos=medium

set -euo pipefail

source "${MICROBEATLAS_ROOT:-${SLURM_SUBMIT_DIR:?}}/slurm/common.sh"
microbeatlas_init_slurm_step "04" "multifasta"

echo "Generando multifasta de OTUs 97% de MicrobeAtlas..."

module load anaconda
conda activate Renv

cd "$PROJECT_ROOT"

Rscript scripts/04_multifasta.R "$@"

echo "Multifasta de OTUs 97% finalizado"

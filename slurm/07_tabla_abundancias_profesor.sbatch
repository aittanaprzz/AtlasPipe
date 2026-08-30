#!/bin/bash

#SBATCH --job-name=tabla_abundancias_profesor
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=80G
#SBATCH --time=12:00:00
#SBATCH --qos=medium

set -euo pipefail

source "${MICROBEATLAS_ROOT:-${SLURM_SUBMIT_DIR:?}}/slurm/common.sh"
microbeatlas_init_slurm_step "07" "tabla_abundancias_profesor"

echo "Generando tabla de abundancias OTU 97% siguiendo las indicaciones del profesor..."

ml anaconda
conda activate Renv

cd "$PROJECT_ROOT"

Rscript scripts/07_tabla_abundancias_profesor.R

echo "Tabla de abundancias OTU 97% finalizada"

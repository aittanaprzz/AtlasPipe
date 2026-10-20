#!/bin/bash

#SBATCH --job-name=reference_fasta
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=80G
#SBATCH --time=12:00:00
#SBATCH --qos=medium

set -euo pipefail

source "${ATLASPIPE_ROOT:-${SLURM_SUBMIT_DIR:?}}/slurm/common.sh"
atlaspipe_init_slurm_step "04" "reference_fasta"

echo "Building MicrobeAtlas OTU97 reference FASTA..."

module load anaconda
conda activate Renv

cd "$PROJECT_ROOT"

Rscript scripts/04_reference_fasta.R "$@"

echo "Reference FASTA finished."

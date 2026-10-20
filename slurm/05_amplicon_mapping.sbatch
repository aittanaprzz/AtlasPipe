#!/bin/bash
  
#SBATCH --job-name=amplicon_mapping
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
#SBATCH --time=24:00:00
#SBATCH --qos=medium

set -euo pipefail

source "${ATLASPIPE_ROOT:-${SLURM_SUBMIT_DIR:?}}/slurm/common.sh"
atlaspipe_init_slurm_step "05" "amplicon_mapping"

echo "Starting amplicon mapping against MicrobeAtlas OTU97..."

ml anaconda
conda activate MicrobeAtlas

cd "$PROJECT_ROOT"

bash scripts/05_amplicon_mapping.sh "$@"

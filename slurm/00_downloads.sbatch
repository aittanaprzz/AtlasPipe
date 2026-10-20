#!/bin/bash

#SBATCH --job-name=downloads
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --ntasks=1                    
#SBATCH --cpus-per-task=1             
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --qos=medium                 

set -euo pipefail

source "${ATLASPIPE_ROOT:-${SLURM_SUBMIT_DIR:?}}/slurm/common.sh"
atlaspipe_init_slurm_step "00" "downloads"
cd "$PROJECT_ROOT"

echo "Starting MicrobeAtlas downloads..."

bash scripts/00_downloads.sh

echo "Downloads finished."

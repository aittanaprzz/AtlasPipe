#!/bin/bash

#SBATCH --job-name=download_microb      
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --ntasks=1                    
#SBATCH --cpus-per-task=1             
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --qos=medium                 

set -euo pipefail

source "${MICROBEATLAS_ROOT:-${SLURM_SUBMIT_DIR:?}}/slurm/common.sh"
microbeatlas_init_slurm_step "00" "data_download"
cd "$PROJECT_ROOT"

echo "Iniciando el trabajo en el clúster Garnatxa..."

bash scripts/00_data_download.sh

echo "Trabajo de garnatxa finalizado"

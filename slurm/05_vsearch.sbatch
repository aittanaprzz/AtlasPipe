#!/bin/bash
  
#SBATCH --job-name=vsearch_global_search
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
#SBATCH --time=24:00:00
#SBATCH --qos=medium

set -euo pipefail

source "${MICROBEATLAS_ROOT:-${SLURM_SUBMIT_DIR:?}}/slurm/common.sh"
microbeatlas_init_slurm_step "05" "vsearch"

echo "Iniciando asignación de ASVs a OTUs 97% con VSEARCH en Garnatxa..."

ml anaconda
conda activate MicrobeAtlas

DATA_DIR="${MICROBEATLAS_DATA_DIR:-$PROJECT_ROOT/data/intermediate}"
ENV_LABEL="${MICROBEATLAS_ENV_LABEL:-aquatic}"
SAMPLES_PER_CC="${MICROBEATLAS_SAMPLES_PER_CC:-1}"

selection_label_from_samples_per_cc() {
    case "$1" in
        all|todos)
            printf "allsamplesperCC"
            ;;
        1)
            printf "1sampelperCC"
            ;;
        ''|*[!0-9]*)
            echo "MICROBEATLAS_SAMPLES_PER_CC debe ser un entero positivo o 'all': $1" >&2
            exit 1
            ;;
        *)
            if [ "$1" -lt 1 ]; then
                echo "MICROBEATLAS_SAMPLES_PER_CC debe ser un entero positivo o 'all': $1" >&2
                exit 1
            fi
            printf "%ssamplesperCC" "$1"
            ;;
    esac
}
SELECTION_LABEL="${MICROBEATLAS_SELECTION_LABEL:-$(selection_label_from_samples_per_cc "$SAMPLES_PER_CC")}"
if [ "$SELECTION_LABEL" = "1sampelperCC" ]; then
    VSEARCH_SUFFIX=""
else
    VSEARCH_SUFFIX="_${ENV_LABEL}_${SELECTION_LABEL}"
fi

IN_DIR="${MICROBEATLAS_VSEARCH_DIR:-$DATA_DIR/05_vsearch}"
DB_FILE="${MICROBEATLAS_OTU_FASTA:-$DATA_DIR/04_multifasta/MicrobeAtlas_selection_${ENV_LABEL}_${SELECTION_LABEL}_otus.97.fasta}"
ASV_FILE="${MICROBEATLAS_ASV_FASTA:-$IN_DIR/ASV.fasta}"
OUT_FILE="${MICROBEATLAS_ASV_OTU_MAP:-$IN_DIR/OTUs97_vsearch${VSEARCH_SUFFIX}.blast6}"
THREADS="${SLURM_CPUS_PER_TASK:-4}"

cd "$IN_DIR"

if [ ! -f "$ASV_FILE" ]; then
    echo "No existe el fichero de ASVs de entrada: $ASV_FILE"
    exit 1
fi

if [ ! -f "$DB_FILE" ]; then
    echo "No existe el multifasta de OTUs 97%: $DB_FILE"
    echo "Error en el paso 04 de generación del multifasta"
    exit 1
fi

vsearch --usearch_global "$ASV_FILE" \
  --db "$DB_FILE" \
  --id 0.97 \
  --blast6out "$OUT_FILE" \
  --strand both \
  --maxrejects 1000 \
  --threads "$THREADS" \
  --output_no_hits \
  --top_hits_only

TOTAL_ASV_COUNT=$(awk 'BEGIN { n = 0 } /^>/ { n++ } END { print n }' "$ASV_FILE")
HITS_COUNT=$(awk 'BEGIN { n = 0 } NF >= 2 && $2 != "*" { n++ } END { print n }' "$OUT_FILE")
NO_HITS_COUNT=$((TOTAL_ASV_COUNT - HITS_COUNT))
echo "ASVs totales: $TOTAL_ASV_COUNT"
echo "ASVs con hit: $HITS_COUNT"
echo "ASVs sin hit: $NO_HITS_COUNT"

echo "Asignación con VSEARCH finalizada: $OUT_FILE"

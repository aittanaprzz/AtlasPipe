#!/bin/bash

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${ATLASPIPE_ROOT:-$(cd "$script_dir/.." && pwd)}"
export ATLASPIPE_ROOT="$PROJECT_ROOT"

DATA_DIR="${ATLASPIPE_DATA_DIR:-$PROJECT_ROOT/data/intermediate}"
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

IN_DIR="${ATLASPIPE_OWN_AMPLICON_MAPPING_DIR:-$DATA_DIR/05_amplicon_mapping}"
DB_FILE="${MICROBEATLAS_OTU_FASTA:-$DATA_DIR/04_reference_fasta/MicrobeAtlas_selection_${ENV_LABEL}_${SELECTION_LABEL}_otus.97.fasta}"
AMPLICON_FILE="${ATLASPIPE_OWN_AMPLICON_FASTA:-$IN_DIR/amplicons.fasta}"
OUT_FILE="${ATLASPIPE_OWN_AMPLICON_OTU97_MAP:-$IN_DIR/OTUs97_vsearch${VSEARCH_SUFFIX}.blast6}"
THREADS="${ATLASPIPE_THREADS:-${SLURM_CPUS_PER_TASK:-4}}"

cd "$IN_DIR"

if [ ! -f "$AMPLICON_FILE" ]; then
    echo "No existe el FASTA de amplicones de entrada: $AMPLICON_FILE"
    exit 1
fi

if [ ! -f "$DB_FILE" ]; then
    echo "No existe el FASTA de referencia de OTUs 97%: $DB_FILE"
    echo "Error en el paso 04 de generación del FASTA de referencia"
    exit 1
fi

vsearch --usearch_global "$AMPLICON_FILE" \
  --db "$DB_FILE" \
  --id 0.97 \
  --blast6out "$OUT_FILE" \
  --strand both \
  --maxrejects 1000 \
  --threads "$THREADS" \
  --output_no_hits \
  --top_hits_only

TOTAL_AMPLICON_COUNT=$(awk 'BEGIN { n = 0 } /^>/ { n++ } END { print n }' "$AMPLICON_FILE")
HITS_COUNT=$(awk 'BEGIN { n = 0 } NF >= 2 && $2 != "*" { n++ } END { print n }' "$OUT_FILE")
NO_HITS_COUNT=$((TOTAL_AMPLICON_COUNT - HITS_COUNT))

echo "Amplicones totales: $TOTAL_AMPLICON_COUNT"
echo "Amplicones con hit OTU97: $HITS_COUNT"
echo "Amplicones sin hit OTU97: $NO_HITS_COUNT"
echo "Mapeo de amplicones finalizado: $OUT_FILE"

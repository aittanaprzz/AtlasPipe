#!/bin/bash

#SBATCH --job-name=biom_subset
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=120G
#SBATCH --time=24:00:00
#SBATCH --qos=medium

set -euo pipefail

source "${MICROBEATLAS_ROOT:-${SLURM_SUBMIT_DIR:?}}/slurm/common.sh"
microbeatlas_init_slurm_step "03" "extract_biome"

echo "Iniciando extracción de subconjuntos BIOM en Garnatxa..."

ml anaconda
conda activate Biom

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

cd "$PROJECT_ROOT"

IN_DIR="${MICROBEATLAS_RAW_DATA_DIR:-$PROJECT_ROOT/data/raw/microbeatlas}"
FILTER_FILE="${MICROBEATLAS_FILTER_FILE:-$DATA_DIR/02_filters_${ENV_LABEL}_${SELECTION_LABEL}/community_cluster_selection_${SELECTION_LABEL}_forbiomsubset.txt}"
OUT_DIR="${MICROBEATLAS_BIOM_SUBSET_DIR:-$DATA_DIR/03_samples_subsets}"
OUT_FILE="${MICROBEATLAS_BIOM_SUBSET_FILE:-$OUT_DIR/samples-otus.97_selection_${ENV_LABEL}_${SELECTION_LABEL}.biom}"
OUT_TSV="${MICROBEATLAS_BIOM_SUBSET_TSV:-$OUT_DIR/samples-otus.97_selection_${ENV_LABEL}_${SELECTION_LABEL}.tsv}"
SUMMARY_FILE="${MICROBEATLAS_BIOM_SUBSET_SUMMARY:-$OUT_DIR/biom_subset_summary_${ENV_LABEL}_${SELECTION_LABEL}.tsv}"
BIOM_TABLE="${MICROBEATLAS_BIOM_TABLE:-$IN_DIR/samples-otus.97.mapped.biom}"

mkdir -p "$OUT_DIR"

if [ ! -f "$BIOM_TABLE" ]; then
    echo "No existe la tabla BIOM de entrada: $BIOM_TABLE"
    echo "Error en el paso 00 de descarga de datos"
    exit 1
fi

if [ ! -f "$FILTER_FILE" ]; then
    echo "No existe el fichero de filtros: $FILTER_FILE"
    echo "Error en el paso 02 de filtrado de muestras"
    exit 1
fi

FILTER_LINES=$(awk 'NF > 0 { n++ } END { print n + 0 }' "$FILTER_FILE")
FILTER_UNIQUE_SAMPLES=$(awk 'NF > 0 { seen[$1] = 1 } END { print length(seen) }' "$FILTER_FILE")
FILTER_DUPLICATE_SAMPLES=$((FILTER_LINES - FILTER_UNIQUE_SAMPLES))

echo "Resumen previo a extracción BIOM:"
echo "  Muestras solicitadas en filtro: $FILTER_LINES"
echo "  Muestras únicas en filtro: $FILTER_UNIQUE_SAMPLES"
echo "  Muestras duplicadas en filtro: $FILTER_DUPLICATE_SAMPLES"

biom subset-table -i "$BIOM_TABLE" -a sample -s "$FILTER_FILE" -o "$OUT_FILE"
biom convert -i "$OUT_FILE" -o "$OUT_TSV" --to-tsv

read -r TSV_SAMPLES TSV_FEATURES_TOTAL TSV_FEATURES_97 TSV_UNMAPPED_ROWS TSV_UNMAPPED_SUM < <(
    awk '
        BEGIN {
            FS = "\t"
            samples = 0
            features_total = 0
            features_97 = 0
            unmapped_rows = 0
            unmapped_sum = 0
        }
        NR == 2 {
            samples = NF - 1
            next
        }
        NR > 2 {
            features_total++
            if ($1 ~ /^97_/) {
                features_97++
            } else if ($1 == "Unmapped") {
                unmapped_rows++
                for (i = 2; i <= NF; i++) {
                    unmapped_sum += $i
                }
            }
        }
        END {
            print samples, features_total, features_97, unmapped_rows, unmapped_sum
        }
    ' "$OUT_TSV"
)

{
    printf "parameter\tvalue\n"
    printf "env_label\t%s\n" "$ENV_LABEL"
    printf "samples_per_cc\t%s\n" "$SAMPLES_PER_CC"
    printf "selection_label\t%s\n" "$SELECTION_LABEL"
    printf "filter_file\t%s\n" "$FILTER_FILE"
    printf "input_biom\t%s\n" "$BIOM_TABLE"
    printf "output_biom\t%s\n" "$OUT_FILE"
    printf "output_tsv\t%s\n" "$OUT_TSV"
    printf "filter_lines\t%s\n" "$FILTER_LINES"
    printf "filter_unique_samples\t%s\n" "$FILTER_UNIQUE_SAMPLES"
    printf "filter_duplicate_samples\t%s\n" "$FILTER_DUPLICATE_SAMPLES"
    printf "tsv_samples\t%s\n" "$TSV_SAMPLES"
    printf "tsv_features_total\t%s\n" "$TSV_FEATURES_TOTAL"
    printf "tsv_features_97\t%s\n" "$TSV_FEATURES_97"
    printf "tsv_unmapped_rows\t%s\n" "$TSV_UNMAPPED_ROWS"
    printf "tsv_unmapped_sum\t%s\n" "$TSV_UNMAPPED_SUM"
} > "$SUMMARY_FILE"

echo "Resumen de extracción BIOM:"
cat "$SUMMARY_FILE"

echo "Extracción de subconjuntos BIOM finalizada."

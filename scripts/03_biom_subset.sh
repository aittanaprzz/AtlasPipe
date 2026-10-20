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

cd "$PROJECT_ROOT"

IN_DIR="${MICROBEATLAS_RAW_DATA_DIR:-$PROJECT_ROOT/data/raw}"
FILTER_FILE="${MICROBEATLAS_FILTER_FILE:-$DATA_DIR/02_filters_${ENV_LABEL}_${SELECTION_LABEL}/community_cluster_selection_${SELECTION_LABEL}_forbiomsubset.txt}"
OUT_DIR="${MICROBEATLAS_BIOM_SUBSET_DIR:-$DATA_DIR/03_samples_subsets}"
OUT_FILE="${MICROBEATLAS_BIOM_SUBSET_FILE:-$OUT_DIR/samples-otus.97_selection_${ENV_LABEL}_${SELECTION_LABEL}.biom}"
OUT_TSV="${MICROBEATLAS_BIOM_SUBSET_TSV:-$OUT_DIR/samples-otus.97_selection_${ENV_LABEL}_${SELECTION_LABEL}.tsv}"
SUMMARY_FILE="${MICROBEATLAS_BIOM_SUBSET_SUMMARY:-$OUT_DIR/biom_subset_summary_${ENV_LABEL}_${SELECTION_LABEL}.tsv}"
BIOM_TABLE="${MICROBEATLAS_BIOM_TABLE:-$IN_DIR/samples-otus.97.metag.minfilter.rod2025companion.biom}"

mkdir -p "$OUT_DIR"

if [ ! -f "$BIOM_TABLE" ]; then
    echo "Input BIOM table does not exist: $BIOM_TABLE"
    echo "Check step 00 downloads or MICROBEATLAS_BIOM_TABLE."
    exit 1
fi

if [ ! -f "$FILTER_FILE" ]; then
    echo "Sample filter file does not exist: $FILTER_FILE"
    echo "Check step 02 sample selection."
    exit 1
fi

FILTER_LINES=$(awk 'NF > 0 { n++ } END { print n + 0 }' "$FILTER_FILE")
FILTER_UNIQUE_SAMPLES=$(awk 'NF > 0 { seen[$1] = 1 } END { print length(seen) }' "$FILTER_FILE")
FILTER_DUPLICATE_SAMPLES=$((FILTER_LINES - FILTER_UNIQUE_SAMPLES))

echo "Pre-subset BIOM summary:"
echo "  Requested sample rows in filter: $FILTER_LINES"
echo "  Unique samples in filter: $FILTER_UNIQUE_SAMPLES"
echo "  Duplicate sample rows in filter: $FILTER_DUPLICATE_SAMPLES"

biom subset-table -i "$BIOM_TABLE" -a sample -s "$FILTER_FILE" -o "$OUT_FILE"
biom convert -i "$OUT_FILE" -o "$OUT_TSV" --to-tsv

read -r TSV_HEADER_UNIQUE_SAMPLES TSV_HEADER_DUPLICATE_SAMPLES TSV_MISSING_FILTER_SAMPLES TSV_EXTRA_TSV_SAMPLES < <(
    awk '
        BEGIN {
            FS = "\t"
        }
        NR == FNR {
            if (NF > 0) {
                want[$1] = 1
            }
            next
        }
        FNR == 2 {
            for (i = 2; i <= NF; i++) {
                got_count[$i]++
                got[$i] = 1
            }
        }
        END {
            duplicate_tsv = 0
            for (sample in got_count) {
                if (got_count[sample] > 1) {
                    duplicate_tsv++
                }
            }
            missing = 0
            for (sample in want) {
                if (!(sample in got)) {
                    missing++
                }
            }
            extra = 0
            for (sample in got) {
                if (!(sample in want)) {
                    extra++
                }
            }
            print length(got), duplicate_tsv, missing, extra
        }
    ' "$FILTER_FILE" "$OUT_TSV"
)

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
            if ($1 ~ /^[ABE]?97_/) {
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
    printf "tsv_header_unique_samples\t%s\n" "$TSV_HEADER_UNIQUE_SAMPLES"
    printf "tsv_header_duplicate_samples\t%s\n" "$TSV_HEADER_DUPLICATE_SAMPLES"
    printf "tsv_missing_sample_selection\t%s\n" "$TSV_MISSING_FILTER_SAMPLES"
    printf "tsv_extra_samples\t%s\n" "$TSV_EXTRA_TSV_SAMPLES"
    printf "tsv_features_total\t%s\n" "$TSV_FEATURES_TOTAL"
    printf "tsv_features_97\t%s\n" "$TSV_FEATURES_97"
    printf "tsv_unmapped_rows\t%s\n" "$TSV_UNMAPPED_ROWS"
    printf "tsv_unmapped_sum\t%s\n" "$TSV_UNMAPPED_SUM"
} > "$SUMMARY_FILE"

if [ "$FILTER_DUPLICATE_SAMPLES" -ne 0 ] ||
   [ "$TSV_HEADER_DUPLICATE_SAMPLES" -ne 0 ] ||
   [ "$TSV_MISSING_FILTER_SAMPLES" -ne 0 ] ||
   [ "$TSV_EXTRA_TSV_SAMPLES" -ne 0 ]; then
    echo "The BIOM subset does not exactly match the selected samples."
    echo "  Duplicate sample rows in filter: $FILTER_DUPLICATE_SAMPLES"
    echo "  Duplicate samples in TSV header: $TSV_HEADER_DUPLICATE_SAMPLES"
    echo "  Requested samples missing from TSV: $TSV_MISSING_FILTER_SAMPLES"
    echo "  Extra samples in TSV: $TSV_EXTRA_TSV_SAMPLES"
    exit 1
fi

echo "BIOM subset summary:"
cat "$SUMMARY_FILE"

echo "BIOM subset extraction finished."

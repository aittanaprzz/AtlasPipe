#!/bin/bash
set -euo pipefail

PROJECT_ROOT="${ATLASPIPE_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$PROJECT_ROOT"

if [ -f "$PROJECT_ROOT/config/env.sh" ]; then
    # shellcheck disable=SC1091
    source "$PROJECT_ROOT/config/env.sh"
else
    echo "No local config found at config/env.sh."
    echo "Create one with: cp config/example_env.sh config/env.sh"
    exit 1
fi

ATLASPIPE_RUN_LABEL="${ATLASPIPE_RUN_LABEL:-default}"
ATLASPIPE_LOG_DIR="${ATLASPIPE_LOG_DIR:-$PROJECT_ROOT/log/$ATLASPIPE_RUN_LABEL}"
export ATLASPIPE_RUN_LABEL ATLASPIPE_LOG_DIR
mkdir -p "$ATLASPIPE_LOG_DIR"

ATLASPIPE_INPUT_DIR="${ATLASPIPE_INPUT_DIR:-$PROJECT_ROOT/data/input}"
MICROBEATLAS_RAW_DATA_DIR="${MICROBEATLAS_RAW_DATA_DIR:-$PROJECT_ROOT/data/raw}"
ATLASPIPE_INTERMEDIATE_DIR="${ATLASPIPE_INTERMEDIATE_DIR:-$PROJECT_ROOT/data/intermediate/$ATLASPIPE_RUN_LABEL}"
ATLASPIPE_RESULTS_DIR="${ATLASPIPE_RESULTS_DIR:-$PROJECT_ROOT/results/$ATLASPIPE_RUN_LABEL}"
ATLASPIPE_DATA_DIR="${ATLASPIPE_DATA_DIR:-$ATLASPIPE_INTERMEDIATE_DIR}"
MICROBEATLAS_BIOM_TABLE="${MICROBEATLAS_BIOM_TABLE:-$MICROBEATLAS_RAW_DATA_DIR/samples-otus.97.metag.minfilter.rod2025companion.biom}"
ATLASPIPE_OWN_AMPLICON_TABLE_RAW="${ATLASPIPE_OWN_AMPLICON_TABLE_RAW:-$ATLASPIPE_INPUT_DIR/final_table_filtered.tsv}"
export ATLASPIPE_INPUT_DIR MICROBEATLAS_RAW_DATA_DIR ATLASPIPE_INTERMEDIATE_DIR
export ATLASPIPE_RESULTS_DIR ATLASPIPE_DATA_DIR MICROBEATLAS_BIOM_TABLE
export ATLASPIPE_OWN_AMPLICON_TABLE_RAW

preflight() {
    local missing=0
    local required_paths=(
        "scripts/00_downloads.sh"
        "scripts/01_prepare_inputs.R"
        "scripts/02_write_biom_sample_ids.py"
        "scripts/02_sample_selection.R"
        "scripts/03_biom_subset.sh"
        "scripts/04_reference_fasta.R"
        "scripts/05_amplicon_mapping.sh"
        "scripts/06_novel_clustering_parse.R"
        "scripts/07_final_table.R"
        "slurm/00_downloads.sbatch"
        "slurm/01_prepare_inputs.sbatch"
        "slurm/02_sample_selection.sbatch"
        "slurm/03_biom_subset.sbatch"
        "slurm/04_reference_fasta.sbatch"
        "slurm/05_amplicon_mapping.sbatch"
        "slurm/06_novel_clustering.sbatch"
        "slurm/07_final_table.sbatch"
    )

    if [ "$ATLASPIPE_RUN_LABEL" = "default" ]; then
        echo "ATLASPIPE_RUN_LABEL cannot be 'default' for a real run." >&2
        missing=1
    fi

    if [ ! -s "$ATLASPIPE_OWN_AMPLICON_TABLE_RAW" ]; then
        echo "User amplicon table does not exist or is empty: $ATLASPIPE_OWN_AMPLICON_TABLE_RAW" >&2
        missing=1
    fi

    if [[ "$(basename "$MICROBEATLAS_BIOM_TABLE")" != *minfilter* ]]; then
        echo "MICROBEATLAS_BIOM_TABLE does not appear to point to a minfilter table:" >&2
        echo "  $MICROBEATLAS_BIOM_TABLE" >&2
        missing=1
    fi

    local default_biom="$MICROBEATLAS_RAW_DATA_DIR/samples-otus.97.metag.minfilter.rod2025companion.biom"
    if [ "$MICROBEATLAS_BIOM_TABLE" != "$default_biom" ] && [ ! -s "$MICROBEATLAS_BIOM_TABLE" ]; then
        echo "MICROBEATLAS_BIOM_TABLE points to a custom path that does not exist or is empty:" >&2
        echo "  $MICROBEATLAS_BIOM_TABLE" >&2
        echo "Fix that path or use the default minfilter BIOM: $default_biom" >&2
        missing=1
    fi

    mkdir -p "$MICROBEATLAS_RAW_DATA_DIR" "$ATLASPIPE_INTERMEDIATE_DIR" "$ATLASPIPE_RESULTS_DIR" "$ATLASPIPE_LOG_DIR"

    for path in "${required_paths[@]}"; do
        if [ ! -e "$PROJECT_ROOT/$path" ]; then
            echo "Missing required pipeline file: $path" >&2
            missing=1
        fi
    done

    if [ -z "${MICROBEATLAS_ENV_IDS:-}" ]; then
        echo "MICROBEATLAS_ENV_IDS is not set; step 02 would select all environments." >&2
        missing=1
    fi

    if [ "$missing" -ne 0 ]; then
        echo "Preflight failed. No SLURM jobs were submitted." >&2
        exit 1
    fi
}

preflight

sbatch_step() {
    local step_id="$1"
    local step_name="$2"
    local dependency="${3:-}"
    local script="$4"
    shift 4

    local args=(
        --parsable
        --export="ALL,ATLASPIPE_ROOT=$PROJECT_ROOT,ATLASPIPE_RUN_LABEL=$ATLASPIPE_RUN_LABEL,ATLASPIPE_LOG_DIR=$ATLASPIPE_LOG_DIR"
        --output="$ATLASPIPE_LOG_DIR/${step_id}_${step_name}_%j.out"
        --error="$ATLASPIPE_LOG_DIR/${step_id}_${step_name}_%j.err"
    )

    if [ -n "$dependency" ]; then
        args+=(--dependency=afterok:"$dependency")
    fi

    sbatch "${args[@]}" "$script" "$@"
}

jid00=$(sbatch_step 00 downloads "" slurm/00_downloads.sbatch)
jid01=$(sbatch_step 01 prepare_inputs "$jid00" slurm/01_prepare_inputs.sbatch)
jid02=$(sbatch_step 02 sample_selection "$jid01" slurm/02_sample_selection.sbatch)
jid03=$(sbatch_step 03 biom_subset "$jid02" slurm/03_biom_subset.sbatch)
jid04=$(sbatch_step 04 reference_fasta "$jid03" slurm/04_reference_fasta.sbatch)
jid05=$(sbatch_step 05 amplicon_mapping "$jid04" slurm/05_amplicon_mapping.sbatch)
jid06=$(sbatch_step 06 novel_clustering "$jid05" slurm/06_novel_clustering.sbatch)
jid07=$(sbatch_step 07 final_table "$jid06" slurm/07_final_table.sbatch)

echo "Submitted SLURM pipeline jobs:"
echo "  00 downloads:              $jid00"
echo "  01 prepare inputs:         $jid01"
echo "  02 sample selection:       $jid02"
echo "  03 BIOM subset:            $jid03"
echo "  04 reference FASTA:        $jid04"
echo "  05 amplicon mapping:       $jid05"
echo "  06 novel clustering:       $jid06"
echo "  07 final table:            $jid07"
echo "Logs:"
echo "  $ATLASPIPE_LOG_DIR"

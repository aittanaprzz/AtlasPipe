#!/bin/bash
set -euo pipefail

PROJECT_ROOT="${MICROBEATLAS_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$PROJECT_ROOT"

if [ -f "$PROJECT_ROOT/config/env.sh" ]; then
    # shellcheck disable=SC1091
    source "$PROJECT_ROOT/config/env.sh"
else
    echo "No local config found at config/env.sh."
    echo "Create one with: cp config/example_env.sh config/env.sh"
    exit 1
fi

MICROBEATLAS_RUN_LABEL="${MICROBEATLAS_RUN_LABEL:-default}"
MICROBEATLAS_LOG_DIR="${MICROBEATLAS_LOG_DIR:-$PROJECT_ROOT/log/$MICROBEATLAS_RUN_LABEL}"
export MICROBEATLAS_RUN_LABEL MICROBEATLAS_LOG_DIR
mkdir -p "$MICROBEATLAS_LOG_DIR"

sbatch_step() {
    local step_id="$1"
    local step_name="$2"
    local dependency="${3:-}"
    local script="$4"
    shift 4

    local args=(
        --parsable
        --export="ALL,MICROBEATLAS_ROOT=$PROJECT_ROOT,MICROBEATLAS_RUN_LABEL=$MICROBEATLAS_RUN_LABEL,MICROBEATLAS_LOG_DIR=$MICROBEATLAS_LOG_DIR"
        --output="$MICROBEATLAS_LOG_DIR/${step_id}_${step_name}_%j.out"
        --error="$MICROBEATLAS_LOG_DIR/${step_id}_${step_name}_%j.err"
    )

    if [ -n "$dependency" ]; then
        args+=(--dependency=afterok:"$dependency")
    fi

    sbatch "${args[@]}" "$script" "$@"
}

jid00=$(sbatch_step 00 data_download "" slurm/00_data_download.sbatch)
jid01=$(sbatch_step 01 prepare_own_asv_inputs "$jid00" slurm/01_prepare_own_asv_inputs.sbatch)
jid02=$(sbatch_step 02 filter_samples "$jid01" slurm/02_filter_samples_1sampelperCC.sbatch)
jid03=$(sbatch_step 03 extract_biome "$jid02" slurm/03_extract_biome.sbatch)
jid04=$(sbatch_step 04 multifasta "$jid03" slurm/04_multifasta.sbatch)
jid05=$(sbatch_step 05 vsearch "$jid04" slurm/05_vsearch.sbatch)
jid06=$(sbatch_step 06 vsearch_clustering "$jid05" slurm/06_vsearch_clustering.sbatch)
jid07=$(sbatch_step 07 tabla_abundancias_profesor "$jid06" slurm/07_tabla_abundancias_profesor.sbatch)

echo "Submitted SLURM pipeline jobs:"
echo "  00 data download:          $jid00"
echo "  01 prepare own ASV input:  $jid01"
echo "  02 filter samples:         $jid02"
echo "  03 extract BIOM subset:    $jid03"
echo "  04 make multifasta:        $jid04"
echo "  05 VSEARCH mapping:        $jid05"
echo "  06 cluster no-hit ASVs:    $jid06"
echo "  07 final abundance table:  $jid07"
echo "Logs:"
echo "  $MICROBEATLAS_LOG_DIR"

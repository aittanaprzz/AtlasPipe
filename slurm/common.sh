#!/bin/bash

atlaspipe_init_slurm_step() {
    local step_id="$1"
    local step_name="$2"
    local caller="${BASH_SOURCE[1]:-${BASH_SOURCE[0]}}"
    local submit_root="${ATLASPIPE_ROOT:-${SLURM_SUBMIT_DIR:-}}"
    local job_id

    if [ -n "$submit_root" ] && [ -d "$submit_root/slurm" ]; then
        PROJECT_ROOT="$(cd "$submit_root" && pwd)"
        SLURM_DIR="$PROJECT_ROOT/slurm"
    else
        SLURM_DIR="$(cd "$(dirname "$caller")" && pwd)"
        PROJECT_ROOT="$(cd "$SLURM_DIR/.." && pwd)"
    fi
    export ATLASPIPE_ROOT="$PROJECT_ROOT"

    if [ -f "$PROJECT_ROOT/config/env.sh" ]; then
        # shellcheck disable=SC1091
        source "$PROJECT_ROOT/config/env.sh"
        PROJECT_ROOT="${ATLASPIPE_ROOT:-$PROJECT_ROOT}"
        export ATLASPIPE_ROOT="$PROJECT_ROOT"
    fi

    ATLASPIPE_RUN_LABEL="${ATLASPIPE_RUN_LABEL:-default}"
    export ATLASPIPE_RUN_LABEL

    ATLASPIPE_LOG_DIR="${ATLASPIPE_LOG_DIR:-$PROJECT_ROOT/log/$ATLASPIPE_RUN_LABEL}"
    export ATLASPIPE_LOG_DIR
    mkdir -p "$ATLASPIPE_LOG_DIR"

    job_id="${SLURM_JOB_ID:-manual_$(date +%Y%m%d_%H%M%S)}"
    ATLASPIPE_STDOUT_LOG="$ATLASPIPE_LOG_DIR/${step_id}_${step_name}_${job_id}.out"
    ATLASPIPE_STDERR_LOG="$ATLASPIPE_LOG_DIR/${step_id}_${step_name}_${job_id}.err"
    export ATLASPIPE_STDOUT_LOG ATLASPIPE_STDERR_LOG

    exec >"$ATLASPIPE_STDOUT_LOG" 2>"$ATLASPIPE_STDERR_LOG"

    echo "AtlasPipe step: ${step_id}_${step_name}"
    echo "Run label: $ATLASPIPE_RUN_LABEL"
    echo "Project root: $PROJECT_ROOT"
    echo "Log directory: $ATLASPIPE_LOG_DIR"
    echo "Stdout log: $ATLASPIPE_STDOUT_LOG"
    echo "Stderr log: $ATLASPIPE_STDERR_LOG"
    echo "Started: $(date -Iseconds)"
}

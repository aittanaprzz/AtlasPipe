#!/bin/bash

microbeatlas_init_slurm_step() {
    local step_id="$1"
    local step_name="$2"
    local caller="${BASH_SOURCE[1]:-${BASH_SOURCE[0]}}"
    local submit_root="${MICROBEATLAS_ROOT:-${SLURM_SUBMIT_DIR:-}}"
    local job_id

    if [ -n "$submit_root" ] && [ -d "$submit_root/slurm" ]; then
        PROJECT_ROOT="$(cd "$submit_root" && pwd)"
        SLURM_DIR="$PROJECT_ROOT/slurm"
    else
        SLURM_DIR="$(cd "$(dirname "$caller")" && pwd)"
        PROJECT_ROOT="$(cd "$SLURM_DIR/.." && pwd)"
    fi
    export MICROBEATLAS_ROOT="$PROJECT_ROOT"

    if [ -f "$PROJECT_ROOT/config/env.sh" ]; then
        # shellcheck disable=SC1091
        source "$PROJECT_ROOT/config/env.sh"
        PROJECT_ROOT="${MICROBEATLAS_ROOT:-$PROJECT_ROOT}"
        export MICROBEATLAS_ROOT="$PROJECT_ROOT"
    fi

    MICROBEATLAS_RUN_LABEL="${MICROBEATLAS_RUN_LABEL:-default}"
    export MICROBEATLAS_RUN_LABEL

    MICROBEATLAS_LOG_DIR="${MICROBEATLAS_LOG_DIR:-$PROJECT_ROOT/log/$MICROBEATLAS_RUN_LABEL}"
    export MICROBEATLAS_LOG_DIR
    mkdir -p "$MICROBEATLAS_LOG_DIR"

    job_id="${SLURM_JOB_ID:-manual_$(date +%Y%m%d_%H%M%S)}"
    MICROBEATLAS_STDOUT_LOG="$MICROBEATLAS_LOG_DIR/${step_id}_${step_name}_${job_id}.out"
    MICROBEATLAS_STDERR_LOG="$MICROBEATLAS_LOG_DIR/${step_id}_${step_name}_${job_id}.err"
    export MICROBEATLAS_STDOUT_LOG MICROBEATLAS_STDERR_LOG

    exec >"$MICROBEATLAS_STDOUT_LOG" 2>"$MICROBEATLAS_STDERR_LOG"

    echo "MicrobeAtlas step: ${step_id}_${step_name}"
    echo "Run label: $MICROBEATLAS_RUN_LABEL"
    echo "Project root: $PROJECT_ROOT"
    echo "Log directory: $MICROBEATLAS_LOG_DIR"
    echo "Stdout log: $MICROBEATLAS_STDOUT_LOG"
    echo "Stderr log: $MICROBEATLAS_STDERR_LOG"
    echo "Started: $(date -Iseconds)"
}

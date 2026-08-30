#!/bin/bash

# Template configuration.
#
# Copy this file to config/env.sh and edit config/env.sh for each local
# installation or analysis run:
#
#   cp config/example_env.sh config/env.sh
#   nano config/env.sh

# Project root.
export MICROBEATLAS_ROOT="/path/to/MicrobeAtlas-pipeline"

# Short name for this test run.
export MICROBEATLAS_RUN_LABEL="aquatic_1sample"

# Main data directories.
export MICROBEATLAS_INPUT_DIR="$MICROBEATLAS_ROOT/data/input"
export MICROBEATLAS_RAW_DATA_DIR="$MICROBEATLAS_ROOT/data/raw/microbeatlas"
export MICROBEATLAS_INTERMEDIATE_DIR="$MICROBEATLAS_ROOT/data/intermediate/$MICROBEATLAS_RUN_LABEL"
export MICROBEATLAS_RESULTS_DIR="$MICROBEATLAS_ROOT/results/$MICROBEATLAS_RUN_LABEL"
export MICROBEATLAS_LOG_DIR="$MICROBEATLAS_ROOT/log/$MICROBEATLAS_RUN_LABEL"

# Backward-compatible alias used internally by several scripts.
export MICROBEATLAS_DATA_DIR="$MICROBEATLAS_INTERMEDIATE_DIR"

# User ASV input table.
export MICROBEATLAS_OWN_ASV_TABLE_RAW="$MICROBEATLAS_INPUT_DIR/final_table_filtered.tsv"

# MicrobeAtlas sample selection.
export MICROBEATLAS_ENV_LABEL="aquatic"
export MICROBEATLAS_SAMPLES_PER_CC="1"

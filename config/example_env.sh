#!/bin/bash

# Template configuration.
#
# Copy this file to config/env.sh and edit config/env.sh for each local
# installation or analysis run:
#
#   cp config/example_env.sh config/env.sh
#   nano config/env.sh

# Project root.
export ATLASPIPE_ROOT="/path/to/AtlasPipe"

# Short name for this test run.
export ATLASPIPE_RUN_LABEL="aquatic_1sample"

# Main data directories.
export ATLASPIPE_INPUT_DIR="$ATLASPIPE_ROOT/data/input"
export MICROBEATLAS_RAW_DATA_DIR="$ATLASPIPE_ROOT/data/raw"
export ATLASPIPE_INTERMEDIATE_DIR="$ATLASPIPE_ROOT/data/intermediate/$ATLASPIPE_RUN_LABEL"
export ATLASPIPE_RESULTS_DIR="$ATLASPIPE_ROOT/results/$ATLASPIPE_RUN_LABEL"
export ATLASPIPE_LOG_DIR="$ATLASPIPE_ROOT/log/$ATLASPIPE_RUN_LABEL"

# Backward-compatible alias used internally by several scripts.
export ATLASPIPE_DATA_DIR="$ATLASPIPE_INTERMEDIATE_DIR"

# User amplicon input table.
export ATLASPIPE_OWN_AMPLICON_TABLE_RAW="$ATLASPIPE_INPUT_DIR/final_table_filtered.tsv"

# MicrobeAtlas BIOM abundance table.
#
# The rod2025companion table is the default because it preserves samples for
# all AMPLICON + robust community clusters represented in the MicrobeAtlas
# community-cluster mapping.
export MICROBEATLAS_BIOM_TABLE="$MICROBEATLAS_RAW_DATA_DIR/samples-otus.97.metag.minfilter.rod2025companion.biom"

# MicrobeAtlas sample selection.
export MICROBEATLAS_ENV_LABEL="aquatic"
export MICROBEATLAS_ENV_IDS="30"
export MICROBEATLAS_SAMPLES_PER_CC="1"

# Optional environment exclusions for step 02.
# Select and exclude environments only by numeric IDs from:
#   Rscript scripts/02_sample_selection.R --list-envs
# Selecting a parent ID also selects its sub-environments.
# export MICROBEATLAS_ENV_IDS="30,45"
# export MICROBEATLAS_EXCLUDE_ENV_IDS="32"

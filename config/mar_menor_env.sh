#!/bin/bash

# To use it:
#   source config/env.sh
#
# This profile uses the Mar Menor input table, the rod2025companion minfilter
# BIOM abundance table, and independent output directories.

# Project root.
export ATLASPIPE_ROOT="/storage/microbiome/home_members/pesai/AtlasPipe"

# Short name for this test run.
export ATLASPIPE_RUN_LABEL="mar_menor_rod2025"

# Main data directories 
export ATLASPIPE_INPUT_DIR="$ATLASPIPE_ROOT/data/input"
export MICROBEATLAS_RAW_DATA_DIR="$ATLASPIPE_ROOT/data/raw"
export ATLASPIPE_INTERMEDIATE_DIR="$ATLASPIPE_ROOT/data/intermediate/$ATLASPIPE_RUN_LABEL"
export ATLASPIPE_RESULTS_DIR="$ATLASPIPE_ROOT/results/$ATLASPIPE_RUN_LABEL"
export ATLASPIPE_LOG_DIR="$ATLASPIPE_ROOT/log/$ATLASPIPE_RUN_LABEL"

# Backward-compatible alias used internally by several scripts.
export ATLASPIPE_DATA_DIR="$ATLASPIPE_INTERMEDIATE_DIR"

# User amplicon input table.
export ATLASPIPE_OWN_AMPLICON_TABLE_RAW="$ATLASPIPE_INPUT_DIR/MarMenos/ASV_count_table_Mar_Menor_microbeatlas_input.tsv"

# MicrobeAtlas BIOM abundance table.
export MICROBEATLAS_BIOM_TABLE="$MICROBEATLAS_RAW_DATA_DIR/samples-otus.97.metag.minfilter.rod2025companion.biom"

# MicrobeAtlas sample selection.
export MICROBEATLAS_ENV_LABEL="aquatic-soil"
# Environment catalog IDs: 30 = aquatic, 51 = soil.
# Parent IDs include their sub-environments.
export MICROBEATLAS_ENV_IDS="30,51"
export MICROBEATLAS_SAMPLES_PER_CC="1"
# The rod2025companion BIOM should preserve all selected AMPLICON + robust
# community clusters, so missing BIOM community clusters should fail loudly.
export MICROBEATLAS_ALLOW_MISSING_BIOM_CC="FALSE"

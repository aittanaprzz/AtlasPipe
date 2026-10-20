#!/usr/bin/env Rscript

  # Prepare user-provided amplicon data for the AtlasPipe pipeline.
  #
  # This script reads a custom amplicon table, validates the required
  # fields, standardizes feature identifiers, and writes the input files
  # required by downstream pipeline steps.
  #
  # Expected input:
  #   A tab-separated table where:
  #     - column 1 is the original amplicon identifier
  #     - column 2 is the nucleotide sequence
  #     - remaining non-metadata columns are sample abundances
  #
  # Outputs:
  #   - Amplicon feature abundance table for step 07.
  #   - Amplicon feature FASTA file for step 05.
  #   - Feature ID mapping table linking original user IDs to internal ext_* IDs.
  #   - Summary file describing the generated inputs.
  #
  # Sample names are preserved exactly as provided in the input table.

library(data.table)

script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
  }
  getwd()
}

PROJECT_ROOT <- Sys.getenv(
  "ATLASPIPE_ROOT",
  normalizePath(file.path(script_dir(), ".."), mustWork = FALSE)
)
INPUT_DIR <- Sys.getenv("ATLASPIPE_INPUT_DIR", file.path(PROJECT_ROOT, "data", "input"))
DATA_DIR <- Sys.getenv("ATLASPIPE_DATA_DIR", file.path(PROJECT_ROOT, "data", "intermediate"))

INPUT_FILE <- Sys.getenv(
  "ATLASPIPE_OWN_AMPLICON_TABLE_RAW",
  file.path(INPUT_DIR, "final_table_filtered.tsv")
)
ABUND_DIR <- Sys.getenv(
  "ATLASPIPE_OWN_ABUND_DIR",
  file.path(DATA_DIR, "01_own_amplicon_inputs")
)
AMPLICON_MAPPING_DIR <- Sys.getenv(
  "ATLASPIPE_OWN_AMPLICON_MAPPING_DIR",
  file.path(DATA_DIR, "05_amplicon_mapping")
)
ABUND_FILE <- Sys.getenv(
  "ATLASPIPE_OWN_ABUND_FILE_OUT",
  file.path(ABUND_DIR, "amplicon_abundances.tsv")
)
AMPLICON_FASTA <- Sys.getenv(
  "ATLASPIPE_OWN_AMPLICON_FASTA",
  file.path(AMPLICON_MAPPING_DIR, "amplicons.fasta")
)
SUMMARY_FILE <- Sys.getenv(
  "ATLASPIPE_OWN_AMPLICON_INPUTS_SUMMARY",
  file.path(ABUND_DIR, "own_amplicon_inputs_summary.tsv")
)
AMPLICON_ID_MAP_FILE <- Sys.getenv(
  "ATLASPIPE_OWN_AMPLICON_ID_MAP",
  file.path(ABUND_DIR, "own_amplicon_id_map.tsv")
)
if (!file.exists(INPUT_FILE)) {
  stop("No existe la tabla de amplicones de entrada: ", INPUT_FILE)
}

dir.create(ABUND_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(AMPLICON_MAPPING_DIR, recursive = TRUE, showWarnings = FALSE)

x <- fread(INPUT_FILE, check.names = FALSE)

if (ncol(x) < 3) {
  stop(
    "La tabla de amplicones debe contener al menos tres columnas: ",
    "ID original, secuencia y una muestra de abundancia. Fichero: ",
    INPUT_FILE
  )
}

amplicon_id_col <- names(x)[1]
sequence_col <- names(x)[2]

metadata_cols <- c(
  amplicon_id_col,
  sequence_col,
  "bim1",
  "Kingdom",
  "Phylum",
  "Class",
  "Order",
  "Family",
  "Genus",
  "Species"
)
sample_cols <- setdiff(names(x), metadata_cols)

if (length(sample_cols) == 0) {
  stop("No se encontraron columnas de muestras en: ", INPUT_FILE)
}

amplicon_ids <- x[[amplicon_id_col]]
amplicon_sequences <- toupper(x[[sequence_col]])

if (any(is.na(amplicon_ids) | trimws(amplicon_ids) == "")) {
  bad_rows <- which(is.na(amplicon_ids) | trimws(amplicon_ids) == "")
  stop("Hay IDs de amplicones vacios en filas. Ejemplos: ", paste(head(bad_rows, 10), collapse = ", "))
}
amplicon_ids <- trimws(as.character(amplicon_ids))

if (anyDuplicated(amplicon_ids)) {
  duplicated_amplicons <- unique(amplicon_ids[duplicated(amplicon_ids)])
  stop("Hay IDs de amplicones duplicados. Ejemplos: ", paste(head(duplicated_amplicons, 10), collapse = ", "))
}

make_external_ids <- function(ids) {
  safe <- gsub("[^A-Za-z0-9_.-]+", "_", ids)
  safe <- gsub("^_+|_+$", "", safe)
  paste0("ext_", safe)
}

external_ids <- make_external_ids(amplicon_ids)
duplicated_external <- unique(external_ids[duplicated(external_ids)])
if (length(duplicated_external) > 0) {
  examples <- vapply(
    head(duplicated_external, 10),
    function(id) paste(amplicon_ids[external_ids == id], collapse = ", "),
    character(1)
  )
  stop(
    "El renombrado interno a ext_* generaria IDs duplicados. ",
    "Cambia los codigos originales para que sean distinguibles tras sustituir caracteres no seguros por '_'. ",
    "Ejemplos: ",
    paste(examples, collapse = " | ")
  )
}

if (any(is.na(amplicon_sequences) | amplicon_sequences == "")) {
  bad_amplicons <- amplicon_ids[is.na(amplicon_sequences) | amplicon_sequences == ""]
  stop("Hay amplicones sin secuencia. Ejemplos: ", paste(head(bad_amplicons, 10), collapse = ", "))
}

if (!all(grepl("^[ACGTUN-]+$", amplicon_sequences))) {
  bad_amplicons <- amplicon_ids[!grepl("^[ACGTUN-]+$", amplicon_sequences)]
  stop("Hay secuencias de amplicones con caracteres inesperados. Ejemplos: ", paste(head(bad_amplicons, 10), collapse = ", "))
}

if (anyDuplicated(sample_cols)) {
  duplicated_samples <- unique(sample_cols[duplicated(sample_cols)])
  stop("Hay nombres de muestra duplicados: ", paste(head(duplicated_samples, 10), collapse = ", "))
}

abund <- x[, c(amplicon_id_col, sample_cols), with = FALSE]
setnames(abund, amplicon_id_col, "AMPLICON_ID")
abund[, AMPLICON_ID := external_ids]

for (col in sample_cols) {
  abund[[col]] <- as.numeric(abund[[col]])
  if (any(is.na(abund[[col]]))) {
    stop("La muestra contiene valores no numericos tras convertir: ", col)
  }
}

amplicon_id_map <- data.table(
  original_amplicon_id = amplicon_ids,
  amplicon_id = external_ids,
  sequence = amplicon_sequences
)
fwrite(amplicon_id_map, AMPLICON_ID_MAP_FILE, sep = "\t")
fwrite(abund, ABUND_FILE, sep = "\t")

fasta_con <- file(AMPLICON_FASTA, open = "w")
on.exit(close(fasta_con), add = TRUE)
for (i in seq_along(external_ids)) {
  writeLines(c(paste0(">", external_ids[[i]]), amplicon_sequences[[i]]), fasta_con)
}

fwrite(
  data.table(
    parameter = c(
      "input_file",
      "abundance_file",
      "amplicon_fasta",
      "amplicon_id_map",
      "n_amplicons",
      "n_original_ids_renamed_to_ext",
      "n_samples",
      "total_abundance",
      "sample_names"
    ),
    value = c(
      INPUT_FILE,
      ABUND_FILE,
      AMPLICON_FASTA,
      AMPLICON_ID_MAP_FILE,
      nrow(abund),
      length(external_ids),
      length(sample_cols),
      sum(as.matrix(abund[, ..sample_cols])),
      paste(sample_cols, collapse = ",")
    )
  ),
  SUMMARY_FILE,
  sep = "\t"
)

cat("Tabla de abundancias de amplicones preparada:", ABUND_FILE, "\n")
cat("FASTA de amplicones preparado:", AMPLICON_FASTA, "\n")
cat("Mapa de IDs de amplicones:", AMPLICON_ID_MAP_FILE, "\n")
cat("Resumen:", SUMMARY_FILE, "\n")

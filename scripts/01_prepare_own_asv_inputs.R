#!/usr/bin/env Rscript
# Prepare own ASV inputs from a DADA-style table.
#
# Input expected:
#   asv_code, asv sequence, sample count columns, optional metadata/taxonomy.
# Outputs:
#   - ASV abundance table for step 07.
#   - ASV FASTA for step 05.
#
# Sample names are kept exactly as provided by the user in the input table.

library(data.table)

script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
  }
  getwd()
}

PROJECT_ROOT <- Sys.getenv(
  "MICROBEATLAS_ROOT",
  normalizePath(file.path(script_dir(), ".."), mustWork = FALSE)
)
INPUT_DIR <- Sys.getenv("MICROBEATLAS_INPUT_DIR", file.path(PROJECT_ROOT, "data", "input"))
DATA_DIR <- Sys.getenv("MICROBEATLAS_DATA_DIR", file.path(PROJECT_ROOT, "data", "intermediate"))

INPUT_FILE <- Sys.getenv(
  "MICROBEATLAS_OWN_ASV_TABLE_RAW",
  file.path(INPUT_DIR, "final_table_filtered.tsv")
)
ABUND_DIR <- Sys.getenv(
  "MICROBEATLAS_OWN_ABUND_DIR",
  file.path(DATA_DIR, "01_own_asv_inputs")
)
VSEARCH_DIR <- Sys.getenv(
  "MICROBEATLAS_VSEARCH_DIR",
  file.path(DATA_DIR, "05_vsearch")
)
ABUND_FILE <- Sys.getenv(
  "MICROBEATLAS_OWN_ABUND_FILE_OUT",
  file.path(ABUND_DIR, "ASV_abundances.tsv")
)
ASV_FASTA <- Sys.getenv(
  "MICROBEATLAS_ASV_FASTA",
  file.path(VSEARCH_DIR, "ASV.fasta")
)
SUMMARY_FILE <- Sys.getenv(
  "MICROBEATLAS_OWN_ASV_INPUTS_SUMMARY",
  file.path(ABUND_DIR, "own_asv_inputs_summary.tsv")
)
if (!file.exists(INPUT_FILE)) {
  stop("No existe la tabla ASV de entrada: ", INPUT_FILE)
}

dir.create(ABUND_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(VSEARCH_DIR, recursive = TRUE, showWarnings = FALSE)

x <- fread(INPUT_FILE, check.names = FALSE)

required <- c("asv_code", "asv")
missing_required <- setdiff(required, names(x))
if (length(missing_required) > 0) {
  stop(
    "La tabla ASV debe contener las columnas requeridas: ",
    paste(required, collapse = ", "),
    ". Faltan: ",
    paste(missing_required, collapse = ", ")
  )
}

metadata_cols <- c(
  "asv_code",
  "asv",
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

asv_ids <- x$asv_code
asv_sequences <- toupper(x$asv)

if (anyDuplicated(asv_ids)) {
  duplicated_asv <- unique(asv_ids[duplicated(asv_ids)])
  stop("Hay ASVs duplicados. Ejemplos: ", paste(head(duplicated_asv, 10), collapse = ", "))
}

if (!all(grepl("^(asv|otu)_", asv_ids))) {
  bad_asv <- asv_ids[!grepl("^(asv|otu)_", asv_ids)]
  stop("Hay IDs que no parecen ASV. Ejemplos: ", paste(head(bad_asv, 10), collapse = ", "))
}

if (any(is.na(asv_sequences) | asv_sequences == "")) {
  bad_asv <- asv_ids[is.na(asv_sequences) | asv_sequences == ""]
  stop("Hay ASVs sin secuencia. Ejemplos: ", paste(head(bad_asv, 10), collapse = ", "))
}

if (!all(grepl("^[ACGTUN-]+$", asv_sequences))) {
  bad_asv <- asv_ids[!grepl("^[ACGTUN-]+$", asv_sequences)]
  stop("Hay secuencias ASV con caracteres inesperados. Ejemplos: ", paste(head(bad_asv, 10), collapse = ", "))
}

if (anyDuplicated(sample_cols)) {
  duplicated_samples <- unique(sample_cols[duplicated(sample_cols)])
  stop("Hay nombres de muestra duplicados: ", paste(head(duplicated_samples, 10), collapse = ", "))
}

abund <- x[, c("asv_code", sample_cols), with = FALSE]
setnames(abund, "asv_code", "ASV")

for (col in sample_cols) {
  abund[[col]] <- as.numeric(abund[[col]])
  if (any(is.na(abund[[col]]))) {
    stop("La muestra contiene valores no numericos tras convertir: ", col)
  }
}

fwrite(abund, ABUND_FILE, sep = "\t")

fasta_con <- file(ASV_FASTA, open = "w")
on.exit(close(fasta_con), add = TRUE)
for (i in seq_along(asv_ids)) {
  writeLines(c(paste0(">", asv_ids[[i]]), asv_sequences[[i]]), fasta_con)
}

fwrite(
  data.table(
    parameter = c(
      "input_file",
      "abundance_file",
      "asv_fasta",
      "n_asvs",
      "n_samples",
      "total_abundance",
      "sample_names"
    ),
    value = c(
      INPUT_FILE,
      ABUND_FILE,
      ASV_FASTA,
      nrow(abund),
      length(sample_cols),
      sum(as.matrix(abund[, ..sample_cols])),
      paste(sample_cols, collapse = ",")
    )
  ),
  SUMMARY_FILE,
  sep = "\t"
)

cat("Tabla ASV preparada:", ABUND_FILE, "\n")
cat("FASTA ASV preparado:", ASV_FASTA, "\n")
cat("Resumen:", SUMMARY_FILE, "\n")

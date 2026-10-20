#!/usr/bin/env Rscript

# Builds a reference FASTA with the OTU sequences present in the filtered
# MicrobeAtlas samples.


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
IN_DIR <- Sys.getenv("ATLASPIPE_DATA_DIR", file.path(PROJECT_ROOT, "data", "intermediate"))
RAW_DATA_DIR <- Sys.getenv("MICROBEATLAS_RAW_DATA_DIR", file.path(PROJECT_ROOT, "data", "raw"))
ENV_LABEL <- Sys.getenv("MICROBEATLAS_ENV_LABEL", "aquatic")

selection_label_from_samples_per_cc <- function(x) {
  x <- tolower(trimws(as.character(x)))
  if (x %in% c("all", "todos")) {
    return("allsamplesperCC")
  }
  n <- suppressWarnings(as.integer(x))
  if (is.na(n) || n < 1L || as.character(n) != x) {
    stop("MICROBEATLAS_SAMPLES_PER_CC debe ser un entero positivo o 'all': ", x)
  }
  if (n == 1L) {
    return("1sampelperCC")
  }
  paste0(n, "samplesperCC")
}
SELECTION_LABEL <- Sys.getenv(
  "MICROBEATLAS_SELECTION_LABEL",
  selection_label_from_samples_per_cc(Sys.getenv("MICROBEATLAS_SAMPLES_PER_CC", "1"))
)
SUMMARY_SUFFIX <- if (SELECTION_LABEL == "1sampelperCC") {
  ""
} else {
  paste0("_", ENV_LABEL, "_", SELECTION_LABEL)
}

ID_OTU_FILE <- file.path(
  IN_DIR,
  "03_samples_subsets",
  paste0("samples-otus.97_selection_", ENV_LABEL, "_", SELECTION_LABEL, ".tsv")
)
ID_OTU_FILE <- Sys.getenv("MICROBEATLAS_ATLAS_TSV", ID_OTU_FILE)
SEQ_FILE <- Sys.getenv(
  "MICROBEATLAS_OTUS_ALLINFO",
  file.path(RAW_DATA_DIR, "otus.97.allinfo")
)

OUT_DIR <- Sys.getenv(
  "MICROBEATLAS_REFERENCE_FASTA_DIR",
  Sys.getenv("MICROBEATLAS_MULTIFASTA_DIR", file.path(IN_DIR, "04_reference_fasta"))
)
OUT_FILE <- file.path(
  OUT_DIR,
  paste0("MicrobeAtlas_selection_", ENV_LABEL, "_", SELECTION_LABEL, "_otus.97.fasta")
)
OUT_FILE <- Sys.getenv("MICROBEATLAS_OTU_FASTA", OUT_FILE)
SUMMARY_FILE <- Sys.getenv(
  "MICROBEATLAS_REFERENCE_FASTA_SUMMARY",
  Sys.getenv(
    "MICROBEATLAS_MULTIFASTA_SUMMARY",
    file.path(dirname(OUT_FILE), paste0("reference_fasta_summary", SUMMARY_SUFFIX, ".tsv"))
  )
)

if (!file.exists(ID_OTU_FILE)) {
  gz_id_otu_file <- paste0(ID_OTU_FILE, ".gz")
  if (file.exists(gz_id_otu_file)) {
    ID_OTU_FILE <- gz_id_otu_file
  } else {
    stop("No existe la tabla de OTUs seleccionados: ", ID_OTU_FILE)
  }
}

if (!file.exists(SEQ_FILE)) {
  stop("No existe el fichero con secuencias de OTUs: ", SEQ_FILE)
}

dir.create(dirname(OUT_FILE), recursive = TRUE, showWarnings = FALSE)


# El TSV generado por biom convert tiene una primera linea de comentario y la
# cabecera real en la segunda linea. Conservamos solo OTUs 97_* con abundancia
# positiva en al menos una muestra; la fila Unmapped no es una secuencia OTU.
# El BIOM rod2025companion prefija algunos OTUs con dominio (A97_, B97_, E97_).
# Mantenemos ese prefijo en el FASTA para no mezclar namespaces, pero buscamos
# las secuencias en otus.97.allinfo usando el ID base 97_*.
is_otu97_id <- function(x) {
  grepl("^[ABE]?97_", x)
}

otu97_base_id <- function(x) {
  sub("^[ABE](97_)", "\\1", x)
}

open_text <- function(path) {
  if (grepl("\\.gz$", path)) {
    return(gzfile(path, open = "rt"))
  }
  file(path, open = "r")
}

selected_lookup <- new.env(hash = TRUE, parent = emptyenv())
selected_by_base <- new.env(hash = TRUE, parent = emptyenv())

otu_con <- open_text(ID_OTU_FILE)
on.exit(close(otu_con), add = TRUE)

first_line <- readLines(otu_con, n = 1, warn = FALSE)
if (length(first_line) == 0) {
  stop("La tabla de OTUs esta vacia: ", ID_OTU_FILE)
}

header_line <- readLines(otu_con, n = 1, warn = FALSE)
if (length(header_line) == 0) {
  stop("La tabla de OTUs no contiene cabecera: ", ID_OTU_FILE)
}

header <- strsplit(header_line, "\t", fixed = TRUE)[[1]]
if (length(header) < 2) {
  stop("La tabla de OTUs no contiene columnas de muestras: ", ID_OTU_FILE)
}

n_samples <- length(header) - 1L
n_rows_total <- 0L
n_rows_97 <- 0L
n_rows_97_positive <- 0L
n_rows_97_zero <- 0L
n_unmapped_ignored <- 0L
n_other_ignored <- 0L
n_selected <- 0L
repeat {
  lines <- readLines(otu_con, n = 5000, warn = FALSE)
  if (length(lines) == 0) {
    break
  }

  for (line in lines) {
    fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(fields) < 2) {
      next
    }

    n_rows_total <- n_rows_total + 1L
    otu_id <- fields[1]
    if (is.na(otu_id) || otu_id == "") {
      n_other_ignored <- n_other_ignored + 1L
      next
    }

    if (otu_id == "Unmapped") {
      n_unmapped_ignored <- n_unmapped_ignored + 1L
      next
    }

    if (!is_otu97_id(otu_id)) {
      n_other_ignored <- n_other_ignored + 1L
      next
    }

    n_rows_97 <- n_rows_97 + 1L
    values <- suppressWarnings(as.numeric(fields[-1]))
    values[is.na(values)] <- 0

    if (any(values > 0)) {
      n_rows_97_positive <- n_rows_97_positive + 1L
      if (!exists(otu_id, envir = selected_lookup, inherits = FALSE)) {
        assign(otu_id, TRUE, envir = selected_lookup)
        base_id <- otu97_base_id(otu_id)
        selected_for_base <- if (exists(base_id, envir = selected_by_base, inherits = FALSE)) {
          get(base_id, envir = selected_by_base, inherits = FALSE)
        } else {
          character()
        }
        assign(base_id, c(selected_for_base, otu_id), envir = selected_by_base)
        n_selected <- n_selected + 1L
      }
    } else {
      n_rows_97_zero <- n_rows_97_zero + 1L
    }
  }
}

selected_otus <- ls(selected_lookup)

if (length(selected_otus) == 0) {
  stop("No se encontraron OTU IDs en: ", ID_OTU_FILE)
}

message("OTUs 97% presentes detectados: ", length(selected_otus))

# otus.97.allinfo no tiene cabecera. El primer campo contiene el identificador
# completo tipo 'MAPv3;90_x;96_y;97_z' y el septimo campo contiene la secuencia.
# Lo recorremos por bloques para evitar cargar todas las secuencias en memoria.
found_lookup <- new.env(hash = TRUE, parent = emptyenv())
n_allinfo_records <- 0L
n_malformed_allinfo_records <- 0L

seq_con <- file(SEQ_FILE, open = "r")
out_con <- file(OUT_FILE, open = "w")
on.exit(close(seq_con), add = TRUE)
on.exit(close(out_con), add = TRUE)

repeat {
  lines <- readLines(seq_con, n = 10000, warn = FALSE)
  if (length(lines) == 0) {
    break
  }

  fields <- strsplit(lines, "\t", fixed = TRUE)

  for (record in fields) {
    if (length(record) < 7) {
      n_malformed_allinfo_records <- n_malformed_allinfo_records + 1L
      next
    }

    n_allinfo_records <- n_allinfo_records + 1L
    otu_parts <- strsplit(record[1], ";", fixed = TRUE)[[1]]
    otu_id <- otu_parts[length(otu_parts)]
    sequence <- record[7]

    if (
      exists(otu_id, envir = selected_by_base, inherits = FALSE) &&
      !is.na(sequence) &&
      sequence != ""
    ) {
      for (selected_id in get(otu_id, envir = selected_by_base, inherits = FALSE)) {
        if (!exists(selected_id, envir = found_lookup, inherits = FALSE)) {
          writeLines(c(paste0(">", selected_id), sequence), out_con)
          assign(selected_id, TRUE, envir = found_lookup)
        }
      }
    }
  }
}

missing_otus <- selected_otus[
  !vapply(
    selected_otus,
    exists,
    logical(1),
    envir = found_lookup,
    inherits = FALSE
  )
]

summary_table <- data.frame(
  parameter = c(
    "input_tsv",
    "input_otus_allinfo",
    "output_fasta",
    "env_label",
    "samples_per_cc",
    "selection_label",
    "n_samples",
    "tsv_rows_total",
    "tsv_rows_97",
    "tsv_rows_97_positive",
    "tsv_rows_97_zero",
    "tsv_unmapped_ignored",
    "tsv_other_ignored",
    "selected_unique_otus",
    "allinfo_records_read",
    "allinfo_malformed_records",
    "sequences_found",
    "sequences_missing",
    "missing_otu_examples"
  ),
  value = c(
    ID_OTU_FILE,
    SEQ_FILE,
    OUT_FILE,
    ENV_LABEL,
    Sys.getenv("MICROBEATLAS_SAMPLES_PER_CC", "1"),
    SELECTION_LABEL,
    n_samples,
    n_rows_total,
    n_rows_97,
    n_rows_97_positive,
    n_rows_97_zero,
    n_unmapped_ignored,
    n_other_ignored,
    length(selected_otus),
    n_allinfo_records,
    n_malformed_allinfo_records,
    length(ls(found_lookup)),
    length(missing_otus),
    paste(head(missing_otus, 20), collapse = ",")
  )
)

write.table(
  summary_table,
  SUMMARY_FILE,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

if (length(ls(found_lookup)) == 0) {
  stop("No se encontro ninguna secuencia para los OTUs seleccionados.")
}

if (length(missing_otus) > 0) {
  warning(
    "Faltan secuencias para ",
    length(missing_otus),
    " OTUs seleccionados; se omitiran del FASTA de referencia. Ejemplos: ",
    paste(head(missing_otus, 10), collapse = ", ")
  )
}

message("Multifasta generado: ", OUT_FILE)
message("Numero de secuencias: ", length(ls(found_lookup)))
message("Reference FASTA summary: ", SUMMARY_FILE)
print(summary_table, row.names = FALSE)

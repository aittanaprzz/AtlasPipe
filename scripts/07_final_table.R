#!/usr/bin/env Rscript

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
DATA_DIR <- Sys.getenv("ATLASPIPE_DATA_DIR", file.path(PROJECT_ROOT, "data", "intermediate"))
RESULTS_DIR <- Sys.getenv("ATLASPIPE_RESULTS_DIR", file.path(PROJECT_ROOT, "results"))
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
VSEARCH_SUFFIX <- if (SELECTION_LABEL == "1sampelperCC") {
  ""
} else {
  paste0("_", ENV_LABEL, "_", SELECTION_LABEL)
}
SUMMARY_SUFFIX <- VSEARCH_SUFFIX
AMPLICON_MAPPING_DIR <- Sys.getenv(
  "ATLASPIPE_OWN_AMPLICON_MAPPING_DIR",
  file.path(DATA_DIR, "05_amplicon_mapping")
)
REFERENCE_FASTA_DIR <- Sys.getenv(
  "MICROBEATLAS_REFERENCE_FASTA_DIR",
  Sys.getenv("MICROBEATLAS_MULTIFASTA_DIR", file.path(DATA_DIR, "04_reference_fasta"))
)

# entradas
ATLAS_FILE <- Sys.getenv(
  "MICROBEATLAS_ATLAS_TSV",
  file.path(DATA_DIR, "03_samples_subsets", paste0("samples-otus.97_selection_", ENV_LABEL, "_", SELECTION_LABEL, ".tsv"))
)
ATLASPIPE_OWN_ABUND_FILE_OUT <- Sys.getenv(
  "ATLASPIPE_OWN_ABUND_FILE_OUT",
  file.path(DATA_DIR, "01_own_amplicon_inputs", "amplicon_abundances.tsv")
)
MAP_FILE <- Sys.getenv(
  "ATLASPIPE_OWN_AMPLICON_OTU97_MAP",
  file.path(AMPLICON_MAPPING_DIR, paste0("OTUs97_vsearch", VSEARCH_SUFFIX, ".blast6"))
)
CLUSTER_UC_FILE <- Sys.getenv(
  "ATLASPIPE_OWN_NO_HITS_UC",
  file.path(AMPLICON_MAPPING_DIR, paste0("amplicon_clusters97", VSEARCH_SUFFIX, ".uc"))
)
NO_HITS_FEATURE_MAP_FILE <- Sys.getenv(
  "ATLASPIPE_OWN_NO_HITS_FEATURE_MAP",
  file.path(AMPLICON_MAPPING_DIR, paste0("amplicon_cluster_map", VSEARCH_SUFFIX, ".tsv"))
)
AMPLICON_FASTA_FILE <- Sys.getenv(
  "ATLASPIPE_OWN_AMPLICON_FASTA",
  file.path(AMPLICON_MAPPING_DIR, "amplicons.fasta")
)
OUT_DIR <- Sys.getenv(
  "ATLASPIPE_FINAL_ABUND_DIR",
  file.path(RESULTS_DIR, "07_final_table")
)
ATLASPIPE_OWN_ALLOW_OTU_INPUT <- toupper(Sys.getenv("ATLASPIPE_OWN_ALLOW_OTU_INPUT", "FALSE")) %in%
  c("TRUE", "T", "1", "YES", "Y")

OUT_LONG_TOTAL <- file.path(OUT_DIR, "contextualized_abundance_table_long.tsv")
OUT_WIDE_TOTAL <- file.path(OUT_DIR, "contextualized_abundance_table.tsv")
OUT_AMPLICON_FEATURE_MAP <- file.path(OUT_DIR, "amplicon_feature_map.tsv")
OUT_AMPLICON_CONTEXT_MAP <- file.path(OUT_DIR, "amplicon_context_map.tsv")
OUT_AMPLICON_TABLE <- file.path(OUT_DIR, "amplicon_abundance_table.tsv")
OUT_AMPLICON_LONG <- file.path(OUT_DIR, "amplicon_abundance_table_long.tsv")
OUT_ATLAS_COVERAGE <- file.path(OUT_DIR, "microbeatlas_mapping_coverage.tsv")
OUT_QUALITY_SUMMARY <- file.path(OUT_DIR, "pipeline_quality_summary.tsv")
OUT_METRICS_SUMMARY <- file.path(OUT_DIR, "pipeline_metrics_summary.tsv")
OUT_REPRODUCIBILITY <- file.path(OUT_DIR, "pipeline_reproducibility.txt")
OUT_SAMPLE_METADATA <- file.path(OUT_DIR, "sample_metadata.tsv")
OUT_SAMPLE_DEPTHS <- file.path(OUT_DIR, "sample_depths.tsv")
ATLAS_SAMPLE_METADATA_FILE <- Sys.getenv(
  "MICROBEATLAS_SAMPLE_METADATA_FILE",
  file.path(
    DATA_DIR,
    paste0("02_filters_", ENV_LABEL, "_", SELECTION_LABEL),
    paste0("sample_metadata_", SELECTION_LABEL, ".tsv")
  )
)
ATLASPIPE_OWN_SAMPLE_METADATA_FILE <- Sys.getenv("ATLASPIPE_OWN_SAMPLE_METADATA", "")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

empty_long <- function() {
  data.table(
    SAMPLE = character(),
    FEATURE_ID = character(),
    FEATURE_TYPE = character(),
    abundance = numeric()
  )
}

empty_wide <- function() {
  data.table(FEATURE_ID = character())
}

empty_amplicon_feature_map <- function() {
  data.table(
    AMPLICON_ID = character(),
    FEATURE_ID = character(),
    FEATURE_TYPE = character(),
    assignment_source = character(),
    pident = numeric(),
    alnlen = numeric(),
    cluster_centroid = character(),
    cluster_size = integer()
  )
}

message_count <- function(label, value) {
  cat(label, value, "\n")
}

is_microbeatlas_otu97_id <- function(x) {
  grepl("^[ABE]?97_", x)
}

feature_type_from_id <- function(feature_id) {
  fifelse(
    is_microbeatlas_otu97_id(feature_id),
    "MicrobeAtlas_OTU97",
    fifelse(
      grepl("^cluster97_", feature_id) | grepl("^cluster_unmapped_", feature_id),
      "Own_cluster97",
      "unknown"
    )
  )
}

is_external_amplicon_id <- function(ids) {
  grepl("^ext_", ids)
}

wide_to_long <- function(x, id_name) {
  setnames(x, 1, id_name)
  long <- melt(
    x,
    id.vars = id_name,
    variable.name = "SAMPLE",
    value.name = "abundance"
  )
  setnames(long, id_name, "FEATURE_ID")
  long[, FEATURE_TYPE := feature_type_from_id(FEATURE_ID)]
  long[, abundance := as.numeric(abundance)]
  long[is.na(abundance), abundance := 0]
  long <- long[abundance > 0]
  long[, .(abundance = sum(abundance)), by = .(SAMPLE, FEATURE_ID, FEATURE_TYPE)]
}

numeric_sample_cols <- function(x, id_col) {
  sample_cols <- setdiff(names(x), id_col)
  x[, (sample_cols) := lapply(.SD, as.numeric), .SDcols = sample_cols]
  x[is.na(x)] <- 0
  sample_cols
}

aggregate_wide_by_otu <- function(x) {
  sample_cols <- numeric_sample_cols(x, "FEATURE_ID")
  x[, lapply(.SD, sum, na.rm = TRUE), by = FEATURE_ID, .SDcols = sample_cols]
}

write_positive_long_from_wide <- function(x, out_file) {
  long <- melt(
    x,
    id.vars = "FEATURE_ID",
    variable.name = "SAMPLE",
    value.name = "abundance"
  )
  long[, FEATURE_TYPE := feature_type_from_id(FEATURE_ID)]
  long[, abundance := as.numeric(abundance)]
  long <- long[!is.na(abundance) & abundance > 0]
  fwrite(long, out_file, sep = "\t")
  long
}

write_positive_long_from_biom_tsv <- function(path, out_file) {
  con <- file(path, open = "r")
  on.exit(close(con), add = TRUE)

  first_line <- readLines(con, n = 1, warn = FALSE)
  if (length(first_line) == 0) {
    stop("La tabla MicrobeAtlas esta vacia: ", path)
  }

  header_line <- readLines(con, n = 1, warn = FALSE)
  if (length(header_line) == 0) {
    stop("La tabla MicrobeAtlas no contiene cabecera: ", path)
  }

  header <- strsplit(header_line, "\t", fixed = TRUE)[[1]]
  samples <- header[-1]
  fwrite(empty_long(), out_file, sep = "\t")

  otu_seen <- new.env(hash = TRUE, parent = emptyenv())
  n_positive <- 0L

  repeat {
    lines <- readLines(con, n = 100, warn = FALSE)
    if (length(lines) == 0) {
      break
    }

    parsed <- strsplit(lines, "\t", fixed = TRUE)
    chunk <- rbindlist(lapply(parsed, function(fields) {
      if (length(fields) < 2) {
        return(empty_long())
      }

      otu <- fields[1]
      if (is.na(otu) || otu == "" || !is_microbeatlas_otu97_id(otu)) {
        return(empty_long())
      }

      values <- suppressWarnings(as.numeric(fields[-1]))
      values[is.na(values)] <- 0
      positive <- which(values > 0)
      if (length(positive) == 0) {
        return(empty_long())
      }

      data.table(
        SAMPLE = samples[positive],
        FEATURE_ID = otu,
        FEATURE_TYPE = "MicrobeAtlas_OTU97",
        abundance = values[positive]
      )
    }), use.names = TRUE)

    if (nrow(chunk) > 0) {
      fwrite(chunk, out_file, sep = "\t", append = TRUE, col.names = FALSE)
      n_positive <- n_positive + nrow(chunk)
      for (otu in unique(chunk$FEATURE_ID)) {
        assign(otu, TRUE, envir = otu_seen)
      }
    }
  }

  list(n_positive = n_positive, n_otus = length(ls(otu_seen)))
}

join_wide_tables <- function(tables) {
  tables <- tables[vapply(tables, nrow, integer(1)) > 0]
  if (length(tables) == 0) {
    return(empty_wide())
  }

  out <- rbindlist(tables, use.names = TRUE, fill = TRUE)
  out[is.na(out)] <- 0
  sample_cols <- numeric_sample_cols(out, "FEATURE_ID")
  out <- out[, lapply(.SD, sum, na.rm = TRUE), by = FEATURE_ID, .SDcols = sample_cols]
  out
}

wide_from_long <- function(x) {
  if (nrow(x) == 0) {
    return(empty_wide())
  }

  wide <- dcast(
    x,
    FEATURE_ID ~ SAMPLE,
    value.var = "abundance",
    fill = 0
  )
  wide
}

append_long_file <- function(x, out_file) {
  if (nrow(x) > 0) {
    fwrite(x, out_file, sep = "\t", append = TRUE, col.names = FALSE)
  }
}

open_biom_tsv_connection <- function(path) {
  if (file.exists(path)) {
    return(file(path, open = "r"))
  }

  gz_path <- paste0(path, ".gz")
  if (file.exists(gz_path)) {
    return(gzfile(gz_path, open = "r"))
  }

  stop("No existe la tabla MicrobeAtlas: ", path, " ni ", gz_path)
}

write_microbeatlas_coverage <- function(path, out_file) {
  con <- open_biom_tsv_connection(path)
  on.exit(close(con), add = TRUE)

  first_line <- readLines(con, n = 1, warn = FALSE)
  if (length(first_line) == 0) {
    stop("La tabla MicrobeAtlas esta vacia: ", path)
  }

  header_line <- readLines(con, n = 1, warn = FALSE)
  if (length(header_line) == 0) {
    stop("La tabla MicrobeAtlas no contiene cabecera: ", path)
  }

  header <- strsplit(header_line, "\t", fixed = TRUE)[[1]]
  samples <- header[-1]
  mapped <- numeric(length(samples))
  unmapped <- numeric(length(samples))
  mapped_rows <- 0L
  unmapped_rows <- 0L
  other_rows <- 0L

  repeat {
    lines <- readLines(con, n = 100, warn = FALSE)
    if (length(lines) == 0) {
      break
    }

    for (line in lines) {
      fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
      if (length(fields) < 2) {
        next
      }

      feature_id <- fields[1]
      values <- suppressWarnings(as.numeric(fields[-1]))
      if (length(values) < length(samples)) {
        values <- c(values, rep(NA_real_, length(samples) - length(values)))
      }
      values <- values[seq_along(samples)]
      values[is.na(values)] <- 0

      if (is_microbeatlas_otu97_id(feature_id)) {
        mapped <- mapped + values
        mapped_rows <- mapped_rows + 1L
      } else if (feature_id == "Unmapped") {
        unmapped <- unmapped + values
        unmapped_rows <- unmapped_rows + 1L
      } else {
        other_rows <- other_rows + 1L
      }
    }
  }

  total <- mapped + unmapped
  coverage <- data.table(
    SAMPLE = samples,
    mapped_otu97_abundance = mapped,
    unmapped_abundance = unmapped,
    total_mapped_plus_unmapped = total,
    mapped_otu97_fraction = fifelse(total > 0, mapped / total, NA_real_),
    unmapped_fraction = fifelse(total > 0, unmapped / total, NA_real_)
  )
  fwrite(coverage, out_file, sep = "\t")

  list(
    n_samples = length(samples),
    mapped_rows = mapped_rows,
    unmapped_rows = unmapped_rows,
    other_rows = other_rows,
    samples_with_unmapped = sum(unmapped > 0),
    zero_mapped_samples = sum(mapped == 0),
    unmapped_fraction_gt_0.9_samples = sum(coverage$unmapped_fraction > 0.9, na.rm = TRUE),
    unmapped_fraction_lt_0.1_samples = sum(coverage$unmapped_fraction < 0.1, na.rm = TRUE),
    total_mapped_abundance = sum(mapped),
    total_unmapped_abundance = sum(unmapped),
    total_mapped_plus_unmapped = sum(total),
    mean_mapped_fraction = mean(coverage$mapped_otu97_fraction, na.rm = TRUE),
    mean_unmapped_fraction = mean(coverage$unmapped_fraction, na.rm = TRUE),
    min_unmapped_fraction = min(coverage$unmapped_fraction, na.rm = TRUE),
    max_unmapped_fraction = max(coverage$unmapped_fraction, na.rm = TRUE)
  )
}

read_summary_value <- function(path, key, default = NA_character_) {
  if (!file.exists(path)) {
    return(default)
  }
  x <- fread(path, fill = TRUE, colClasses = "character")
  if (!all(c("parameter", "value") %in% names(x))) {
    setnames(x, 1:2, c("parameter", "value"))
  }
  value <- x[parameter == key]$value
  if (length(value) == 0) {
    return(default)
  }
  value[[1]]
}

write_sample_depths_from_wide <- function(path, out_file) {
  con <- file(path, open = "r")
  on.exit(close(con), add = TRUE)

  header_line <- readLines(con, n = 1, warn = FALSE)
  if (length(header_line) == 0) {
    stop("La tabla final ancha esta vacia: ", path)
  }

  header <- strsplit(header_line, "\t", fixed = TRUE)[[1]]
  if (length(header) < 3) {
    stop("La tabla final ancha no contiene columnas de muestras: ", path)
  }

  samples <- header[-c(1, 2)]
  depths <- numeric(length(samples))

  repeat {
    lines <- readLines(con, n = 50, warn = FALSE)
    if (length(lines) == 0) {
      break
    }

    for (line in lines) {
      fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
      values <- suppressWarnings(as.numeric(fields[-c(1, 2)]))
      if (length(values) < length(samples)) {
        values <- c(values, rep(NA_real_, length(samples) - length(values)))
      }
      values <- values[seq_along(samples)]
      values[is.na(values)] <- 0
      depths <- depths + values
    }
  }

  out <- data.table(
    SAMPLE = samples,
    final_feature_depth = depths,
    zero_final_feature_depth = depths == 0
  )
  fwrite(out, out_file, sep = "\t")
  out
}

write_final_sample_metadata <- function(sample_depths, atlas_metadata_file, own_wide, out_file, own_metadata_file = "") {
  out <- copy(sample_depths)
  out[, SAMPLE_SOURCE := fifelse(SAMPLE %in% setdiff(names(own_wide), "FEATURE_ID"), "Own", "MicrobeAtlas")]
  out[, `:=`(
    ComClustID = NA_character_,
    Env_cons = NA_character_,
    MICROBEATLAS_ENV_ROOT = NA_character_,
    SeqTech_cons = NA_character_,
    Robust_projection = NA_character_,
    Keyw_cons = NA_character_,
    Proj_cons = NA_character_
  )]

  if (file.exists(atlas_metadata_file)) {
    atlas_meta <- fread(atlas_metadata_file, colClasses = "character")
    if ("SampleID" %in% names(atlas_meta)) {
      setnames(atlas_meta, "SampleID", "SAMPLE")
    }
    atlas_cols <- intersect(
      c("SAMPLE", "ComClustID", "Env_cons", "MICROBEATLAS_ENV_ROOT", "SeqTech_cons", "Robust_projection", "Keyw_cons", "Proj_cons"),
      names(atlas_meta)
    )
    atlas_meta <- unique(atlas_meta[, ..atlas_cols], by = "SAMPLE")
    out <- merge(out, atlas_meta, by = "SAMPLE", all.x = TRUE, suffixes = c("", ".atlas"))
    for (col in setdiff(atlas_cols, "SAMPLE")) {
      atlas_col <- paste0(col, ".atlas")
      if (atlas_col %in% names(out)) {
        out[!is.na(get(atlas_col)) & get(atlas_col) != "", (col) := get(atlas_col)]
        out[, (atlas_col) := NULL]
      }
    }
  }

  if (own_metadata_file != "" && file.exists(own_metadata_file)) {
    own_meta <- fread(own_metadata_file, colClasses = "character")
    if ("SampleID" %in% names(own_meta) && !("SAMPLE" %in% names(own_meta))) {
      setnames(own_meta, "SampleID", "SAMPLE")
    }
    if (!("SAMPLE" %in% names(own_meta))) {
      warning("ATLASPIPE_OWN_SAMPLE_METADATA no contiene columna SAMPLE/SampleID: ", own_metadata_file)
    } else {
      own_meta <- unique(own_meta, by = "SAMPLE")
      own_extra_cols <- setdiff(names(own_meta), names(out))
      out <- merge(out, own_meta[, c("SAMPLE", own_extra_cols), with = FALSE], by = "SAMPLE", all.x = TRUE)
    }
  }

  setcolorder(out, c("SAMPLE", "SAMPLE_SOURCE", setdiff(names(out), c("SAMPLE", "SAMPLE_SOURCE"))))
  setorder(out, SAMPLE_SOURCE, SAMPLE)
  fwrite(out, out_file, sep = "\t")
  out
}

write_amplicon_outputs <- function(own_abund_file, amplicon_feature_map, out_wide, out_long, out_map) {
  own <- fread(own_abund_file, check.names = FALSE)
  if (ncol(own) < 2) {
    stop("La tabla propia preparada no contiene columnas de muestras: ", own_abund_file)
  }

  setnames(own, 1, "OWN_FEATURE_ID")
  own <- own[!is.na(OWN_FEATURE_ID) & OWN_FEATURE_ID != ""]
  duplicated_own <- unique(own$OWN_FEATURE_ID[duplicated(own$OWN_FEATURE_ID)])
  if (length(duplicated_own) > 0) {
    stop("La tabla propia preparada contiene IDs duplicados. Ejemplos: ", paste(head(duplicated_own, 10), collapse = ", "))
  }

  sample_cols <- setdiff(names(own), "OWN_FEATURE_ID")
  own[, (sample_cols) := lapply(.SD, as.numeric), .SDcols = sample_cols]
  for (col in sample_cols) {
    if (any(is.na(own[[col]]))) {
      stop("La tabla propia preparada contiene valores no numericos en la muestra: ", col)
    }
  }

  context <- copy(amplicon_feature_map)
  if (nrow(context) > 0) {
    setnames(
      context,
      old = c("AMPLICON_ID", "FEATURE_ID", "FEATURE_TYPE"),
      new = c("OWN_FEATURE_ID", "CONTEXTUALIZED_FEATURE_ID", "CONTEXTUALIZED_FEATURE_TYPE")
    )
  } else {
    context <- data.table(
      OWN_FEATURE_ID = character(),
      CONTEXTUALIZED_FEATURE_ID = character(),
      CONTEXTUALIZED_FEATURE_TYPE = character(),
      assignment_source = character(),
      pident = numeric(),
      alnlen = numeric(),
      cluster_centroid = character(),
      cluster_size = integer()
    )
  }

  context_cols <- c(
    "OWN_FEATURE_ID",
    "CONTEXTUALIZED_FEATURE_ID",
    "CONTEXTUALIZED_FEATURE_TYPE",
    "assignment_source",
    "pident",
    "alnlen",
    "cluster_centroid",
    "cluster_size"
  )
  missing_context_cols <- setdiff(context_cols, names(context))
  for (col in missing_context_cols) {
    context[, (col) := NA]
  }
  context <- unique(context[, ..context_cols], by = "OWN_FEATURE_ID")

  direct_otu_ids <- setdiff(own[is_microbeatlas_otu97_id(OWN_FEATURE_ID)]$OWN_FEATURE_ID, context$OWN_FEATURE_ID)
  if (length(direct_otu_ids) > 0) {
    context <- rbindlist(
      list(
        context,
        data.table(
          OWN_FEATURE_ID = direct_otu_ids,
          CONTEXTUALIZED_FEATURE_ID = direct_otu_ids,
          CONTEXTUALIZED_FEATURE_TYPE = "MicrobeAtlas_OTU97",
          assignment_source = "own_otu97_input",
          pident = NA_real_,
          alnlen = NA_real_,
          cluster_centroid = NA_character_,
          cluster_size = NA_integer_
        )
      ),
      use.names = TRUE,
      fill = TRUE
    )
  }

  missing_context <- setdiff(own$OWN_FEATURE_ID, context$OWN_FEATURE_ID)
  if (length(missing_context) > 0) {
    stop(
      "Hay features propias sin asignacion contextual. Ejecuta los pasos 05/06 antes del paso 07. Ejemplos: ",
      paste(head(missing_context, 10), collapse = ", ")
    )
  }

  context[, assignment_type := fifelse(
    assignment_source == "own_otu97_input",
    "own_otu97_input",
    fifelse(
      CONTEXTUALIZED_FEATURE_TYPE == "MicrobeAtlas_OTU97",
      "microbeatlas_otu97_hit",
      "own_nohit_novel_cluster97"
    )
  )]
  context <- context[
    ,
    .(
      OWN_FEATURE_ID,
      CONTEXTUALIZED_FEATURE_ID,
      CONTEXTUALIZED_FEATURE_TYPE,
      assignment_type,
      assignment_source,
      pident,
      alnlen,
      cluster_centroid,
      cluster_size
    )
  ]
  setorder(context, OWN_FEATURE_ID)
  fwrite(context, out_map, sep = "\t")

  out <- merge(context, own, by = "OWN_FEATURE_ID", all.x = TRUE)
  setcolorder(out, c(names(context), sample_cols))
  fwrite(out, out_wide, sep = "\t")

  long <- melt(
    out,
    id.vars = names(context),
    measure.vars = sample_cols,
    variable.name = "SAMPLE",
    value.name = "abundance"
  )
  long[, abundance := as.numeric(abundance)]
  long <- long[!is.na(abundance) & abundance > 0]
  fwrite(long, out_long, sep = "\t")

  list(
    n_own_features = nrow(out),
    n_own_samples = length(sample_cols),
    n_positive_records = nrow(long)
  )
}

write_pipeline_summaries <- function(
  quality_file,
  metrics_file,
  reproducibility_file,
  amplicon_identity_stats,
  own_tables,
  atlas_stats,
  wide_stats,
  coverage_stats,
  sample_depth_stats,
  own_wide,
  own_otu_long,
  amplicon_otu_long,
  amplicon_long_unmapped
) {
  amplicon_feature_map <- own_tables$amplicon_feature_map
  own_amplicon_total <- nrow(amplicon_feature_map)
  own_amplicon_mapped <- sum(amplicon_feature_map$FEATURE_TYPE == "MicrobeAtlas_OTU97")
  own_amplicon_nohit <- sum(amplicon_feature_map$FEATURE_TYPE == "Own_cluster97")
  own_samples <- setdiff(names(own_wide), "FEATURE_ID")

  final_own_novel_features <- uniqueN(
    amplicon_long_unmapped[FEATURE_TYPE == "Own_cluster97"]$FEATURE_ID
  )

  quality <- data.table(
    metric = c(
      "own_table_mode",
      "own_amplicon_reference_hit_rate",
      "own_amplicon_novel_fraction",
      "microbeatlas_non_otu97_rows_discarded",
      "microbeatlas_unmapped_rows",
      "microbeatlas_zero_mapped_samples",
      "microbeatlas_unmapped_fraction_gt_0.9_samples",
      "microbeatlas_unmapped_fraction_lt_0.1_samples",
      "microbeatlas_mean_mapped_fraction",
      "microbeatlas_mean_unmapped_fraction",
      "microbeatlas_min_unmapped_fraction",
      "microbeatlas_max_unmapped_fraction",
      "own_amplicon_table_fasta_same_ids",
      "own_amplicon_vsearch_queries_in_table",
      "own_amplicon_cluster_uc_ids_in_table",
      "final_unmapped_feature_rows",
      "final_zero_depth_samples",
      "final_min_sample_depth",
      "final_median_sample_depth",
      "final_max_sample_depth"
    ),
    value = c(
      own_tables$mode,
      ifelse(own_amplicon_total > 0, own_amplicon_mapped / own_amplicon_total, NA_real_),
      ifelse(own_amplicon_total > 0, own_amplicon_nohit / own_amplicon_total, NA_real_),
      wide_stats$atlas_unmapped_count,
      coverage_stats$unmapped_rows,
      coverage_stats$zero_mapped_samples,
      coverage_stats$unmapped_fraction_gt_0.9_samples,
      coverage_stats$unmapped_fraction_lt_0.1_samples,
      coverage_stats$mean_mapped_fraction,
      coverage_stats$mean_unmapped_fraction,
      coverage_stats$min_unmapped_fraction,
      coverage_stats$max_unmapped_fraction,
      amplicon_identity_stats$table_fasta_same_ids,
      amplicon_identity_stats$vsearch_queries_in_table,
      amplicon_identity_stats$cluster_uc_amplicons_in_table,
      0,
      sample_depth_stats$zero_depth_samples,
      sample_depth_stats$min_depth,
      sample_depth_stats$median_depth,
      sample_depth_stats$max_depth
    )
  )
  fwrite(quality, quality_file, sep = "\t")

  metrics <- data.table(
    metric = c(
      "own_samples",
      "own_amplicon_total",
      "own_amplicon_fasta_total",
      "own_amplicon_vsearch_query_total",
      "own_amplicon_cluster_uc_total",
      "own_amplicon_mapped_to_microbeatlas",
      "own_amplicon_nohit_total",
      "own_amplicon_novel_cluster_amplicons",
      "own_amplicon_novel_cluster_features",
      "microbeatlas_samples",
      "microbeatlas_otu97_features",
      "microbeatlas_samples_with_unmapped",
      "microbeatlas_total_mapped_abundance",
      "microbeatlas_total_unmapped_abundance",
      "microbeatlas_total_mapped_plus_unmapped",
      "final_samples",
      "final_features_total",
      "final_microbeatlas_otu97_features",
      "final_own_novel_cluster_features"
    ),
    value = c(
      length(own_samples),
      own_amplicon_total,
      amplicon_identity_stats$fasta_amplicons,
      amplicon_identity_stats$vsearch_query_amplicons,
      amplicon_identity_stats$cluster_uc_amplicons,
      own_amplicon_mapped,
      own_amplicon_nohit,
      own_amplicon_nohit,
      final_own_novel_features,
      coverage_stats$n_samples,
      atlas_stats$n_otus,
      coverage_stats$samples_with_unmapped,
      coverage_stats$total_mapped_abundance,
      coverage_stats$total_unmapped_abundance,
      coverage_stats$total_mapped_plus_unmapped,
      wide_stats$n_samples,
      wide_stats$atlas_otu_count + wide_stats$own_only_otu_count,
      atlas_stats$n_otus,
      final_own_novel_features
    )
  )
  fwrite(metrics, metrics_file, sep = "\t")

  filter_summary <- file.path(DATA_DIR, paste0("02_filters_", ENV_LABEL, "_", SELECTION_LABEL), "filter_summary.tsv")
  biom_summary <- file.path(
    DATA_DIR,
    "03_samples_subsets",
    paste0("biom_subset_summary_", ENV_LABEL, "_", SELECTION_LABEL, ".tsv")
  )
  reference_fasta_summary <- file.path(REFERENCE_FASTA_DIR, paste0("reference_fasta_summary", SUMMARY_SUFFIX, ".tsv"))
  own_summary <- file.path(DATA_DIR, "01_own_amplicon_inputs", "own_amplicon_inputs_summary.tsv")
  cluster_summary <- file.path(AMPLICON_MAPPING_DIR, paste0("cluster_summary97", VSEARCH_SUFFIX, ".tsv"))
  input_biom <- read_summary_value(biom_summary, "input_biom")

  reproducibility <- c(
    "AtlasPipe pipeline reproducibility",
    paste0("generated_at\t", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("project_root\t", PROJECT_ROOT),
    paste0("data_dir\t", DATA_DIR),
    paste0("raw_data_dir\t", dirname(input_biom)),
    paste0("input_biom\t", input_biom),
    "",
    "[configuration]",
    paste0("own_table_mode\t", own_tables$mode),
    paste0("env_id\t", read_summary_value(filter_summary, "env_id")),
    paste0("exclude_env_id\t", read_summary_value(filter_summary, "exclude_env_id")),
    paste0("selected_envs\t", read_summary_value(filter_summary, "selected_envs")),
    paste0("excluded_envs\t", read_summary_value(filter_summary, "excluded_envs")),
    paste0("seqtech\t", read_summary_value(filter_summary, "seqtech")),
    paste0("robust_only\t", read_summary_value(filter_summary, "robust_only")),
    paste0("seed\t", read_summary_value(filter_summary, "seed")),
    paste0("seed_source\t", read_summary_value(filter_summary, "seed_source")),
    paste0("samples_per_cc\t", read_summary_value(filter_summary, "samples_per_cc")),
    paste0("selection_label\t", SELECTION_LABEL),
    paste0("allow_own_otu_input\t", ATLASPIPE_OWN_ALLOW_OTU_INPUT),
    "",
    "[scripts_and_conda_envs]",
    paste0("step_01_script\t", file.path(script_dir(), "01_prepare_inputs.R")),
    "step_01_conda_env\tRenv",
    paste0("step_02_script\t", file.path(script_dir(), "02_sample_selection.R")),
    "step_02_conda_env\tRenv",
    paste0("step_03_script\t", file.path(PROJECT_ROOT, "slurm", "03_biom_subset.sbatch")),
    "step_03_conda_env\tBiom",
    paste0("step_04_script\t", file.path(script_dir(), "04_reference_fasta.R")),
    "step_04_conda_env\tRenv",
    paste0("step_05_script\t", file.path(PROJECT_ROOT, "slurm", "05_amplicon_mapping.sbatch")),
    "step_05_conda_env\tMicrobeAtlas",
    paste0("step_06_script\t", file.path(PROJECT_ROOT, "slurm", "06_novel_clustering.sbatch")),
    "step_06_conda_env\tMicrobeAtlas",
    paste0("step_07_script\t", file.path(script_dir(), "07_final_table.R")),
    "step_07_conda_env\tRenv",
    "",
    "[inputs]",
    paste0("own_amplicon_raw_table\t", read_summary_value(own_summary, "input_file", ATLASPIPE_OWN_ABUND_FILE_OUT)),
    paste0("own_abundance_table\t", ATLASPIPE_OWN_ABUND_FILE_OUT),
    paste0("atlas_tsv\t", ATLAS_FILE),
    paste0("amplicon_fasta\t", AMPLICON_FASTA_FILE),
    paste0("amplicon_otu_map\t", MAP_FILE),
    paste0("no_hits_uc\t", CLUSTER_UC_FILE),
    paste0("no_hits_feature_map\t", NO_HITS_FEATURE_MAP_FILE),
    paste0("otu_allinfo\t", read_summary_value(reference_fasta_summary, "input_otus_allinfo")),
    "",
    "[outputs]",
    paste0("final_wide_table\t", OUT_WIDE_TOTAL),
    paste0("final_long_table\t", OUT_LONG_TOTAL),
    paste0("amplicon_feature_map\t", OUT_AMPLICON_FEATURE_MAP),
    paste0("amplicon_context_map\t", OUT_AMPLICON_CONTEXT_MAP),
    paste0("amplicon_abundance_table\t", OUT_AMPLICON_TABLE),
    paste0("amplicon_abundance_table_long\t", OUT_AMPLICON_LONG),
    paste0("microbeatlas_mapping_coverage\t", OUT_ATLAS_COVERAGE),
    paste0("sample_metadata\t", OUT_SAMPLE_METADATA),
    paste0("sample_depths\t", OUT_SAMPLE_DEPTHS),
    paste0("pipeline_quality_summary\t", quality_file),
    paste0("pipeline_metrics_summary\t", metrics_file),
    paste0("pipeline_reproducibility\t", reproducibility_file),
    "",
    "[step_summaries_kept]",
    paste0("own_amplicon_inputs_summary\t", own_summary),
    paste0("filter_summary\t", filter_summary),
    paste0("biom_subset_summary\t", biom_summary),
    paste0("reference_fasta_summary\t", reference_fasta_summary),
    paste0("no_hits_clustering_summary\t", cluster_summary)
  )
  writeLines(reproducibility, reproducibility_file)
}

write_final_wide_from_biom_tsv <- function(path, own_wide, out_file) {
  con <- file(path, open = "r")
  on.exit(close(con), add = TRUE)

  first_line <- readLines(con, n = 1, warn = FALSE)
  if (length(first_line) == 0) {
    stop("La tabla MicrobeAtlas esta vacia: ", path)
  }

  header_line <- readLines(con, n = 1, warn = FALSE)
  if (length(header_line) == 0) {
    stop("La tabla MicrobeAtlas no contiene cabecera: ", path)
  }

  header <- strsplit(header_line, "\t", fixed = TRUE)[[1]]
  atlas_samples <- header[-1]
  own_samples <- setdiff(names(own_wide), "FEATURE_ID")
  sample_name_collisions <- intersect(atlas_samples, own_samples)
  if (length(sample_name_collisions) > 0) {
    stop(
      "Hay nombres de muestras propias que coinciden con muestras de MicrobeAtlas. ",
      "Renombra las muestras propias antes de construir la tabla final. Ejemplos: ",
      paste(head(sample_name_collisions, 10), collapse = ", ")
    )
  }

  own_rows <- split(own_wide, own_wide$FEATURE_ID)
  own_seen <- new.env(hash = TRUE, parent = emptyenv())
  zero_own <- rep("0", length(own_samples))
  zero_atlas <- rep("0", length(atlas_samples))

  out_con <- file(out_file, open = "w")
  on.exit(close(out_con), add = TRUE)
  writeLines(
    paste(c("FEATURE_ID", "FEATURE_TYPE", atlas_samples, own_samples), collapse = "\t"),
    out_con
  )

  atlas_otu_count <- 0L
  atlas_unmapped_count <- 0L

  repeat {
    lines <- readLines(con, n = 100, warn = FALSE)
    if (length(lines) == 0) {
      break
    }

    for (line in lines) {
      fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
      if (length(fields) < 2) {
        next
      }

      otu <- fields[1]
      if (is.na(otu) || otu == "" || !is_microbeatlas_otu97_id(otu)) {
        atlas_unmapped_count <- atlas_unmapped_count + 1L
        next
      }

      atlas_otu_count <- atlas_otu_count + 1L
      own_values <- zero_own
      if (otu %in% names(own_rows)) {
        own_values <- as.character(unlist(own_rows[[otu]][, ..own_samples], use.names = FALSE))
        assign(otu, TRUE, envir = own_seen)
      }

      writeLines(
        paste(c(otu, "MicrobeAtlas_OTU97", fields[-1], own_values), collapse = "\t"),
        out_con
      )
    }
  }

  missing_own_otus <- setdiff(names(own_rows), ls(own_seen))
  if (length(missing_own_otus) > 0) {
    for (otu in missing_own_otus) {
      own_values <- as.character(unlist(own_rows[[otu]][, ..own_samples], use.names = FALSE))
      writeLines(
        paste(c(otu, feature_type_from_id(otu), zero_atlas, own_values), collapse = "\t"),
        out_con
      )
    }
  }

  list(
    atlas_otu_count = atlas_otu_count,
    atlas_unmapped_count = atlas_unmapped_count,
    own_only_otu_count = length(missing_own_otus),
    n_samples = length(atlas_samples) + length(own_samples)
  )
}

read_amplicon_to_otu <- function(path) {
  x <- fread(
    path,
    select = 1:4,
    fill = TRUE,
    col.names = c("AMPLICON_ID", "FEATURE_ID", "pident", "alnlen")
  )
  x <- x[!is.na(AMPLICON_ID) & AMPLICON_ID != ""]
  x <- x[!is.na(FEATURE_ID) & FEATURE_ID != "" & FEATURE_ID != "*"]
  x <- x[is_microbeatlas_otu97_id(FEATURE_ID)]
  x[, FEATURE_TYPE := "MicrobeAtlas_OTU97"]
  x[, pident := as.numeric(pident)]
  x[, alnlen := as.numeric(alnlen)]
  setorder(x, AMPLICON_ID, -pident, -alnlen, FEATURE_ID)
  unique(x, by = "AMPLICON_ID")
}

read_fasta_ids <- function(path) {
  if (!file.exists(path)) {
    stop("No existe el FASTA de amplicones propios: ", path)
  }

  con <- file(path, open = "r")
  on.exit(close(con), add = TRUE)

  ids <- character()
  repeat {
    lines <- readLines(con, n = 10000, warn = FALSE)
    if (length(lines) == 0) {
      break
    }
    headers <- lines[grepl("^>", lines)]
    headers <- sub("^>", "", headers)
    headers <- sub("[ \t].*$", "", headers)
    ids <- c(ids, headers)
  }

  ids[!is.na(ids) & ids != ""]
}

read_vsearch_query_ids <- function(path) {
  if (!file.exists(path)) {
    stop("No existe la salida VSEARCH AMPLICON_ID -> OTU97: ", path)
  }

  x <- fread(path, select = 1, fill = TRUE, col.names = "AMPLICON_ID")
  unique(x[!is.na(AMPLICON_ID) & AMPLICON_ID != ""]$AMPLICON_ID)
}

read_uc_query_ids <- function(path) {
  if (!file.exists(path)) {
    return(character())
  }

  uc <- fread(path, header = FALSE, fill = TRUE)
  if (ncol(uc) < 9) {
    return(character())
  }

  unique(uc[V1 %in% c("S", "H") & !is.na(V9) & V9 != ""]$V9)
}

validate_own_amplicon_identity <- function(own_path, fasta_path, vsearch_path, cluster_uc_path) {
  own <- fread(own_path, select = 1, colClasses = "character", check.names = FALSE)
  setnames(own, 1, "AMPLICON_ID")
  table_ids <- own[!is.na(AMPLICON_ID) & AMPLICON_ID != ""]$AMPLICON_ID
  fasta_ids <- read_fasta_ids(fasta_path)
  vsearch_ids <- read_vsearch_query_ids(vsearch_path)
  uc_ids <- read_uc_query_ids(cluster_uc_path)

  duplicated_table <- unique(table_ids[duplicated(table_ids)])
  if (length(duplicated_table) > 0) {
    stop("La tabla propia contiene amplicones duplicados. Ejemplos: ", paste(head(duplicated_table, 10), collapse = ", "))
  }

  duplicated_fasta <- unique(fasta_ids[duplicated(fasta_ids)])
  if (length(duplicated_fasta) > 0) {
    stop("amplicons.fasta contiene IDs duplicados. Ejemplos: ", paste(head(duplicated_fasta, 10), collapse = ", "))
  }

  missing_in_fasta <- setdiff(table_ids, fasta_ids)
  extra_in_fasta <- setdiff(fasta_ids, table_ids)
  if (length(missing_in_fasta) > 0 || length(extra_in_fasta) > 0) {
    stop(
      "Incoherencia tabla AMPLICON_ID vs amplicons.fasta. ",
      "amplicones de tabla ausentes en FASTA: ", length(missing_in_fasta),
      if (length(missing_in_fasta) > 0) paste0(" (ejemplos: ", paste(head(missing_in_fasta, 10), collapse = ", "), ")") else "",
      ". amplicones extra en FASTA: ", length(extra_in_fasta),
      if (length(extra_in_fasta) > 0) paste0(" (ejemplos: ", paste(head(extra_in_fasta, 10), collapse = ", "), ")") else "",
      ". Esto sugiere cambio de nombres o entradas mezcladas."
    )
  }

  vsearch_not_in_table <- setdiff(vsearch_ids, table_ids)
  if (length(vsearch_not_in_table) > 0) {
    stop(
      "La salida VSEARCH contiene amplicones que no existen en la tabla propia. Ejemplos: ",
      paste(head(vsearch_not_in_table, 10), collapse = ", ")
    )
  }

  uc_not_in_table <- setdiff(uc_ids, table_ids)
  if (length(uc_not_in_table) > 0) {
    stop(
      "El clustering .uc contiene amplicones que no existen en la tabla propia. Ejemplos: ",
      paste(head(uc_not_in_table, 10), collapse = ", ")
    )
  }

  list(
    table_amplicons = length(table_ids),
    fasta_amplicons = length(fasta_ids),
    vsearch_query_amplicons = length(vsearch_ids),
    cluster_uc_amplicons = length(uc_ids),
    table_fasta_same_ids = identical(sort(table_ids), sort(fasta_ids)),
    vsearch_queries_in_table = length(vsearch_not_in_table) == 0,
    cluster_uc_amplicons_in_table = length(uc_not_in_table) == 0
  )
}

read_amplicon_cluster_feature_map <- function(path) {
  if (!file.exists(path)) {
    stop(
      "No existe el mapa formal AMPLICON_ID no-hit -> feature propia del paso 06: ", path,
      ". Ejecuta 06_novel_clustering.sbatch para generarlo."
    )
  }

  x <- fread(path, colClasses = "character")
  required_cols <- c(
    "AMPLICON_ID",
    "FEATURE_ID",
    "FEATURE_TYPE",
    "assignment_source",
    "cluster_centroid",
    "cluster_size"
  )
  missing_cols <- setdiff(required_cols, names(x))
  if (length(missing_cols) > 0) {
    stop(
      "El mapa formal de clusters propios no contiene columnas requeridas: ",
      paste(missing_cols, collapse = ", "),
      ". Fichero: ",
      path
    )
  }

  x <- x[!is.na(AMPLICON_ID) & AMPLICON_ID != ""]

  duplicated_amplicon <- unique(x$AMPLICON_ID[duplicated(x$AMPLICON_ID)])
  if (length(duplicated_amplicon) > 0) {
    stop(
      "El mapa formal de clusters propios contiene amplicones duplicados. Ejemplos: ",
      paste(head(duplicated_amplicon, 10), collapse = ", ")
    )
  }

  invalid_type <- x[!FEATURE_TYPE %in% c("Own_cluster97")]
  if (nrow(invalid_type) > 0) {
    stop(
      "El mapa formal de clusters propios contiene FEATURE_TYPE no validos. Ejemplos: ",
      paste(head(unique(invalid_type$FEATURE_TYPE), 10), collapse = ", ")
    )
  }

  x[, cluster_size := as.integer(cluster_size)]
  invalid_size <- x[is.na(cluster_size) | cluster_size < 1]
  if (nrow(invalid_size) > 0) {
    stop(
      "El mapa formal de clusters propios contiene cluster_size invalido. Ejemplos AMPLICON_ID: ",
      paste(head(invalid_size$AMPLICON_ID, 10), collapse = ", ")
    )
  }

  unique(
    x[, .(AMPLICON_ID, FEATURE_ID, FEATURE_TYPE, assignment_source, cluster_centroid, cluster_size)],
    by = "AMPLICON_ID"
  )
}

make_own_otu_long <- function(
  path,
  map_path,
  cluster_feature_map_path,
  allow_own_otu_input = FALSE
) {
  own <- fread(path, check.names = FALSE)
  setnames(own, 1, "FEATURE")
  ids <- own$FEATURE
  valid_ids <- ids[!is.na(ids) & ids != ""]
  is_otu_id <- is_microbeatlas_otu97_id(valid_ids)
  is_amplicon_id <- is_external_amplicon_id(valid_ids)
  dropped_ids <- valid_ids[!is_otu_id & !is_amplicon_id]

  if (any(is_otu_id) && !allow_own_otu_input) {
    stop(
      "La tabla de abundancias propia contiene IDs 97_*/A97_*/B97_*/E97_ y parece estar ya colapsada a OTU97. ",
      "Este paso espera amplicones para mapearlos con VSEARCH y clusterizar los no-hit. ",
      "Proporciona la tabla de amplicones original en ATLASPIPE_OWN_ABUND_FILE o, si quieres aceptar IDs OTU97 propios, ",
      "ejecuta con ATLASPIPE_OWN_ALLOW_OTU_INPUT=TRUE."
    )
  }

  if (length(dropped_ids) > 0) {
    warning(
      "Se descartan ",
      length(unique(dropped_ids)),
      " features propias no reconocidas como OTU97 ni feature local ext_*. Ejemplos: ",
      paste(head(unique(dropped_ids), 10), collapse = ", ")
    )
    own <- own[FEATURE %in% valid_ids[is_otu_id | is_amplicon_id]]
    valid_ids <- own$FEATURE
    is_otu_id <- is_microbeatlas_otu97_id(valid_ids)
    is_amplicon_id <- is_external_amplicon_id(valid_ids)
  }

  if (length(valid_ids) == 0) {
    stop("La tabla de abundancias propia no contiene features reconocibles: ", path)
  }

  own_otu_wide <- empty_wide()
  own_otu_long <- empty_long()
  if (any(is_otu_id)) {
    own_otu <- copy(own[is_microbeatlas_otu97_id(FEATURE)])
    setnames(own_otu, 1, "FEATURE_ID")
    own_otu_wide <- aggregate_wide_by_otu(own_otu)
    own_otu_long <- wide_to_long(copy(own_otu_wide), "OTU")
  }

  own_amplicon <- copy(own[is_external_amplicon_id(FEATURE)])
  if (nrow(own_amplicon) == 0) {
    return(list(
      own_otu_long = own_otu_long,
      amplicon_otu_long = empty_long(),
      amplicon_unmapped_long = empty_long(),
      own_otu_wide = own_otu_wide,
      amplicon_otu_wide = empty_wide(),
      amplicon_unmapped_wide = empty_wide(),
      amplicon_feature_map = empty_amplicon_feature_map(),
      mode = "otu"
    ))
  }

  setnames(own_amplicon, 1, "AMPLICON_ID")
  own_amplicon_ids <- unique(own_amplicon$AMPLICON_ID)
  own_long <- melt(
    own_amplicon,
    id.vars = "AMPLICON_ID",
    variable.name = "SAMPLE",
    value.name = "abundance"
  )
  own_long[, abundance := as.numeric(abundance)]
  own_long[is.na(abundance), abundance := 0]
  own_long <- own_long[abundance > 0]

  amplicon_to_otu <- read_amplicon_to_otu(map_path)
  cluster_map <- read_amplicon_cluster_feature_map(cluster_feature_map_path)

  amplicon_mapped <- merge(own_long, amplicon_to_otu, by = "AMPLICON_ID")
  amplicon_otu_long <- amplicon_mapped[
    ,
    .(abundance = sum(abundance)),
    by = .(SAMPLE, FEATURE_ID, FEATURE_TYPE)
  ]
  mapped_feature_map <- unique(amplicon_to_otu[
    AMPLICON_ID %in% own_amplicon_ids,
    .(
      AMPLICON_ID,
      FEATURE_ID,
      FEATURE_TYPE,
      assignment_source = "MicrobeAtlas_vsearch",
      pident,
      alnlen,
      cluster_centroid = NA_character_,
      cluster_size = NA_integer_
    )
  ])

  directly_mapped <- unique(amplicon_to_otu$AMPLICON_ID)
  own_unmapped_ids <- setdiff(own_amplicon_ids, directly_mapped)
  own_unmapped <- own_long[!AMPLICON_ID %in% directly_mapped]
  amplicon_clustered <- merge(own_unmapped, cluster_map, by = "AMPLICON_ID")
  missing_cluster <- setdiff(own_unmapped_ids, unique(cluster_map$AMPLICON_ID))
  if (length(missing_cluster) > 0) {
    stop(
      "Hay ",
      length(missing_cluster),
      " amplicones sin OTU97 que no aparecen en el mapa formal del paso 06. Ejemplos: ",
      paste(head(missing_cluster, 10), collapse = ", ")
    )
  }

  amplicon_unmapped_long <- amplicon_clustered[
    ,
    .(abundance = sum(abundance)),
    by = .(SAMPLE, FEATURE_ID, FEATURE_TYPE)
  ]
  unmapped_feature_map <- unique(cluster_map[
    AMPLICON_ID %in% own_unmapped_ids,
    .(
      AMPLICON_ID,
      FEATURE_ID,
      FEATURE_TYPE,
      assignment_source,
      pident = NA_real_,
      alnlen = NA_real_,
      cluster_centroid,
      cluster_size
    )
  ])

  list(
    own_otu_long = own_otu_long,
    amplicon_otu_long = amplicon_otu_long,
    amplicon_unmapped_long = amplicon_unmapped_long,
    own_otu_wide = own_otu_wide,
    amplicon_otu_wide = wide_from_long(amplicon_otu_long),
    amplicon_unmapped_wide = wide_from_long(amplicon_unmapped_long),
    amplicon_feature_map = rbindlist(
      list(mapped_feature_map, unmapped_feature_map),
      use.names = TRUE,
      fill = TRUE
    ),
    mode = if (any(is_otu_id)) "mixed_otu_amplicon" else "amplicon"
  )
}

amplicon_identity_stats <- validate_own_amplicon_identity(
  ATLASPIPE_OWN_ABUND_FILE_OUT,
  AMPLICON_FASTA_FILE,
  MAP_FILE,
  CLUSTER_UC_FILE
)

own_tables <- make_own_otu_long(
  ATLASPIPE_OWN_ABUND_FILE_OUT,
  MAP_FILE,
  NO_HITS_FEATURE_MAP_FILE,
  allow_own_otu_input = ATLASPIPE_OWN_ALLOW_OTU_INPUT
)
own_otu_long <- own_tables$own_otu_long
amplicon_otu_long <- own_tables$amplicon_otu_long
amplicon_long_unmapped <- own_tables$amplicon_unmapped_long
own_wide <- join_wide_tables(
  list(
    own_tables$own_otu_wide,
    own_tables$amplicon_otu_wide,
    own_tables$amplicon_unmapped_wide
  )
)

# tabla larga dispersa final
atlas_stats <- write_positive_long_from_biom_tsv(ATLAS_FILE, OUT_LONG_TOTAL)
append_long_file(own_otu_long, OUT_LONG_TOTAL)
append_long_file(amplicon_otu_long, OUT_LONG_TOTAL)
append_long_file(amplicon_long_unmapped, OUT_LONG_TOTAL)

# formato ancho
wide_stats <- write_final_wide_from_biom_tsv(ATLAS_FILE, own_wide, OUT_WIDE_TOTAL)

fwrite(own_tables$amplicon_feature_map, OUT_AMPLICON_FEATURE_MAP, sep = "\t")
own_feature_stats <- write_amplicon_outputs(
  ATLASPIPE_OWN_ABUND_FILE_OUT,
  own_tables$amplicon_feature_map,
  OUT_AMPLICON_TABLE,
  OUT_AMPLICON_LONG,
  OUT_AMPLICON_CONTEXT_MAP
)
coverage_stats <- write_microbeatlas_coverage(ATLAS_FILE, OUT_ATLAS_COVERAGE)
sample_depths <- write_sample_depths_from_wide(OUT_WIDE_TOTAL, OUT_SAMPLE_DEPTHS)
write_final_sample_metadata(
  sample_depths,
  ATLAS_SAMPLE_METADATA_FILE,
  own_wide,
  OUT_SAMPLE_METADATA,
  ATLASPIPE_OWN_SAMPLE_METADATA_FILE
)
nonzero_depths <- sample_depths[final_feature_depth > 0]$final_feature_depth
sample_depth_stats <- list(
  zero_depth_samples = sum(sample_depths$zero_final_feature_depth),
  min_depth = if (length(nonzero_depths) > 0) min(nonzero_depths) else 0,
  median_depth = if (length(nonzero_depths) > 0) median(nonzero_depths) else 0,
  max_depth = if (length(nonzero_depths) > 0) max(nonzero_depths) else 0
)
write_pipeline_summaries(
  OUT_QUALITY_SUMMARY,
  OUT_METRICS_SUMMARY,
  OUT_REPRODUCIBILITY,
  amplicon_identity_stats,
  own_tables,
  atlas_stats,
  wide_stats,
  coverage_stats,
  sample_depth_stats,
  own_wide,
  own_otu_long,
  amplicon_otu_long,
  amplicon_long_unmapped
)

cat("Guardado tabla final ancha:", OUT_WIDE_TOTAL, "\n")
cat("Guardado tabla final larga:", OUT_LONG_TOTAL, "\n")
cat("Mapa AMPLICON_ID -> feature:", OUT_AMPLICON_FEATURE_MAP, "\n")
cat("Mapa amplicon -> contexto:", OUT_AMPLICON_CONTEXT_MAP, "\n")
cat("Tabla de amplicones preservada:", OUT_AMPLICON_TABLE, "\n")
cat("Tabla de amplicones preservada larga:", OUT_AMPLICON_LONG, "\n")
cat("Cobertura MicrobeAtlas OTU97/Unmapped:", OUT_ATLAS_COVERAGE, "\n")
cat("Metadatos de muestras:", OUT_SAMPLE_METADATA, "\n")
cat("Profundidades finales por muestra:", OUT_SAMPLE_DEPTHS, "\n")
cat("Resumen de calidad del pipeline:", OUT_QUALITY_SUMMARY, "\n")
cat("Resumen de métricas del pipeline:", OUT_METRICS_SUMMARY, "\n")
cat("Reproducibilidad del pipeline:", OUT_REPRODUCIBILITY, "\n")
cat("Modo de tabla propia:", own_tables$mode, "\n")
message_count("Features MicrobeAtlas OTU97:", atlas_stats$n_otus)
message_count("Features MicrobeAtlas descartadas no OTU97:", wide_stats$atlas_unmapped_count)
message_count("Filas Unmapped MicrobeAtlas:", coverage_stats$unmapped_rows)
message_count("Muestras MicrobeAtlas con Unmapped:", coverage_stats$samples_with_unmapped)
message_count("Media unmapped_fraction MicrobeAtlas:", coverage_stats$mean_unmapped_fraction)
message_count("Max unmapped_fraction MicrobeAtlas:", coverage_stats$max_unmapped_fraction)
message_count("Features propias OTU97 directas:", uniqueN(own_otu_long$FEATURE_ID))
message_count("Features AMPLICON_ID mapeadas a OTU97:", uniqueN(amplicon_otu_long$FEATURE_ID))
message_count(
  "Features AMPLICON_ID no mapeadas en clusters novel:",
  uniqueN(amplicon_long_unmapped[FEATURE_TYPE == "Own_cluster97"]$FEATURE_ID)
)
message_count("Filas finales OTU/AMPLICON_ID:", wide_stats$atlas_otu_count + wide_stats$own_only_otu_count)
message_count("Muestras:", wide_stats$n_samples)
message_count("Features propias preservadas:", own_feature_stats$n_own_features)

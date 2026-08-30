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
  "MICROBEATLAS_ROOT",
  normalizePath(file.path(script_dir(), ".."), mustWork = FALSE)
)
DATA_DIR <- Sys.getenv("MICROBEATLAS_DATA_DIR", file.path(PROJECT_ROOT, "data", "intermediate"))
RESULTS_DIR <- Sys.getenv("MICROBEATLAS_RESULTS_DIR", file.path(PROJECT_ROOT, "results"))
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

# entradas
ATLAS_FILE <- Sys.getenv(
  "MICROBEATLAS_ATLAS_TSV",
  file.path(DATA_DIR, "03_samples_subsets", paste0("samples-otus.97_selection_", ENV_LABEL, "_", SELECTION_LABEL, ".tsv"))
)
OWN_ABUND_FILE <- Sys.getenv(
  "MICROBEATLAS_OWN_ABUND_FILE",
  file.path(DATA_DIR, "01_own_asv_inputs", "ASV_abundances.tsv")
)
MAP_FILE <- Sys.getenv(
  "MICROBEATLAS_ASV_OTU_MAP",
  file.path(DATA_DIR, "05_vsearch", paste0("OTUs97_vsearch", VSEARCH_SUFFIX, ".blast6"))
)
CLUSTER_UC_FILE <- Sys.getenv(
  "MICROBEATLAS_NO_HITS_UC",
  file.path(DATA_DIR, "05_vsearch", paste0("ASVs_no_hits_clustering97", VSEARCH_SUFFIX, ".uc"))
)
NO_HITS_FEATURE_MAP_FILE <- Sys.getenv(
  "MICROBEATLAS_NO_HITS_FEATURE_MAP",
  file.path(DATA_DIR, "05_vsearch", paste0("ASVs_no_hits_feature_map", VSEARCH_SUFFIX, ".tsv"))
)
ASV_FASTA_FILE <- Sys.getenv(
  "MICROBEATLAS_ASV_FASTA",
  file.path(DATA_DIR, "05_vsearch", "ASV.fasta")
)
OUT_DIR <- Sys.getenv(
  "MICROBEATLAS_FINAL_ABUND_DIR",
  file.path(RESULTS_DIR, "07_tabla_abundancias_profesor")
)
ALLOW_OWN_OTU_INPUT <- toupper(Sys.getenv("MICROBEATLAS_ALLOW_OWN_OTU_INPUT", "FALSE")) %in%
  c("TRUE", "T", "1", "YES", "Y")

OUT_LONG_TOTAL <- file.path(OUT_DIR, "tabla_abundancias_features_contextualizadas_largo.tsv")
OUT_WIDE_TOTAL <- file.path(OUT_DIR, "tabla_abundancias_features_contextualizadas.tsv")
OUT_ASV_FEATURE_MAP <- file.path(OUT_DIR, "asv_feature_map.tsv")
OUT_ATLAS_COVERAGE <- file.path(OUT_DIR, "microbeatlas_mapping_coverage.tsv")
OUT_QUALITY_SUMMARY <- file.path(OUT_DIR, "pipeline_quality_summary.tsv")
OUT_METRICS_SUMMARY <- file.path(OUT_DIR, "pipeline_metrics_summary.tsv")
OUT_REPRODUCIBILITY <- file.path(OUT_DIR, "pipeline_reproducibility.txt")

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

empty_asv_feature_map <- function() {
  data.table(
    ASV = character(),
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

feature_type_from_id <- function(feature_id) {
  fifelse(
    grepl("^97_", feature_id),
    "MicrobeAtlas_OTU97",
    fifelse(
      grepl("^ASVcluster97_", feature_id),
      "Own_ASVcluster97",
      fifelse(
        grepl("^ASVsingleton_", feature_id) | grepl("^ASVunmapped_", feature_id),
        "Own_ASVsingleton",
        "unknown"
      )
    )
  )
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
      if (is.na(otu) || otu == "" || !grepl("^97_", otu)) {
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

      if (grepl("^97_", feature_id)) {
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

write_pipeline_summaries <- function(
  quality_file,
  metrics_file,
  reproducibility_file,
  asv_identity_stats,
  own_tables,
  atlas_stats,
  wide_stats,
  coverage_stats,
  own_wide,
  own_otu_long,
  asv_otu_long,
  asv_long_unmapped
) {
  asv_feature_map <- own_tables$asv_feature_map
  own_asv_total <- nrow(asv_feature_map)
  own_asv_mapped <- sum(asv_feature_map$FEATURE_TYPE == "MicrobeAtlas_OTU97")
  own_asv_clustered <- sum(asv_feature_map$FEATURE_TYPE == "Own_ASVcluster97")
  own_asv_singleton <- sum(asv_feature_map$FEATURE_TYPE == "Own_ASVsingleton")
  own_asv_nohit <- own_asv_clustered + own_asv_singleton
  own_samples <- setdiff(names(own_wide), "FEATURE_ID")

  final_own_cluster_features <- uniqueN(
    asv_long_unmapped[FEATURE_TYPE == "Own_ASVcluster97"]$FEATURE_ID
  )
  final_own_singleton_features <- uniqueN(
    asv_long_unmapped[FEATURE_TYPE == "Own_ASVsingleton"]$FEATURE_ID
  )

  quality <- data.table(
    metric = c(
      "own_table_mode",
      "own_asv_reference_hit_rate",
      "own_asv_novel_fraction",
      "microbeatlas_non_otu97_rows_discarded",
      "microbeatlas_unmapped_rows",
      "microbeatlas_zero_mapped_samples",
      "microbeatlas_unmapped_fraction_gt_0.9_samples",
      "microbeatlas_unmapped_fraction_lt_0.1_samples",
      "microbeatlas_mean_mapped_fraction",
      "microbeatlas_mean_unmapped_fraction",
      "microbeatlas_min_unmapped_fraction",
      "microbeatlas_max_unmapped_fraction",
      "own_asv_table_fasta_same_ids",
      "own_asv_vsearch_queries_in_table",
      "own_asv_cluster_uc_ids_in_table",
      "final_unmapped_feature_rows"
    ),
    value = c(
      own_tables$mode,
      ifelse(own_asv_total > 0, own_asv_mapped / own_asv_total, NA_real_),
      ifelse(own_asv_total > 0, own_asv_nohit / own_asv_total, NA_real_),
      wide_stats$atlas_unmapped_count,
      coverage_stats$unmapped_rows,
      coverage_stats$zero_mapped_samples,
      coverage_stats$unmapped_fraction_gt_0.9_samples,
      coverage_stats$unmapped_fraction_lt_0.1_samples,
      coverage_stats$mean_mapped_fraction,
      coverage_stats$mean_unmapped_fraction,
      coverage_stats$min_unmapped_fraction,
      coverage_stats$max_unmapped_fraction,
      asv_identity_stats$table_fasta_same_ids,
      asv_identity_stats$vsearch_queries_in_table,
      asv_identity_stats$cluster_uc_asvs_in_table,
      0
    )
  )
  fwrite(quality, quality_file, sep = "\t")

  metrics <- data.table(
    metric = c(
      "own_samples",
      "own_asv_total",
      "own_asv_fasta_total",
      "own_asv_vsearch_query_total",
      "own_asv_cluster_uc_total",
      "own_asv_mapped_to_microbeatlas",
      "own_asv_nohit_total",
      "own_asv_clustered_asvs",
      "own_asv_singleton_asvs",
      "own_asv_cluster_features",
      "own_asv_singleton_features",
      "microbeatlas_samples",
      "microbeatlas_otu97_features",
      "microbeatlas_samples_with_unmapped",
      "microbeatlas_total_mapped_abundance",
      "microbeatlas_total_unmapped_abundance",
      "microbeatlas_total_mapped_plus_unmapped",
      "final_samples",
      "final_features_total",
      "final_microbeatlas_otu97_features",
      "final_own_cluster_features",
      "final_own_singleton_features"
    ),
    value = c(
      length(own_samples),
      own_asv_total,
      asv_identity_stats$fasta_asvs,
      asv_identity_stats$vsearch_query_asvs,
      asv_identity_stats$cluster_uc_asvs,
      own_asv_mapped,
      own_asv_nohit,
      own_asv_clustered,
      own_asv_singleton,
      final_own_cluster_features,
      final_own_singleton_features,
      coverage_stats$n_samples,
      atlas_stats$n_otus,
      coverage_stats$samples_with_unmapped,
      coverage_stats$total_mapped_abundance,
      coverage_stats$total_unmapped_abundance,
      coverage_stats$total_mapped_plus_unmapped,
      wide_stats$n_samples,
      wide_stats$atlas_otu_count + wide_stats$own_only_otu_count,
      atlas_stats$n_otus,
      final_own_cluster_features,
      final_own_singleton_features
    )
  )
  fwrite(metrics, metrics_file, sep = "\t")

  filter_summary <- file.path(DATA_DIR, paste0("02_filters_", ENV_LABEL, "_", SELECTION_LABEL), "filter_summary.tsv")
  biom_summary <- file.path(
    DATA_DIR,
    "03_samples_subsets",
    paste0("biom_subset_summary_", ENV_LABEL, "_", SELECTION_LABEL, ".tsv")
  )
  multifasta_summary <- file.path(DATA_DIR, "04_multifasta", paste0("multifasta_summary", SUMMARY_SUFFIX, ".tsv"))
  own_summary <- file.path(DATA_DIR, "01_own_asv_inputs", "own_asv_inputs_summary.tsv")
  cluster_summary <- file.path(DATA_DIR, "05_vsearch", paste0("ASVs_no_hits_clustering97_summary", VSEARCH_SUFFIX, ".tsv"))
  input_biom <- read_summary_value(biom_summary, "input_biom")

  reproducibility <- c(
    "MicrobeAtlas pipeline reproducibility",
    paste0("generated_at\t", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("project_root\t", PROJECT_ROOT),
    paste0("data_dir\t", DATA_DIR),
    paste0("raw_data_dir\t", dirname(input_biom)),
    paste0("input_biom\t", input_biom),
    "",
    "[configuration]",
    paste0("own_table_mode\t", own_tables$mode),
    paste0("env_argument\t", read_summary_value(filter_summary, "env_argument")),
    paste0("env_prefix\t", read_summary_value(filter_summary, "env_prefix")),
    paste0("seqtech\t", read_summary_value(filter_summary, "seqtech")),
    paste0("robust_only\t", read_summary_value(filter_summary, "robust_only")),
    paste0("seed\t", read_summary_value(filter_summary, "seed")),
    paste0("seed_source\t", read_summary_value(filter_summary, "seed_source")),
    paste0("samples_per_cc\t", read_summary_value(filter_summary, "samples_per_cc")),
    paste0("selection_label\t", SELECTION_LABEL),
    paste0("allow_own_otu_input\t", ALLOW_OWN_OTU_INPUT),
    "",
    "[scripts_and_conda_envs]",
    paste0("step_01_script\t", file.path(script_dir(), "01_prepare_own_asv_inputs.R")),
    "step_01_conda_env\tRenv",
    paste0("step_02_script\t", file.path(script_dir(), "02_filter_samples_1sampelperCC.R")),
    "step_02_conda_env\tRenv",
    paste0("step_03_script\t", file.path(PROJECT_ROOT, "slurm", "03_extract_biome.sbatch")),
    "step_03_conda_env\tBiom",
    paste0("step_04_script\t", file.path(script_dir(), "04_multifasta.R")),
    "step_04_conda_env\tRenv",
    paste0("step_05_script\t", file.path(PROJECT_ROOT, "slurm", "05_vsearch.sbatch")),
    "step_05_conda_env\tMicrobeAtlas",
    paste0("step_06_script\t", file.path(PROJECT_ROOT, "slurm", "06_vsearch_clustering.sbatch")),
    "step_06_conda_env\tMicrobeAtlas",
    paste0("step_07_script\t", file.path(script_dir(), "07_tabla_abundancias_profesor.R")),
    "step_07_conda_env\tRenv",
    "",
    "[inputs]",
    paste0("own_asv_raw_table\t", read_summary_value(own_summary, "input_file", OWN_ABUND_FILE)),
    paste0("own_abundance_table\t", OWN_ABUND_FILE),
    paste0("atlas_tsv\t", ATLAS_FILE),
    paste0("asv_fasta\t", ASV_FASTA_FILE),
    paste0("asv_otu_map\t", MAP_FILE),
    paste0("no_hits_uc\t", CLUSTER_UC_FILE),
    paste0("no_hits_feature_map\t", NO_HITS_FEATURE_MAP_FILE),
    paste0("otu_allinfo\t", read_summary_value(multifasta_summary, "input_otus_allinfo")),
    "",
    "[outputs]",
    paste0("final_wide_table\t", OUT_WIDE_TOTAL),
    paste0("final_long_table\t", OUT_LONG_TOTAL),
    paste0("asv_feature_map\t", OUT_ASV_FEATURE_MAP),
    paste0("microbeatlas_mapping_coverage\t", OUT_ATLAS_COVERAGE),
    paste0("pipeline_quality_summary\t", quality_file),
    paste0("pipeline_metrics_summary\t", metrics_file),
    paste0("pipeline_reproducibility\t", reproducibility_file),
    "",
    "[step_summaries_kept]",
    paste0("own_asv_inputs_summary\t", own_summary),
    paste0("filter_summary\t", filter_summary),
    paste0("biom_subset_summary\t", biom_summary),
    paste0("multifasta_summary\t", multifasta_summary),
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
      if (is.na(otu) || otu == "" || !grepl("^97_", otu)) {
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

read_asv_to_otu <- function(path) {
  x <- fread(
    path,
    select = 1:4,
    fill = TRUE,
    col.names = c("ASV", "FEATURE_ID", "pident", "alnlen")
  )
  x <- x[!is.na(ASV) & ASV != ""]
  x <- x[!is.na(FEATURE_ID) & FEATURE_ID != "" & FEATURE_ID != "*"]
  x <- x[grepl("^97_", FEATURE_ID)]
  x[, FEATURE_TYPE := "MicrobeAtlas_OTU97"]
  x[, pident := as.numeric(pident)]
  x[, alnlen := as.numeric(alnlen)]
  setorder(x, ASV, -pident, -alnlen, FEATURE_ID)
  unique(x, by = "ASV")
}

read_fasta_ids <- function(path) {
  if (!file.exists(path)) {
    stop("No existe el FASTA de ASVs propios: ", path)
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
    stop("No existe la salida VSEARCH ASV -> OTU97: ", path)
  }

  x <- fread(path, select = 1, fill = TRUE, col.names = "ASV")
  unique(x[!is.na(ASV) & ASV != ""]$ASV)
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

validate_own_asv_identity <- function(own_path, fasta_path, vsearch_path, cluster_uc_path) {
  own <- fread(own_path, select = 1, colClasses = "character", check.names = FALSE)
  setnames(own, 1, "ASV")
  table_ids <- own[!is.na(ASV) & ASV != ""]$ASV
  fasta_ids <- read_fasta_ids(fasta_path)
  vsearch_ids <- read_vsearch_query_ids(vsearch_path)
  uc_ids <- read_uc_query_ids(cluster_uc_path)

  duplicated_table <- unique(table_ids[duplicated(table_ids)])
  if (length(duplicated_table) > 0) {
    stop("La tabla propia contiene ASVs duplicados. Ejemplos: ", paste(head(duplicated_table, 10), collapse = ", "))
  }

  duplicated_fasta <- unique(fasta_ids[duplicated(fasta_ids)])
  if (length(duplicated_fasta) > 0) {
    stop("ASV.fasta contiene IDs duplicados. Ejemplos: ", paste(head(duplicated_fasta, 10), collapse = ", "))
  }

  missing_in_fasta <- setdiff(table_ids, fasta_ids)
  extra_in_fasta <- setdiff(fasta_ids, table_ids)
  if (length(missing_in_fasta) > 0 || length(extra_in_fasta) > 0) {
    stop(
      "Incoherencia tabla ASV vs ASV.fasta. ",
      "ASVs de tabla ausentes en FASTA: ", length(missing_in_fasta),
      if (length(missing_in_fasta) > 0) paste0(" (ejemplos: ", paste(head(missing_in_fasta, 10), collapse = ", "), ")") else "",
      ". ASVs extra en FASTA: ", length(extra_in_fasta),
      if (length(extra_in_fasta) > 0) paste0(" (ejemplos: ", paste(head(extra_in_fasta, 10), collapse = ", "), ")") else "",
      ". Esto sugiere cambio de nombres o entradas mezcladas."
    )
  }

  vsearch_not_in_table <- setdiff(vsearch_ids, table_ids)
  if (length(vsearch_not_in_table) > 0) {
    stop(
      "La salida VSEARCH contiene ASVs que no existen en la tabla propia. Ejemplos: ",
      paste(head(vsearch_not_in_table, 10), collapse = ", ")
    )
  }

  uc_not_in_table <- setdiff(uc_ids, table_ids)
  if (length(uc_not_in_table) > 0) {
    stop(
      "El clustering .uc contiene ASVs que no existen en la tabla propia. Ejemplos: ",
      paste(head(uc_not_in_table, 10), collapse = ", ")
    )
  }

  list(
    table_asvs = length(table_ids),
    fasta_asvs = length(fasta_ids),
    vsearch_query_asvs = length(vsearch_ids),
    cluster_uc_asvs = length(uc_ids),
    table_fasta_same_ids = identical(sort(table_ids), sort(fasta_ids)),
    vsearch_queries_in_table = length(vsearch_not_in_table) == 0,
    cluster_uc_asvs_in_table = length(uc_not_in_table) == 0
  )
}

read_asv_cluster_feature_map <- function(path) {
  if (!file.exists(path)) {
    stop(
      "No existe el mapa formal ASV no-hit -> feature propia del paso 06: ", path,
      ". Ejecuta 06_vsearch_clustering.sbatch para generarlo."
    )
  }

  x <- fread(path, colClasses = "character")
  required_cols <- c(
    "ASV",
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

  x <- x[!is.na(ASV) & ASV != ""]
  duplicated_asv <- unique(x$ASV[duplicated(x$ASV)])
  if (length(duplicated_asv) > 0) {
    stop(
      "El mapa formal de clusters propios contiene ASVs duplicados. Ejemplos: ",
      paste(head(duplicated_asv, 10), collapse = ", ")
    )
  }

  invalid_type <- x[!FEATURE_TYPE %in% c("Own_ASVcluster97", "Own_ASVsingleton")]
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
      "El mapa formal de clusters propios contiene cluster_size invalido. Ejemplos ASV: ",
      paste(head(invalid_size$ASV, 10), collapse = ", ")
    )
  }

  inconsistent_cluster <- x[FEATURE_TYPE == "Own_ASVcluster97" & cluster_size < 2]
  inconsistent_singleton <- x[FEATURE_TYPE == "Own_ASVsingleton" & cluster_size != 1]
  if (nrow(inconsistent_cluster) > 0 || nrow(inconsistent_singleton) > 0) {
    stop(
      "El mapa formal de clusters propios no respeta la definicion cluster_size >= 2 vs singleton == 1. ",
      "Clusters invalidos: ", nrow(inconsistent_cluster),
      ". Singletons invalidos: ", nrow(inconsistent_singleton)
    )
  }

  unique(
    x[, .(ASV, FEATURE_ID, FEATURE_TYPE, assignment_source, cluster_centroid, cluster_size)],
    by = "ASV"
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
  is_otu_id <- grepl("^97_", valid_ids)
  is_asv_id <- grepl("^(asv|otu)_", valid_ids)
  dropped_ids <- valid_ids[!is_otu_id & !is_asv_id]

  if (any(is_otu_id) && !allow_own_otu_input) {
    stop(
      "La tabla de abundancias propia contiene IDs 97_* y parece estar ya colapsada a OTU97. ",
      "Este paso espera ASVs de DADA para mapearlos con VSEARCH y clusterizar los no-hit. ",
      "Proporciona la tabla ASV original en MICROBEATLAS_OWN_ABUND_FILE o, solo para compatibilidad, ",
      "ejecuta con MICROBEATLAS_ALLOW_OWN_OTU_INPUT=TRUE."
    )
  }

  if (length(dropped_ids) > 0) {
    warning(
      "Se descartan ",
      length(unique(dropped_ids)),
      " features propias no reconocidas como OTU97 ni ASV/OTU local. Ejemplos: ",
      paste(head(unique(dropped_ids), 10), collapse = ", ")
    )
    own <- own[FEATURE %in% valid_ids[is_otu_id | is_asv_id]]
    valid_ids <- own$FEATURE
    is_otu_id <- grepl("^97_", valid_ids)
    is_asv_id <- grepl("^(asv|otu)_", valid_ids)
  }

  if (length(valid_ids) == 0) {
    stop("La tabla de abundancias propia no contiene features reconocibles: ", path)
  }

  own_otu_wide <- empty_wide()
  own_otu_long <- empty_long()
  if (any(is_otu_id)) {
    own_otu <- copy(own[grepl("^97_", FEATURE)])
    setnames(own_otu, 1, "FEATURE_ID")
    own_otu_wide <- aggregate_wide_by_otu(own_otu)
    own_otu_long <- wide_to_long(copy(own_otu_wide), "OTU")
  }

  own_asv <- copy(own[grepl("^(asv|otu)_", FEATURE)])
  if (nrow(own_asv) == 0) {
    return(list(
      own_otu_long = own_otu_long,
      asv_otu_long = empty_long(),
      asv_unmapped_long = empty_long(),
      own_otu_wide = own_otu_wide,
      asv_otu_wide = empty_wide(),
      asv_unmapped_wide = empty_wide(),
      asv_feature_map = empty_asv_feature_map(),
      mode = "otu"
    ))
  }

  setnames(own_asv, 1, "ASV")
  own_long <- melt(
    own_asv,
    id.vars = "ASV",
    variable.name = "SAMPLE",
    value.name = "abundance"
  )
  own_long[, abundance := as.numeric(abundance)]
  own_long[is.na(abundance), abundance := 0]
  own_long <- own_long[abundance > 0]

  asv_to_otu <- read_asv_to_otu(map_path)
  cluster_map <- read_asv_cluster_feature_map(cluster_feature_map_path)

  asv_mapped <- merge(own_long, asv_to_otu, by = "ASV")
  asv_otu_long <- asv_mapped[
    ,
    .(abundance = sum(abundance)),
    by = .(SAMPLE, FEATURE_ID, FEATURE_TYPE)
  ]
  mapped_feature_map <- unique(asv_to_otu[
    ASV %in% unique(own_long$ASV),
    .(
      ASV,
      FEATURE_ID,
      FEATURE_TYPE,
      assignment_source = "MicrobeAtlas_vsearch",
      pident,
      alnlen,
      cluster_centroid = NA_character_,
      cluster_size = NA_integer_
    )
  ])

  directly_mapped <- unique(asv_to_otu$ASV)
  own_unmapped <- own_long[!ASV %in% directly_mapped]
  asv_clustered <- merge(own_unmapped, cluster_map, by = "ASV")
  missing_cluster <- setdiff(unique(own_unmapped$ASV), unique(cluster_map$ASV))
  if (length(missing_cluster) > 0) {
    stop(
      "Hay ",
      length(missing_cluster),
      " ASVs sin OTU97 que no aparecen en el mapa formal del paso 06. Ejemplos: ",
      paste(head(missing_cluster, 10), collapse = ", ")
    )
  }

  asv_unmapped_long <- asv_clustered[
    ,
    .(abundance = sum(abundance)),
    by = .(SAMPLE, FEATURE_ID, FEATURE_TYPE)
  ]
  unmapped_feature_map <- unique(cluster_map[
    ASV %in% unique(own_unmapped$ASV),
    .(
      ASV,
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
    asv_otu_long = asv_otu_long,
    asv_unmapped_long = asv_unmapped_long,
    own_otu_wide = own_otu_wide,
    asv_otu_wide = wide_from_long(asv_otu_long),
    asv_unmapped_wide = wide_from_long(asv_unmapped_long),
    asv_feature_map = rbindlist(
      list(mapped_feature_map, unmapped_feature_map),
      use.names = TRUE,
      fill = TRUE
    ),
    mode = if (any(is_otu_id)) "mixed_otu_asv" else "asv"
  )
}

asv_identity_stats <- validate_own_asv_identity(
  OWN_ABUND_FILE,
  ASV_FASTA_FILE,
  MAP_FILE,
  CLUSTER_UC_FILE
)

own_tables <- make_own_otu_long(
  OWN_ABUND_FILE,
  MAP_FILE,
  NO_HITS_FEATURE_MAP_FILE,
  allow_own_otu_input = ALLOW_OWN_OTU_INPUT
)
own_otu_long <- own_tables$own_otu_long
asv_otu_long <- own_tables$asv_otu_long
asv_long_unmapped <- own_tables$asv_unmapped_long
own_wide <- join_wide_tables(
  list(
    own_tables$own_otu_wide,
    own_tables$asv_otu_wide,
    own_tables$asv_unmapped_wide
  )
)

# tabla larga dispersa final
atlas_stats <- write_positive_long_from_biom_tsv(ATLAS_FILE, OUT_LONG_TOTAL)
append_long_file(own_otu_long, OUT_LONG_TOTAL)
append_long_file(asv_otu_long, OUT_LONG_TOTAL)
append_long_file(asv_long_unmapped, OUT_LONG_TOTAL)

# formato ancho
wide_stats <- write_final_wide_from_biom_tsv(ATLAS_FILE, own_wide, OUT_WIDE_TOTAL)

fwrite(own_tables$asv_feature_map, OUT_ASV_FEATURE_MAP, sep = "\t")
coverage_stats <- write_microbeatlas_coverage(ATLAS_FILE, OUT_ATLAS_COVERAGE)
write_pipeline_summaries(
  OUT_QUALITY_SUMMARY,
  OUT_METRICS_SUMMARY,
  OUT_REPRODUCIBILITY,
  asv_identity_stats,
  own_tables,
  atlas_stats,
  wide_stats,
  coverage_stats,
  own_wide,
  own_otu_long,
  asv_otu_long,
  asv_long_unmapped
)

cat("Guardado tabla final ancha:", OUT_WIDE_TOTAL, "\n")
cat("Guardado tabla final larga:", OUT_LONG_TOTAL, "\n")
cat("Mapa ASV -> feature:", OUT_ASV_FEATURE_MAP, "\n")
cat("Cobertura MicrobeAtlas OTU97/Unmapped:", OUT_ATLAS_COVERAGE, "\n")
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
message_count("Features ASV mapeadas a OTU97:", uniqueN(asv_otu_long$FEATURE_ID))
message_count(
  "Features ASV no mapeadas clusterizadas:",
  uniqueN(asv_long_unmapped[FEATURE_TYPE == "Own_ASVcluster97"]$FEATURE_ID)
)
message_count(
  "Features ASV no mapeadas singleton:",
  uniqueN(asv_long_unmapped[FEATURE_TYPE == "Own_ASVsingleton"]$FEATURE_ID)
)
message_count("Filas finales OTU/ASV:", wide_stats$atlas_otu_count + wide_stats$own_only_otu_count)
message_count("Muestras:", wide_stats$n_samples)

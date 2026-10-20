#!/usr/bin/env Rscript
# Alternative sample filter for MicrobeAtlas.
#
# Filters community clusters by environment/seqtech/robust projection and
# selects a reproducible random number of samples per community cluster.

# TRAER EL SCRIPT CATÁLOGO
# AÑADIR PASO DE OBSERVAR LAS MUESTRAS

library(tidyverse)
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
IN_DIR <- Sys.getenv("MICROBEATLAS_RAW_DATA_DIR", file.path(PROJECT_ROOT, "data", "raw"))
BIOM_SAMPLE_IDS_FILE <- Sys.getenv(
  "MICROBEATLAS_BIOM_SAMPLE_IDS",
  file.path(DATA_DIR, "02_filters", "biom_sample_ids.txt")
)
ENV_CATALOG_FILE <- Sys.getenv(
  "MICROBEATLAS_ENV_CATALOG_FILE",
  file.path(DATA_DIR, "02_filters", "environment_catalog.tsv")
)

AUTOMATIC_SEED <- 123L

defaults <- list(
  env_id = "",
  exclude_env_id = "",
  seqtech = "AMPLICON",
  robust_only = "TRUE",
  samples_per_cc = Sys.getenv("MICROBEATLAS_SAMPLES_PER_CC", "1"),
  allow_missing_biom_cc = Sys.getenv("MICROBEATLAS_ALLOW_MISSING_BIOM_CC", "FALSE"),
  list_envs = FALSE
)

is_true <- function(x) {
  tolower(trimws(as.character(x))) == "true"
}

split_csv <- function(x) {
  x <- trimws(unlist(strsplit(x, ",", fixed = TRUE)))
  x[x != ""]
}

parse_args <- function(args) {
  opts <- defaults
  i <- 1L

  while (i <= length(args)) {
    arg <- args[[i]]

    if (arg %in% c("--help", "-h")) {
      cat("Uso: Rscript 02_sample_selection.R [--list-envs] [--env-id 1,3] [--exclude-env-id 2,4] [--seqtech AMPLICON] [--robust-only TRUE|FALSE] [--samples-per-cc 1|2|3|all]\n")
      quit(status = 0)
    }

    if (arg == "--list-envs") {
      opts$list_envs <- TRUE
      i <- i + 1L
      next
    }

    if (grepl("^--[^=]+=", arg)) {
      key <- sub("^--", "", sub("=.*$", "", arg))
      value <- sub("^[^=]+=", "", arg)
    } else if (grepl("^--", arg)) {
      key <- sub("^--", "", arg)
      if (i == length(args) || grepl("^--", args[[i + 1L]])) {
        stop("Falta valor para argumento: ", arg)
      }
      if (key %in% c("env-id", "exclude-env-id")) {
        j <- i + 1L
        values <- character()
        while (j <= length(args) && !grepl("^--", args[[j]])) {
          values <- c(values, args[[j]])
          j <- j + 1L
        }
        value <- paste(values, collapse = ",")
        i <- j - 1L
      } else {
        value <- args[[i + 1L]]
        i <- i + 1L
      }
    } else {
      stop("Argumento no reconocido: ", arg)
    }

    if (key == "env-id") {
      opts$env_id <- value
    } else if (key == "exclude-env-id") {
      opts$exclude_env_id <- value
    } else if (key == "seqtech") {
      opts$seqtech <- value
    } else if (key == "robust-only") {
      opts$robust_only <- value
    } else if (key == "samples-per-cc") {
      opts$samples_per_cc <- value
    } else if (key == "seed") {
      warning("--seed ya no es un parametro de usuario; se ignora y se usa la seed automatica interna: ", AUTOMATIC_SEED)
    } else {
      stop("Argumento no reconocido: --", key)
    }

    i <- i + 1L
  }

  opts
}

parse_samples_per_cc <- function(x) {
  x <- tolower(trimws(as.character(x)))
  if (x %in% c("all", "todos")) {
    return(NA_integer_)
  }

  out <- suppressWarnings(as.integer(x))
  if (is.na(out) || out < 1L || as.character(out) != x) {
    stop("--samples-per-cc debe ser un entero positivo o 'all': ", x)
  }
  out
}

selection_label <- function(samples_per_cc) {
  if (is.na(samples_per_cc)) {
    return("allsamplesperCC")
  }
  if (samples_per_cc == 1L) {
    return("1sampelperCC")
  }
  paste0(samples_per_cc, "samplesperCC")
}

sanitize_label <- function(x) {
  x <- tolower(trimws(x))
  x <- gsub("[^a-z0-9]+", "-", x)
  x <- gsub("^-+|-+$", "", x)
  ifelse(x == "", "selection", x)
}

default_env_label <- function(opts) {
  env_id <- split_csv(opts$env_id)
  exclude_env_id <- split_csv(opts$exclude_env_id)
  if (length(env_id) == 0) {
    return("all")
  }
  label <- paste0("env", paste(env_id, collapse = "-"))
  if (length(exclude_env_id) > 0) {
    label <- paste0(label, "_exclude", paste(exclude_env_id, collapse = "-"))
  }
  sanitize_label(label)
}

resolve_env_ids <- function(refs, env_catalog, option_name) {
  refs <- split_csv(refs)
  if (length(refs) == 0) {
    return(character())
  }

  invalid_refs <- refs[!grepl("^[0-9]+$", refs)]
  if (length(invalid_refs) > 0) {
    stop(
      option_name,
      " solo acepta numeros del catalogo de ambientes. Valores invalidos: ",
      paste(invalid_refs, collapse = ", "),
      ". Usa --list-envs para ver los numeros disponibles."
    )
  }

  idx <- as.integer(refs)
  invalid_idx <- idx[idx < 1L | idx > nrow(env_catalog)]
  if (length(invalid_idx) > 0) {
    stop(
      option_name,
      " fuera de rango: ",
      paste(invalid_idx, collapse = ", "),
      ". Mira el catalogo en: ",
      ENV_CATALOG_FILE
    )
  }

  unique(env_catalog$Env_cons[idx])
}

matches_env_tree <- function(env_values, roots) {
  if (length(roots) == 0) {
    return(rep(FALSE, length(env_values)))
  }
  Reduce(`|`, lapply(roots, function(root) {
    env_values == root | startsWith(env_values, paste0(root, ";"))
  }))
}

filter_base_metadata <- function(metadata, opts) {
  out <- metadata %>%
    filter(!is.na(SeqTech_cons), SeqTech_cons == opts$seqtech) %>%
    filter(!is.na(Env_cons), Env_cons != "")

  if (is_true(opts$robust_only)) {
    out <- out %>% filter(is_true(Robust_projection))
  }

  out
}

select_environments <- function(metadata, opts, env_catalog) {
  selected_roots <- resolve_env_ids(opts$env_id, env_catalog, "--env-id")
  excluded_roots <- resolve_env_ids(opts$exclude_env_id, env_catalog, "--exclude-env-id")

  if (length(selected_roots) == 0) {
    keep <- rep(TRUE, nrow(metadata))
  } else {
    keep <- matches_env_tree(metadata$Env_cons, selected_roots)
  }

  if (length(excluded_roots) > 0) {
    keep <- keep & !matches_env_tree(metadata$Env_cons, excluded_roots)
  }

  metadata[keep, ]
}

make_env_summary <- function(metadata, mapping) {
  if (nrow(metadata) == 0) {
    return(data.table(
      idx = integer(),
      Env_cons = character(),
      N_ComClusters = integer(),
      N_Samples = integer()
    ))
  }

  metadata_dt <- as.data.table(metadata)
  mapping_dt <- as.data.table(mapping)

  clusters_per_env <- unique(metadata_dt[, .(ComClustID, Env_cons)])[
    ,
    .(N_ComClusters = uniqueN(ComClustID)),
    by = Env_cons
  ]
  merged <- merge(
    mapping_dt,
    metadata_dt[, .(ComClustID, Env_cons)],
    by = "ComClustID",
    allow.cartesian = TRUE
  )
  samples_per_env <- unique(merged[, .(SampleID, Env_cons)])[
    ,
    .(N_Samples = uniqueN(SampleID)),
    by = Env_cons
  ]
  summary <- merge(clusters_per_env, samples_per_env, by = "Env_cons", all.x = TRUE)
  summary[is.na(N_Samples), N_Samples := 0L]
  setorder(summary, Env_cons)
  summary[, idx := .I]
  setcolorder(summary, c("idx", "Env_cons", "N_ComClusters", "N_Samples"))
  summary
}

read_biom_sample_ids <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }

  ids <- fread(path, header = FALSE, colClasses = "character")[[1]]
  ids <- ids[!is.na(ids) & ids != ""]
  unique(ids)
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))
samples_per_cc <- parse_samples_per_cc(opts$samples_per_cc)
selection_name <- selection_label(samples_per_cc)
env_label <- Sys.getenv("MICROBEATLAS_ENV_LABEL", default_env_label(opts))
OUT_DIR <- Sys.getenv(
  "MICROBEATLAS_FILTER_DIR",
  file.path(DATA_DIR, paste0("02_filters_", env_label, "_", selection_name))
)
SUMMARY_FILE <- Sys.getenv(
  "MICROBEATLAS_FILTER_SUMMARY_FILE",
  file.path(OUT_DIR, "filter_summary.tsv")
)
SELECTION_FILE <- Sys.getenv(
  "MICROBEATLAS_FILTER_SELECTION_FILE",
  file.path(OUT_DIR, paste0("community_cluster_selection_", selection_name, ".tsv"))
)
BIOM_SUBSET_IDS_FILE <- Sys.getenv(
  "MICROBEATLAS_FILTER_BIOM_SUBSET_IDS_FILE",
  file.path(OUT_DIR, paste0("community_cluster_selection_", selection_name, "_forbiomsubset.txt"))
)
SAMPLE_METADATA_FILE <- Sys.getenv(
  "MICROBEATLAS_SAMPLE_METADATA_FILE",
  file.path(OUT_DIR, paste0("sample_metadata_", selection_name, ".tsv"))
)
SELECTED_ENV_COUNTS_FILE <- Sys.getenv(
  "MICROBEATLAS_SELECTED_ENV_COUNTS_FILE",
  file.path(OUT_DIR, paste0("selected_environment_counts_", selection_name, ".tsv"))
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(ENV_CATALOG_FILE), recursive = TRUE, showWarnings = FALSE)

set.seed(AUTOMATIC_SEED)

cc_metadata <- fread(file.path(IN_DIR, "community_cluster_metadata.tsv"), sep = "\t")
cc_mapping <- fread(file.path(IN_DIR, "community_cluster_mapping.tsv"), sep = "\t")

cc_metadata_base <- filter_base_metadata(cc_metadata, opts)
env_catalog <- make_env_summary(cc_metadata_base, cc_mapping)
fwrite(env_catalog, file = ENV_CATALOG_FILE, sep = "\t")

if (opts$list_envs) {
  fwrite(env_catalog, sep = "\t")
  quit(status = 0)
}

cc_metadata_filtered <- select_environments(cc_metadata_base, opts, env_catalog) %>%
  as.data.frame()

cc_mapping_filtered <- cc_metadata_filtered %>%
  inner_join(cc_mapping, by = "ComClustID") %>%
  filter(!is.na(SampleID), SampleID != "")

if (nrow(cc_mapping_filtered) == 0) {
  stop("No hay muestras tras filtrar CC robustos, AMPLICON y ambientes seleccionados")
}

biom_sample_ids <- read_biom_sample_ids(BIOM_SAMPLE_IDS_FILE)
if (is.null(biom_sample_ids)) {
  message("No se encontro lista de muestras BIOM; se muestreara sin interseccion previa: ", BIOM_SAMPLE_IDS_FILE)
  cc_mapping_for_sampling <- cc_mapping_filtered
} else {
  cc_mapping_for_sampling <- cc_mapping_filtered %>%
    filter(SampleID %in% biom_sample_ids)
}

missing_cc <- setdiff(unique(cc_mapping_filtered$ComClustID), unique(cc_mapping_for_sampling$ComClustID))
if (length(missing_cc) > 0) {
  missing_message <- paste0(
    "Hay ",
    length(missing_cc),
    " community clusters sin ninguna muestra presente en el BIOM. Ejemplos: ",
    paste(head(missing_cc, 10), collapse = ", ")
  )

  if (!is_true(opts$allow_missing_biom_cc)) {
    stop(missing_message)
  }

  warning(
    missing_message,
    ". Se continuara solo con los community clusters disponibles en el BIOM porque ",
    "MICROBEATLAS_ALLOW_MISSING_BIOM_CC=TRUE."
  )
}

if (nrow(cc_mapping_for_sampling) == 0) {
  stop("No quedan muestras disponibles tras intersectar la seleccion con el BIOM")
}

if (is.na(samples_per_cc)) {
  cc_mapping_sampled <- cc_mapping_for_sampling %>%
    distinct()
} else {
  cc_mapping_sampled <- cc_mapping_for_sampling %>%
    group_by(ComClustID) %>%
    group_modify(~ slice_sample(.x, n = min(nrow(.x), samples_per_cc))) %>%
    ungroup()
}

clusters_below_requested_samples <- if (is.na(samples_per_cc)) {
  0L
} else {
  cc_mapping_for_sampling %>%
    count(ComClustID, name = "n_available_samples") %>%
    filter(n_available_samples < samples_per_cc) %>%
    nrow()
}

fwrite(
  cc_mapping_sampled,
  file = SELECTION_FILE,
  sep = "\t"
)
fwrite(
  cc_mapping_sampled %>% ungroup() %>% distinct(SampleID),
  file = BIOM_SUBSET_IDS_FILE,
  sep = "\t",
  col.names = FALSE
)

sample_metadata <- cc_mapping_sampled %>%
  ungroup() %>%
  mutate(
    MICROBEATLAS_ENV_ROOT = sub(";.*$", "", Env_cons),
    SAMPLE_SOURCE = "MicrobeAtlas"
  ) %>%
  select(
    SampleID,
    SAMPLE_SOURCE,
    ComClustID,
    Env_cons,
    MICROBEATLAS_ENV_ROOT,
    SeqTech_cons,
    Robust_projection,
    Keyw_cons,
    Proj_cons
  ) %>%
  distinct()
fwrite(sample_metadata, file = SAMPLE_METADATA_FILE, sep = "\t")

selected_env_counts <- as.data.table(sample_metadata)[
  ,
  .(
    n_samples = uniqueN(SampleID),
    n_community_clusters = uniqueN(ComClustID)
  ),
  by = .(MICROBEATLAS_ENV_ROOT, Env_cons)
]
setorder(selected_env_counts, MICROBEATLAS_ENV_ROOT, Env_cons)
fwrite(selected_env_counts, file = SELECTED_ENV_COUNTS_FILE, sep = "\t")

selected_envs <- sort(unique(cc_mapping_filtered$Env_cons))
excluded_roots <- resolve_env_ids(opts$exclude_env_id, env_catalog, "--exclude-env-id")
excluded_envs <- sort(unique(env_catalog$Env_cons[matches_env_tree(env_catalog$Env_cons, excluded_roots)]))
summary_table <- data.table(
  parameter = c(
    "env_id",
    "exclude_env_id",
    "seqtech",
    "robust_only",
    "allow_missing_biom_cc",
    "seed",
    "seed_source",
    "samples_per_cc",
    "selection_label",
    "environment_catalog_file",
    "biom_sample_ids_file",
    "selection_file",
    "biom_subset_ids_file",
    "sample_metadata_file",
    "selected_environment_counts_file",
    "biom_sample_ids_available",
    "selected_envs",
    "excluded_envs",
    "n_selected_envs",
    "n_community_clusters",
    "n_samples_all",
    "n_samples_available_for_sampling",
    "n_community_clusters_missing_from_biom",
    "missing_from_biom_examples",
    "n_community_clusters_below_requested_samples",
    "n_sampled_rows",
    "n_sampled_unique_samples"
  ),
  value = c(
    opts$env_id,
    opts$exclude_env_id,
    opts$seqtech,
    opts$robust_only,
    opts$allow_missing_biom_cc,
    AUTOMATIC_SEED,
    "automatic_internal",
    ifelse(is.na(samples_per_cc), "all", samples_per_cc),
    selection_name,
    ENV_CATALOG_FILE,
    BIOM_SAMPLE_IDS_FILE,
    SELECTION_FILE,
    BIOM_SUBSET_IDS_FILE,
    SAMPLE_METADATA_FILE,
    SELECTED_ENV_COUNTS_FILE,
    !is.null(biom_sample_ids),
    paste(selected_envs, collapse = ","),
    paste(excluded_envs, collapse = ","),
    length(selected_envs),
    uniqueN(cc_mapping_filtered$ComClustID),
    uniqueN(cc_mapping_filtered$SampleID),
    uniqueN(cc_mapping_for_sampling$SampleID),
    length(missing_cc),
    paste(head(missing_cc, 50), collapse = ","),
    clusters_below_requested_samples,
    nrow(cc_mapping_sampled),
    uniqueN(cc_mapping_sampled$SampleID)
  )
)
fwrite(summary_table, file = SUMMARY_FILE, sep = "\t")

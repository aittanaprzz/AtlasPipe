#!/usr/bin/env Rscript
# Alternative sample filter for MicrobeAtlas.
#
# Filters community clusters by environment/seqtech/robust projection and
# selects a reproducible random number of samples per community cluster.

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
  "MICROBEATLAS_ROOT",
  normalizePath(file.path(script_dir(), ".."), mustWork = FALSE)
)
DATA_DIR <- Sys.getenv("MICROBEATLAS_DATA_DIR", file.path(PROJECT_ROOT, "data", "intermediate"))
IN_DIR <- Sys.getenv("MICROBEATLAS_RAW_DATA_DIR", file.path(PROJECT_ROOT, "data", "raw", "microbeatlas"))
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
  env = "all",
  env_id = "",
  env_prefix = "",
  env_file = "",
  exclude_env = "",
  exclude_prefix = "",
  seqtech = "AMPLICON",
  robust_only = "TRUE",
  samples_per_cc = Sys.getenv("MICROBEATLAS_SAMPLES_PER_CC", "1"),
  list_envs = FALSE
)

is_true <- function(x) {
  tolower(trimws(as.character(x))) == "true"
}

split_csv <- function(x) {
  x <- trimws(unlist(strsplit(x, ",", fixed = TRUE)))
  x[x != ""]
}

read_env_file <- function(path) {
  if (path == "") {
    return(character())
  }
  if (!file.exists(path)) {
    stop("No existe el fichero de ambientes: ", path)
  }

  envs <- trimws(readLines(path, warn = FALSE))
  envs <- envs[envs != "" & !grepl("^#", envs)]
  unique(envs)
}

parse_args <- function(args) {
  opts <- defaults
  i <- 1L

  while (i <= length(args)) {
    arg <- args[[i]]

    if (arg %in% c("--help", "-h")) {
      cat("Uso: Rscript 02_filter_samples_1sampelperCC.R [--list-envs] [--env all|env1,env2] [--env-id 1,3] [--env-prefix prefix] [--env-file file] [--exclude-env 1 2|env1,env2] [--exclude-prefix prefix] [--seqtech AMPLICON] [--robust-only TRUE|FALSE] [--samples-per-cc 1|2|3|all]\n")
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
      if (key %in% c("exclude-env", "exclude-prefix")) {
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

    if (key == "env") {
      opts$env <- value
    } else if (key == "env-id") {
      opts$env_id <- value
    } else if (key == "env-prefix") {
      opts$env_prefix <- value
    } else if (key == "env-file") {
      opts$env_file <- value
    } else if (key == "exclude-env") {
      opts$exclude_env <- value
    } else if (key == "exclude-prefix") {
      opts$exclude_prefix <- value
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
  env_prefix <- split_csv(opts$env_prefix)
  if (length(env_prefix) > 0) {
    return(paste(sanitize_label(env_prefix), collapse = "-"))
  }

  env_id <- split_csv(opts$env_id)
  env_file <- trimws(opts$env_file)
  env <- split_csv(opts$env)
  if (length(env_id) > 0 || env_file != "") {
    return("customenv")
  }
  if (length(env) == 0 || identical(tolower(env), "all")) {
    return("all")
  }
  paste(sanitize_label(env), collapse = "-")
}

resolve_env_refs <- function(refs, env_catalog, option_name) {
  refs <- split_csv(refs)
  if (length(refs) == 0) {
    return(character())
  }

  is_idx <- grepl("^[0-9]+$", refs)
  if (all(is_idx)) {
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
    return(unique(env_catalog$Env_cons[idx]))
  }

  if (any(is_idx)) {
    stop(option_name, " no puede mezclar indices numericos y nombres de ambientes en el mismo argumento")
  }

  invalid_env <- setdiff(refs, env_catalog$Env_cons)
  if (length(invalid_env) > 0) {
    stop(
      option_name,
      " contiene ambientes no disponibles tras los filtros base: ",
      paste(invalid_env, collapse = ", "),
      ". Mira el catalogo en: ",
      ENV_CATALOG_FILE
    )
  }

  unique(refs)
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
  env_exact <- split_csv(opts$env)
  env_prefix <- split_csv(opts$env_prefix)
  env_from_file <- read_env_file(opts$env_file)

  use_all_envs <- length(env_exact) == 0 || identical(tolower(env_exact), "all")
  if (
    use_all_envs &&
    length(split_csv(opts$env_id)) == 0 &&
    length(env_prefix) == 0 &&
    length(env_from_file) == 0
  ) {
    keep <- rep(TRUE, nrow(metadata))
  } else {
    env_exact <- env_exact[tolower(env_exact) != "all"]
    env_exact <- unique(c(env_exact, env_from_file))
    env_exact <- unique(c(env_exact, resolve_env_refs(opts$env_id, env_catalog, "--env-id")))

    keep <- rep(FALSE, nrow(metadata))
    if (length(env_exact) > 0) {
      keep <- keep | metadata$Env_cons %in% env_exact
    }
    if (length(env_prefix) > 0) {
      keep <- keep | Reduce(`|`, lapply(env_prefix, function(prefix) {
        startsWith(metadata$Env_cons, prefix)
      }))
    }
  }

  exclude_env <- resolve_env_refs(opts$exclude_env, env_catalog, "--exclude-env")
  exclude_prefix <- split_csv(opts$exclude_prefix)
  if (length(exclude_env) > 0) {
    keep <- keep & !(metadata$Env_cons %in% exclude_env)
  }
  if (length(exclude_prefix) > 0) {
    drop_by_prefix <- Reduce(`|`, lapply(exclude_prefix, function(prefix) {
      startsWith(metadata$Env_cons, prefix)
    }))
    keep <- keep & !drop_by_prefix
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
  stop(
    "Hay ",
    length(missing_cc),
    " community clusters sin ninguna muestra presente en el BIOM. Ejemplos: ",
    paste(head(missing_cc, 10), collapse = ", ")
  )
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

selected_envs <- sort(unique(cc_mapping_filtered$Env_cons))
excluded_by_env <- resolve_env_refs(opts$exclude_env, env_catalog, "--exclude-env")
exclude_prefix <- split_csv(opts$exclude_prefix)
excluded_by_prefix <- character()
if (length(exclude_prefix) > 0) {
  excluded_by_prefix <- unique(unlist(lapply(exclude_prefix, function(prefix) {
    env_catalog$Env_cons[startsWith(env_catalog$Env_cons, prefix)]
  })))
}
excluded_envs <- sort(unique(c(excluded_by_env, excluded_by_prefix)))
summary_table <- data.table(
  parameter = c(
    "env_argument",
    "env_id",
    "env_prefix",
    "env_file",
    "exclude_env",
    "exclude_prefix",
    "seqtech",
    "robust_only",
    "seed",
    "seed_source",
    "samples_per_cc",
    "selection_label",
    "environment_catalog_file",
    "biom_sample_ids_file",
    "selection_file",
    "biom_subset_ids_file",
    "biom_sample_ids_available",
    "selected_envs",
    "excluded_envs",
    "n_selected_envs",
    "n_community_clusters",
    "n_samples_all",
    "n_samples_available_for_sampling",
    "n_community_clusters_below_requested_samples",
    "n_sampled_rows",
    "n_sampled_unique_samples"
  ),
  value = c(
    opts$env,
    opts$env_id,
    opts$env_prefix,
    opts$env_file,
    opts$exclude_env,
    opts$exclude_prefix,
    opts$seqtech,
    opts$robust_only,
    AUTOMATIC_SEED,
    "automatic_internal",
    ifelse(is.na(samples_per_cc), "all", samples_per_cc),
    selection_name,
    ENV_CATALOG_FILE,
    BIOM_SAMPLE_IDS_FILE,
    SELECTION_FILE,
    BIOM_SUBSET_IDS_FILE,
    !is.null(biom_sample_ids),
    paste(selected_envs, collapse = ","),
    paste(excluded_envs, collapse = ","),
    length(selected_envs),
    uniqueN(cc_mapping_filtered$ComClustID),
    uniqueN(cc_mapping_filtered$SampleID),
    uniqueN(cc_mapping_for_sampling$SampleID),
    clusters_below_requested_samples,
    nrow(cc_mapping_sampled),
    uniqueN(cc_mapping_sampled$SampleID)
  )
)
fwrite(summary_table, file = SUMMARY_FILE, sep = "\t")

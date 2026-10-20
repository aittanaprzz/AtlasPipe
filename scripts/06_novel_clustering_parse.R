#!/usr/bin/env Rscript

library(data.table)

script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
  }
  getwd()
}

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!grepl("^--", key)) {
      stop("Argumento inesperado: ", key)
    }
    if (i == length(args)) {
      stop("Falta valor para argumento: ", key)
    }
    out[[sub("^--", "", key)]] <- args[[i + 1L]]
    i <- i + 2L
  }
  out
}

read_fasta_ids <- function(path) {
  if (!file.exists(path)) {
    stop("No existe el FASTA de amplicones no mapeados: ", path)
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

empty_feature_map <- function() {
  data.table(
    AMPLICON_ID = character(),
    FEATURE_ID = character(),
    FEATURE_TYPE = character(),
    assignment_source = character(),
    cluster_id = character(),
    cluster_centroid = character(),
    cluster_size = integer(),
    is_centroid = logical()
  )
}

empty_cluster_sizes <- function() {
  data.table(
    cluster_id = character(),
    FEATURE_ID = character(),
    FEATURE_TYPE = character(),
    cluster_centroid = character(),
    cluster_size = integer()
  )
}

build_cluster_tables <- function(uc_path, nohit_fasta_path) {
  nohit_ids <- read_fasta_ids(nohit_fasta_path)
  duplicated_nohit <- unique(nohit_ids[duplicated(nohit_ids)])
  if (length(duplicated_nohit) > 0) {
    stop(
      "amplicones duplicados en amplicons_no_hits.fasta. Ejemplos: ",
      paste(head(duplicated_nohit, 10), collapse = ", ")
    )
  }

  if (length(nohit_ids) == 0) {
    return(list(
      feature_map = empty_feature_map(),
      cluster_sizes = empty_cluster_sizes()
    ))
  }

  if (!file.exists(uc_path)) {
    stop("No existe el fichero .uc de clustering de amplicones no mapeados: ", uc_path)
  }

  uc <- fread(uc_path, header = FALSE, fill = TRUE, colClasses = "character")
  if (ncol(uc) < 10) {
    stop("El fichero .uc no tiene las 10 columnas esperadas de VSEARCH: ", uc_path)
  }

  seeds <- uc[
    V1 == "S" & !is.na(V2) & V2 != "" & !is.na(V9) & V9 != "",
    .(cluster_id = V2, AMPLICON_ID = V9, cluster_centroid = V9)
  ]
  hits <- uc[
    V1 == "H" & !is.na(V2) & V2 != "" & !is.na(V9) & V9 != "" & !is.na(V10) & V10 != "*",
    .(cluster_id = V2, AMPLICON_ID = V9, cluster_centroid = V10)
  ]

  cluster_map <- rbindlist(list(seeds, hits), use.names = TRUE)
  if (nrow(cluster_map) == 0) {
    stop("El fichero .uc no contiene filas S/H validas: ", uc_path)
  }

  duplicated_assignments <- unique(cluster_map$AMPLICON_ID[duplicated(cluster_map$AMPLICON_ID)])
  if (length(duplicated_assignments) > 0) {
    stop(
      "amplicones con mas de una asignacion en .uc. Ejemplos: ",
      paste(head(duplicated_assignments, 10), collapse = ", ")
    )
  }

  missing_in_uc <- setdiff(nohit_ids, cluster_map$AMPLICON_ID)
  extra_in_uc <- setdiff(cluster_map$AMPLICON_ID, nohit_ids)
  if (length(missing_in_uc) > 0 || length(extra_in_uc) > 0) {
    stop(
      "Incoherencia amplicons_no_hits.fasta vs .uc. ",
      "amplicones del FASTA ausentes en .uc: ", length(missing_in_uc),
      if (length(missing_in_uc) > 0) paste0(" (ejemplos: ", paste(head(missing_in_uc, 10), collapse = ", "), ")") else "",
      ". amplicones extra en .uc: ", length(extra_in_uc),
      if (length(extra_in_uc) > 0) paste0(" (ejemplos: ", paste(head(extra_in_uc, 10), collapse = ", "), ")") else ""
    )
  }

  cluster_sizes <- cluster_map[
    ,
    .(
      cluster_centroid = cluster_centroid[1],
      cluster_size = uniqueN(AMPLICON_ID)
    ),
    by = cluster_id
  ]
  cluster_sizes[, FEATURE_TYPE := "Own_cluster97"]
  cluster_sizes[, FEATURE_ID := paste0("cluster97_", cluster_centroid)]
  setcolorder(cluster_sizes, c("cluster_id", "FEATURE_ID", "FEATURE_TYPE", "cluster_centroid", "cluster_size"))
  cluster_sizes[, cluster_order := as.integer(cluster_id)]
  setorder(cluster_sizes, cluster_order)
  cluster_sizes[, cluster_order := NULL]

  feature_map <- merge(cluster_map, cluster_sizes, by = c("cluster_id", "cluster_centroid"))
  feature_map[, is_centroid := AMPLICON_ID == cluster_centroid]
  feature_map[, assignment_source := "own_nohit_97cluster"]
  feature_map <- feature_map[
    ,
    .(AMPLICON_ID, FEATURE_ID, FEATURE_TYPE, assignment_source, cluster_id, cluster_centroid, cluster_size, is_centroid)
  ]
  feature_map[, cluster_order := as.integer(cluster_id)]
  feature_map[, non_centroid_order := !is_centroid]
  setorder(feature_map, cluster_order, non_centroid_order, AMPLICON_ID)
  feature_map[, c("cluster_order", "non_centroid_order") := NULL]

  list(
    feature_map = feature_map,
    cluster_sizes = cluster_sizes
  )
}

write_summary <- function(tables, out_file) {
  sizes <- tables$cluster_sizes$cluster_size
  nohit_amplicons_total <- nrow(tables$feature_map)
  own_novel_clusters_total <- sum(tables$cluster_sizes$FEATURE_TYPE == "Own_cluster97")
  amplicons_in_own_clusters97 <- sum(tables$cluster_sizes[FEATURE_TYPE == "Own_cluster97"]$cluster_size)

  summary <- data.table(
    metric = c(
      "nohit_amplicons_total",
      "own_novel_clusters_total",
      "own_novel_features_total",
      "amplicons_in_own_novel_clusters",
      "smallest_cluster_size",
      "median_cluster_size",
      "largest_cluster_size"
    ),
    value = c(
      nohit_amplicons_total,
      own_novel_clusters_total,
      own_novel_clusters_total,
      amplicons_in_own_clusters97,
      if (length(sizes) > 0) min(sizes) else 0,
      if (length(sizes) > 0) median(sizes) else 0,
      if (length(sizes) > 0) max(sizes) else 0
    )
  )
  fwrite(summary, out_file, sep = "\t")
}

args <- parse_args()
PROJECT_ROOT <- Sys.getenv(
  "ATLASPIPE_ROOT",
  normalizePath(file.path(script_dir(), ".."), mustWork = FALSE)
)
DATA_DIR <- Sys.getenv("ATLASPIPE_DATA_DIR", file.path(PROJECT_ROOT, "data", "intermediate"))
IN_DIR <- Sys.getenv(
  "ATLASPIPE_OWN_AMPLICON_MAPPING_DIR",
  file.path(DATA_DIR, "05_amplicon_mapping")
)

uc_file <- args[["uc"]]
if (is.null(uc_file)) {
  uc_file <- Sys.getenv("ATLASPIPE_OWN_NO_HITS_UC", file.path(IN_DIR, "amplicon_clusters97.uc"))
}
nohit_fasta <- args[["nohit-fasta"]]
if (is.null(nohit_fasta)) {
  nohit_fasta <- Sys.getenv("ATLASPIPE_OWN_NO_HITS_FASTA", file.path(IN_DIR, "amplicons_no_hits.fasta"))
}
feature_map_file <- args[["feature-map"]]
if (is.null(feature_map_file)) {
  feature_map_file <- Sys.getenv("ATLASPIPE_OWN_NO_HITS_FEATURE_MAP", file.path(IN_DIR, "amplicon_cluster_map.tsv"))
}
cluster_sizes_file <- args[["cluster-sizes"]]
if (is.null(cluster_sizes_file)) {
  cluster_sizes_file <- Sys.getenv("ATLASPIPE_OWN_NO_HITS_CLUSTER_SIZES", file.path(IN_DIR, "cluster_sizes.tsv"))
}
summary_file <- args[["summary"]]
if (is.null(summary_file)) {
  summary_file <- Sys.getenv("ATLASPIPE_OWN_NO_HITS_CLUSTER_SUMMARY", file.path(IN_DIR, "cluster_summary97.tsv"))
}

dir.create(dirname(feature_map_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(cluster_sizes_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(summary_file), recursive = TRUE, showWarnings = FALSE)

tables <- build_cluster_tables(uc_file, nohit_fasta)
fwrite(tables$feature_map, feature_map_file, sep = "\t")
fwrite(tables$cluster_sizes, cluster_sizes_file, sep = "\t")
write_summary(tables, summary_file)

cat("Mapa AMPLICON_ID no-hit -> feature propia:", feature_map_file, "\n")
cat("Tamanos de clusters propios:", cluster_sizes_file, "\n")
cat("Resumen clustering amplicones no-hit:", summary_file, "\n")

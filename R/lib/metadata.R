# Shared metadata/config loading + validation.
#
# `dataset_metadata.yml` describes ONE dataset x ONE matrix, resolved one of
# three ways (see R/lib/matrices.R::load_input_matrix() for the precedence):
#   1. `dataset.matrix_path` -- an already-built matrix, used as-is.
#   2. `dataset.preprocessing_script` -- a script YOU write that reads the
#      raw parquet trio, does arbitrary manipulation, and writes the result
#      out (see R/README.md's "Preprocessing script contract"). Any
#      sample/feature filtering belongs here, not as a post-hoc id list --
#      if you need more than a bare long-to-wide pivot, write a script.
#   3. neither given -- the framework reads `dataset.expression_path`
#      (long format) and pivots it wide with no other transformation.
# `cluster_config.yml` describes SLURM resource settings per method, kept
# separate so the same dataset metadata can be pointed at a different
# cluster.
#
# There is no "bases" concept here: this framework always runs on exactly
# one matrix per config. If you want to run, say, Day 0 and Day 2 samples
# separately, build/point at two matrices and write two dataset_metadata.yml
# files (one per run) -- see R/README.md.
#
# Both configs are validated eagerly (before any slurm_apply()/slurm_call()
# is constructed) so a typo'd path or empty grid fails fast on the login
# node instead of surfacing as a cryptic error deep in a batch array.

library(yaml)

#' Read and validate a dataset_metadata.yml file.
read_dataset_metadata <- function(path) {
  meta <- yaml::read_yaml(path)
  validate_dataset_metadata(meta)
  meta
}

#' Read and validate a cluster_config.yml file.
read_cluster_config <- function(path) {
  cfg <- yaml::read_yaml(path)
  validate_cluster_config(cfg)
  cfg
}

validate_dataset_metadata <- function(meta) {
  stopifnot(
    "dataset_metadata.yml must have a top-level `dataset` block" = !is.null(meta$dataset),
    "dataset block must have an `id`" = !is.null(meta$dataset$id)
  )

  matrix_path <- meta$dataset$matrix_path
  script_path <- meta$dataset$preprocessing_script

  if (!is.null(matrix_path)) {
    if (!file.exists(matrix_path)) {
      stop("dataset$matrix_path does not exist: ", matrix_path)
    }
  } else if (!is.null(script_path)) {
    if (!file.exists(script_path)) {
      stop("dataset$preprocessing_script does not exist: ", script_path)
    }
  } else {
    if (is.null(meta$dataset$expression_path) || !file.exists(meta$dataset$expression_path)) {
      stop("None of dataset$matrix_path, dataset$preprocessing_script are set, so ",
           "dataset$expression_path must point to an existing long-format ",
           "expression parquet file to pivot by default. See R/README.md's ",
           "\"Preparing your input matrix\".")
    }
  }

  if (!is.null(meta$methods)) {
    for (nm in names(meta$methods)) {
      m <- meta$methods[[nm]]
      if (isTRUE(m$enabled) && nm %in% c("pca", "nmf", "cogaps")) {
        if (is.null(m$rank_range) || length(m$rank_range) != 2) {
          stop("methods$", nm, "$rank_range must be a two-element [min, max]")
        }
      }
      if (isTRUE(m$enabled) && !is.null(m$masking_cv) && isTRUE(m$masking_cv$enabled)) {
        if (is.null(m$masking_cv$mask_frac) || is.null(m$masking_cv$mask_seed)) {
          stop("methods$", nm, "$masking_cv is enabled but missing mask_frac/mask_seed")
        }
      }
    }
  }

  invisible(meta)
}

validate_cluster_config <- function(cfg) {
  required <- c("mem", "cpus_per_task", "time", "container", "libPaths",
                "sh_template", "rscript_path")
  check_entry <- function(entry, label) {
    missing <- setdiff(required, names(entry))
    if (length(missing) > 0) {
      stop("cluster_config entry '", label, "' is missing required field(s): ",
           paste(missing, collapse = ", "))
    }
  }
  for (nm in names(cfg)) {
    entry <- cfg[[nm]]
    # `network` nests one level deeper (per-backend), everything else is flat.
    if (nm == "network") {
      for (backend in names(entry)) check_entry(entry[[backend]], paste0("network.", backend))
    } else {
      check_entry(entry, nm)
    }
  }
  invisible(cfg)
}

#' Expand a [min, max] rank_range list from YAML into an integer vector.
rank_seq <- function(rank_range) seq.int(rank_range[[1]], rank_range[[2]])

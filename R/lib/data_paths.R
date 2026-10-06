# Machine-independent dataset paths. Dataset configs give data paths relative
# to a data root: `SEPSIS_DATA_ROOT`, else `data_root` in config/pipeline.yml;
# absolute paths are used as written. Read configs through read_dataset_yaml()
# (docs: Data, "Data root").

#' The data root on this machine
#'
#' @param pipeline_path Path to `config/pipeline.yml`.
#' @return The `SEPSIS_DATA_ROOT` environment variable if set, else `data_root`
#'   from `config/pipeline.yml`, else `NULL`.
default_data_root <- function(pipeline_path = "config/pipeline.yml") {
  env <- Sys.getenv("SEPSIS_DATA_ROOT")
  if (nzchar(env)) return(env)
  if (!file.exists(pipeline_path) && requireNamespace("here", quietly = TRUE)) {
    pipeline_path <- here::here(pipeline_path)
  }
  if (file.exists(pipeline_path)) return(yaml::read_yaml(pipeline_path)$data_root)
  NULL
}

#' Resolve a dataset's relative data paths
#'
#' Prefixes `data_root` onto the relative `expression_path`,
#' `sample_metadata_path` and `feature_metadata_path` of a config's `dataset:`
#' block; absolute and `~` paths are left as written.
#'
#' @param ds A config's `dataset:` block (named list).
#' @param data_root Data root directory, or `NULL` to leave paths unchanged.
#' @return `ds` with its data paths resolved.
resolve_data_paths <- function(ds, data_root = default_data_root()) {
  if (is.null(data_root)) return(ds)
  for (key in c("expression_path", "sample_metadata_path", "feature_metadata_path")) {
    p <- ds[[key]]
    if (!is.null(p) && !startsWith(p, "/") && !startsWith(p, "~")) {
      ds[[key]] <- file.path(data_root, p)
    }
  }
  ds
}

#' Read a dataset config with its data paths resolved
#'
#' The one way configs should be read (`read_dataset_metadata()` adds
#' validation on top).
#'
#' @param path Path to `config/<dataset>_config.yml`.
#' @param data_root Data root directory (`default_data_root()`).
#' @return The parsed config, with `resolve_data_paths()` applied to `dataset:`.
read_dataset_yaml <- function(path, data_root = default_data_root()) {
  cfg <- yaml::read_yaml(path)
  cfg$dataset <- resolve_data_paths(cfg$dataset, data_root)
  cfg
}

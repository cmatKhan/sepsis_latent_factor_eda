# Reading and validating dataset configs (config/<dataset>_config.yml).
# A method runs when it appears under `methods:`; each method block holds
# overrides of that method's registry `defaults`. Validation runs when
# _targets.R is read, so a bad config fails before any job is submitted
# (docs: Data, "Dataset configs").

library(yaml)

if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

#' Read and validate a dataset config
#'
#' @param path Path to `config/<dataset>_config.yml`.
#' @param manifest Method manifest from discover_method_registries(), or
#'   `NULL` to skip the method checks in validate_dataset_metadata().
#' @param data_root Data root that relative data paths resolve against
#'   (default_data_root()).
#' @return The parsed config, data paths resolved, after validation.
read_dataset_metadata <- function(path, manifest = NULL, data_root = default_data_root()) {
  meta <- read_dataset_yaml(path, data_root)
  validate_dataset_metadata(meta, manifest)
  meta
}

#' Validate a parsed dataset config
#'
#' Checks that the matrix can be built (`matrix_path` or
#' `preprocessing_script` exists, or `expression_path` is set), that at least
#' one method is configured, that every configured method was discovered, and
#' that no `methods:` key is a YAML boolean token read as `TRUE`/`FALSE`
#' (docs: Data, "YAML keys that are booleans"). Errors on the first problem.
#'
#' @param meta A parsed config.
#' @param manifest Method manifest, or `NULL` to skip the method checks.
#' @return `meta`, invisibly.
validate_dataset_metadata <- function(meta, manifest = NULL) {
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

  if (is.null(meta$methods) || length(meta$methods) == 0) {
    stop("dataset_metadata.yml must have at least one method under `methods:` ",
         "(presence = enabled -- there is no `enabled:` flag)")
  }

  # Catches a real R `yaml` package footgun (confirmed, 2.3.12): an
  # UNQUOTED key that happens to be a YAML 1.1 legacy boolean token (y/Y/
  # yes/Yes/YES/n/N/no/No/NO/true/True/TRUE/false/False/FALSE/on/On/ON/
  # off/Off/OFF) parses as a literal TRUE/FALSE, not the string you wrote --
  # e.g. a real argument name `n` would silently become the key `FALSE`,
  # which downstream becomes the column name "FALSE." (R's make.names()
  # escaping the reserved word) and breaks the tool call with "unused
  # argument". Quote any such key in YAML (`"n": [...]`) -- see
  # docs/data.qmd ("YAML keys that are booleans"). Applied recursively (e.g. CoGAPS's
  # nested `params`/`distributed_params`/`run`), since there's no longer a
  # fixed wrapper depth to stop at.
  check_bool_key_typo_recursive <- function(block, label) {
    for (nm in names(block)) {
      if (nm %in% c("TRUE", "FALSE")) {
        stop(label, " has a key literally named `", nm, "` -- almost certainly an unquoted ",
             "YAML boolean-token key (y/n/yes/no/on/off/true/false, in any case) that got ",
             "coerced instead of staying a literal string. Quote it in the YAML, e.g. `\"n\": ",
             "[...]` instead of `n: [...]`. See R/README.md's config reference.")
      }
      val <- block[[nm]]
      if (is.list(val) && !is.null(names(val))) check_bool_key_typo_recursive(val, paste0(label, ".", nm))
    }
  }

  # Only checks "does this method exist" when a manifest is supplied (i.e.
  # from _targets.R, which sources R/methods/*.R before
  # calling read_dataset_metadata()) -- meaningless without it. Legacy
  # callers (manifest = NULL) skip that half; the boolean-token-key check
  # always runs regardless.
  for (nm in names(meta$methods)) {
    if (nm == "network") {
      for (backend in names(meta$methods$network)) {
        if (!is.null(manifest) && is.null(manifest[[backend]]$registry)) {
          stop("methods.network.", backend, " is configured but no such method was discovered under R/methods/")
        }
        check_bool_key_typo_recursive(meta$methods$network[[backend]], paste0("methods.network.", backend))
      }
    } else {
      if (!is.null(manifest) && is.null(manifest[[nm]]$registry)) {
        stop("methods.", nm, " is configured but no such method was discovered under R/methods/")
      }
      check_bool_key_typo_recursive(meta$methods[[nm]], paste0("methods.", nm))
    }
  }

  # Tensor methods (CP/Tucker) used to need dataset:-level subject/
  # timepoint columns here, driven by each method's registry
  # (`requires_subject_timepoint`, see R/lib/method_registry.R) rather
  # than a hardcoded name list -- this is why `manifest` is threaded
  # through to this function at all. Removed along with CP/Tucker
  # 2026-09-29 (see app.R's FACTORIZATION_METHODS header) -- no current
  # method sets `requires_subject_timepoint`, so this loop is a no-op
  # today, but left in place (rather than deleted) since it's generic,
  # registry-driven infrastructure a future method could still use.
  if (!is.null(manifest)) {
    check_subject_timepoint <- function(nm) {
      reg <- manifest[[nm]]$registry
      if (is.null(reg) || !isTRUE(reg$requires_subject_timepoint)) return(invisible())
      if (is.null(meta$dataset$subject_id_col) || is.null(meta$dataset$timepoint_col)) {
        stop("methods.", nm, " is configured but dataset.subject_id_col and ",
             "dataset.timepoint_col are not both set -- required to build the genes x ",
             "subjects x timepoints tensor. See R/README.md's config reference.")
      }
      if (is.null(meta$dataset$sample_metadata_path)) {
        stop("methods.", nm, " is configured but dataset.sample_metadata_path ",
             "is not set -- required to look up each sample's subject/timepoint")
      }
    }
    for (nm in setdiff(names(meta$methods), "network")) check_subject_timepoint(nm)
  }

  invisible(meta)
}

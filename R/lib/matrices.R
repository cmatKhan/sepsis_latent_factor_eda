# Input-matrix loading. The final matrix is resolved one of three ways, in
# this precedence order (see load_input_matrix() below):
#   1. `dataset.matrix_path` is set      -> read that RDS matrix as-is.
#   2. `dataset.preprocessing_script`    -> run your script (arbitrary
#      is set                               manipulation of the raw parquet
#                                           trio -- see R/README.md's
#                                           "Preprocessing script contract").
#   3. neither is set                    -> read `dataset.expression_path`
#                                           (long format) and pivot it wide,
#                                           with no other transformation.
# Deciding normalization/filtering is your call in all three cases except
# (3), which is deliberately a no-op transform (long -> wide only) for
# datasets that don't need anything fancier. There is deliberately no
# post-hoc sample/feature id filtering step here: if a matrix needs
# restricting to a sample or feature subset, do it inside a
# preprocessing_script (mode 2) or bake it into the matrix you point
# matrix_path at (mode 1) -- not via a separate id list applied after the
# fact.

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Run a user-supplied preprocessing script against the raw parquet trio.
#'
#' Contract (see R/README.md's "Preprocessing script contract"): the script
#' is sourced with `expression_path`, `sample_metadata_path`,
#' `feature_metadata_path`, `sample_id_col`, `feature_id_col`, `value_col`,
#' `ensembl_col`, and `tmp_dir` already defined in its evaluation
#' environment. It must write its final feature x sample matrix to
#' `file.path(tmp_dir, "matrix.rds")`. Any sample/feature filtering the
#' matrix needs belongs inside the script itself.
#'
#' `tmp_dir` lives under R's per-session temp directory, which is not
#' cleaned up until this R process exits -- long enough for
#' `submit_job_family()` to serialize the resulting matrix into each job's
#' input `.RData` later in the same script run. Don't delete it yourself.
run_preprocessing_script <- function(ds) {
  tmp_dir <- tempfile("dataset_preprocessing_")
  dir.create(tmp_dir, recursive = TRUE)

  env <- new.env(parent = globalenv())
  env$expression_path       <- ds$expression_path
  env$sample_metadata_path  <- ds$sample_metadata_path
  env$feature_metadata_path <- ds$feature_metadata_path
  env$sample_id_col         <- ds$sample_id_col %||% "sample_id"
  env$feature_id_col        <- ds$feature_id_col %||% "feature_id"
  env$value_col             <- ds$value_col %||% "value"
  env$ensembl_col           <- ds$ensembl_col %||% "ensembl"
  env$tmp_dir               <- tmp_dir

  source(ds$preprocessing_script, local = env)

  out_path <- file.path(tmp_dir, "matrix.rds")
  if (!file.exists(out_path)) {
    stop("preprocessing_script '", ds$preprocessing_script, "' did not write ",
         "the required output file ", out_path, " -- see R/README.md's ",
         "\"Preprocessing script contract\".")
  }
  readRDS(out_path)
}

#' Default (no preprocessing_script given) matrix builder: read the raw
#' long-format expression parquet and pivot it wide. No normalization, no
#' filtering -- just long -> wide.
pivot_expression_long <- function(ds) {
  if (!requireNamespace("arrow", quietly = TRUE)) {
    stop("Package 'arrow' is required to read parquet expression data")
  }
  sample_col  <- ds$sample_id_col %||% "sample_id"
  feature_col <- ds$feature_id_col %||% "feature_id"
  value_col   <- ds$value_col %||% "value"

  expr_long <- arrow::read_parquet(ds$expression_path)
  wide <- tidyr::pivot_wider(
    expr_long[, c(feature_col, sample_col, value_col)],
    id_cols = dplyr::all_of(feature_col),
    names_from = dplyr::all_of(sample_col),
    values_from = dplyr::all_of(value_col)
  )

  mat <- as.matrix(wide[, setdiff(names(wide), feature_col)])
  rownames(mat) <- wide[[feature_col]]
  mat
}

#' Load the dataset's single input matrix (feature x sample), resolved via
#' matrix_path -> preprocessing_script -> default long-to-wide pivot (see
#' file header).
load_input_matrix <- function(dataset_meta) {
  ds <- dataset_meta$dataset

  if (!is.null(ds$matrix_path)) {
    readRDS(ds$matrix_path)
  } else if (!is.null(ds$preprocessing_script)) {
    run_preprocessing_script(ds)
  } else {
    pivot_expression_long(ds)
  }
}

#' Shift a matrix's rows so every entry is non-negative -- required by
#' NMF/CoGAPS, not a dataset-specific normalization choice, so this stays a
#' generic transform applied inside the method's setup script rather than
#' baked into your input matrix.
shift_nonneg <- function(mat) mat - matrixStats::rowMins(mat)

#' Draw one fixed held-out mask (linear indices) for masking-CV, shared
#' across every parameter value scored so the comparison is fair.
make_mask_index <- function(mat, mask_frac, mask_seed) {
  set.seed(mask_seed)
  sample(length(mat), round(mask_frac * length(mat)))
}

# Building a dataset's input matrix: a prebuilt `matrix_path`, a
# `preprocessing_script`, or a plain long-to-wide pivot of `expression_path`,
# first match wins (docs: Data, "Building a matrix").

#' Null default
#'
#' `a` unless it is `NULL`, then `b` (base R's `%||%`, defined here for older R).
#'
#' @param a A value, possibly `NULL`.
#' @param b The fallback.
#' @return `a`, or `b` when `a` is `NULL`.
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Build a matrix with a dataset's preprocessing script
#'
#' Sources the script with `expression_path`, `sample_metadata_path`,
#' `feature_metadata_path`, `sample_id_col`, `feature_id_col`, `value_col`,
#' `ensembl_col` and `tmp_dir` defined; the script must write its feature x
#' sample matrix to `file.path(tmp_dir, "matrix.rds")` (docs: Data,
#' "Preprocessing-script contract").
#'
#' @param ds A resolved config's `dataset:` block.
#' @return The feature x sample matrix the script wrote.
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

#' Pivot long-format expression wide
#'
#' The default matrix builder: reads `expression_path` and pivots it to feature
#' x sample, with no normalization or filtering. Also used directly by the
#' differential-expression scripts to get raw values.
#'
#' @param ds A resolved config's `dataset:` block.
#' @return Feature x sample numeric matrix.
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

#' Load a dataset's input matrix
#'
#' Uses `matrix_path`, else `preprocessing_script`, else
#' pivot_expression_long() (the `mat_<dataset>` target).
#'
#' @param dataset_meta A config from read_dataset_metadata().
#' @return Feature x sample numeric matrix.
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

#' Shift a matrix to be non-negative
#'
#' Subtracts each row's minimum, as NMF and CoGAPS require (the
#' `mat_nn_<dataset>` target). An algorithm requirement, not a normalization.
#'
#' @param mat Feature x sample matrix.
#' @return The shifted matrix, every entry >= 0.
shift_nonneg <- function(mat) mat - matrixStats::rowMins(mat)

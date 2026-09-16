# Shared tensor-building helper for CP and Tucker decomposition (see
# R/methods/cp.R, R/methods/tucker.R) -- both need the exact
# same genes x subjects x timepoints array, built once and cached (see
# R/create_slurm_bundle.R's get_tensor()).
#
# rTensor's base cp()/tucker() require a COMPLETE (dense) array -- no
# native missing-entry handling -- so a subject missing any timepoint
# can't be included. Per your choice, incomplete subjects are dropped
# automatically (never silently: always reported).

library(rTensor)

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Reshape `mat` (features x samples, this framework's standard matrix)
#' into a features x subjects x timepoints rTensor `Tensor` object, using
#' `dataset.subject_id_col`/`dataset.timepoint_col` (new dataset:-level
#' config fields) to look up each sample's (subject, timepoint) via
#' `dataset.sample_metadata_path`.
#'
#' @param mat features x samples matrix (as returned by load_input_matrix())
#' @param dataset_meta the full dataset_metadata list (needs
#'   dataset$sample_metadata_path, dataset$sample_id_col,
#'   dataset$subject_id_col, dataset$timepoint_col)
#' @return an rTensor::Tensor (3-way), with dimnames set to list(gene =
#'   rownames(mat), subject = <kept subject ids>, time = <timepoint
#'   levels, sorted>)
build_tensor <- function(mat, dataset_meta) {
  ds <- dataset_meta$dataset
  subject_col <- ds$subject_id_col
  time_col <- ds$timepoint_col
  sample_col <- ds$sample_id_col %||% "sample_id"
  if (is.null(subject_col) || is.null(time_col)) {
    stop("dataset.subject_id_col and dataset.timepoint_col must both be set ",
         "to use cp/tucker (tensor methods) -- see R/README.md's config reference")
  }

  sm <- read_sample_metadata(ds$sample_metadata_path)
  for (col in c(sample_col, subject_col, time_col)) {
    if (!(col %in% names(sm))) stop("column '", col, "' not found in ", ds$sample_metadata_path)
  }
  sm <- sm[sm[[sample_col]] %in% colnames(mat), , drop = FALSE]

  time_levels <- sort(unique(as.character(sm[[time_col]])))
  subjects_all <- unique(as.character(sm[[subject_col]]))

  # keep only subjects with exactly one sample at EVERY timepoint level
  lookup <- list()
  complete <- character(0)
  incomplete <- list()
  for (subj in subjects_all) {
    rows <- sm[as.character(sm[[subject_col]]) == subj, , drop = FALSE]
    have <- as.character(rows[[time_col]])
    missing_tp <- setdiff(time_levels, have)
    dup_tp <- have[duplicated(have)]
    if (length(missing_tp) > 0 || length(dup_tp) > 0) {
      incomplete[[subj]] <- list(missing = missing_tp, duplicated = unique(dup_tp))
      next
    }
    complete <- c(complete, subj)
    for (tp in time_levels) {
      lookup[[paste(subj, tp, sep = "\r")]] <- rows[[sample_col]][have == tp][1]
    }
  }

  if (length(incomplete) > 0) {
    detail <- vapply(names(incomplete), function(s) {
      bits <- character(0)
      if (length(incomplete[[s]]$missing) > 0) bits <- c(bits, paste("missing", paste(incomplete[[s]]$missing, collapse = ",")))
      if (length(incomplete[[s]]$duplicated) > 0) bits <- c(bits, paste("duplicated", paste(incomplete[[s]]$duplicated, collapse = ",")))
      paste0(s, " (", paste(bits, collapse = "; "), ")")
    }, character(1))
    message("build_tensor(): dropping ", length(incomplete), " subject(s) without a complete, ",
            "unique time course: ", paste(detail, collapse = "; "))
  }
  if (length(complete) == 0) {
    stop("build_tensor(): no subject has a complete time course across all ",
         length(time_levels), " timepoint levels (", paste(time_levels, collapse = ", "), ")")
  }
  message("build_tensor(): keeping ", length(complete), " of ", length(subjects_all),
          " subjects x ", length(time_levels), " timepoints (", nrow(mat), " features)")

  arr <- array(NA_real_, dim = c(nrow(mat), length(complete), length(time_levels)),
               dimnames = list(gene = rownames(mat), subject = complete, time = time_levels))
  for (si in seq_along(complete)) {
    for (ti in seq_along(time_levels)) {
      sample_id <- lookup[[paste(complete[si], time_levels[ti], sep = "\r")]]
      arr[, si, ti] <- mat[, sample_id]
    }
  }
  rTensor::as.tensor(arr)
}

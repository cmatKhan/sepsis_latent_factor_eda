# projectR projection of one representative fit's loadings onto another
# dataset's matrix, run as one dynamic branch of projectr_<dataset>
# (row built by projection_row(), R/targets/projections.R).
# `include_intercept` is decided by the caller: off within a dataset family
# (config/dataset_families.yml), on across families.
#
# Both sides are remapped to Ensembl gene ids first -- the canonical
# cross-dataset id space (R/lib/ingest/symbol_mapping.R). For PCA, the fit
# is rebuilt as a `prcomp` object (rotation, the source matrix's gene means
# as `center`, per-component score SDs as `sdev`) and projected with the
# projectR fork's `center_by_loadings = TRUE`.

#' Project one fit's loadings onto another dataset's matrix
#'
#' Both sides are remapped to Ensembl ids first. A PCA fit is rebuilt as a
#' `prcomp` object (rotation, the source matrix's gene means as `center`,
#' score SDs as `sdev`) and projected with the projectR fork's
#' `center_by_loadings = TRUE`; other methods use plain loadings.
#'
#' @param source_fit_key The source fit's key.
#' @param source_method Its method.
#' @param source_dataset_id Its dataset.
#' @param loadings_file Absolute path to its loadings.
#' @param target_dataset_id The dataset projected onto.
#' @param target_matrix_file Absolute path to that dataset's matrix.
#' @param projection_type `"within_dataset"` or `"cross_dataset"`.
#' @param include_intercept Whether projectR fits an intercept (cross-family only).
#' @param source_ensembl_map,target_ensembl_map The two datasets' Ensembl maps.
#' @param source_mat_file,source_scores_file The source's matrix and scores
#'   (PCA only).
#' @return List with the job's identifiers, `n_genes_matched` and `result`
#'   (projectR's output, or `NULL` on failure).
run_projectr_job <- function(source_fit_key, source_method, source_dataset_id, loadings_file,
                              target_dataset_id, target_matrix_file,
                              projection_type, include_intercept,
                              source_ensembl_map, target_ensembl_map,
                              source_mat_file = NULL, source_scores_file = NULL) {
  library(projectR)
  loadings <- as.matrix(readRDS(loadings_file))
  target_mat <- as.matrix(readRDS(target_matrix_file))

  loadings_ens <- remap_to_ensembl(loadings, source_ensembl_map)
  target_ens   <- remap_to_ensembl(target_mat, target_ensembl_map)

  if (source_method == "pca" && !is.null(source_mat_file) && !is.null(source_scores_file)) {
    # Reconstruct the native prcomp object projectR's dispatch needs for
    # center_by_loadings/pvar -- sdev/center are cheaply recoverable from
    # the source matrix + scores (see the method-object audit earlier in
    # this project's history), so nothing extra needed to have been stored.
    source_mat <- as.matrix(readRDS(source_mat_file))
    scores <- as.matrix(readRDS(source_scores_file))
    # `center` must be remapped to the SAME id space as `loadings_ens`
    # (confirmed directly, 2026-09-16): projectR matches data/loadings by
    # rownames first, then looks up `loadings$center` by those SAME
    # (now-remapped) names for center_by_loadings -- a center vector still
    # indexed by the original probe/GeneID rownames fails with "gene(s) in
    # the matched data are absent from loadings$center" even though the
    # loadings themselves matched fine. Reuse remap_to_ensembl() (built for
    # matrices) via a 1-column matrix rather than duplicating its
    # probe-collapse logic.
    center_raw <- matrix(rowMeans(source_mat), ncol = 1, dimnames = list(rownames(source_mat), "center"))
    center_ens <- remap_to_ensembl(center_raw, source_ensembl_map)[, 1]
    pcfit <- structure(list(rotation = loadings_ens, center = center_ens[rownames(loadings_ens)],
                             sdev = apply(scores, 2, sd)), class = "prcomp")
    out <- tryCatch(projectR(data = target_ens, loadings = pcfit, full = TRUE, center_by_loadings = TRUE),
                     error = function(e) { message("projectR (pca) failed: ", conditionMessage(e)); NULL })
  } else {
    out <- tryCatch(projectR(data = target_ens, loadings = loadings_ens, full = TRUE,
                              include_intercept = include_intercept),
                     error = function(e) { message("projectR failed: ", conditionMessage(e)); NULL })
  }

  # Actual data/loadings intersection size -- NOT nrow(loadings_ens), which
  # is just the source's own total remapped gene count (confirmed
  # directly, 2026-09-16: projectR logs e.g. "200 row names matched
  # between data and loadings" internally, which can be smaller than
  # nrow(loadings_ens) whenever the target has fewer/different genes).
  n_genes_matched <- length(intersect(rownames(loadings_ens), rownames(target_ens)))

  list(source_fit_key = source_fit_key, source_dataset_id = source_dataset_id,
       target_dataset_id = target_dataset_id, method = source_method,
       projection_type = projection_type, include_intercept = include_intercept,
       n_genes_matched = n_genes_matched, result = out)
}

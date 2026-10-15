# projectR projection of one representative fit's loadings onto another
# dataset's matrix, run inside one dynamic branch of projectr_<dataset>
# (projection_rows(), R/targets/projections.R), which converts the target
# matrix to Ensembl ids once for all of its fits.
# `include_intercept` is decided by the caller: off within a dataset family
# (config/dataset_families.yml), on across families.
#
# Ensembl gene ids are the cross-dataset id space (R/lib/ingest/
# symbol_mapping.R). For PCA, the fit is rebuilt as a `prcomp` object
# (rotation, the source matrix's gene means as `center`, per-component score
# SDs as `sdev`) and projected with the projectR fork's
# `center_by_loadings = TRUE`.

#' Project one fit's loadings onto another dataset's matrix
#'
#' The loadings are converted to Ensembl ids with the source dataset's map;
#' the target arrives converted. A PCA fit is rebuilt as a `prcomp` object
#' (rotation, the source matrix's gene means as `center`, score SDs as
#' `sdev`) and projected with the projectR fork's `center_by_loadings = TRUE`;
#' other methods use plain loadings.
#'
#' @param source_fit_key The source fit's key.
#' @param source_method Its method.
#' @param source_dataset_id Its dataset.
#' @param loadings_file Absolute path to its loadings.
#' @param target_dataset_id The dataset projected onto.
#' @param target_ens That dataset's matrix, rows already Ensembl ids.
#' @param projection_type `"within_dataset"` or `"cross_dataset"`.
#' @param include_intercept Whether projectR fits an intercept (cross-family only).
#' @param source_ensembl_map The source dataset's Ensembl map.
#' @param source_mat_file,source_scores_file The source's matrix and scores
#'   (PCA only).
#' @return List with the job's identifiers, `n_genes_matched` and `result`
#'   (projectR's output, or `NULL` on failure).
run_projectr_job <- function(source_fit_key, source_method, source_dataset_id, loadings_file,
                             target_dataset_id, target_ens, projection_type, include_intercept,
                             source_ensembl_map, source_mat_file = NULL, source_scores_file = NULL) {
  library(projectR)
  loadings_ens <- remap_to_ensembl(as.matrix(readRDS(loadings_file)), source_ensembl_map)

  if (source_method == "pca" && !is.null(source_mat_file) && !is.null(source_scores_file)) {
    # sdev and center are recoverable from the source matrix and scores, so
    # nothing extra was stored at fit time. `center` must be in the same id
    # space as the loadings: projectR looks it up by the matched (Ensembl)
    # row names.
    source_mat <- as.matrix(readRDS(source_mat_file))
    scores <- as.matrix(readRDS(source_scores_file))
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

  # The genes both sides share (projectR uses this intersection), not the
  # source's own gene count.
  n_genes_matched <- length(intersect(rownames(loadings_ens), rownames(target_ens)))

  list(source_fit_key = source_fit_key, source_dataset_id = source_dataset_id,
       target_dataset_id = target_dataset_id, method = source_method,
       projection_type = projection_type, include_intercept = include_intercept,
       n_genes_matched = n_genes_matched, result = out)
}

# Phase-2 slurm ARRAY job (cheap relative to fgsea -- plain array
# parallelism, 1 cpu/task): one task per (source fit, target dataset)
# pair from build_projectr_pairs(). `include_intercept` flows straight
# through from the jobs_df row -- fixed at grid-build time by which of
# projectr_within_grid/projectr_cross_grid the row belongs to, never
# decided here.

run_projectr_job <- function(source_fit_id, source_method, source_dataset_id, loadings_file,
                              target_dataset_id, target_matrix_file,
                              projection_type, include_intercept,
                              source_mat_file = NULL, source_scores_file = NULL,
                              symbol_maps = list()) {
  library(projectR)
  loadings <- as.matrix(readRDS(loadings_file))
  target_mat <- as.matrix(readRDS(target_matrix_file))

  loadings_sym <- remap_to_symbol(loadings, symbol_maps[[source_dataset_id]])
  target_sym   <- remap_to_symbol(target_mat, symbol_maps[[target_dataset_id]])

  if (source_method == "pca" && !is.null(source_mat_file) && !is.null(source_scores_file)) {
    # Reconstruct the native prcomp object projectR's dispatch needs for
    # center_by_loadings/pvar -- sdev/center are cheaply recoverable from
    # the source matrix + scores (see the method-object audit earlier in
    # this project's history), so nothing extra needed to have been stored.
    source_mat <- as.matrix(readRDS(source_mat_file))
    scores <- as.matrix(readRDS(source_scores_file))
    # `center` must be remapped to the SAME symbol space as `loadings_sym`
    # (confirmed directly, 2026-09-16): projectR matches data/loadings by
    # rownames first, then looks up `loadings$center` by those SAME
    # (now-symbol) names for center_by_loadings -- a center vector still
    # indexed by the original probe/GeneID rownames fails with "gene(s) in
    # the matched data are absent from loadings$center" even though the
    # loadings themselves matched fine. Reuse remap_to_symbol() (built for
    # matrices) via a 1-column matrix rather than duplicating its
    # probe-collapse logic.
    center_raw <- matrix(rowMeans(source_mat), ncol = 1, dimnames = list(rownames(source_mat), "center"))
    center_sym <- remap_to_symbol(center_raw, symbol_maps[[source_dataset_id]])[, 1]
    pcfit <- structure(list(rotation = loadings_sym, center = center_sym[rownames(loadings_sym)],
                             sdev = apply(scores, 2, sd)), class = "prcomp")
    out <- tryCatch(projectR(data = target_sym, loadings = pcfit, full = TRUE, center_by_loadings = TRUE),
                     error = function(e) { message("projectR (pca) failed: ", conditionMessage(e)); NULL })
  } else {
    out <- tryCatch(projectR(data = target_sym, loadings = loadings_sym, full = TRUE,
                              include_intercept = include_intercept),
                     error = function(e) { message("projectR failed: ", conditionMessage(e)); NULL })
  }

  # Actual data/loadings intersection size -- NOT nrow(loadings_sym), which
  # is just the source's own total remapped gene count (confirmed
  # directly, 2026-09-16: projectR logs e.g. "200 row names matched
  # between data and loadings" internally, which can be smaller than
  # nrow(loadings_sym) whenever the target has fewer/different genes).
  n_genes_matched <- length(intersect(rownames(loadings_sym), rownames(target_sym)))

  list(source_fit_id = source_fit_id, source_dataset_id = source_dataset_id,
       target_dataset_id = target_dataset_id, method = source_method,
       projection_type = projection_type, include_intercept = include_intercept,
       n_genes_matched = n_genes_matched, result = out)
}

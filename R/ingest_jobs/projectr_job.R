# Phase-2 slurm ARRAY job (cheap relative to fgsea -- plain array
# parallelism, 1 cpu/task): one task per (source fit, target dataset)
# pair from build_projectr_pairs(). `include_intercept` flows straight
# through from the jobs_df row -- fixed at grid-build time by which of
# projectr_within_grid/projectr_cross_grid the row belongs to, never
# decided here.
#
# Remapping to Ensembl (not gene symbol) is the canonical cross-dataset
# identifier space here -- see R/lib/ingest/symbol_mapping.R's header for
# why (symbol is display-only).
#
# `ensembl_maps` is read as a FREE VARIABLE, NOT a function parameter --
# this is an ARRAY job (slurm_apply(), jobs_df = one row per source-fit x
# target-dataset pair), and `ensembl_maps` is baked in via
# `global_objects` only, never as a jobs_df column (there's exactly ONE
# shared copy for the whole grid, not a per-row value). An EARLIER version
# of this function declared it as a formal parameter instead (named
# `symbol_maps`, with a `= list()` default) -- since jobs_df never
# supplies a column of that name, `do.call(f, params_row)` never
# overrides that default, so the function ALWAYS silently received an
# empty list (confirmed directly), never the real baked-in map. Every
# projectr_job task was therefore comparing raw, un-remapped native
# platform ids the whole time -- almost certainly the actual cause of
# "0 row names matched between data and loadings" on cross-dataset
# pairs using different platforms, not (only) a symbol-column-detection
# issue. Compare `pathways` in run_fgsea_job() (R/ingest_jobs/
# fgsea_job.R), which was already done correctly this way.
run_projectr_job <- function(source_fit_id, source_method, source_dataset_id, loadings_file,
                              target_dataset_id, target_matrix_file,
                              projection_type, include_intercept,
                              source_mat_file = NULL, source_scores_file = NULL) {
  library(projectR)
  loadings <- as.matrix(readRDS(loadings_file))
  target_mat <- as.matrix(readRDS(target_matrix_file))

  loadings_ens <- remap_to_ensembl(loadings, ensembl_maps[[source_dataset_id]])
  target_ens   <- remap_to_ensembl(target_mat, ensembl_maps[[target_dataset_id]])

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
    center_ens <- remap_to_ensembl(center_raw, ensembl_maps[[source_dataset_id]])[, 1]
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

  list(source_fit_id = source_fit_id, source_dataset_id = source_dataset_id,
       target_dataset_id = target_dataset_id, method = source_method,
       projection_type = projection_type, include_intercept = include_intercept,
       n_genes_matched = n_genes_matched, result = out)
}

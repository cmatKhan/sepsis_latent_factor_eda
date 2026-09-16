# Phase-2 slurm ARRAY job: one task per (dataset, method, representative
# fit) -- parallelized across the array (nodes) plus BiocParallel workers
# within each task (cpus_per_task). NOT run across seeds -- jobs_df is
# built from representative_fit_ids() (see R/lib/ingest/redundancy.R),
# which already collapses to one best-seed fit per rank (or PCA's single
# max-rank fit). `pathways` (msigdbr, built once) and `symbol_maps` (one
# per dataset) are baked in as global objects at grid-build time -- no DB
# access from compute nodes.
#
# CoGAPS fits additionally get a fora()-based overrepresentation test
# against their already-computed patternMarkers() gene sets (free reuse
# of R/lib/ingest/redundancy.R's output) -- CoGAPS::getPatternGeneSet()
# itself isn't called directly since it needs the CogapsResult's own
# (platform-native) rownames; we already have the symbol-mapped loadings
# and marker gene lists in hand, so fora() is called directly instead.

run_fgsea_job <- function(dataset_id, method, fit_id, loadings_file,
                          symbol_map = NULL, cogaps_marker_genes = NULL, cpus_per_task = 8) {
  library(fgsea); library(BiocParallel)
  register(MulticoreParam(cpus_per_task))

  L <- as.matrix(readRDS(loadings_file))
  L_sym <- remap_to_symbol(L, symbol_map)

  gsea_results <- lapply(seq_len(ncol(L_sym)), function(fi) {
    ranks <- sort(L_sym[, fi], decreasing = TRUE)
    list(factor_index = fi,
         result = tryCatch(fgsea(pathways = pathways, stats = ranks, minSize = 10, maxSize = 500),
                            error = function(e) NULL))
  })

  fora_results <- NULL
  if (method == "cogaps" && !is.null(cogaps_marker_genes)) {
    fora_results <- lapply(cogaps_marker_genes, function(genes) {
      tryCatch(fora(pathways = pathways, genes = genes, universe = rownames(L_sym),
                     minSize = 10, maxSize = 500),
               error = function(e) NULL)
    })
  }

  list(dataset_id = dataset_id, method = method, fit_id = fit_id,
       gsea = gsea_results, fora = fora_results)
}

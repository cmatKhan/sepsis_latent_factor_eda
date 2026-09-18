# Phase-2 slurm ARRAY job: one task per (dataset, method, representative
# fit) -- parallelized across the array (nodes) plus BiocParallel workers
# within each task (cpus_per_task). NOT run across seeds -- jobs_df is
# built from representative_fit_ids() (see R/lib/ingest/redundancy.R),
# which already collapses to one best-seed fit per rank (or PCA's single
# max-rank fit). `pathways` (msigdbr, Ensembl-gene-keyed, built once) and
# `ensembl_map` (one per dataset) are baked in as global objects /
# per-row params at grid-build time -- no DB access from compute nodes.
# Remapping to Ensembl (not gene symbol) is the canonical cross-dataset
# identifier space here -- see R/lib/ingest/symbol_mapping.R's header for
# why.
#
# CoGAPS fits additionally get a fora()-based overrepresentation test
# against their already-computed patternMarkers() gene sets (free reuse
# of R/lib/ingest/redundancy.R's output) -- CoGAPS::getPatternGeneSet()
# itself isn't called directly since it needs the CogapsResult's own
# (platform-native) rownames; we already have the Ensembl-remapped
# loadings and marker gene lists in hand, so fora() is called directly
# instead.

run_fgsea_job <- function(dataset_id, method, fit_id, loadings_file,
                          ensembl_map = NULL, cogaps_marker_genes = NULL, cpus_per_task = 8) {
  library(fgsea); library(BiocParallel)
  register(MulticoreParam(cpus_per_task))

  L <- as.matrix(readRDS(loadings_file))
  L_ens <- remap_to_ensembl(L, ensembl_map)

  gsea_results <- lapply(seq_len(ncol(L_ens)), function(fi) {
    ranks <- sort(L_ens[, fi], decreasing = TRUE)
    list(factor_index = fi,
         result = tryCatch(fgsea(pathways = pathways, stats = ranks, minSize = 10, maxSize = 500),
                            error = function(e) NULL))
  })

  fora_results <- NULL
  if (method == "cogaps" && !is.null(cogaps_marker_genes)) {
    fora_results <- lapply(cogaps_marker_genes, function(genes) {
      # cogaps_marker_genes come from pattern_markers (R/lib/ingest/
      # redundancy.R), which is computed against the CogapsResult's OWN
      # native rownames -- NEVER remapped. Must go through the SAME
      # ensembl_map as the loadings/universe below, or `genes` and
      # `universe` are silently in mismatched id spaces (genes would
      # almost never actually be found in the universe).
      genes_ens <- if (!is.null(ensembl_map)) {
        mapped <- ensembl_map[genes]
        unique(mapped[!is.na(mapped) & nzchar(mapped)])
      } else {
        genes
      }
      tryCatch(fora(pathways = pathways, genes = genes_ens, universe = rownames(L_ens),
                     minSize = 10, maxSize = 500),
               error = function(e) NULL)
    })
  }

  list(dataset_id = dataset_id, method = method, fit_id = fit_id,
       gsea = gsea_results, fora = fora_results)
}

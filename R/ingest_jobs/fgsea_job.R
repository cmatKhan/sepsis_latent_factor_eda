# Phase-2 slurm ARRAY job: one task per (dataset, method, representative
# fit) -- parallelized across the array (nodes) plus BiocParallel workers
# within each task (cpus_per_task). NOT run across seeds -- jobs_df is
# built from representative_fit_ids() (see R/lib/ingest/redundancy.R),
# which already collapses to one best-seed fit per rank (or PCA's single
# max-rank fit). `pathways` (Hallmark only, msigdbr, Ensembl-gene-keyed,
# built once -- unchanged since before 2026-09-19) and `ensembl_map` (one
# per dataset) are baked in as global objects / per-row params at
# grid-build time -- no DB access from compute nodes. Remapping to Ensembl
# (not gene symbol) is the canonical cross-dataset identifier space here --
# see R/lib/ingest/symbol_mapping.R's header for why.
#
# CoGAPS fits additionally get a fora()-based overrepresentation test
# against their already-computed patternMarkers() gene sets (free reuse
# of R/lib/ingest/redundancy.R's output) -- CoGAPS::getPatternGeneSet()
# itself isn't called directly since it needs the CogapsResult's own
# (platform-native) rownames; we already have the Ensembl-remapped
# loadings and marker gene lists in hand, so fora() is called directly
# instead.
#
# NEW (2026-09-19): `run_local_enrichment()` below adds a local,
# fully-parallel replacement for the retired gprofiler_grid job family,
# which used to make one live gprofiler2::gost() HTTP request per (fit,
# factor, direction) combination against g:Profiler's public API --
# 47,517 such combinations for this project's grid, strictly serial (API
# rate limits), and NOT parallelizable the way this array job is.
# gprofiler2 itself has no local/offline mode (confirmed: every function
# that queries annotation data goes through `gprofiler_request()`, a
# plain HTTP POST -- see gprofiler2's own docs/vignette -- and g:Profiler's
# FAQ explicitly says they don't offer a self-hosted/Docker instance for
# high query volumes). `pathways_by_source` (built alongside `pathways` in
# R/create_ingest_slurm_bundle.R, ALSO baked in as a global object) is a
# named list of msigdbr collections -- HALLMARK/GO:BP/GO:MF/KEGG/REAC/WP --
# in the SAME Ensembl-gene id space as `pathways`/the remapped loadings,
# letting `fgsea()` (preranked GSEA) and `fora()` (hypergeometric ORA)
# reproduce gprofiler's GSEA/ORA modes entirely locally:
#   - GSEA: same preranked approach as the existing Hallmark-only pass
#     below, just repeated per collection (skipping HALLMARK there since
#     the existing `gsea_results` pass already covers it under the
#     pre-existing `query_type = "fgsea"` rows -- see
#     R/ingest_enrichment_results.R).
#   - ORA: gprofiler's ORA mode tested the TOP-100 genes by loading value
#     (both positive- and negative-loading ends for pca/ica/spca, whose
#     loadings are signed; just the positive end for nmf/cogaps/spca's
#     non-negative loadings) against each source -- reproduced here via
#     the same top-100 selection `gprofiler_job.R` used to do, now
#     followed by a local `fora()` call per collection instead of a
#     remote gost() call.
# The app's own per-factor, on-demand gprofiler2::gost() calls (app/app.R)
# are a SEPARATE, still-live code path, deliberately kept -- see
# R/ingest_enrichment_results.R's header for how both write into the same
# `enrichment_cache` schema.

run_fgsea_job <- function(dataset_id, method, fit_id, loadings_file,
                          ensembl_map = NULL, cogaps_marker_genes = NULL, cpus_per_task = 8) {
  library(fgsea); library(BiocParallel)
  register(MulticoreParam(cpus_per_task))

  L <- as.matrix(readRDS(loadings_file))
  L_ens <- remap_to_ensembl(L, ensembl_map)
  n_factors <- ncol(L_ens)

  gsea_results <- lapply(seq_len(n_factors), function(fi) {
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

  # ---- local ORA/GSEA replacing gprofiler_grid (see this file's header) ----
  local_gsea_results <- list()
  local_ora_results <- list()
  if (exists("pathways_by_source", inherits = TRUE) && length(pathways_by_source) > 0) {
    gsea_sources <- setdiff(names(pathways_by_source), "HALLMARK")   # HALLMARK GSEA already covered above
    dirs <- if (method %in% c("pca", "ica", "spca")) c("pos", "neg") else "pos"

    for (fi in seq_len(n_factors)) {
      ranks <- sort(L_ens[, fi], decreasing = TRUE)
      for (src in gsea_sources) {
        res <- tryCatch(fgsea(pathways = pathways_by_source[[src]], stats = ranks, minSize = 10, maxSize = 500),
                         error = function(e) NULL)
        local_gsea_results[[length(local_gsea_results) + 1]] <- list(factor_index = fi, source = src, result = res)
      }

      v <- L_ens[, fi]
      for (dir_i in dirs) {
        genes <- if (dir_i == "neg") names(sort(v))[seq_len(min(100, length(v)))]
                 else names(sort(v, decreasing = TRUE))[seq_len(min(100, length(v)))]
        for (src in names(pathways_by_source)) {
          res <- tryCatch(fora(pathways = pathways_by_source[[src]], genes = genes, universe = rownames(L_ens),
                                minSize = 10, maxSize = 500),
                           error = function(e) NULL)
          local_ora_results[[length(local_ora_results) + 1]] <-
            list(factor_index = fi, direction = dir_i, source = src, n_genes = length(genes), result = res)
        }
      }
    }
  }

  list(dataset_id = dataset_id, method = method, fit_id = fit_id,
       gsea = gsea_results, fora = fora_results,
       local_gsea = local_gsea_results, local_ora = local_ora_results)
}

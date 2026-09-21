# Phase-2 slurm ARRAY job: one task per (dataset, method, representative
# fit) -- parallelized across the array (nodes), plus `cpus_per_task`
# BiocParallel workers WITHIN each task, dispatched across that fit's
# flattened factor x source(/direction) work list (see run_fgsea_job()'s
# own doc for why it's structured that way rather than relying on each
# individual fgsea()/fora() call's own internal parallelism). NOT run
# across seeds -- jobs_df is
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
# The app's own per-factor, on-demand queries (app/app.R) used to be a
# SEPARATE gprofiler2::gost() code path, then briefly a local fora()/
# fgsea() on-demand equivalent, then retired ENTIRELY -- the app has no
# live-compute enrichment path at all now, it's a pure read-only viewer of
# whatever this job family (+ R/ingest_jobs/wgcna_ora_job.R for WGCNA) has
# already computed. See R/ingest_enrichment_results.R's header for how
# both write into the same `enrichment_cache` schema.

# Logging: one line per PHASE per fit (start, Hallmark GSEA, cogaps ORA
# when applicable, local GSEA, local ORA, done-with-elapsed) -- NOT
# per-factor or per-source, which would be hundreds of lines/fit (this
# runs as a batched array task -- see submit_job_family()'s `max_array_size`
# doc -- so a single task's log can cover many fits' worth of these lines).
# `n_sig()` below just counts padj < 0.05 rows across a list of fgsea()/
# fora() result data.frames, reused by every phase's summary line.
n_sig <- function(results_list) {
  sum(vapply(results_list, function(r) {
    res <- if (is.list(r) && "result" %in% names(r)) r$result else r
    if (is.null(res) || nrow(res) == 0) return(0L)
    sum(res$padj < 0.05, na.rm = TRUE)
  }, integer(1)))
}

#' fgsea::collapsePathways()'s pathway-redundancy reduction on one fgsea()
#' call's significant subset -- run HERE (ingest time), not later on demand,
#' because it needs the full fgsea result object (specifically `ES`) plus
#' the exact `pathways`/`stats` used, none of which are cheaply
#' reconstructable from enrichment_cache's summary columns alone (see
#' R/lib/ingest/db.R's enrichment_cache.is_main_pathway comment). Same
#' NA-safe padj filter as R/ingest_enrichment_results.R's `sig <- res[...]`
#' lines (fgsea legitimately returns NA padj for "unbalanced" pathways --
#' `res$padj < 0.05` alone would keep those as spurious all-NA rows).
#' Returns a character vector of "main" (non-redundant) pathway names, or
#' NULL when there's nothing meaningful to collapse (0 or 1 significant
#' pathways) or collapsePathways() itself errors -- never fails the
#' enclosing fgsea job over this.
main_pathways_for <- function(res, pathways, ranks) {
  if (is.null(res) || nrow(res) == 0) return(NULL)
  sig <- res[!is.na(res$padj) & res$padj < 0.05, ]
  if (nrow(sig) == 0) return(NULL)
  if (nrow(sig) == 1) return(sig$pathway)
  tryCatch(collapsePathways(sig, pathways, ranks)$mainPathways, error = function(e) NULL)
}

#' Registering MulticoreParam(cpus_per_task) alone is NOT enough to
#' actually use `cpus_per_task` cores: `fora()` doesn't accept a BPPARAM
#' argument at all (confirmed directly against the installed fgsea package
#' -- its body has no BiocParallel dispatch, it's always single-threaded),
#' and fgseaMultilevel()'s own internal parallelism only splits work
#' WITHIN one fgsea() call, which is a weak effect for anything smaller
#' than the biggest collection (GO:BP) -- the outer loop that calls
#' fgsea()/fora() hundreds of times per fit (once per factor, or per
#' factor x source, or per factor x direction x source) was, and remains
#' here, an ordinary sequential R for-loop, so those hundreds of small
#' single-threaded calls never overlap. Confirmed on a real cluster run:
#' 23.46% CPU efficiency on an 8-core task (seff: 01:03:58 CPU time over
#' 00:34:05 wall-clock). The fix is to parallelize the OUTER loop instead
#' -- flatten factor/source/direction combinations into one job list and
#' dispatch it with bplapply(), forcing each INNER fgsea() call to
#' BPPARAM = SerialParam() so it doesn't also try to fork out to
#' `cpus_per_task` workers from inside an already-forked worker (nested
#' MulticoreParam forking oversubscribes cores and can be slower than no
#' parallelism at all, or error outright in restricted environments).
#' fora() needs no such override -- it was never going to parallelize
#' either way.
run_fgsea_job <- function(dataset_id, method, fit_id, loadings_file,
                          ensembl_map = NULL, cogaps_marker_genes = NULL, cpus_per_task = 8) {
  library(fgsea); library(BiocParallel)
  bp <- MulticoreParam(cpus_per_task)
  register(bp)
  t0 <- Sys.time()
  tag <- sprintf("[fgsea_grid] %s/%s fit=%d", dataset_id, method, fit_id)

  L <- as.matrix(readRDS(loadings_file))
  L_ens <- remap_to_ensembl(L, ensembl_map)
  n_factors <- ncol(L_ens)
  message(sprintf("%s: starting -- %d factor(s), %d genes post-remap", tag, n_factors, nrow(L_ens)))

  # Hallmark GSEA: one fgsea() call per factor, dispatched across
  # `cpus_per_task` workers instead of relying on each individual call's
  # own (weak, for a 50-pathway collection) internal parallelism.
  gsea_results <- bplapply(seq_len(n_factors), function(fi) {
    ranks <- sort(L_ens[, fi], decreasing = TRUE)
    res <- tryCatch(fgsea(pathways = pathways, stats = ranks, minSize = 10, maxSize = 500,
                           BPPARAM = SerialParam()),
                     error = function(e) NULL)
    list(factor_index = fi, result = res, main_pathways = main_pathways_for(res, pathways, ranks))
  }, BPPARAM = bp)
  message(sprintf("%s: Hallmark GSEA done -- %d significant term-hit(s) across %d factor(s)",
                   tag, n_sig(gsea_results), n_factors))

  fora_results <- NULL
  if (method == "cogaps" && !is.null(cogaps_marker_genes)) {
    fora_results <- bplapply(cogaps_marker_genes, function(genes) {
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
    }, BPPARAM = bp)
    message(sprintf("%s: CoGAPS marker-gene ORA done -- %d significant term-hit(s) across %d pattern(s)",
                     tag, n_sig(fora_results), length(fora_results)))
  }

  # ---- local ORA/GSEA replacing gprofiler_grid (see this file's header) ----
  # Both flattened into ONE combined (factor x source[/direction]) work
  # list each, dispatched as a single bplapply() across `cpus_per_task`
  # workers -- this is where most of a fit's wall-clock time goes (up to
  # n_factors x 5 GSEA calls + n_factors x 2 x 6 ORA calls), so it's where
  # actually using all requested cores matters most.
  local_gsea_results <- list()
  local_ora_results <- list()
  if (exists("pathways_by_source", inherits = TRUE) && length(pathways_by_source) > 0) {
    gsea_sources <- setdiff(names(pathways_by_source), "HALLMARK")   # HALLMARK GSEA already covered above
    dirs <- if (method %in% c("pca", "ica", "spca")) c("pos", "neg") else "pos"

    gsea_jobs <- expand.grid(fi = seq_len(n_factors), src = gsea_sources, stringsAsFactors = FALSE)
    local_gsea_results <- bplapply(seq_len(nrow(gsea_jobs)), function(i) {
      fi <- gsea_jobs$fi[i]; src <- gsea_jobs$src[i]
      ranks <- sort(L_ens[, fi], decreasing = TRUE)
      res <- tryCatch(fgsea(pathways = pathways_by_source[[src]], stats = ranks, minSize = 10, maxSize = 500,
                             BPPARAM = SerialParam()),
                       error = function(e) NULL)
      list(factor_index = fi, source = src, result = res,
           main_pathways = main_pathways_for(res, pathways_by_source[[src]], ranks))
    }, BPPARAM = bp)

    ora_jobs <- expand.grid(fi = seq_len(n_factors), dir_i = dirs, src = names(pathways_by_source),
                             stringsAsFactors = FALSE)
    local_ora_results <- bplapply(seq_len(nrow(ora_jobs)), function(i) {
      fi <- ora_jobs$fi[i]; dir_i <- ora_jobs$dir_i[i]; src <- ora_jobs$src[i]
      v <- L_ens[, fi]
      genes <- if (dir_i == "neg") names(sort(v))[seq_len(min(100, length(v)))]
               else names(sort(v, decreasing = TRUE))[seq_len(min(100, length(v)))]
      res <- tryCatch(fora(pathways = pathways_by_source[[src]], genes = genes, universe = rownames(L_ens),
                            minSize = 10, maxSize = 500),
                       error = function(e) NULL)
      list(factor_index = fi, direction = dir_i, source = src, n_genes = length(genes), result = res)
    }, BPPARAM = bp)

    message(sprintf("%s: local GSEA done -- %d source(s) x %d factor(s), %d significant term-hit(s)",
                     tag, length(gsea_sources), n_factors, n_sig(local_gsea_results)))
    message(sprintf("%s: local ORA done -- %d source(s) x %d direction(s) x %d factor(s), %d significant term-hit(s)",
                     tag, length(pathways_by_source), length(dirs), n_factors, n_sig(local_ora_results)))
  } else {
    message(sprintf("%s: no `pathways_by_source` global found -- skipping local GSEA/ORA", tag))
  }

  message(sprintf("%s: done in %.1fs", tag, as.numeric(difftime(Sys.time(), t0, units = "secs"))))

  list(dataset_id = dataset_id, method = method, fit_id = fit_id,
       gsea = gsea_results, fora = fora_results,
       local_gsea = local_gsea_results, local_ora = local_ora_results)
}

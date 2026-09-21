# Phase-2 slurm ARRAY job: one task per WGCNA fit -- computes local fora()
# ORA (same msigdbr collections as fgsea_grid's local ORA) for every
# module's full gene-membership list. WGCNA modules have no continuous
# ranking (no loadings to sort by), so there's no GSEA equivalent here --
# ORA only.
#
# Added when the app's live-compute enrichment path was retired
# (app/app.R's old run_wgcna_enrich_all_modules_core(), which called
# gprofiler2::gost() and later a local fora() call on demand) -- WGCNA was
# never part of fgsea_grid's scope (representative_fit_ids() has no
# ranking concept to select "representative" fits by there, and modules
# aren't loadings), so it gets its own job family here instead, following
# the exact same "compute in isolation, write your own results_*.RDS,
# ingest separately on the login node" pattern fgsea_grid/projectr_*_grid
# use (see R/ingest_jobs/fgsea_job.R's header for why: this is a slurm
# ARRAY job, and SQLite only tolerates one writer at a time -- see
# R/ingest_jobs/ingest_core_job.R's header).
#
# `module_genes` (named list, module number [as character] -> native gene
# ids, module 0 "unassigned" excluded -- see below) and `universe_genes`
# (every gene assigned to ANY module, INCLUDING module 0 -- it's still
# real tested network genes) are baked in as per-row jobs_df values at
# grid-build time (R/create_ingest_slurm_bundle.R) -- no DB access from
# compute nodes, same convention as fgsea_grid's cogaps_marker_genes.
# `pathways_by_source` is a shared global_object, the same msigdbr
# collections fgsea_grid uses.
#
# Module 0 is never itself queried as an ORA gene set: ingest_dataset.R
# never gives it a `factors` row either (`mod_ids <- setdiff(sort(unique(
# ext$modules$module)), 0L)`), so there'd be nowhere in enrichment_cache
# to attach a module-0 result even if computed.

run_wgcna_ora_job <- function(dataset_id, fit_id, module_genes, universe_genes,
                               ensembl_map = NULL, cpus_per_task = 4) {
  library(fgsea); library(BiocParallel)
  bp <- MulticoreParam(cpus_per_task)
  register(bp)
  t0 <- Sys.time()
  tag <- sprintf("[wgcna_ora_grid] %s fit=%d", dataset_id, fit_id)

  remap <- function(genes) {
    if (is.null(ensembl_map)) return(genes)
    mapped <- ensembl_map[genes]
    unique(mapped[!is.na(mapped) & nzchar(mapped)])
  }
  universe_ens <- remap(universe_genes)
  modules <- names(module_genes)
  message(sprintf("%s: starting -- %d module(s), %d universe genes post-remap",
                   tag, length(modules), length(universe_ens)))

  # Flattened (module x source) work list dispatched as one bplapply() --
  # same rationale as fgsea_job.R's restructuring: fora() has no BPPARAM
  # support at all, so the only way to actually use `cpus_per_task` cores
  # is to parallelize ACROSS these many small calls, not rely on any one
  # of them being multi-threaded internally.
  jobs <- expand.grid(module = modules, src = names(pathways_by_source), stringsAsFactors = FALSE)
  results <- bplapply(seq_len(nrow(jobs)), function(i) {
    mod <- jobs$module[i]; src <- jobs$src[i]
    genes_ens <- remap(module_genes[[mod]])
    res <- tryCatch(fora(pathways = pathways_by_source[[src]], genes = genes_ens, universe = universe_ens,
                          minSize = 10, maxSize = 500),
                     error = function(e) NULL)
    list(module = as.integer(mod), source = src, n_genes = length(genes_ens), result = res)
  }, BPPARAM = bp)

  n_sig_total <- sum(vapply(results, function(r) {
    if (is.null(r$result) || nrow(r$result) == 0) return(0L)
    sum(!is.na(r$result$padj) & r$result$padj < 0.05)
  }, integer(1)))
  message(sprintf("%s: done in %.1fs -- %d module(s) x %d source(s), %d significant term-hit(s)",
                   tag, as.numeric(difftime(Sys.time(), t0, units = "secs")),
                   length(modules), length(pathways_by_source), n_sig_total))

  list(dataset_id = dataset_id, fit_id = fit_id, ora = results)
}

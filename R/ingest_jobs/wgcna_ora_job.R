# Module ORA for one WGCNA fit, run as one dynamic branch of
# wgcna_ora_<dataset> (rows built by wgcna_ora_rows(),
# R/targets/enrichment.R): fgsea::fora() of every module's gene list
# against each MSigDB collection. WGCNA modules have no continuous ranking,
# so there is no GSEA equivalent -- ORA only.
#
# `module_genes`: module number (as character) -> native gene ids, module 0
# ("unassigned") excluded -- it has no `factors` row to attach results to.
# `universe_genes`: every gene in the network, module 0 included. Both are
# remapped to Ensembl ids before testing.

#' Module ORA for one WGCNA fit
#'
#' fgsea::fora() of every module's genes against each MSigDB collection, over
#' the network's genes (module 0 included in the universe, not tested).
#' Modules have no continuous ranking, so there is no GSEA here.
#'
#' @param dataset_id Dataset id (for log lines).
#' @param fit_key The fit's key.
#' @param module_genes Named list (by module number) of native gene ids.
#' @param universe_genes Every gene in the network.
#' @param pathways_by_source fetch_msigdb_pathways() output.
#' @param ensembl_map The dataset's Ensembl map.
#' @param cpus_per_task Worker CPUs.
#' @return `list(dataset_id, fit_key, ora)`, one ORA result per (module, collection).
run_wgcna_ora_job <- function(dataset_id, fit_key, module_genes, universe_genes, pathways_by_source,
                               ensembl_map = NULL, cpus_per_task = 4) {
  library(fgsea); library(BiocParallel)
  bp <- MulticoreParam(cpus_per_task)
  register(bp)
  t0 <- Sys.time()
  tag <- sprintf("[wgcna_ora] %s %s", dataset_id, fit_key)

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

  list(dataset_id = dataset_id, fit_key = fit_key, ora = results)
}

# Enrichment for one representative fit, run as one fgsea_<method>_<dataset>
# branch (R/targets/factories.R):
#   - Hallmark preranked GSEA per factor, with collapsePathways() marking the
#     non-redundant "main" pathways;
#   - for CoGAPS, ORA of each pattern's marker genes;
#   - GSEA per factor against the other MSigDB collections (GO:BP, GO:MF,
#     KEGG, REACTOME, WikiPathways), and ORA of each factor's top 100 genes
#     (both ends for sign-ambiguous methods) against every collection.
# Loadings are remapped to Ensembl gene ids first. The result is trimmed to
# significant terms by trim_enrichment_result() and turned into DB rows by
# enrichment_rows() (R/targets/enrichment.R).

#' Count significant terms across results (for log lines)
#'
#' @param results_list List of fgsea()/fora() tables, or lists with a `result` table.
#' @return Number of rows with padj < 0.05.
n_sig <- function(results_list) {
  sum(vapply(results_list, function(r) {
    res <- if (is.list(r) && "result" %in% names(r)) r$result else r
    if (is.null(res) || nrow(res) == 0) return(0L)
    sum(res$padj < 0.05, na.rm = TRUE)
  }, integer(1)))
}

#' Non-redundant main pathways of one fgsea() result
#'
#' fgsea::collapsePathways() on the significant subset. Run here because it
#' needs the full result, pathways and ranking.
#'
#' @param res An fgsea() table.
#' @param pathways The gene sets tested.
#' @param ranks The ranking statistic.
#' @return Main pathway names, or `NULL` when nothing is significant or the
#'   collapse fails.
main_pathways_for <- function(res, pathways, ranks) {
  if (is.null(res) || nrow(res) == 0) return(NULL)
  sig <- res[!is.na(res$padj) & res$padj < 0.05, ]
  if (nrow(sig) == 0) return(NULL)
  if (nrow(sig) == 1) return(sig$pathway)
  tryCatch(collapsePathways(sig, pathways, ranks)$mainPathways, error = function(e) NULL)
}

#' Enrichment for one representative fit
#'
#' Parallelized across factors and collections with BiocParallel (the many
#' small fgsea()/fora() calls are the cost; each inner call runs serially to
#' avoid nested forking).
#'
#' @param dataset_id Dataset id (for log lines).
#' @param method Method name.
#' @param fit_key The fit's key.
#' @param loadings_file Absolute path to its loadings.
#' @param pathways_by_source fetch_msigdb_pathways() output.
#' @param ensembl_map The dataset's Ensembl map.
#' @param cogaps_marker_genes CoGAPS marker genes per pattern, or `NULL`.
#' @param cpus_per_task Worker CPUs (BiocParallel workers).
#' @param both_directions Test both ends of each factor in ORA
#'   (sign-ambiguous methods).
#' @return `list(dataset_id, method, fit_key, gsea, fora, local_gsea, local_ora)`
#'   of fgsea()/fora() results.
run_fgsea_job <- function(dataset_id, method, fit_key, loadings_file, pathways_by_source,
                          ensembl_map = NULL, cogaps_marker_genes = NULL, cpus_per_task = 8,
                          both_directions = FALSE) {
  pathways <- pathways_by_source$HALLMARK
  library(fgsea); library(BiocParallel)
  bp <- MulticoreParam(cpus_per_task)
  register(bp)
  t0 <- Sys.time()
  tag <- sprintf("[fgsea] %s/%s %s", dataset_id, method, fit_key)

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
  if (length(pathways_by_source) > 0) {
    gsea_sources <- setdiff(names(pathways_by_source), "HALLMARK")   # HALLMARK GSEA already covered above
    # both ends of a factor for sign-ambiguous methods (<method>_ingest$sign_ambiguous)
    dirs <- if (both_directions) c("pos", "neg") else "pos"

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
    message(sprintf("%s: empty `pathways_by_source` -- skipping local GSEA/ORA", tag))
  }

  message(sprintf("%s: done in %.1fs", tag, as.numeric(difftime(Sys.time(), t0, units = "secs"))))

  list(dataset_id = dataset_id, method = method, fit_key = fit_key,
       gsea = gsea_results, fora = fora_results,
       local_gsea = local_gsea_results, local_ora = local_ora_results)
}

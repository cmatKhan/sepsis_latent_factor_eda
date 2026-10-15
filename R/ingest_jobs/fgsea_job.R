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
#' All GSEA calls (every factor x collection, Hallmark included) run as one
#' BiocParallel round and all ORA calls (CoGAPS marker genes, each factor's top
#' genes) as a second, each call its own task with the largest collections
#' first, so the workers stay busy until the end; each inner call runs
#' serially to avoid nested forking.
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
#'   of fgsea()/fora() results: `gsea` is Hallmark, `local_gsea` the other
#'   collections.
run_fgsea_job <- function(dataset_id, method, fit_key, loadings_file, pathways_by_source,
                          ensembl_map = NULL, cogaps_marker_genes = NULL, cpus_per_task = 8,
                          both_directions = FALSE) {
  library(fgsea); library(BiocParallel)
  t0 <- Sys.time()
  tag <- sprintf("[fgsea] %s/%s %s", dataset_id, method, fit_key)

  L <- as.matrix(readRDS(loadings_file))
  L_ens <- remap_to_ensembl(L, ensembl_map)
  n_factors <- ncol(L_ens)
  message(sprintf("%s: starting -- %d factor(s), %d genes post-remap", tag, n_factors, nrow(L_ens)))

  # One task per call (`tasks`): BiocParallel's default splits the list into
  # one fixed chunk per worker, which leaves workers idle behind a chunk of
  # large collections. Largest collections first, so the small ones fill in.
  run_round <- function(items, f) {
    if (length(items) == 0) return(list())
    size <- vapply(items, function(it) length(pathways_by_source[[it$src]]), integer(1))
    ord <- order(-size)
    out <- bplapply(items[ord], f, BPPARAM = MulticoreParam(cpus_per_task, tasks = length(items)))
    out[order(ord)]
  }
  # both ends of a factor for sign-ambiguous methods (<method>_ingest$sign_ambiguous)
  dirs <- if (both_directions) c("pos", "neg") else "pos"
  grid_items <- function(...) {
    g <- expand.grid(..., stringsAsFactors = FALSE)
    lapply(seq_len(nrow(g)), function(i) as.list(g[i, , drop = FALSE]))
  }

  # ---- GSEA: every factor x collection ---------------------------------------
  gsea_items <- grid_items(fi = seq_len(n_factors), src = names(pathways_by_source))
  gsea_out <- run_round(gsea_items, function(it) {
    ranks <- sort(L_ens[, it$fi], decreasing = TRUE)
    res <- tryCatch(fgsea(pathways = pathways_by_source[[it$src]], stats = ranks, minSize = 10, maxSize = 500,
                          BPPARAM = SerialParam()),
                    error = function(e) NULL)
    list(factor_index = it$fi, source = it$src, result = res,
         main_pathways = main_pathways_for(res, pathways_by_source[[it$src]], ranks))
  })
  is_hallmark <- vapply(gsea_out, function(x) x$source == "HALLMARK", logical(1))
  gsea_results <- lapply(gsea_out[is_hallmark], function(x) x[c("factor_index", "result", "main_pathways")])
  local_gsea_results <- gsea_out[!is_hallmark]

  # ---- ORA: CoGAPS marker genes (Hallmark), each factor's top 100 genes ------
  markers <- if (method == "cogaps" && !is.null(cogaps_marker_genes)) cogaps_marker_genes else list()
  ora_items <- c(
    lapply(seq_along(markers), function(i) list(kind = "markers", i = i, src = "HALLMARK")),
    lapply(grid_items(fi = seq_len(n_factors), dir_i = dirs, src = names(pathways_by_source)),
           function(it) c(list(kind = "top"), it)))
  ora_out <- run_round(ora_items, function(it) {
    if (it$kind == "markers") {
      # The marker genes are in the fit's native ids: map them through the
      # same map as the loadings, or they never meet the universe.
      genes <- markers[[it$i]]
      genes_ens <- if (!is.null(ensembl_map)) {
        mapped <- ensembl_map[genes]
        unique(mapped[!is.na(mapped) & nzchar(mapped)])
      } else {
        genes
      }
      return(tryCatch(fora(pathways = pathways_by_source$HALLMARK, genes = genes_ens, universe = rownames(L_ens),
                           minSize = 10, maxSize = 500),
                      error = function(e) NULL))
    }
    v <- L_ens[, it$fi]
    genes <- if (it$dir_i == "neg") names(sort(v))[seq_len(min(100, length(v)))]
             else names(sort(v, decreasing = TRUE))[seq_len(min(100, length(v)))]
    res <- tryCatch(fora(pathways = pathways_by_source[[it$src]], genes = genes, universe = rownames(L_ens),
                         minSize = 10, maxSize = 500),
                    error = function(e) NULL)
    list(factor_index = it$fi, direction = it$dir_i, source = it$src, n_genes = length(genes), result = res)
  })
  is_marker <- vapply(ora_items, function(it) it$kind == "markers", logical(1))
  fora_results <- if (length(markers)) stats::setNames(ora_out[is_marker], names(markers)) else NULL
  local_ora_results <- ora_out[!is_marker]

  message(sprintf("%s: GSEA done -- %d collection(s) x %d factor(s), %d significant term-hit(s)",
                  tag, length(pathways_by_source), n_factors, n_sig(gsea_out)))
  message(sprintf("%s: ORA done -- %d call(s), %d significant term-hit(s)", tag, length(ora_out), n_sig(ora_out)))
  message(sprintf("%s: done in %.1fs", tag, as.numeric(difftime(Sys.time(), t0, units = "secs"))))

  list(dataset_id = dataset_id, method = method, fit_key = fit_key,
       gsea = gsea_results, fora = fora_results,
       local_gsea = local_gsea_results, local_ora = local_ora_results)
}

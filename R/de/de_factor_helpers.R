# Correlating DE genes with the factor-analysis pipeline's latent factors
# (PCA/NMF/CoGAPS/sPCA/ICA), for the SAME dataset -- two methods:
#   (a) overlap_test(): a binary contingency-table test (Fisher's exact)
#       between "genes in a factor" and "genes significant in a DE
#       contrast."
#   (b) factor_gsea_test(): a GSEA-style test -- rank every gene by a DE
#       contrast's signed t-statistic, ask whether a factor's own gene set
#       is enriched at either end of that ranking (fgsea::fgsea(), same
#       engine as R/de/de_fgsea_helpers.R, just with a factor's gene set
#       as the "pathway" instead of an MSigDB collection).
#
# Unlike R/de/de_fgsea_helpers.R (which needs everything in Ensembl to
# match MSigDB), this file works entirely in each dataset's OWN native
# gene-id space -- no Ensembl remap needed, just reconciling the DE
# topTable's convention against that SAME dataset's factor loadings'
# convention (see reconcile_loadings_to_de_ids()'s header).
#
# Needs R/db/connect.R (open_stability_db()/resolve_artifact(), the
# TWO-ARG version -- not app/R/db_helpers.R's single-arg app-only shadow)
# and R/lib/ingest/symbol_mapping.R (remap_to_symbol()) sourced first.

suppressMessages({
  library(fgsea)
})

#' The fit used for a (dataset, method) pair
#'
#' The representative fit whose rank is nearest `target_rank` -- a fixed rank
#' rather than model selection, for this first pass.
#'
#' @param con Stability DB connection.
#' @param dataset_id Dataset id.
#' @param method Method name.
#' @param target_rank Rank (K for sPCA) to aim for.
#' @return One-row data frame (fit_id, rank, loadings_file), or `NULL`.
pick_fit_at_rank <- function(con, dataset_id, method, target_rank = 10) {
  fits <- DBI::dbGetQuery(con,
    "SELECT f.fit_id, f.rank, a.path AS loadings_file FROM fits f
     JOIN fit_artifacts a ON a.fit_id = f.fit_id AND a.kind = 'loadings'
     WHERE f.dataset_id = ? AND f.method = ? AND f.status = 'ok' AND f.bootstrap = 0",
    params = list(dataset_id, method))
  if (nrow(fits) == 0) return(NULL)
  rep_ids <- representative_fit_ids(con, dataset_id, method)
  fits <- fits[fits$fit_id %in% rep_ids, ]
  if (nrow(fits) == 0) return(NULL)
  fits[which.min(abs(fits$rank - target_rank)), ]
}

#' A factor's gene set
#'
#' Non-zero genes for sPCA; the top `n` by |weight| for dense methods.
#'
#' @param loadings_col One loadings column (named numeric, names = genes).
#' @param method Method name.
#' @param n Genes to keep for dense methods.
#' @return Character vector of genes.
factor_gene_set <- function(loadings_col, method, n = 100) {
  if (identical(method, "spca")) {
    names(loadings_col)[loadings_col != 0]
  } else {
    names(sort(abs(loadings_col), decreasing = TRUE))[seq_len(min(n, length(loadings_col)))]
  }
}

#' Put loadings in the DE results' gene ids
#'
#' RNA-seq DE uses the native feature ids, so loadings already match; array DE
#' collapsed to gene symbols, so probe-keyed loadings are remapped with
#' remap_to_symbol().
#'
#' @param loadings Genes x factors matrix.
#' @param ds_meta The dataset's `dataset:` block.
#' @param feature_meta Feature-metadata data frame.
#' @param platform `"array"` or `"rnaseq"`.
#' @return The loadings with DE-compatible row names.
reconcile_loadings_to_de_ids <- function(loadings, ds_meta, feature_meta, platform = c("array", "rnaseq")) {
  platform <- match.arg(platform)
  if (platform == "rnaseq") return(loadings)
  symbol_map <- build_symbol_map(list(dataset = ds_meta), fm = feature_meta)
  remap_to_symbol(loadings, symbol_map)
}

#' Fisher's exact test of a factor's genes against DE genes
#'
#' Over `universe`, the genes both sides tested.
#'
#' @param factor_genes The factor's gene set.
#' @param de_sig_genes A contrast's significant genes.
#' @param universe Genes tested by both.
#' @return One-row data frame (counts, odds ratio, p-value).
overlap_test <- function(factor_genes, de_sig_genes, universe) {
  factor_genes <- intersect(factor_genes, universe)
  de_sig_genes <- intersect(de_sig_genes, universe)
  in_factor <- universe %in% factor_genes
  in_de     <- universe %in% de_sig_genes
  tab <- table(factor(in_factor, c(FALSE, TRUE)), factor(in_de, c(FALSE, TRUE)))
  ft <- stats::fisher.test(tab)
  data.frame(n_factor_genes = length(factor_genes), n_de_genes = length(de_sig_genes),
             n_overlap = sum(in_factor & in_de), universe_size = length(universe),
             odds_ratio = unname(ft$estimate), p_value = ft$p.value)
}

#' Are factors' gene sets enriched at either end of a DE ranking?
#'
#' fgsea() with each factor's gene set as a pathway; `minSize` 5 and no
#' `maxSize`, since a large set (e.g. a loosely sparse sPCA factor) is the
#' question here, not something to exclude.
#'
#' @param rank_vector One contrast's signed `t` statistics (named by gene,
#'   sorted decreasing).
#' @param factor_gene_sets Named list, one gene set per factor.
#' @return fgsea() result table.
factor_gsea_test <- function(rank_vector, factor_gene_sets) {
  tryCatch(
    fgsea::fgsea(pathways = factor_gene_sets, stats = rank_vector, minSize = 5, maxSize = Inf),
    error = function(e) NULL
  )
}

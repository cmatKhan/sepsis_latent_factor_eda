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
# Needs R/lib/ingest/db.R (open_stability_db()/resolve_artifact(), the
# TWO-ARG version -- not app/R/db_helpers.R's single-arg app-only shadow)
# and R/lib/ingest/symbol_mapping.R (remap_to_symbol()) sourced first.

suppressMessages({
  library(fgsea)
})

#' The one fit used for a (dataset, method) pair in this analysis: a
#' FIXED default rank (10 for pca/nmf/cogaps/ica, K=10 for spca),
#' deliberately sidestepping model-selection entirely for this first pass
#' (see R/de/README.md) -- falls back to the closest available rank if 10
#' isn't in that dataset's grid. Uses representative_fit_ids() (R/lib/
#' ingest/redundancy.R) to collapse the seed/para sweep first, so this is
#' never picking an arbitrary single seed.
pick_fit_at_rank <- function(con, dataset_id, method, target_rank = 10) {
  fits <- DBI::dbGetQuery(con,
    "SELECT fit_id, rank, loadings_file FROM fits
     WHERE dataset_id = ? AND method = ? AND status = 'ok' AND loadings_file IS NOT NULL
       AND (bootstrap IS NULL OR bootstrap = 0)",
    params = list(dataset_id, method))
  if (nrow(fits) == 0) return(NULL)
  rep_ids <- representative_fit_ids(con, dataset_id, method)
  fits <- fits[fits$fit_id %in% rep_ids, ]
  if (nrow(fits) == 0) return(NULL)
  fits[which.min(abs(fits$rank - target_rank)), ]
}

#' A factor's gene set: non-zero-weight genes for sPCA (genuinely sparse),
#' top-`n` by |weight| for every other (dense) method -- the two rules the
#' user specified directly. `loadings_col` is one column of a loadings
#' matrix (named numeric vector, names = genes).
factor_gene_set <- function(loadings_col, method, n = 100) {
  if (identical(method, "spca")) {
    names(loadings_col)[loadings_col != 0]
  } else {
    names(sort(abs(loadings_col), decreasing = TRUE))[seq_len(min(n, length(loadings_col)))]
  }
}

#' Reconcile a loadings matrix's rownames to the SAME gene-id convention
#' R/de/<DATASET>_de.R's topTable uses for that dataset (see
#' R/de/README.md's per-dataset table):
#'   - RNA-seq datasets: DE never collapses (R/de/de_helpers.R's
#'     collapse_to_symbol() is array-only) -- loadings are already fit on
#'     the same native feature_id DE uses, no remap needed.
#'   - Array datasets: DE collapsed to gene SYMBOL, but loadings are still
#'     probe-keyed -- remapped via R/lib/ingest/symbol_mapping.R's
#'     remap_to_symbol() (built for app display, reused here for a real
#'     computational join; same many-probes-to-one-symbol collapse rule
#'     as remap_to_ensembl(), see that file's .remap_ids()).
reconcile_loadings_to_de_ids <- function(loadings, ds_meta, feature_meta, platform = c("array", "rnaseq")) {
  platform <- match.arg(platform)
  if (platform == "rnaseq") return(loadings)
  symbol_map <- build_symbol_map(list(dataset = ds_meta), fm = feature_meta)
  remap_to_symbol(loadings, symbol_map)
}

#' Fisher's exact test on the 2x2 table (in factor & DE-sig, in factor &
#' not, not-in-factor & DE-sig, not-in-factor & not) over `universe` (the
#' genes actually tested by BOTH sides -- required for a valid test, not
#' "every gene in the genome"). Returns a one-row data.frame.
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

#' GSEA-style test: is a factor's gene set enriched at either end of a DE
#' contrast's ranking? `rank_vector` = one contrast's signed `t` stat
#' (named by gene, sorted decreasing); `factor_gene_sets` = named list,
#' one entry per factor (see factor_gene_set()). minSize is lower than
#' R/de/de_fgsea_helpers.R's MSigDB-collection default (5 vs 10) since
#' factor gene sets are far smaller/user-defined; no maxSize cap -- unlike
#' testing against MSigDB collections, here a large gene set (e.g. a
#' loosely-sparse sPCA factor) is exactly what's being asked about, not
#' something to exclude for being "too big a pathway."
factor_gsea_test <- function(rank_vector, factor_gene_sets) {
  tryCatch(
    fgsea::fgsea(pathways = factor_gene_sets, stats = rank_vector, minSize = 5, maxSize = Inf),
    error = function(e) NULL
  )
}

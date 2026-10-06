# "Too many patterns?" diagnostics for one representative fit per rank: how
# distinct a fit's own factors are from each other (within-fit cosine), and
# which genes mark each factor. CoGAPS uses CoGAPS::patternMarkers() on its
# raw result; other methods use a generic analog (each gene assigned to the
# factor it loads most strongly on, after z-scoring across factors). A factor
# with no marker genes suggests too many patterns. Not run for PCA, whose
# components are orthogonal by construction.

#' Cosine similarity among one fit's own factors
#'
#' @param loadings Genes x factors matrix.
#' @return K x K signed cosine matrix (diagonal 1).
within_fit_similarity <- function(loadings) {
  prep <- prepare_loadings(loadings)
  pair_similarities(prep, prep)$cosine   # K x K, diagonal == 1
}

#' Marker genes by strongest relative loading
#'
#' Each gene is assigned to the factor where its row-z-scored loading is
#' largest (z-scoring across factors keeps a gene's overall expression level
#' from deciding it). Genes with identical loadings everywhere (common in
#' sparse fits) go to the first factor.
#'
#' @param loadings Genes x factors matrix with gene row names.
#' @return Data frame (gene, factor_index, score), one row per gene.
generic_pattern_markers <- function(loadings) {
  z <- t(scale(t(loadings)))             # row (gene) z-score across factors
  # A gene with (near-)identical loadings across every factor -- common for
  # sparse methods like sPCA, especially at low K, where most loadings are
  # exactly 0 -- gets zero column variance from scale(), producing an
  # all-NaN row here. which.max() on an all-NaN vector returns integer(0),
  # which breaks apply()'s result simplification (it silently returns a
  # list instead of an atomic vector, and everything downstream indexing
  # into it fails with "invalid subscript type 'list'"). Treat such a gene
  # as tied across all factors (z = 0 everywhere) rather than erroring --
  # which.max() then deterministically picks the first factor.
  z[!is.finite(z)] <- 0
  best <- apply(z, 1, which.max)
  data.frame(gene = rownames(loadings), factor_index = best,
             score = z[cbind(seq_len(nrow(z)), best)], stringsAsFactors = FALSE)
}

#' CoGAPS marker genes
#'
#' CoGAPS::patternMarkers() in the generic_pattern_markers() shape. `score` is
#' the rank within the pattern, not comparable to other methods' scores.
#'
#' @param cogaps_result A `CogapsResult`.
#' @return Data frame (gene, factor_index, score).
cogaps_pattern_markers <- function(cogaps_result) {
  pm <- CoGAPS::patternMarkers(cogaps_result)
  out <- do.call(rbind, lapply(seq_along(pm$PatternMarkers), function(k) {
    genes <- pm$PatternMarkers[[k]]
    if (length(genes) == 0) return(NULL)
    data.frame(gene = genes, factor_index = k, score = seq_along(genes), stringsAsFactors = FALSE)
  }))
  if (is.null(out)) out <- data.frame(gene = character(0), factor_index = integer(0), score = numeric(0))
  out
}

#' Within-fit redundancy summary
#'
#' @param loadings Genes x factors matrix.
#' @param markers Marker-gene data frame (gene, factor_index, score).
#' @return `list(max_offdiag, median_offdiag, n_factors_with_no_markers,
#'   matrix)`: the largest and median off-diagonal |cosine|, the number of
#'   factors without a marker gene, and the full K x K matrix.
redundancy_summary <- function(loadings, markers) {
  sim <- within_fit_similarity(loadings)
  diag(sim) <- NA
  list(
    max_offdiag = max(abs(sim), na.rm = TRUE),
    median_offdiag = median(abs(sim), na.rm = TRUE),
    n_factors_with_no_markers = length(setdiff(seq_len(ncol(loadings)), unique(markers$factor_index))),
    matrix = sim
  )
}

#' Redundancy diagnostic for one fit
#'
#' @param method Method name (CoGAPS uses its raw result for markers).
#' @param loadings_path Absolute path to the fit's loadings artifact.
#' @param raw_result_path Absolute path to the raw CoGAPS result, or `NA`.
#' @return `list(markers, summary)`, or `NULL` if the loadings can't be read.
run_redundancy_for_fit <- function(method, loadings_path, raw_result_path = NA_character_) {
  if (is.na(loadings_path) || !file.exists(loadings_path)) return(NULL)
  L <- as.matrix(readRDS(loadings_path))

  markers <- if (method == "cogaps" && !is.na(raw_result_path) && file.exists(raw_result_path)) {
    cogaps_raw <- tryCatch(readRDS(raw_result_path), error = function(e) NULL)
    if (is.null(cogaps_raw)) generic_pattern_markers(L) else cogaps_pattern_markers(cogaps_raw)
  } else {
    generic_pattern_markers(L)
  }
  list(markers = markers, summary = redundancy_summary(L, markers))
}

#' A method's representative fits in the DB
#'
#' @param con DB connection.
#' @param dataset_id Dataset id.
#' @param method Method name.
#' @return Integer `fit_id`s with `representative = 1`, ordered by rank.
representative_fit_ids <- function(con, dataset_id, method) {
  DBI::dbGetQuery(con,
    "SELECT fit_id FROM fits WHERE dataset_id = ? AND method = ? AND representative = 1
     ORDER BY rank", params = list(dataset_id, method))$fit_id
}


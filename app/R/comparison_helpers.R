# Cross-method comparison of fits within one dataset: top-gene-set overlap
# (Jaccard) and sample-space agreement. Loadings similarity itself lives in
# R/lib/ingest/similarity.R (prepare_loadings(), hungarian_match_abs()).

#' Top `n` gene ids of one loading vector. `sign`: "both" = largest |loading|;
#' "positive"/"negative" = largest / most negative values.
top_genes <- function(loadings_vector, n, sign = c("both", "positive", "negative")) {
  sign <- match.arg(sign)
  v <- loadings_vector
  ord <- switch(sign,
    both     = order(abs(v), decreasing = TRUE),
    positive = order(v, decreasing = TRUE),
    negative = order(v, decreasing = FALSE))
  names(v)[ord[seq_len(min(n, length(v)))]]
}

gene_sets_from_loadings <- function(L, n, sign = "both") {
  sets <- lapply(seq_len(ncol(L)), function(j) top_genes(L[, j], n, sign))
  names(sets) <- colnames(L) %||% paste0("factor_", seq_len(ncol(L)))
  sets
}

#' WGCNA's analogue: one gene set per module (full membership; module 0, the
#' unassigned genes, is excluded). Names are the module numbers.
gene_sets_from_modules <- function(con, fit_id) {
  d <- DBI::dbGetQuery(con, "SELECT gene, module FROM gene_modules WHERE fit_id = ? AND module != 0",
                       params = list(fit_id))
  if (nrow(d) == 0) return(list())
  split(d$gene, d$module)
}

#' Gene sets for one fit: top-`n` loading genes per factor, or module
#' membership for WGCNA. `n_factors` keeps the first n factors (PCA's nested
#' fit at a smaller rank). NULL if the fit has nothing to compare.
fit_gene_sets <- function(con, method, fit_id, n = 50, n_factors = NULL) {
  if (method == "wgcna") {
    s <- gene_sets_from_modules(con, fit_id)
    return(if (length(s) == 0) NULL else s)
  }
  L <- load_loadings(con, fit_id, n_factors)
  if (is.null(L)) return(NULL)
  gene_sets_from_loadings(L, n)
}

jaccard <- function(a, b) {
  u <- union(a, b)
  if (length(u) == 0) return(NA_real_)
  length(intersect(a, b)) / length(u)
}

jaccard_matrix <- function(sets_a, sets_b) {
  m <- matrix(NA_real_, length(sets_a), length(sets_b), dimnames = list(names(sets_a), names(sets_b)))
  for (i in seq_along(sets_a)) for (j in seq_along(sets_b)) m[i, j] <- jaccard(sets_a[[i]], sets_b[[j]])
  m
}

#' Gene-set agreement between two fits: Hungarian-match the factors on |Jaccard|
#' and report the matched values. Returns list(matrix, matches, mean) where
#' `matches` has a_label/b_label/jaccard rows, or NULL if either side is empty.
gene_set_agreement <- function(sets_a, sets_b) {
  if (is.null(sets_a) || is.null(sets_b)) return(NULL)
  jm <- jaccard_matrix(sets_a, sets_b)
  idx <- hungarian_match_abs(jm)
  matches <- data.frame(a = rownames(jm)[idx[, "a"]], b = colnames(jm)[idx[, "b"]],
                        jaccard = jm[idx], stringsAsFactors = FALSE)
  list(matrix = jm, matches = matches[order(-matches$jaccard), ], mean = mean(matches$jaccard, na.rm = TRUE),
       idx = idx)
}

#' Cluster each fit's samples (hierarchical, `k` clusters) and compare the two
#' assignments on shared samples by Adjusted Rand Index. NULL if fewer than
#' max(k, 3) samples are shared or either fit has no scores.
cluster_and_ari <- function(con, fit_a, fit_b, n_a = NULL, n_b = NULL, k = 4) {
  La <- load_scores(con, fit_a, n_a); Lb <- load_scores(con, fit_b, n_b)
  if (is.null(La) || is.null(Lb)) return(NULL)
  shared <- intersect(rownames(La), rownames(Lb))
  if (length(shared) < max(k, 3)) return(NULL)
  cl_a <- cutree(hclust(dist(La[shared, , drop = FALSE])), k = k)
  cl_b <- cutree(hclust(dist(Lb[shared, , drop = FALSE])), k = k)
  list(ari = mclust::adjustedRandIndex(cl_a, cl_b), n_shared = length(shared))
}

#' Method x method agreement for one dataset's optimal fits (`optimal` is
#' optimal_fits() output). Returns list(methods, gene_jaccard, sample_ari,
#' pairs) where the matrices are symmetric with diagonal 1 and `pairs` holds
#' each pair's gene_set_agreement() for drill-down.
method_pairwise_agreement <- function(con, optimal, top_n = 50) {
  methods <- optimal$method
  n <- length(methods)
  gj <- sa <- matrix(NA_real_, n, n, dimnames = list(methods, methods))
  diag(gj) <- 1; diag(sa) <- 1
  pairs <- list()
  sets <- lapply(seq_len(n), function(i) fit_gene_sets(con, methods[i], optimal$fit_id[i], top_n, optimal$n_factors[i]))
  for (i in seq_len(n)) for (j in seq_len(n)) {
    if (i >= j) next
    ag <- gene_set_agreement(sets[[i]], sets[[j]])
    if (!is.null(ag)) {
      gj[i, j] <- gj[j, i] <- ag$mean
      pairs[[paste(methods[i], methods[j], sep = "|")]] <- ag
    }
    ca <- tryCatch(cluster_and_ari(con, optimal$fit_id[i], optimal$fit_id[j],
                                   if (methods[i] == "pca") optimal$n_factors[i],
                                   if (methods[j] == "pca") optimal$n_factors[j]),
                   error = function(e) NULL)
    if (!is.null(ca)) sa[i, j] <- sa[j, i] <- ca$ari
  }
  list(methods = methods, gene_jaccard = gj, sample_ari = sa, pairs = pairs)
}

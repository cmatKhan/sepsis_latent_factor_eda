# Factor similarity: cosine, Pearson and Spearman between the columns of two
# loading matrices (genes x factors), and 1-to-1 Hungarian matching. Each
# metric is a crossprod of column-normalized matrices, so prepare_loadings()
# transforms each fit once and every pair comparison is three t(A) %*% B.

#' Normalize a loading matrix for similarity
#'
#' @param mat Genes x factors matrix with gene row names.
#' @return `list(cosine, pearson, spearman)`: the matrix with unit-norm
#'   columns, centered then unit-norm, and ranked, centered then unit-norm.
prepare_loadings <- function(mat) {
  stopifnot(!is.null(rownames(mat)))
  unit_cols <- function(m) {
    nrm <- sqrt(colSums(m^2))
    nrm[nrm == 0] <- 1
    sweep(m, 2, nrm, "/")
  }
  centered <- sweep(mat, 2, colMeans(mat), "-")
  ranked   <- apply(mat, 2, rank)
  rownames(ranked) <- rownames(mat)
  ranked_c <- sweep(ranked, 2, colMeans(ranked), "-")
  list(
    cosine   = unit_cols(mat),
    pearson  = unit_cols(centered),
    spearman = unit_cols(ranked_c)
  )
}

#' Similarity between the factors of two fits
#'
#' Computed on the genes both fits share. Normalization happened on each fit's
#' full gene set, which is exact when both share all genes (the usual case: same
#' dataset, same matrix); otherwise the shared subset is renormalized.
#'
#' @param prep_a,prep_b prepare_loadings() outputs.
#' @return `list(cosine, pearson, spearman)` of factors_a x factors_b
#'   matrices, or `NULL` when fewer than 2 genes are shared.
pair_similarities <- function(prep_a, prep_b) {
  shared <- intersect(rownames(prep_a$cosine), rownames(prep_b$cosine))
  if (length(shared) < 2) return(NULL)
  full_a <- nrow(prep_a$cosine); full_b <- nrow(prep_b$cosine)
  if (length(shared) == full_a && length(shared) == full_b &&
      identical(rownames(prep_a$cosine), rownames(prep_b$cosine))) {
    list(
      cosine   = crossprod(prep_a$cosine,   prep_b$cosine),
      pearson  = crossprod(prep_a$pearson,  prep_b$pearson),
      spearman = crossprod(prep_a$spearman, prep_b$spearman)
    )
  } else {
    # gene sets differ: renormalize on the shared subset so each metric is
    # a true cosine/correlation over the genes actually compared
    sub <- function(m) {
      m <- m[shared, , drop = FALSE]
      nrm <- sqrt(colSums(m^2)); nrm[nrm == 0] <- 1
      sweep(m, 2, nrm, "/")
    }
    list(
      cosine   = crossprod(sub(prep_a$cosine),   sub(prep_b$cosine)),
      pearson  = crossprod(sub(prep_a$pearson),  sub(prep_b$pearson)),
      spearman = crossprod(sub(prep_a$spearman), sub(prep_b$spearman))
    )
  }
}

#' Optimal 1-to-1 factor matching
#'
#' Maximizes total (signed) similarity with `clue::solve_LSAP()`.
#'
#' @param sim Factors_a x factors_b similarity matrix.
#' @return Two-column integer matrix (`a`, `b`) of matched row and column
#'   indices, `min(nrow, ncol)` rows.
hungarian_match <- function(sim) {
  transposed <- FALSE
  if (nrow(sim) > ncol(sim)) {
    sim <- t(sim)
    transposed <- TRUE
  }
  cost <- max(sim) - sim
  assignment <- clue::solve_LSAP(cost)
  a <- seq_len(nrow(sim))
  b <- as.integer(assignment)
  out <- if (transposed) cbind(a = b, b = a) else cbind(a = a, b = b)
  out
}

#' Optimal 1-to-1 factor matching on |similarity|
#'
#' For sign-ambiguous methods (PCA, sPCA, ICA), where a component can flip sign
#' between fits: a strongly negative similarity is as good a match as a
#' strongly positive one. Also used across methods (e.g. an NMF factor matching
#' the negative tail of a PCA component) and on Jaccard matrices.
#'
#' @param sim Factors_a x factors_b similarity matrix.
#' @return Two-column integer matrix (`a`, `b`) of matched indices.
hungarian_match_abs <- function(sim) {
  transposed <- FALSE
  if (nrow(sim) > ncol(sim)) { sim <- t(sim); transposed <- TRUE }
  a_sim <- abs(sim)
  cost <- max(a_sim) - a_sim
  assignment <- clue::solve_LSAP(cost)
  a <- seq_len(nrow(sim))
  b <- as.integer(assignment)
  if (transposed) cbind(a = b, b = a) else cbind(a = a, b = b)
}

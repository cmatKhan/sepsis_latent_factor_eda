# Factor-similarity machinery: cosine / Pearson / Spearman between the
# columns of two loading matrices (genes x factors), plus 1-to-1 Hungarian
# matching. All three metrics are computed and stored side by side (user
# decision -- space is not a concern yet).
#
# Efficiency note: each metric reduces to a crossprod of appropriately
# transformed column-unit-norm matrices, so `prepare_loadings()` transforms
# each fit's loadings ONCE (cosine: L2-normalized columns; Pearson:
# column-centered then normalized; Spearman: column-ranked, centered,
# normalized) and every pair comparison is then just three t(A) %*% B.

#' Transform one loading matrix into the three normalized variants.
#' Rows must be named by gene.
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

#' Similarity matrices (factors_a x factors_b) between two prepared fits,
#' computed on the intersection of their gene sets. Note: normalization
#' happened on the FULL gene set at prepare time; when the two fits share
#' all genes (the usual case here -- same dataset, same input matrix) this
#' is exact. Fits with disjoint genes return NULL.
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

#' Optimal 1-to-1 assignment maximizing total similarity, via
#' clue::solve_LSAP (which minimizes cost and requires a non-negative
#' matrix with nrow <= ncol -- handled by shifting and transposing).
#' Returns a two-column matrix (a = row factor index, b = col factor
#' index) with min(nrow, ncol) rows.
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

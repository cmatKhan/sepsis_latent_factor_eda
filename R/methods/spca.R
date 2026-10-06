# Sparse PCA (elasticnet::spca()), deterministic given its parameters.
#
# Sparsity: `sparse = "varnum"`, under which `para` is the fraction of genes
# allowed a non-zero loading per component, converted to a count from the
# dataset's gene number -- portable across datasets, unlike the raw L1
# penalty of `sparse = "penalty"`, whose effective sparsity varied with gene
# count (and was sensitive to `lambda`).
#
# Search (Guerra-Urzola, Van Deun, Vera & Sijtsma 2021, Psychometrika
# 86(4):893-919, https://pmc.ncbi.nlm.nih.gov/articles/PMC8636462/): fix K,
# choose the sparsity maximizing the Index of Sparseness
# IS = PEV_sparse * PEV_pca * PS. A coarse `para_coarse` grid at every K, then
# spca_refine_grid() around each K's peak, run as its own targets stage;
# spca_select_fits() picks the representative para per K. `lambda` is one
# value: under varnum it changed total PEV by at most 5.7e-4 on EARLI.
#
# A fit still running after `max_hours` stops itself and returns
# `killed_not_converged = TRUE` (high sparsity and K can run all `max.iter`
# iterations without converging). Fits the matrix `mat`.
#
# The fit targets depend on this file's registry and fit function: editing
# their code refits the method (comments don't count). Docs: Methods.

#' Fit sPCA for one (K, para)
#'
#' @param K Number of components.
#' @param para Fraction of genes non-zero per component (`sparse = "varnum"`),
#'   or the L1 penalty (`sparse = "penalty"`).
#' @param type,sparse,use.corr,lambda,max.iter,eps.conv Passed to elasticnet::spca().
#' @param max_hours Wall-time limit; a longer fit is abandoned.
#' @return `list(rank, killed_not_converged = FALSE, elapsed, mse, loadings,
#'   scores, center, pev, var.all, spca)` -- scores on the centered matrix,
#'   `spca` the rest of spca()'s output; or `list(rank, killed_not_converged =
#'   TRUE, elapsed)` when stopped.
run_spca_job <- function(K, para, type = "predictor", sparse = "varnum",
                          use.corr = FALSE, lambda = 1e-6, max.iter = 200, eps.conv = 1e-3,
                          max_hours = 8) {
  library(elasticnet)
  # `para` is a fraction of genes when sparse="varnum" (the default here
  # and in every dataset config) -- convert to an absolute nonzero-
  # loadings target using this dataset's REAL gene count, not a value
  # calibrated elsewhere. Left as a raw penalty (previous behavior) if
  # `sparse="penalty"` is ever explicitly requested again.
  para_vec <- if (identical(sparse, "varnum")) {
    rep(max(1L, round(para * nrow(mat))), K)
  } else {
    rep(para, K)
  }

  # Stopping condition: a fit still running after `max_hours` is abandoned
  # and reported as killed_not_converged = TRUE with no other results.
  # At high sparsity and high K, spca()'s alternating algorithm appears not
  # to converge and runs all max.iter iterations (ANEMONES para 0.05 jumped
  # from <15 min at K 8 to 2.9 h at K 9 and >8 h at K 10-12, 2026-10-04);
  # elasticnet doesn't report its iteration count, so elapsed time is the
  # only signal. setTimeLimit() interrupts at R's next interrupt check,
  # which spca()'s R-level loops reach constantly.
  t0 <- proc.time()[["elapsed"]]
  setTimeLimit(elapsed = max_hours * 3600, transient = TRUE)
  on.exit(setTimeLimit(elapsed = Inf), add = TRUE)
  fit <- tryCatch(
    spca(t(mat), K = K, para = para_vec, type = type, sparse = sparse,
         use.corr = use.corr, lambda = lambda, max.iter = max.iter, eps.conv = eps.conv),
    error = function(e) {
      if (grepl("elapsed time limit", conditionMessage(e))) NULL else stop(e)
    }
  )
  setTimeLimit(elapsed = Inf)
  elapsed <- proc.time()[["elapsed"]] - t0
  if (is.null(fit)) {
    return(list(rank = K, killed_not_converged = TRUE, elapsed = elapsed))
  }

  loadings <- fit$loadings
  rownames(loadings) <- rownames(mat)
  colnames(loadings) <- paste0("Component_", seq_len(ncol(loadings)))
  # spca() centers the data internally (scale(x, center = TRUE, ...)), so
  # scores are the centered samples projected onto the loadings.
  center <- rowMeans(mat)
  scores <- t(mat - center) %*% loadings

  list(rank = K, killed_not_converged = FALSE, elapsed = elapsed,
       # mse = 1 - CUMULATIVE PEV (sum across all K components): elasticnet's
       # `pev` is per-component (each entry is that ONE component's
       # incremental, non-overlapping contribution via the QR-based
       # adjustment for non-orthogonality; see ?elasticnet::spca), so its
       # sum is the variance explained by the whole K-component solution.
       mse = 1 - sum(fit$pev), loadings = loadings, scores = scores, center = center,
       # pev's full per-component curve -- elasticnet's own adjusted-variance
       # bookkeeping, NOT recomputable by naively projecting mat onto
       # loadings (sPCA's components aren't orthogonal, so that double-counts
       # shared variance)
       pev = fit$pev, var.all = fit$var.all,
       # the rest of spca()'s return value (call, type, K, para, lambda,
       # sparse, vn), so ingestion never needs to refit
       spca = unclass(fit)[setdiff(names(fit), c("loadings", "pev", "var.all"))])
}

#' Index of Sparseness for every fit
#'
#' Fits stopped at `max_hours` have no loadings and are left out.
#'
#' @param fits List of `list(params, result)`, both stages.
#' @param pca_fit The dataset's PCA result (its full `sdev`).
#' @return Data frame (K, para, PEV_sparse, PEV_pca, PS, IS, selected), one row
#'   per fit; `selected` marks the IS maximum within each K.
spca_select_fits <- function(fits, pca_fit) {
  pev_pca <- cumsum(pca_fit$sdev^2) / sum(pca_fit$sdev^2)
  # fits stopped by run_spca_job()'s max_hours have no loadings/PEV and
  # can't be scored; they're absent from the IS table
  fits <- Filter(function(f) !isTRUE(f$result$killed_not_converged), fits)
  tab <- do.call(rbind, lapply(fits, function(f) {
    L <- as.matrix(f$result$loadings)
    data.frame(K = f$params$K, para = f$params$para,
               PEV_sparse = sum(f$result$pev), PEV_pca = pev_pca[f$params$K],
               PS = mean(L == 0))
  }))
  tab$IS <- tab$PEV_sparse * tab$PEV_pca * tab$PS
  tab <- tab[order(tab$K, tab$para), ]
  tab$selected <- stats::ave(tab$IS, tab$K, FUN = function(x) seq_along(x) == which.max(x)) == 1
  rownames(tab) <- NULL
  tab
}

#' Second-stage para values around each K's IS peak
#'
#' The midpoints to each neighbouring coarse value, plus one step outward
#' (+0.1 above, or half below) when the peak is at the grid's edge; only
#' values not already fit.
#'
#' @param fits The coarse stage's fits.
#' @param pca_fit The dataset's PCA result.
#' @param p Resolved sPCA parameters.
#' @return Grid rows shaped like spca_build_grid()'s, or `NULL`.
spca_refine_grid <- function(fits, pca_fit, p) {
  tab <- spca_select_fits(fits, pca_fit)
  rows <- lapply(split(tab, tab$K), function(t) {
    tried <- sort(unique(t$para))
    best <- t$para[which.max(t$IS)]
    i <- match(best, tried)
    new <- c(if (i > 1) (tried[i - 1] + best) / 2 else best / 2,
             if (i < length(tried)) (tried[i + 1] + best) / 2 else min(1, best + 0.1))
    new <- setdiff(round(new, 4), round(tried, 4))
    if (length(new) == 0) return(NULL)
    spca_build_grid(utils::modifyList(p, list(K = t$K[1], para_coarse = new)))
  })
  do.call(rbind, rows)
}

#' sPCA's coarse grid
#'
#' Crosses `K` x `para_coarse`; errors if `lambda` has more than one value.
#'
#' @param p Resolved sPCA parameters.
#' @return Data frame, one row per fit.
spca_build_grid <- function(p) {
  if (length(p$lambda) != 1) {
    stop("spca's `lambda` must be a single value, not swept -- see R/README.md's ",
         "\"sPCA\" section (no practical effect under sparse=\"varnum\").")
  }
  expand.grid(K = p$K, para = p$para_coarse, sparse = p$sparse, use.corr = p$use.corr,
              lambda = p$lambda, max.iter = p$max.iter, eps.conv = p$eps.conv,
              max_hours = p$max_hours, stringsAsFactors = FALSE)
}

spca_registry <- list(
  needs_nonneg = FALSE,
  global_object = "mat",
  jobname = "spca_grid",
  fn = run_spca_job,
  pkgs = "elasticnet",
  # para_coarse: the first-stage grid (fraction of genes nonzero per
  # component under sparse="varnum"); spca_refine_grid() adds points
  # around each K's IS peak. Reaches down to 0.02 because the paper's
  # gene-expression IS peak was at 97% sparsity, and low-para fits are
  # the cheapest. See the header and docs/methods.qmd, "sPCA: choosing the sparsity".
  defaults = list(K = 2:12, para_coarse = c(0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5),
                   sparse = "varnum", lambda = 1e-6,
                   use.corr = FALSE, max.iter = 200, eps.conv = 1e-3,
                   max_hours = 8),
  build_grid = spca_build_grid,
  refine_grid = spca_refine_grid,
  select_fits = spca_select_fits
)

# Ingest contract (R/targets/ingest.R): kept apart from spca_registry, which
# the fit targets depend on, so editing it never refits. `mse` stays NULL:
# sPCA's natural summary is variance explained (pev_sparse in fit_metrics),
# not a reconstruction error. The IS selection metrics come from
# selection_spca_<dataset> (spca_select_fits()).
spca_ingest <- list(
  family = "param_grid", sign_ambiguous = TRUE, has_loadings = TRUE,
  resource_params = character(0),
  extract = function(result, params) {
    L <- result$loadings
    nonzero <- if (!is.null(L)) colSums(as.matrix(L) != 0) else NULL
    fit_record_from(
      rank = params$K, loadings = L, scores = result$scores,
      diag = c(list(pev = result$pev, var_all = result$var.all, n_nonzero = nonzero,
                    center = result$center), result$spca),
      metrics = c(pev_sparse = if (!is.null(result$pev)) sum(result$pev) else NA_real_),
      factor_n_genes = nonzero)
  }
)

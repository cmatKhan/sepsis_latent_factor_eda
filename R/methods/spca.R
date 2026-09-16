# sPCA (elasticnet::spca()) method registry.
#
# Deterministic given (K, para) -- no seed argument.
#
# `para` is spca()'s per-component sparsity penalty (a length-K vector),
# but here it's a single SCALAR applied uniformly across every component
# (`para_vec <- rep(para, K)` inside the job function below) -- lets `K`
# sweep freely over an arbitrary rank range instead of requiring a
# hand-curated penalty vector per K value. `K` and `para` are independent,
# CROSSED dimensions (`expand.grid()` below handles this the same way it
# handles any other pair of swept arguments -- no special-casing needed).
#
# Operates on the single matrix `mat` (set as a global object by the
# orchestrator script) -- there is no basis argument.

run_spca_job <- function(K, para, type = "predictor", sparse = "penalty",
                          use.corr = FALSE, lambda = 1e-6, max.iter = 200, eps.conv = 1e-3) {
  library(elasticnet)
  para_vec <- rep(para, K)
  fit <- spca(t(mat), K = K, para = para_vec, type = type, sparse = sparse,
              use.corr = use.corr, lambda = lambda, max.iter = max.iter, eps.conv = eps.conv)
  loadings <- fit$loadings
  rownames(loadings) <- rownames(mat)
  colnames(loadings) <- paste0("Component_", seq_len(ncol(loadings)))
  scores <- t(mat) %*% loadings

  list(rank = K, mse = 1 - fit$pev[length(fit$pev)], loadings = loadings, scores = scores,
       # pev's full per-component curve -- elasticnet's own adjusted-variance
       # bookkeeping, NOT recomputable by naively projecting mat onto
       # loadings (sPCA's components aren't orthogonal, so that double-counts
       # shared variance)
       pev = fit$pev, var.all = fit$var.all)
}

spca_registry <- list(
  needs_nonneg = FALSE,
  global_object = "mat",
  jobname = "spca_grid",
  fn = run_spca_job,
  pkgs = "elasticnet",
  defaults = list(K = 2:20, para = c(0.05, 0.1, 0.2), sparse = "penalty",
                   use.corr = FALSE, lambda = 1e-6, max.iter = 200, eps.conv = 1e-3),
  build_grid = function(p) {
    expand.grid(K = p$K, para = p$para, sparse = p$sparse, use.corr = p$use.corr,
                lambda = p$lambda, max.iter = p$max.iter, eps.conv = p$eps.conv,
                stringsAsFactors = FALSE)
  }
)

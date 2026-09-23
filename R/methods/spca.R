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

  list(rank = K,
       # mse = 1 - CUMULATIVE PEV (sum across all K components), not the
       # last component's own individual contribution -- elasticnet's own
       # `pev` vector is per-component (each entry is that ONE component's
       # incremental, non-overlapping contribution via the QR-based
       # adjustment for non-orthogonality; see ?elasticnet::spca), so
       # summing it gives the total variance explained by the whole
       # K-component solution together. Bug fixed 2026-09-22: this used
       # to be `1 - fit$pev[length(fit$pev)]`, silently reporting mse
       # near 1 (implying ~0% variance explained) for every fit, since a
       # single late component's own incremental contribution is always
       # small on real data -- confirmed on real data the true cumulative
       # value was ~40x larger (0.058 vs 0.0013 for one real K=10 fit).
       mse = 1 - sum(fit$pev), loadings = loadings, scores = scores,
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
  # para: widened 2026-09-22 -- real-data testing (ANEMONES, K=10, full
  # 8000-gene matrix) showed the old c(0.05, 0.1, 0.2) grid was already
  # fully saturated (~94% sparsity, matching production) at para=0.001,
  # 50x smaller than the old minimum; the true transition to lower
  # sparsity lies somewhere below that, unlocated locally (each full-scale
  # fit costs ~9min single-threaded) -- this fallback now matches every
  # real dataset config's grid (config/*.yml) so a new config that omits
  # `para` gets the same wide sweep, not the old saturated one.
  defaults = list(K = 2:20, para = c(1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 0.01, 0.05, 0.1, 0.2),
                   sparse = "penalty",
                   use.corr = FALSE, lambda = 1e-6, max.iter = 200, eps.conv = 1e-3),
  build_grid = function(p) {
    expand.grid(K = p$K, para = p$para, sparse = p$sparse, use.corr = p$use.corr,
                lambda = p$lambda, max.iter = p$max.iter, eps.conv = p$eps.conv,
                stringsAsFactors = FALSE)
  }
)

# sPCA (elasticnet::spca()) method registry.
#
# Deterministic given (K, para, lambda) -- no seed argument.
#
# `para`/`sparse`/`lambda` history (revised 2026-09-29): this used to run
# `sparse="penalty"` with `para` as spca()'s raw L1 penalty and `lambda`
# (the elastic-net ridge/quadratic term) fixed at the package default
# (1e-6) always. Real-data testing that day found two problems with that:
#   1. A fixed penalty's EFFECTIVE sparsity scales strongly with gene
#      count -- the exact same para value went from 0% sparse (no effect
#      at all) at p=300 genes to ~50% sparse at p=1000, on the same real
#      dataset. At this project's real scale (~8000 genes), the grid's
#      whole range was very plausibly sitting in an already-fully-
#      saturated regime (matching the ~94% sparsity seen at para=0.001,
#      itself 50x smaller than the OLDER grid's minimum) -- a single fixed
#      penalty grid can't be calibrated consistently across datasets with
#      different gene counts, or even reliably within one.
#   2. `lambda` is not safe to hold near-zero: at a fixed para, raising
#      lambda from 1e-6 to 1 flipped a ~90%-sparse real fit to fully
#      dense. para and lambda are not independent knobs, and elastic net's
#      whole reason for existing over plain LASSO is the ridge term
#      helping in exactly the p >> n regime this project always runs in
#      (thousands of genes, tens-to-hundreds of samples) -- Zou & Hastie
#      2006.
# Fix: `sparse="varnum"` instead -- `para` is now a FRACTION of genes
# (0-1], converted to an absolute nonzero-loadings-per-component target at
# RUN TIME using this dataset's own real gene count (`nrow(mat)`, below).
# This is portable across datasets by construction (a fraction always
# means the same thing regardless of gene count) and directly
# interpretable, unlike an abstract penalty scale. `lambda` is now a real
# swept dimension (config/*.yml's `methods.spca.lambda`), not a fixed
# default.
#
# `K`, `para`, `lambda` are independent, CROSSED dimensions
# (`expand.grid()` below handles this the same way it handles any other
# swept arguments -- no special-casing needed).
#
# Operates on the single matrix `mat` (set as a global object by the
# orchestrator script) -- there is no basis argument.

run_spca_job <- function(K, para, type = "predictor", sparse = "varnum",
                          use.corr = FALSE, lambda = 1e-6, max.iter = 200, eps.conv = 1e-3) {
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
  # para/sparse/lambda: revised 2026-09-29, see run_spca_job()'s header for
  # the real-data evidence -- para is now a fraction of genes (0-1] under
  # sparse="varnum", not an absolute penalty; lambda is a real swept
  # dimension, not a fixed near-zero default. This fallback matches every
  # real dataset config's grid (config/*.yml) so a new config that omits
  # these fields gets the same sweep, not the old penalty-based one.
  defaults = list(K = 2:20, para = seq(0.1, 1.0, by = 0.1),
                   sparse = "varnum", lambda = c(1e-6, 0.1, 1.0),
                   use.corr = FALSE, max.iter = 200, eps.conv = 1e-3),
  build_grid = function(p) {
    expand.grid(K = p$K, para = p$para, sparse = p$sparse, use.corr = p$use.corr,
                lambda = p$lambda, max.iter = p$max.iter, eps.conv = p$eps.conv,
                stringsAsFactors = FALSE)
  }
)

# ICA (fastICA::fastICA()) method registry.
#
# fastICA's fixed-point iteration starts from a random unmixing matrix --
# so, per design, this genuinely needs a seed sweep to assess cross-seed
# factor stability, same as NMF/CoGAPS.
#
# `seed` has no fastICA() equivalent -- passed to `set.seed()` before the
# call, same idiom as R/methods/nmf.R. fastICA's own `method` argument is
# fixed to `"C"` (the faster of its two implementations), not configurable.
#
# `alpha` is deliberately kept a STATIC (length-1) value in `defaults`
# below, never crossed with `n.comp`/`seed` -- letting it sweep would (a)
# multiply job count for no rank-selection benefit, and (b) break the
# app's `same_rank` cross-seed-stability comparison, which only compares
# `n.comp`, not `alpha` (see R/lib/ingest/pairs.R). The `stopifnot()` in
# `build_grid` below enforces this explicitly instead of relying on
# nobody ever configuring a multi-value `alpha` override.
#
# PCA-PREWHITENING: fastICA()'s own whitening step forms the full genes x
# genes covariance matrix and takes its SVD -- an O(genes^3) cost,
# confirmed prohibitively slow at real gene counts. The fix: reduce to a
# small number of principal components first (via prcomp()) and run
# fastICA() on that reduced (samples x n_pcs) matrix instead of the full
# (samples x genes) one -- a standard, well-precedented approximation for
# ICA-on-genomics. `n_pcs` is `n.comp` plus a small margin (5 extra PCs,
# capped at samples/genes - 1).
#
# `fit$A` (n.comp x n_pcs, the mixing matrix in PC space) is mapped back to
# gene space as `pcfit$rotation %*% t(fit$A)` (genes x n.comp) -- `fit$S`
# (samples x n.comp) is the sample score matrix directly.

run_ica_seed_sweep_job <- function(n.comp, alpha = 1, seed, fun = "logcosh",
                                    alg.typ = "parallel", maxit = 200, tol = 1e-4, ...) {
  library(fastICA)
  set.seed(seed)
  n_pcs <- min(n.comp + 5, ncol(mat) - 1, nrow(mat) - 1)
  pcfit <- prcomp(t(mat), rank. = n_pcs, center = TRUE, scale. = FALSE)
  fit <- fastICA(pcfit$x, n.comp = n.comp, alpha = alpha, fun = fun,
                  alg.typ = alg.typ, method = "C", maxit = maxit, tol = tol, ...)

  loadings <- pcfit$rotation %*% t(fit$A)
  rownames(loadings) <- rownames(mat)
  colnames(loadings) <- paste0("Component_", seq_len(ncol(loadings)))
  scores <- fit$S
  rownames(scores) <- colnames(mat)
  colnames(scores) <- colnames(loadings)

  recon_pc <- fit$S %*% fit$A
  recon <- sweep(recon_pc %*% t(pcfit$rotation), 2, pcfit$center, "+")   # samples x genes
  mse <- mean((t(mat) - recon)^2)

  list(rank = n.comp, seed = seed, mse = mse, loadings = loadings, scores = scores,
       # W (unmixing matrix): fastICA's own output has no self-reported
       # convergence/quality diagnostic -- W's orthonormality is the cheap
       # real check (a genuinely converged W should be close to
       # orthonormal). K (whitening matrix): the other half of fastICA's
       # own return value, kept for the same reason. pcfit$sdev: the
       # prewhitening PCA step's full spectrum -- the only way to verify
       # how much variance the n_pcs-dimensional reduction (this file's
       # header comment calls it "a standard, well-precedented
       # approximation") actually retained for a given dataset/rank,
       # currently unverifiable after the fact.
       W = fit$W, K = fit$K, prewhiten_sdev = pcfit$sdev)
}

ica_registry <- list(
  needs_nonneg = FALSE,
  global_object = "mat",
  jobname = "ica_grid",
  fn = run_ica_seed_sweep_job,
  pkgs = "fastICA",
  defaults = list(n.comp = c(5, 7, 10, 15, 20), alpha = 1,
                   seed = c(42, 123, 456, 7, 99, 2024, 8675309, 271828, 31415, 90210),
                   fun = "logcosh", alg.typ = "parallel", maxit = 200, tol = 1e-4),
  build_grid = function(p) {
    if (length(p$alpha) != 1) {
      stop("ica's `alpha` must stay a single static value, not swept -- see this file's header ",
           "for why (breaks the app's same_rank cross-seed comparison). Got length ", length(p$alpha))
    }
    expand.grid(n.comp = p$n.comp, alpha = p$alpha, seed = p$seed, fun = p$fun,
                alg.typ = p$alg.typ, maxit = p$maxit, tol = p$tol, stringsAsFactors = FALSE)
  }
)

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
#
# ICASSO (Himberg, Hyvärinen & Esposito 2004, "Validating the independent
# components of neuroimaging time series via clustering and
# visualization," NeuroImage): the standard rigor-check for ICA stability
# -- run ICA many times under two kinds of randomization (reinitialization
# AND bootstrap-resampling the data), then cluster every recovered
# component across all runs by |correlation| and score each cluster's
# robustness via Iq = intra-cluster similarity - inter-cluster similarity.
# `bootstrap` below adds the second randomization axis on top of the
# seed's own reinitialization -- crossed with n.comp/seed in build_grid(),
# doubling job count per n.comp (10 reinit-only + 10 bootstrap+reinit,
# matching Himberg et al.'s typical ~20-50 total ICASSO run count). The
# clustering/Iq computation itself lives in R/lib/ingest/icasso.R, run at
# ingest time over every ok fit (seed-sweep and bootstrap alike) at a
# given rank -- see that file's header.
#
# When `bootstrap = TRUE`, the job resamples `mat`'s SAMPLE COLUMNS with
# replacement (using this job's own `seed`, so it's reproducible) before
# the usual prewhiten -> fastICA pipeline. The resulting `scores` would
# have duplicate sample ids (not a real 1:1 sample correspondence), so
# it's dropped (NULL) rather than stored malformed -- a bootstrap fit
# contributes its LOADINGS to the ICASSO similarity/clustering
# computation only, never to the app's per-sample drill-down views. Its
# `mse` is reconstruction error against the resampled data, not a claim
# about the real dataset -- app-side queries exclude bootstrap fits from
# every reconstruction-quality/representative-fit comparison (see
# R/lib/ingest/db.R's `fits.bootstrap` column comment).
#
# Per-component excess kurtosis (Lee & Batzoglou 2003: real
# biological-process ICA components are expected to be highly
# non-Gaussian/"super-Gaussian" -- positive excess kurtosis -- unlike
# Gaussian noise) is computed directly on fastICA's raw source output
# `fit$S`, available for bootstrap fits too (unlike deriving it post hoc
# from `scores_file`, which bootstrap fits don't have) -- stored via
# R/lib/ingest/extract.R's `diag$kurtosis` into the `ica_component_kurtosis`
# table.

run_ica_seed_sweep_job <- function(n.comp, alpha = 1, seed, fun = "logcosh",
                                    alg.typ = "parallel", maxit = 200, tol = 1e-4,
                                    bootstrap = FALSE, ...) {
  library(fastICA)
  set.seed(seed)

  mat_use <- mat
  if (isTRUE(bootstrap)) {
    boot_idx <- sample(ncol(mat), ncol(mat), replace = TRUE)
    mat_use <- mat[, boot_idx, drop = FALSE]
  }

  n_pcs <- min(n.comp + 5, ncol(mat_use) - 1, nrow(mat_use) - 1)
  pcfit <- prcomp(t(mat_use), rank. = n_pcs, center = TRUE, scale. = FALSE)
  fit <- fastICA(pcfit$x, n.comp = n.comp, alpha = alpha, fun = fun,
                  alg.typ = alg.typ, method = "C", maxit = maxit, tol = tol, ...)

  loadings <- pcfit$rotation %*% t(fit$A)
  rownames(loadings) <- rownames(mat)
  colnames(loadings) <- paste0("Component_", seq_len(ncol(loadings)))

  # Non-Gaussianity of fastICA's own raw source estimate -- see this
  # file's header. Computed on fit$S directly (before any PCA-rotation
  # mapping), matching the formula R/backfill_diagnostics.R's
  # backfill_ica_kurtosis() uses on scores_file for older fits.
  mu <- colMeans(fit$S)
  centered_S <- sweep(fit$S, 2, mu, "-")
  sigma <- sqrt(colMeans(centered_S^2))
  ok <- sigma > 0
  kurtosis <- rep(NA_real_, ncol(fit$S))
  kurtosis[ok] <- colMeans(centered_S[, ok, drop = FALSE]^4) / sigma[ok]^4

  if (isTRUE(bootstrap)) {
    scores <- NULL
  } else {
    scores <- fit$S
    rownames(scores) <- colnames(mat)
    colnames(scores) <- colnames(loadings)
  }

  recon_pc <- fit$S %*% fit$A
  recon <- sweep(recon_pc %*% t(pcfit$rotation), 2, pcfit$center, "+")   # samples x genes
  mse <- mean((t(mat_use) - recon)^2)

  list(rank = n.comp, seed = seed, bootstrap = as.integer(isTRUE(bootstrap)),
       mse = mse, loadings = loadings, scores = scores, kurtosis = kurtosis,
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
                   bootstrap = c(FALSE, TRUE),
                   fun = "logcosh", alg.typ = "parallel", maxit = 200, tol = 1e-4),
  build_grid = function(p) {
    if (length(p$alpha) != 1) {
      stop("ica's `alpha` must stay a single static value, not swept -- see this file's header ",
           "for why (breaks the app's same_rank cross-seed comparison). Got length ", length(p$alpha))
    }
    expand.grid(n.comp = p$n.comp, alpha = p$alpha, seed = p$seed, bootstrap = p$bootstrap,
                fun = p$fun, alg.typ = p$alg.typ, maxit = p$maxit, tol = p$tol,
                stringsAsFactors = FALSE)
  }
)

# NMF (NNLM::nnmf) method registry.
#
# `k` is nnmf()'s real rank argument. `seed` has no nnmf() equivalent --
# this framework calls `set.seed(seed)` before `nnmf()` itself. The job
# function takes `...` so any nnmf() argument not explicitly named as a
# formal here (alpha, beta, method, loss, max.iter, ...) still passes
# through if you add it to `defaults`/`build_grid` below and wire it into
# the `expand.grid()` call -- see R/README.md's "Adding a new method".
#
# Operates on the single matrix `mat_nn` (non-negative-shifted, set as a
# global object by the orchestrator script) -- there is no basis argument.

run_nmf_seed_sweep_job <- function(k, seed, max.iter = 10000, verbose = 0L, ...) {
  library(NNLM)
  set.seed(seed)
  fit <- nnmf(mat_nn, k = k, max.iter = max.iter, verbose = verbose, ...)
  rownames(fit$W) <- rownames(mat_nn)
  colnames(fit$W) <- paste0("Pattern_", seq_len(ncol(fit$W)))
  colnames(fit$H) <- colnames(mat_nn)
  rownames(fit$H) <- colnames(fit$W)

  recon <- fit$W %*% fit$H
  mse <- mean((mat_nn - recon)^2)
  list(rank = k, seed = seed, mse = mse, W = fit$W, H = fit$H,
       # convergence diagnostics -- NOT recoverable from W/H alone, needed
       # to tell whether nnmf() actually converged or just hit max.iter
       n.iteration = fit$n.iteration, target.loss = fit$target.loss,
       average.epochs = fit$average.epochs)
}

# 10-seed default reused by ica/cogaps below -- not shared code on
# purpose (each method script is self-contained), just a coincidentally
# common choice across every real dataset config so far.
nmf_registry <- list(
  needs_nonneg = TRUE,
  global_object = "mat_nn",
  jobname = "nmf_grid",
  fn = run_nmf_seed_sweep_job,
  pkgs = "NNLM",
  defaults = list(k = 5:20, seed = c(42, 123, 456, 7, 99, 2024, 8675309, 271828, 31415, 90210)),
  build_grid = function(p) expand.grid(k = p$k, seed = p$seed, stringsAsFactors = FALSE)
)

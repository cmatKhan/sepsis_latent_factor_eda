# NMF (NNLM::nnmf()). `k` is the rank. `seed` has no nnmf() argument: the fit
# calls set.seed() first. `n.threads` (OpenMP) defaults to the nmf
# controller's `cpus_per_task` and is left out of a fit's identity. Other
# nnmf() arguments pass through `...` once added to `defaults` and
# build_grid(). Fits the non-negative shift `mat_nn`.
#
# The fit targets depend on this file's registry and fit function: editing
# their code refits the method (comments don't count). Docs: Methods.

#' Fit NMF for one (k, seed)
#'
#' @param k Rank.
#' @param seed Random seed (set.seed() before fitting).
#' @param max.iter,verbose,n.threads,... Passed to NNLM::nnmf().
#' @return `list(rank, seed, mse, W, H, n.iteration, target.loss,
#'   average.epochs, nnmf, elapsed)`: W is genes x patterns, H patterns x
#'   samples, `nnmf` the rest of nnmf()'s output (per-iteration traces,
#'   run time); or `list(rank, seed, mse = NA, error, elapsed)` on failure.
run_nmf_seed_sweep_job <- function(k, seed, max.iter = 10000, verbose = 0L, n.threads = 1L, ...) capture_fit(base = list(rank = k, seed = seed, mse = NA_real_), {
  library(NNLM)
  set.seed(seed)
  fit <- nnmf(mat_nn, k = k, max.iter = max.iter, verbose = verbose, n.threads = n.threads, ...)
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
       average.epochs = fit$average.epochs,
       # the rest of nnmf()'s return value (per-iteration mse/mkl/target-loss
       # traces, run.time, options, call) -- small, and kept so ingestion
       # can use any of it without refitting. W/H are above.
       nnmf = unclass(fit)[setdiff(names(fit), c("W", "H"))])
})

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
  resource_defaults = function(p, slurm_cfg) {
    if (is.null(p[["n.threads"]])) p$n.threads <- slurm_cfg$cpus_per_task
    p
  },
  build_grid = function(p) expand.grid(k = p$k, seed = p$seed, n.threads = p$n.threads, stringsAsFactors = FALSE)
)

# Ingest contract (R/targets/ingest.R): kept apart from nmf_registry, which
# the fit targets depend on, so editing it never refits.
nmf_ingest <- list(
  family = "seed_sweep", sign_ambiguous = FALSE, has_loadings = TRUE,
  resource_params = "n.threads",
  extract = function(result, params) {
    extra <- result$nnmf %||% result[c("n.iteration", "target.loss", "average.epochs")]
    last <- function(x) if (length(x)) as.numeric(utils::tail(x, 1)) else NA_real_
    fit_record_from(
      rank = params$k, seed = params$seed, mse = result$mse,
      loadings = result$W, scores = if (!is.null(result$H)) t(result$H),
      diag = extra,
      metrics = c(n_iteration = last(extra$n.iteration), target_loss = last(extra$target.loss),
                  average_epochs = last(extra$average.epochs),
                  run_seconds = if (!is.null(extra$run.time)) unname(extra$run.time[["elapsed"]]) else NA_real_))
  }
)

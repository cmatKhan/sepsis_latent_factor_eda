# PCA (prcomp()). Deterministic, and its components are nested: the rank-k
# solution is the first k columns of any higher-rank one. So PCA is fit once
# per dataset at the largest configured `rank` (a list's maximum is used);
# `mse_by_rank` gives every truncation's reconstruction error, and choosing
# the first n components is left to ingestion and the app. `rank` stands for
# prcomp()'s `rank.`. Fits the matrix `mat`.
#
# The fit targets depend on this file's registry and fit function: editing
# their code refits the method (comments don't count). Docs: Methods.

#' Fit PCA at one rank
#'
#' @param rank Number of components kept.
#' @return `list(rank, seed = NA, mse, mse_by_rank, rotation, scores, sdev,
#'   center, elapsed)`; `rotation` is genes x components, `sdev` the full
#'   spectrum.
run_pca_seed_sweep_job <- function(rank) capture_fit(base = list(rank = rank, mse = NA_real_), {
  fit <- prcomp(t(mat), rank. = rank, center = TRUE, scale. = FALSE)

  recon <- fit$x %*% t(fit$rotation)
  recon <- sweep(recon, 2, fit$center, "+")
  mse <- mean((t(mat) - recon)^2)

  # Reconstruction MSE of the rank-k truncation, k = 1..rank, from the full
  # spectrum: the centered data's squared singular values are
  # sdev^2 * (n - 1), and dropping components k+1.. leaves exactly their
  # sum as residual sum of squares (Eckart-Young). mse_by_rank[rank] == mse.
  n <- ncol(mat)
  ss <- fit$sdev^2 * (n - 1)
  mse_by_rank <- vapply(seq_len(rank), function(k) sum(ss[-seq_len(k)]), numeric(1)) / length(mat)

  list(rank = rank, seed = NA_integer_, mse = mse, mse_by_rank = mse_by_rank,
       rotation = fit$rotation, scores = fit$x,
       # sdev: prcomp()'s FULL per-component spectrum, regardless of
       # rank. truncation -- the only way to get real per-component/
       # cumulative proportion-of-variance-explained (what
       # summary.prcomp()/screeplot() report), never recoverable after
       # the fact from the rank-truncated rotation/scores alone. center:
       # needed alongside sdev/rotation to reconstruct a real `prcomp`-
       # classed object for projectR's PCA-mode dispatch (see
       # R/ingest_jobs/projectr_job.R's pca branch, which currently
       # recomputes an equivalent by hand from the cached matrix instead).
       sdev = fit$sdev, center = fit$center)
})

pca_registry <- list(
  needs_nonneg = FALSE,
  global_object = "mat",
  jobname = "pca_grid",
  fn = run_pca_seed_sweep_job,
  pkgs = character(0),
  defaults = list(rank = 20),
  # One fit at the largest configured rank (see header).
  build_grid = function(p) data.frame(rank = max(p$rank))
)

# Ingest contract (R/targets/ingest.R): kept apart from pca_registry, which
# the fit targets depend on, so editing it never refits.
pca_ingest <- list(
  family = "seed_sweep", sign_ambiguous = TRUE, has_loadings = TRUE,
  resource_params = character(0),
  extract = function(result, params) {
    fit_record_from(
      rank = params$rank, mse = result$mse,
      loadings = result$rotation, scores = result$scores,
      diag = list(sdev = result$sdev, center = result$center, mse_by_rank = result$mse_by_rank))
  }
)

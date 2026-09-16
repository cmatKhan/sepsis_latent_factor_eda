# PCA method registry.
#
# PCA (via prcomp/svd) is deterministic -- there is no seed-stability
# question in the NMF/CoGAPS sense. It's swept purely by rank so
# rank-vs-reconstruction-error curves are directly comparable to the
# other methods' output.
#
# `defaults$rank` below is a documented exception to "config key == tool
# argument name": prcomp()'s real argument is `rank.` (trailing dot -- a
# base-R naming quirk to avoid colliding with the `rank()` function, not a
# meaningful concept worth reproducing verbatim in the config).
#
# Operates on the single matrix `mat` (set as a global object by the
# orchestrator script) -- there is no basis argument.

run_pca_seed_sweep_job <- function(rank) {
  fit <- prcomp(t(mat), rank. = rank, center = TRUE, scale. = FALSE)

  recon <- fit$x %*% t(fit$rotation)
  recon <- sweep(recon, 2, fit$center, "+")
  mse <- mean((t(mat) - recon)^2)

  list(rank = rank, seed = NA_integer_, mse = mse, rotation = fit$rotation, scores = fit$x)
}

pca_registry <- list(
  needs_nonneg = FALSE,
  global_object = "mat",
  jobname = "pca_grid",
  fn = run_pca_seed_sweep_job,
  pkgs = character(0),
  defaults = list(rank = 2:20),
  build_grid = function(p) expand.grid(rank = p$rank, stringsAsFactors = FALSE)
)

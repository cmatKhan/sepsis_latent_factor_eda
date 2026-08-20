# PCA method registry.
#
# PCA (via prcomp/svd) is deterministic -- there is no seed-stability
# question in the NMF/CoGAPS sense. Its "seed_sweep" family therefore
# degenerates to a single fit per rank: included anyway so the same
# orchestrator shape works uniformly across methods, and so
# rank-vs-reconstruction-error curves are directly comparable to the other
# methods' full-grid output. Masking-CV is retained as PCA's genuine
# rank-selection design (held-out-entry reconstruction via a rank-k SVD
# approximation).
#
# Operates on the single matrix `mat` (set as a global object by the
# orchestrator script) -- there is no basis argument.

pca_stability_designs <- c("seed_sweep", "masking_cv")

run_pca_seed_sweep_job <- function(rank) {
  fit <- prcomp(t(mat), rank. = rank, center = TRUE, scale. = FALSE)

  recon <- fit$x %*% t(fit$rotation)
  recon <- sweep(recon, 2, fit$center, "+")
  mse <- mean((t(mat) - recon)^2)

  list(rank = rank, seed = NA_integer_, mse = mse, rotation = fit$rotation, scores = fit$x)
}

run_pca_masking_cv_job <- function(rank) {
  mat_masked <- mat
  mat_masked[mask_idx] <- NA

  # Simple iterative low-rank imputation (Gabriel/EM-style): fill NAs with
  # row means, refit rank-k SVD, refill NAs from the reconstruction, repeat.
  # This is the standard "PCA masking-CV" idiom (same spirit as
  # missMDA::imputePCA) without adding a new package dependency.
  row_means <- rowMeans(mat_masked, na.rm = TRUE)
  filled <- mat_masked
  filled[mask_idx] <- row_means[row(mat_masked)[mask_idx]]

  for (iter in seq_len(25)) {
    fit <- prcomp(t(filled), rank. = rank, center = TRUE, scale. = FALSE)
    recon <- fit$x %*% t(fit$rotation)
    recon <- sweep(recon, 2, fit$center, "+")
    recon <- t(recon)
    filled[mask_idx] <- recon[mask_idx]
  }

  mse <- mean((mat[mask_idx] - filled[mask_idx])^2)
  list(rank = rank, mse = mse)
}

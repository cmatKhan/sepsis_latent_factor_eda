# NMF (NNLM::nnmf) method registry -- mechanics unchanged from
# R/legacy/nmf_rslurm.R and R/legacy/nmf_maskcv_rslurm.R, parameterized by
# metadata rather than hardcoded literals.
#
# Operates on the single matrix `mat_nn` (non-negative-shifted, set as a
# global object by the orchestrator script) -- there is no basis argument.

nmf_stability_designs <- c("seed_sweep", "masking_cv")

run_nmf_seed_sweep_job <- function(rank, seed) {
  library(NNLM)
  set.seed(seed)
  fit <- nnmf(mat_nn, k = rank, max.iter = 10000, verbose = 0L)
  rownames(fit$W) <- rownames(mat_nn)
  colnames(fit$W) <- paste0("Pattern_", seq_len(ncol(fit$W)))
  colnames(fit$H) <- colnames(mat_nn)
  rownames(fit$H) <- colnames(fit$W)

  recon <- fit$W %*% fit$H
  mse <- mean((mat_nn - recon)^2)
  list(rank = rank, seed = seed, mse = mse, W = fit$W, H = fit$H)
}

run_nmf_masking_cv_job <- function(rank) {
  library(NNLM)
  mat_masked <- mat_nn
  mat_masked[mask_idx] <- NA

  fit <- nnmf(mat_masked, k = rank, max.iter = 10000, verbose = 0L)
  recon <- fit$W %*% fit$H
  mse <- mean((mat_nn[mask_idx] - recon[mask_idx])^2)

  list(rank = rank, mse = mse)
}

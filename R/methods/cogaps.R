# CoGAPS method registry -- mechanics unchanged from
# R/legacy/cogaps_rslurm.R and R/legacy/cogaps_maskcv_rslurm.R, parameterized
# by metadata rather than hardcoded literals.
#
# Operates on the single matrix `mat_nn` (non-negative-shifted, set as a
# global object by the orchestrator script) -- there is no basis argument.
#
# Masking is done via CoGAPS's `uncertainty` matrix, not NA (CoGAPS has no
# documented native missing-data handling): kept cells get the textbook
# default "10% of signal, floored" uncertainty, held-out cells get an
# enormous uncertainty so they contribute ~nothing to the fit.

cogaps_stability_designs <- c("seed_sweep", "masking_cv")

run_cogaps_seed_sweep_job <- function(rank, seed, n_iterations, nsets) {
  library(CoGAPS)
  params <- CogapsParams(
    nPatterns   = rank,
    nIterations = n_iterations,
    seed        = seed,
    distributed = "genome-wide"
  )
  params <- setDistributedParams(params, nSets = nsets)
  result <- tryCatch(
    CoGAPS(mat_nn, params = params, nThreads = 1),
    error = function(e) NULL
  )
  if (is.null(result)) {
    return(list(rank = rank, seed = seed, mse = NA_real_, result = NULL))
  }

  recon <- result@featureLoadings %*% t(result@sampleFactors)
  mse   <- mean((mat_nn - recon)^2)
  list(rank = rank, seed = seed, mse = mse, result = result)
}

run_cogaps_masking_cv_job <- function(rank, alpha, maskcv_seed, n_iterations, nsets) {
  library(CoGAPS)
  unc <- uncertainty
  unc[mask_idx] <- 1e4 * max(mat_nn)

  params <- CogapsParams(
    nPatterns   = rank,
    nIterations = n_iterations,
    seed        = maskcv_seed,
    alphaA      = alpha,
    alphaP      = alpha,
    distributed = "genome-wide"
  )
  params <- setDistributedParams(params, nSets = nsets)
  result <- tryCatch(
    CoGAPS(mat_nn, params = params, uncertainty = unc, nThreads = 1),
    error = function(e) NULL
  )
  if (is.null(result)) return(list(rank = rank, alpha = alpha, mse = NA_real_))

  recon <- result@featureLoadings %*% t(result@sampleFactors)
  mse   <- mean((mat_nn[mask_idx] - recon[mask_idx])^2)
  list(rank = rank, alpha = alpha, mse = mse)
}

#' Default-uncertainty ("10% of signal, floored") matrix -- shared setup
#' needed before submitting the masking-CV family (held-out cells get
#' overwritten to an enormous value inside the job function itself).
build_uncertainty <- function(mat_nn) pmax(0.1 * mat_nn, 0.1 * mean(mat_nn))

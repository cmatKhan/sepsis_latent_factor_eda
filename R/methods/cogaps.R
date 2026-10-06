# CoGAPS. Its arguments span three call sites, so the config has three
# sub-blocks: `params` -> CogapsParams() (nPatterns, seed, nIterations,
# distributed, ...), `distributed_params` -> setDistributedParams() (nSets,
# cut, ...), `run` -> CoGAPS() itself (nThreads, ...). `distributed_params`
# and `run` pass through as-is; a new CogapsParams argument must be added to
# `defaults` and cogaps_build_grid(). nSets and nThreads default to the
# cogaps controller's `cpus_per_task` and are left out of a fit's identity.
# Fits the non-negative shift `mat_nn`.
#
# The fit targets depend on this file's registry and fit function: editing
# their code refits the method (comments don't count). Docs: Methods.

#' Fit CoGAPS for one grid row
#'
#' @param params Named list for CogapsParams() (nPatterns, seed, ...).
#' @param distributed_params Named list for setDistributedParams().
#' @param run Named list of extra CoGAPS() arguments.
#' @return `list(rank, seed, mse, result, elapsed)` with `result` the
#'   `CogapsResult`; on failure `result = NULL` plus `error`.
run_cogaps_seed_sweep_job <- function(params, distributed_params, run) capture_fit(
    base = list(rank = params$nPatterns, seed = params$seed, mse = NA_real_, result = NULL), {
  library(CoGAPS)
  p <- do.call(CogapsParams, params)
  if (length(distributed_params) > 0) p <- do.call(setDistributedParams, c(list(p), distributed_params))

  # data = quote(mat_nn), not the matrix itself: CoGAPS prints
  # deparse(substitute(data)) to label the run, and with the value spliced
  # into the call that is the whole matrix -- minutes to hours of log
  # output per fit. The symbol is evaluated here, where mat_nn is bound.
  result <- do.call(CoGAPS, c(list(data = quote(mat_nn), params = p), run))

  recon <- result@featureLoadings %*% t(result@sampleFactors)
  mse   <- mean((mat_nn - recon)^2)
  list(rank = params$nPatterns, seed = params$seed, mse = mse, result = result)
})

#' Default nSets and nThreads from the worker's CPUs
#'
#' @param p Resolved CoGAPS parameters.
#' @param slurm_cfg `list(cpus_per_task)`.
#' @return `p` with `distributed_params$nSets` and `run$nThreads` filled in
#'   where unset.
cogaps_resource_defaults <- function(p, slurm_cfg) {
  p$distributed_params <- p$distributed_params %||% list()
  p$run <- p$run %||% list()
  if (is.null(p$distributed_params[["nSets"]])) p$distributed_params$nSets <- slurm_cfg$cpus_per_task
  if (is.null(p$run[["nThreads"]])) p$run$nThreads <- slurm_cfg$cpus_per_task
  p
}

#' CoGAPS's parameter grid
#'
#' Crosses `nPatterns` x `seed`; each row carries its own `params` list plus
#' the shared `distributed_params` and `run`.
#'
#' @param p Resolved CoGAPS parameters.
#' @return Data frame with list-columns `params`, `distributed_params`, `run`.
cogaps_build_grid <- function(p) {
  base <- expand.grid(nPatterns = p$params$nPatterns, seed = p$params$seed, stringsAsFactors = FALSE)
  n <- nrow(base)
  params_col <- lapply(seq_len(n), function(i) {
    list(nPatterns = base$nPatterns[i], seed = base$seed[i],
         nIterations = p$params$nIterations, distributed = p$params$distributed)
  })
  data.frame(
    params = I(params_col),
    distributed_params = I(rep(list(p$distributed_params %||% list()), n)),
    run = I(rep(list(p$run %||% list()), n))
  )
}

cogaps_registry <- list(
  needs_nonneg = TRUE,
  global_object = "mat_nn",
  jobname = "cogaps_grid",
  fn = run_cogaps_seed_sweep_job,
  pkgs = "CoGAPS",
  defaults = list(
    params = list(
      nPatterns = 5:10,
      seed = c(42, 123, 456, 7, 99, 2024, 8675309, 271828, 31415, 90210),
      nIterations = 15000, distributed = "genome-wide"
    ),
    distributed_params = list(), run = list()
  ),
  resource_defaults = cogaps_resource_defaults,
  build_grid = cogaps_build_grid
)

# Ingest contract (R/targets/ingest.R): kept apart from cogaps_registry,
# which the fit targets depend on, so editing it never refits. nSets and
# nThreads only set parallelism, so they're not part of a fit's identity.
cogaps_ingest <- list(
  family = "seed_sweep", sign_ambiguous = FALSE, has_loadings = TRUE,
  resource_params = c("distributed_params.nSets", "run.nThreads"),
  extract = function(result, params) {
    res <- result$result
    if (is.null(res)) return(fit_failed("CoGAPS returned no result", rank = params$nPatterns, seed = params$seed))
    chisq <- tryCatch(as.numeric(res@metadata$meanChiSq), error = function(e) NA_real_)
    fit_record_from(
      rank = params$nPatterns, seed = params$seed, mse = result$mse,
      loadings = res@featureLoadings, scores = res@sampleFactors,
      diag = list(loading_sd = res@loadingStdDev, factor_sd = res@factorStdDev),
      raw = res,
      metrics = c(mean_chisq = if (length(chisq) == 1) chisq else NA_real_))
  }
)

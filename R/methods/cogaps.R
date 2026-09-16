# CoGAPS method registry.
#
# CoGAPS genuinely spans THREE call sites with distinct (and in one
# documented case, colliding) argument namespaces:
#   - CogapsParams(...)            -- ~25 slots (nPatterns, nIterations,
#                                      seed, alphaA, alphaP, distributed,
#                                      nSets, checkpointInFile, ...)
#   - setDistributedParams(params, nSets=, cut=, minNS=, maxNS=)
#   - CoGAPS(data, params=, nThreads=, uncertainty=, checkpointOutFile=, ...)
# The job function below do.call()s each resolved sub-block into the
# matching function. `distributed_params`/`run` are carried through as
# OPAQUE named lists (never inspected key-by-key by build_grid), so any
# current or future setDistributedParams()/CoGAPS() argument is still
# configurable there with zero code changes. CogapsParams arguments are
# NOT similarly generic: `cogaps_build_grid()` below constructs each row's
# `params` list by explicitly naming every key it forwards (nPatterns,
# seed, nIterations, distributed) -- adding ANY new CogapsParams argument
# (swept or static) means adding it to `defaults` AND to
# `cogaps_build_grid()`, matching how every other R/methods/*.R file's
# build_grid() works (it can't go under `distributed_params`/`run`
# instead -- those are the wrong call site).
#
# CogapsParams()'s own arguments live directly under `cogaps.params:` in
# the dataset config -- `distributed_params`/`run` are separate sibling
# keys, one level of nesting each, matching setDistributedParams()/
# CoGAPS()'s own distinct argument sets. (Unlike every other method here,
# this nesting is intrinsic to CoGAPS's three call sites, not a schema
# wrapper -- there's no flatter form available.)
#
# `defaults` sweeps `nPatterns` x `seed` (crossed, like any other pair of
# arguments).
#
# Operates on the single matrix `mat_nn` (non-negative-shifted, set as a
# global object by the orchestrator script) -- there is no basis argument.

run_cogaps_seed_sweep_job <- function(params, distributed_params, run) {
  library(CoGAPS)
  p <- do.call(CogapsParams, params)
  if (length(distributed_params) > 0) p <- do.call(setDistributedParams, c(list(p), distributed_params))

  result <- tryCatch(
    do.call(CoGAPS, c(list(data = mat_nn, params = p), run)),
    error = function(e) NULL
  )
  if (is.null(result)) {
    return(list(rank = params$nPatterns, seed = params$seed, mse = NA_real_, result = NULL))
  }

  recon <- result@featureLoadings %*% t(result@sampleFactors)
  mse   <- mean((mat_nn - recon)^2)
  list(rank = params$nPatterns, seed = params$seed, mse = mse, result = result)
}

#' `nSets`/`nThreads` default to the method's slurm.cogaps.cpus_per_task
#' unless set explicitly (under `distributed_params`/`run` in the config).
#' Applied to the resolved (defaults + config overrides) parameter list,
#' before build_grid() runs.
cogaps_resource_defaults <- function(p, slurm_cfg) {
  p$distributed_params <- p$distributed_params %||% list()
  p$run <- p$run %||% list()
  if (is.null(p$distributed_params[["nSets"]])) p$distributed_params$nSets <- slurm_cfg$cpus_per_task
  if (is.null(p$run[["nThreads"]])) p$run$nThreads <- slurm_cfg$cpus_per_task
  p
}

#' Crosses nPatterns x seed (the only two CogapsParams keys ever swept);
#' nIterations/distributed are carried through unchanged on every row.
#' Every key referenced here is one `defaults` below actually declares --
#' a stray/typo'd config key under `params:` is silently ignored rather
#' than automatically becoming a real CogapsParams argument, matching how
#' every other R/methods/*.R file's build_grid() behaves.
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

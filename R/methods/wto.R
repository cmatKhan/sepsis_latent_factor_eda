# wTO method registry -- generalizes R/legacy/rose_wto_rslurm.R's
# wTO.Complete family.
#
# Operates on the single matrix `mat` (set as a global object by the
# orchestrator script) -- there is no basis argument.
#
# wTO has no held-out-entry masking-CV analogue (it's a correlation-network
# method, not a reconstructive factorization), so only a `param_grid` family
# is built, sweeping the bootstrap replicate count `n` (the "structural"
# parameter -- more replicates -> a more stable/precise edge-significance
# estimate) and, optionally, an outer seed controlling wTO's internal
# bootstrap RNG stream (its own stability-across-randomness question).

wto_stability_designs <- c("param_grid")

run_wto_param_job <- function(n, delta, seed = NA_integer_, wto_k) {
  library(wTO)
  if (!is.na(seed)) set.seed(seed)

  res <- wTO.Complete(
    k = wto_k, n = n,
    Data              = as.data.frame(mat),
    Overlap           = rownames(mat),
    method            = "p",
    method_resampling = "Bootstrap",
    pvalmethod        = "BH",
    plot              = FALSE
  )
  list(n = n, delta = delta, seed = seed, result = res)
}

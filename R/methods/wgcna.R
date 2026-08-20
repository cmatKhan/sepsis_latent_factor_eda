# WGCNA method registry -- generalizes R/legacy/wgcna_rslurm.R's per-basis
# network + module detection family.
#
# Operates on the single matrix `mat` (set as a global object by the
# orchestrator script) -- there is no basis argument.
#
# power_grid plays the "structural" parameter role here (higher power ->
# a sparser, more stringent adjacency -- WGCNA's analogue of increasing
# rank/neighborhood size). Unlike the tutorial convention of auto-picking a
# single power via pickSoftThreshold, this sweeps `power_grid` directly so
# the resulting module sets can be compared across powers (the whole point
# of this stability-characterization pipeline). No seed-sweep or
# masking-CV family: blockwiseModules's module detection isn't RNG-driven
# by default, and there's no natural held-out-entry reconstruction concept
# for a correlation network. Module preservation (cross-matrix) is left as
# a documented future extension (see R/legacy/wgcna_rslurm.R family 3), not
# built here.

wgcna_stability_designs <- c("param_grid")

run_wgcna_param_job <- function(power, min_module_size, merge_cut_height,
                                 network_type, wgcna_cpu_per_task) {
  library(WGCNA)
  enableWGCNAThreads(nThreads = wgcna_cpu_per_task)
  options(stringsAsFactors = FALSE)

  datExpr <- t(mat)

  gsg <- goodSamplesGenes(datExpr, verbose = 0)
  if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]

  net <- blockwiseModules(
    datExpr,
    power = power, networkType = network_type, TOMType = network_type,
    minModuleSize = min_module_size, mergeCutHeight = merge_cut_height,
    numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = FALSE,
    nThreads = wgcna_cpu_per_task, verbose = 0
  )

  list(power = power, net = net, genes = colnames(datExpr))
}

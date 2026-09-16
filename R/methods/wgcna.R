# WGCNA method registry -- generalizes R/legacy/wgcna_rslurm.R's per-basis
# network + module detection family.
#
# `power` plays the "structural" parameter role here (higher power -> a
# sparser, more stringent adjacency -- WGCNA's analogue of increasing
# rank/neighborhood size). Unlike the tutorial convention of auto-picking a
# single power via pickSoftThreshold, this sweeps `power` directly so the
# resulting module sets can be compared across powers.
#
# `nThreads` defaults to `slurm.network.wgcna.cpus_per_task` (via
# `resource_defaults` below) unless set explicitly in `defaults`/config.
#
# Operates on the single matrix `mat` (set as a global object by the
# orchestrator script) -- there is no basis argument.
#
# `network = TRUE` below tells R/lib/method_registry.R's discovery (and
# R/create_slurm_bundle.R) that this method's config/slurm entries nest
# under methods$network$wgcna / slurm$network$wgcna instead of the flat
# methods$wgcna / slurm$wgcna every other method uses.

run_wgcna_param_job <- function(power, minModuleSize, mergeCutHeight, networkType,
                                 nThreads = 1, ...) {
  library(WGCNA)
  enableWGCNAThreads(nThreads = nThreads)
  options(stringsAsFactors = FALSE)

  datExpr <- t(mat)

  gsg <- goodSamplesGenes(datExpr, verbose = 0)
  if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]

  net <- blockwiseModules(
    datExpr,
    power = power, networkType = networkType, TOMType = networkType,
    minModuleSize = minModuleSize, mergeCutHeight = mergeCutHeight,
    numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = FALSE,
    nThreads = nThreads, verbose = 0, ...
  )

  list(power = power, net = net, genes = colnames(datExpr), samples = rownames(datExpr))
}

wgcna_registry <- list(
  needs_nonneg = FALSE,
  network = TRUE,
  global_object = "mat",
  jobname = "wgcna_grid",
  fn = run_wgcna_param_job,
  pkgs = "WGCNA",
  defaults = list(power = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20),
                   minModuleSize = 30, mergeCutHeight = 0.25, networkType = "signed"),
  resource_defaults = function(p, slurm_cfg) {
    if (is.null(p[["nThreads"]])) p$nThreads <- slurm_cfg$cpus_per_task
    p
  },
  build_grid = function(p) {
    expand.grid(power = p$power, minModuleSize = p$minModuleSize,
                mergeCutHeight = p$mergeCutHeight, networkType = p$networkType,
                nThreads = p$nThreads, stringsAsFactors = FALSE)
  }
)

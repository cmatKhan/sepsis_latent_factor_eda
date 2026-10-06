# WGCNA (blockwiseModules()). `power` is swept rather than auto-picked with
# pickSoftThreshold(), so module sets can be compared across powers.
# `maxBlockSize = 50000` keeps each network in one block over all genes:
# the default (5000) would split ~7,500 genes into separately built blocks
# whose genes are never compared. `nThreads` defaults to the wgcna
# controller's `cpus_per_task` (WGCNA needs at least 2) and is left out of a
# fit's identity. `network = TRUE`: the config nests under
# `methods.network.wgcna`. Fits the matrix `mat`.
#
# The fit targets depend on this file's registry and fit function: editing
# their code refits the method (comments don't count). Docs: Methods.

#' Build one WGCNA network
#'
#' @param power Soft-thresholding power.
#' @param minModuleSize,mergeCutHeight,networkType,maxBlockSize,nThreads,...
#'   Passed to WGCNA::blockwiseModules().
#' @return `list(power, net, genes, samples, elapsed)`: `net` is
#'   blockwiseModules()'s output (module colors, eigengenes, gene trees),
#'   after the goodSamplesGenes() filter.
run_wgcna_param_job <- function(power, minModuleSize, mergeCutHeight, networkType,
                                 maxBlockSize = 50000, nThreads = 1, ...) capture_fit(base = list(power = power), {
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
    maxBlockSize = maxBlockSize, nThreads = nThreads, verbose = 0, ...
  )

  list(power = power, net = net, genes = colnames(datExpr), samples = rownames(datExpr))
})

wgcna_registry <- list(
  needs_nonneg = FALSE,
  network = TRUE,
  global_object = "mat",
  jobname = "wgcna_grid",
  fn = run_wgcna_param_job,
  pkgs = "WGCNA",
  defaults = list(power = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20),
                   minModuleSize = 30, mergeCutHeight = 0.25, networkType = "signed",
                   maxBlockSize = 50000),
  resource_defaults = function(p, slurm_cfg) {
    if (is.null(p[["nThreads"]])) p$nThreads <- slurm_cfg$cpus_per_task
    p
  },
  build_grid = function(p) {
    expand.grid(power = p$power, minModuleSize = p$minModuleSize,
                mergeCutHeight = p$mergeCutHeight, networkType = p$networkType,
                maxBlockSize = p$maxBlockSize, nThreads = p$nThreads, stringsAsFactors = FALSE)
  }
)

# Ingest contract (R/targets/ingest.R): kept apart from wgcna_registry,
# which the fit targets depend on, so editing it never refits. Modules are
# WGCNA's factors (module 0 = unassigned genes, not a factor). With one
# block (maxBlockSize above the gene count) net$dendrograms[[1]] is the gene
# tree the modules were cut from.
wgcna_ingest <- list(
  family = "param_grid", sign_ambiguous = FALSE, has_loadings = FALSE,
  resource_params = "nThreads",
  extract = function(result, params) {
    net <- result$net
    if (is.null(net) || is.null(net$colors)) {
      return(fit_failed("blockwiseModules returned no module colors", power = params$power))
    }
    cols <- net$colors
    genes <- if (!is.null(names(cols))) names(cols) else result$genes
    modules <- data.frame(gene = as.character(genes), module = as.integer(cols))
    ids <- sort(setdiff(unique(modules$module), 0L))
    scores <- NULL
    if (!is.null(net$MEs) && !is.null(result$samples)) {
      scores <- as.matrix(net$MEs)
      rownames(scores) <- result$samples
    }
    dendro <- NULL
    if (length(net$dendrograms) == 1) {
      hc <- net$dendrograms[[1]]
      dendro <- list(merge = hc$merge, height = hc$height, order = hc$order,
                     labels = as.character(genes[net$blockGenes[[1]]]))
    }
    rec <- fit_record_from(scores = scores, diag = list(n_blocks = length(net$dendrograms)),
                           raw = net, dendrogram = dendro, modules = modules,
                           metrics = c(n_blocks = length(net$dendrograms),
                                       n_unassigned = sum(modules$module == 0L)))
    rec$n_factors <- length(ids)
    rec$factors <- data.frame(factor_index = ids, label = paste0("ME", ids),
                              n_genes = as.integer(table(modules$module)[as.character(ids)]),
                              kurtosis = NA_real_, excess_kurtosis = NA_real_)
    rec
  }
)

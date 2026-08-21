# Per-method extractors: map one job family's result element (one
# `results_<i>.RDS`, unwrapped) + its params row -> a `fits` row and any
# artifacts (loading matrix, module assignments, edge table).
#
# Result shapes (verified against real GSE110487 outputs):
#   pca_grid:      list(rank, seed=NA, mse, rotation, scores)
#   pca_maskcv:    list(rank, mse)
#   nmf_grid:      list(rank, seed, mse, W, H)
#   nmf_maskcv:    list(rank, mse)
#   cogaps_grid:   list(rank, seed, mse, result = CogapsResult | NULL)
#   cogaps_maskcv: list(rank, alpha, mse)   -- mse NA when the fit failed
#   wgcna_grid:    list(power, net = blockwiseModules output, genes,
#                       samples)   -- `samples` added alongside the
#                       feature/sample-metadata drill-down work; results
#                       ingested before that change lack it, so eigengene
#                       scores are simply NULL for those fits (see below)
#   wto_grid:      list(n, delta, seed, result = wTO.Complete output)
#
# In addition to the feature-loadings matrix (rows = features), pca/nmf/
# cogaps/wgcna results also yield a sample-level "scores" matrix (rows =
# samples, columns = factor_index for pca/nmf/cogaps or module number for
# wgcna) used for the sample/metadata drill-down views. Always oriented
# samples x columns regardless of the underlying method's native shape.

#' jobname -> (method, family) mapping. Generalizes by suffix so new
#' methods following the same naming convention need no change here.
classify_jobname <- function(jobname) {
  method <- sub("_(grid|maskcv)$", "", jobname)
  family <- if (grepl("_maskcv$", jobname)) {
    "maskcv"
  } else if (method %in% c("wgcna", "wto")) {
    "param_grid"
  } else {
    "seed_sweep"
  }
  list(method = method, family = family)
}

#' Standardize a wTO result into an edge data.frame with columns
#' node1, node2, wto, pval, padj. Written defensively since wTO.Complete's
#' return shape varies slightly by version/options.
standardize_wto_edges <- function(result) {
  candidates <- if (is.data.frame(result)) list(result) else Filter(is.data.frame, result)
  edge_df <- NULL
  for (cand in candidates) {
    nms <- tolower(names(cand))
    if (all(c("node.1", "node.2") %in% nms) || all(c("node1", "node2") %in% nms)) {
      edge_df <- cand
      break
    }
  }
  if (is.null(edge_df)) stop("could not locate an edge table (Node.1/Node.2 columns) in wTO result")
  nms <- tolower(names(edge_df))
  pick <- function(...) {
    for (cand in c(...)) {
      i <- which(nms == cand)
      if (length(i) == 1) return(edge_df[[i]])
    }
    rep(NA_real_, nrow(edge_df))
  }
  data.frame(
    node1 = as.character(pick("node.1", "node1")),
    node2 = as.character(pick("node.2", "node2")),
    wto   = as.numeric(pick("wto", "wto_sign", "wto.abs")),
    pval  = as.numeric(pick("pval", "p.value", "p")),
    padj  = as.numeric(pick("pval.adj", "padj", "p.adj", "pval.fdr")),
    stringsAsFactors = FALSE
  )
}

#' Extract one result element into fit metadata + artifacts. `params_row`
#' is the corresponding row of the bundle's params.RDS (authoritative for
#' the swept parameters even if the result is malformed/failed).
extract_result <- function(jobname, result, params_row) {
  cls <- classify_jobname(jobname)
  method <- cls$method
  family <- cls$family

  fit <- list(
    rank = NA_integer_, seed = NA_integer_, alpha = NA_real_,
    power = NA_integer_, n_boot = NA_integer_, delta = NA_real_,
    mse = NA_real_, n_factors = NA_integer_, status = "ok"
  )
  grab <- function(name, from = params_row) {
    if (!is.null(from[[name]])) from[[name]] else NULL
  }
  if (!is.null(grab("rank")))  fit$rank   <- as.integer(params_row$rank)
  if (!is.null(grab("seed")))  fit$seed   <- suppressWarnings(as.integer(params_row$seed))
  if (!is.null(grab("alpha"))) fit$alpha  <- as.numeric(params_row$alpha)
  if (!is.null(grab("power"))) fit$power  <- as.integer(params_row$power)
  if (!is.null(grab("n")))     fit$n_boot <- as.integer(params_row$n)
  if (!is.null(grab("delta"))) fit$delta  <- as.numeric(params_row$delta)

  loadings <- NULL
  modules  <- NULL
  edges    <- NULL
  scores   <- NULL

  if (is.null(result)) {
    fit$status <- "missing"
    return(list(fit = fit, method = method, family = family,
                loadings = NULL, modules = NULL, edges = NULL, scores = NULL))
  }

  if (!is.null(result$mse)) fit$mse <- as.numeric(result$mse)

  if (family == "maskcv") {
    if (is.na(fit$mse)) fit$status <- "failed"
    return(list(fit = fit, method = method, family = family,
                loadings = NULL, modules = NULL, edges = NULL, scores = NULL))
  }

  if (method == "pca") {
    loadings <- result$rotation
    scores   <- result$scores
  } else if (method == "nmf") {
    loadings <- result$W
    if (!is.null(result$H)) scores <- t(result$H)  # patterns x samples -> samples x patterns
  } else if (method == "cogaps") {
    if (is.null(result$result)) {
      fit$status <- "failed"
    } else {
      loadings <- result$result@featureLoadings
      scores   <- result$result@sampleFactors
    }
  } else if (method == "wgcna") {
    if (is.null(result$net) || is.null(result$net$colors)) {
      fit$status <- "failed"
    } else {
      cols <- result$net$colors
      genes <- if (!is.null(names(cols))) names(cols) else result$genes
      modules <- data.frame(gene = as.character(genes),
                            module = as.integer(cols),
                            stringsAsFactors = FALSE)
      fit$n_factors <- length(setdiff(unique(modules$module), 0L))  # module 0 = WGCNA's "unassigned"
      # module eigengenes: only available for results produced after
      # R/methods/wgcna.R started returning `samples` (see comment above)
      if (!is.null(result$net$MEs) && !is.null(result$samples)) {
        me <- as.matrix(result$net$MEs)
        rownames(me) <- result$samples
        scores <- me
      }
    }
  } else if (method == "wto") {
    if (is.null(result$result)) {
      fit$status <- "failed"
    } else {
      edges <- standardize_wto_edges(result$result)
    }
  } else {
    stop("no extractor for method: ", method)
  }

  if (!is.null(loadings)) {
    loadings <- as.matrix(loadings)
    if (is.null(rownames(loadings))) {
      fit$status <- "failed"
      loadings <- NULL
      scores <- NULL
    } else {
      fit$n_factors <- ncol(loadings)
    }
  }

  # scores must be a proper samples x columns matrix with sample-id
  # rownames, or it's not usable for the metadata drill-down views --
  # dropped silently (never fails the fit; loadings/modules are what
  # determine fit status) rather than stored malformed
  if (!is.null(scores)) {
    scores <- as.matrix(scores)
    if (is.null(rownames(scores))) scores <- NULL
  }

  list(fit = fit, method = method, family = family,
       loadings = loadings, modules = modules, edges = edges, scores = scores)
}

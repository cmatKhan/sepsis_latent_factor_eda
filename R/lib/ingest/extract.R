# Per-method extractors: map one job family's result element (one
# `results_<i>.RDS`, unwrapped) + its params row -> a `fits` row and any
# artifacts (loading matrix, module assignments).
#
# Result shapes (verified against real GSE110487 outputs, or against
# synthetic data for the newer methods -- see each R/methods/*.R):
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
#   spca_grid:     list(rank, mse, loadings, scores)
#   spca_maskcv:   list(rank, mse)
#   ica_grid:      list(rank, seed, mse, loadings, scores)
#   ica_maskcv:    list(rank, mse)
#   cp_grid:       list(rank, mse, converged, loadings, scores, time_loadings)
#   tucker_grid:   list(rank_genes, rank_subjects, rank_time, mse, converged,
#                       loadings, scores, time_loadings, core)
#
# params.RDS column shapes (post config-driven-method-registry refactor --
# see R/create_slurm_bundle.R / R/methods/*.R). These are the AUTHORITATIVE
# source for rank/seed/alpha/power (read below), since they're correct
# even for a missing/failed task the result list can't speak for:
#   pca:    rank                              (flat; unchanged)
#   nmf:    k, seed                           (flat; was `rank`, now nnmf()'s real arg name)
#   cogaps: params (list-column: nPatterns, seed, alphaA, ...)  -- see
#           R/methods/cogaps.R's three-call-site sub-block design
#   wgcna:  power                             (flat; unchanged)
#   spca:   K, para                           (flat; `para` -> fit$alpha, see below)
#   ica:    n.comp, alpha, seed               (flat; fastICA's own `alpha`, unrelated to CoGAPS/sPCA's use of the same fits.alpha column)
#   cp:     num_components                    (flat)
#   tucker: rank_genes, rank_subjects, rank_time  (flat, three separate columns)
#
# In addition to the feature-loadings matrix (rows = features), pca/nmf/
# cogaps/wgcna/spca/cp/tucker results also yield a sample-level "scores"
# matrix (rows = samples for pca/nmf/cogaps/spca, module number for wgcna,
# or SUBJECTS -- not samples -- for cp/tucker; see R/methods/cp.R)
# used for the sample/metadata drill-down views. cp/tucker additionally
# yield a third, time-mode matrix (rows = timepoint levels).

#' jobname -> (method, family) mapping. Generalizes by suffix so new
#' methods following the same naming convention need no change here.
#' `PARAM_GRID_METHODS` are deterministic, parameter-only sweeps with no
#' genuine cross-seed stability question (unlike pca/nmf/cogaps's
#' "seed_sweep", which does) -- see each R/methods/<name>.R.
PARAM_GRID_METHODS <- c("wgcna", "spca", "cp", "tucker")

classify_jobname <- function(jobname) {
  method <- sub("_(grid|maskcv)$", "", jobname)
  family <- if (grepl("_maskcv$", jobname)) {
    "maskcv"
  } else if (method %in% PARAM_GRID_METHODS) {
    "param_grid"
  } else {
    "seed_sweep"
  }
  list(method = method, family = family)
}

#' Extract one result element into fit metadata + artifacts. `params_row`
#' is the corresponding row of the bundle's params.RDS (authoritative for
#' the swept parameters even if the result is malformed/failed).
extract_result <- function(jobname, result, params_row) {
  cls <- classify_jobname(jobname)
  method <- cls$method
  family <- cls$family

  fit <- list(
    rank = NA_integer_, seed = NA_integer_, alpha = NA_real_, power = NA_integer_,
    rank_genes = NA_integer_, rank_subjects = NA_integer_, rank_time = NA_integer_,
    mse = NA_real_, n_factors = NA_integer_, status = "ok"
  )

  # method-specific: params.RDS's column shape differs by method (see file
  # header) -- CoGAPS's `params` is a list-column (one named list per row,
  # holding whatever CogapsParams(...) arguments the config set).
  if (method == "pca") {
    if (!is.null(params_row$rank)) fit$rank <- as.integer(params_row$rank)
  } else if (method == "nmf") {
    if (!is.null(params_row$k))    fit$rank <- as.integer(params_row$k)
    if (!is.null(params_row$seed)) fit$seed <- suppressWarnings(as.integer(params_row$seed))
  } else if (method == "cogaps") {
    p <- params_row$params
    if (is.list(p) && is.null(names(p))) p <- p[[1]]   # unwrap the 1-row list-column
    if (!is.null(p$nPatterns)) fit$rank  <- as.integer(p$nPatterns)
    if (!is.null(p$seed))      fit$seed  <- suppressWarnings(as.integer(p$seed))
    if (!is.null(p$alphaA))    fit$alpha <- as.numeric(p$alphaA)
  } else if (method == "wgcna") {
    if (!is.null(params_row$power)) fit$power <- as.integer(params_row$power)
  } else if (method == "spca") {
    if (!is.null(params_row$K))    fit$rank  <- as.integer(params_row$K)
    # sPCA has no real seed -- `alpha` (otherwise CoGAPS-only) doubles as
    # storage for the per-fit `para` sparsity penalty, so the app can
    # label/disambiguate multiple fits at the same rank (see
    # app/R/db_helpers.R::fits_at_rank()).
    if (!is.null(params_row$para)) fit$alpha <- as.numeric(params_row$para)
  } else if (method == "ica") {
    if (!is.null(params_row$n.comp)) fit$rank  <- as.integer(params_row$n.comp)
    if (!is.null(params_row$seed))   fit$seed  <- suppressWarnings(as.integer(params_row$seed))
    if (!is.null(params_row$alpha))  fit$alpha <- as.numeric(params_row$alpha)
  } else if (method == "cp") {
    if (!is.null(params_row$num_components)) fit$rank <- as.integer(params_row$num_components)
  } else if (method == "tucker") {
    if (!is.null(params_row$rank_genes))    fit$rank_genes    <- as.integer(params_row$rank_genes)
    if (!is.null(params_row$rank_subjects)) fit$rank_subjects <- as.integer(params_row$rank_subjects)
    if (!is.null(params_row$rank_time))     fit$rank_time     <- as.integer(params_row$rank_time)
  }

  loadings <- NULL
  modules  <- NULL
  scores   <- NULL
  time_loadings <- NULL

  if (is.null(result)) {
    fit$status <- "missing"
    return(list(fit = fit, method = method, family = family,
                loadings = NULL, modules = NULL, scores = NULL, time_loadings = NULL))
  }

  if (!is.null(result$mse)) fit$mse <- as.numeric(result$mse)

  if (family == "maskcv") {
    if (is.na(fit$mse)) fit$status <- "failed"
    return(list(fit = fit, method = method, family = family,
                loadings = NULL, modules = NULL, scores = NULL, time_loadings = NULL))
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
  } else if (method == "spca") {
    loadings <- result$loadings
    scores   <- result$scores
  } else if (method == "ica") {
    loadings <- result$loadings
    scores   <- result$scores
  } else if (method %in% c("cp", "tucker")) {
    if (isFALSE(result$converged) && is.null(result$loadings)) {
      fit$status <- "failed"
    } else {
      loadings      <- result$loadings
      scores        <- result$scores        # SUBJECT-mode, not sample-mode -- see file header
      time_loadings <- result$time_loadings
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
      time_loadings <- NULL
    } else {
      fit$n_factors <- ncol(loadings)
    }
  }

  # scores/time_loadings must be proper matrices with id rownames, or
  # they're not usable for the drill-down views -- dropped silently
  # (never fails the fit; loadings/modules are what determine fit status)
  # rather than stored malformed
  if (!is.null(scores)) {
    scores <- as.matrix(scores)
    if (is.null(rownames(scores))) scores <- NULL
  }
  if (!is.null(time_loadings)) {
    time_loadings <- as.matrix(time_loadings)
    if (is.null(rownames(time_loadings))) time_loadings <- NULL
  }

  list(fit = fit, method = method, family = family,
       loadings = loadings, modules = modules, scores = scores, time_loadings = time_loadings)
}

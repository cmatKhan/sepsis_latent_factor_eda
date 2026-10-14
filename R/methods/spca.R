# Sparse PCA (elasticnet::spca()), deterministic given its parameters.
#
# Sparsity: `sparse = "varnum"`, under which `para` is the fraction of genes
# allowed a non-zero loading per component, converted to a count from the
# dataset's gene number -- portable across datasets, unlike the raw L1
# penalty of `sparse = "penalty"`, whose effective sparsity varied with gene
# count (and was sensitive to `lambda`).
#
# Choosing the sparsity (Guerra-Urzola, Van Deun, Vera & Sijtsma 2021,
# Psychometrika 86(4):893-919): fix K, take the para maximizing the Index of
# Sparseness IS = PEV_sparse * PEV_pca * PS. Two targets stages:
#   1. spca_build_grid(): a scan at one K (`K_scan`, default 5) over
#      `para_scan` (0.2, 0.3, 0.4, 0.5);
#   2. spca_refine_grid(): the scan's IS winner p*, then every K at
#      p* - para_step, p*, p* + para_step (para_step 0.05), minus what the
#      scan already fit -- about 36 fits per dataset.
# spca_select_fits() then marks the IS maximum within each K over both
# stages. Why one para per dataset: on the earlier full grid (30 datasets,
# K 2-12, para 0.02-0.5; 2,968 fits, 1,451 CPU-hours) every IS winner was at
# para 0.25-0.50, the winner moved by at most 0.10 across K within a dataset
# but differed between datasets, and a scan at K = 5 plus +/- 0.05 for the
# other Ks landed within 0.05 of every K's winner at >= 80% of Ks in all 30
# datasets (mean IS loss 0.004%, worst 0.63%), for about half the compute.
# para <= 0.2 never won. `K_scan`, `para_scan` and `para_step` can be set
# per dataset in its `methods: spca:` block, like `K`. `lambda` is one
# value: under varnum it changed total PEV by at most 5.7e-4 on EARLI.
#
# Convergence: spca() is run with `trace = TRUE`, which prints "iterations
# k" every 10 iterations; the captured trace gives the iterations reached
# (to within 10) and whether `max.iter` was hit without converging. A fit
# still running after `max_hours` stops itself (`killed_not_converged`),
# keeping the iterations it reached. Fits the matrix `mat`.
#
# The fit targets depend on this file's registry and fit function: editing
# their code refits the method (comments don't count).

#' Iterations reached, from spca()'s trace output
#'
#' @param trace Captured output lines of spca(trace = TRUE).
#' @return The last "iterations k" value printed (a multiple of 10), or 0.
spca_trace_iterations <- function(trace) {
  k <- suppressWarnings(as.integer(sub("^\\s*iterations\\s+([0-9]+).*$", "\\1",
                                       grep("iterations", trace, value = TRUE))))
  if (length(k) == 0 || all(is.na(k))) 0L else max(k, na.rm = TRUE)
}

#' Fit sPCA for one (K, para)
#'
#' @param K Number of components.
#' @param para Fraction of genes non-zero per component (`sparse = "varnum"`),
#'   or the L1 penalty (`sparse = "penalty"`).
#' @param type,sparse,use.corr,lambda,max.iter,eps.conv Passed to elasticnet::spca().
#' @param max_hours Wall-time limit; a longer fit is abandoned.
#' @return `list(rank, killed_not_converged, iterations, hit_max_iter, trace,
#'   mse, loadings, scores, center, pev, var.all, spca, elapsed)` -- scores on
#'   the centered matrix, `spca` the rest of spca()'s output, `iterations` the
#'   last value spca()'s trace printed (a lower bound, within 10), and
#'   `hit_max_iter` whether it reached `max.iter` without converging. A fit
#'   stopped at `max_hours` returns `rank`, `killed_not_converged = TRUE`,
#'   `iterations`, `trace` and `elapsed`; an error returns `rank`, `error` and
#'   `elapsed` (capture_fit()).
run_spca_job <- function(K, para, type = "predictor", sparse = "varnum",
                          use.corr = FALSE, lambda = 1e-6, max.iter = 200, eps.conv = 1e-3,
                          max_hours = 8) capture_fit(base = list(rank = K), {
  library(elasticnet)
  # `para` is a fraction of genes when sparse="varnum" (the default here
  # and in every dataset config): converted to a count of non-zero loadings
  # from this dataset's gene number. Used as a raw penalty if
  # `sparse = "penalty"` is requested.
  para_vec <- if (identical(sparse, "varnum")) {
    rep(max(1L, round(para * nrow(mat))), K)
  } else {
    rep(para, K)
  }

  # A fit still running after `max_hours` is abandoned and reported as
  # killed_not_converged = TRUE. At high sparsity and high K spca()'s
  # alternating algorithm can run all max.iter iterations without
  # converging (ANEMONES para 0.05 went from <15 min at K 8 to >8 h at K
  # 10-12, 2026-10-04). setTimeLimit() interrupts at R's next interrupt
  # check, which spca()'s R-level loops reach constantly. The trace is
  # written to a file so the iterations reached survive an interrupt.
  trace_file <- tempfile("spca_trace_")
  on.exit(unlink(trace_file), add = TRUE)
  con <- file(trace_file, open = "wt")
  sink(con)
  setTimeLimit(elapsed = max_hours * 3600, transient = TRUE)
  fit <- tryCatch(
    spca(t(mat), K = K, para = para_vec, type = type, sparse = sparse,
         use.corr = use.corr, lambda = lambda, max.iter = max.iter, trace = TRUE,
         eps.conv = eps.conv),
    error = function(e) {
      if (grepl("elapsed time limit", conditionMessage(e))) NULL else e
    },
    finally = {
      setTimeLimit(elapsed = Inf)
      sink()
      close(con)
    }
  )
  trace <- readLines(trace_file, warn = FALSE)
  iterations <- spca_trace_iterations(trace)
  if (inherits(fit, "error")) stop(fit)
  if (is.null(fit)) {
    return(list(rank = K, killed_not_converged = TRUE, iterations = iterations, trace = trace))
  }

  loadings <- fit$loadings
  rownames(loadings) <- rownames(mat)
  colnames(loadings) <- paste0("Component_", seq_len(ncol(loadings)))
  # spca() centers the data internally (scale(x, center = TRUE, ...)), so
  # scores are the centered samples projected onto the loadings.
  center <- rowMeans(mat)
  scores <- t(mat - center) %*% loadings

  list(rank = K, killed_not_converged = FALSE, iterations = iterations,
       hit_max_iter = iterations >= max.iter, trace = trace,
       # mse = 1 - CUMULATIVE PEV (sum across all K components): elasticnet's
       # `pev` is per-component (each entry is that ONE component's
       # incremental, non-overlapping contribution via the QR-based
       # adjustment for non-orthogonality; see ?elasticnet::spca), so its
       # sum is the variance explained by the whole K-component solution.
       mse = 1 - sum(fit$pev), loadings = loadings, scores = scores, center = center,
       # pev's full per-component curve -- elasticnet's own adjusted-variance
       # bookkeeping, NOT recomputable by naively projecting mat onto
       # loadings (sPCA's components aren't orthogonal, so that double-counts
       # shared variance)
       pev = fit$pev, var.all = fit$var.all,
       # the rest of spca()'s return value (call, type, K, para, lambda,
       # sparse, vn), so ingestion never needs to refit
       spca = unclass(fit)[setdiff(names(fit), c("loadings", "pev", "var.all"))])
})

#' Index of Sparseness for every fit
#'
#' Fits without loadings (stopped at `max_hours`, or errored) are left out;
#' fits that hit `max.iter` stay in, flagged.
#'
#' @param fits List of `list(params, result)`, both stages.
#' @param pca_fit The dataset's PCA result (its full `sdev`).
#' @param p Resolved sPCA parameters, or `NULL` (then `at_scan_edge` is `NA`).
#' @return Data frame (K, para, PEV_sparse, PEV_pca, PS, IS, hit_max_iter,
#'   selected, at_scan_edge), one row per fit; `selected` marks the IS maximum
#'   within each K; `at_scan_edge` (same for every row) is whether the scan's
#'   winner sat at the edge of `para_scan`.
spca_select_fits <- function(fits, pca_fit, p = NULL) {
  pev_pca <- cumsum(pca_fit$sdev^2) / sum(pca_fit$sdev^2)
  fits <- Filter(function(f) !is.null(f$result$loadings), fits)
  if (length(fits) == 0) return(NULL)
  tab <- do.call(rbind, lapply(fits, function(f) {
    L <- as.matrix(f$result$loadings)
    data.frame(K = f$params$K, para = f$params$para,
               PEV_sparse = sum(f$result$pev), PEV_pca = pev_pca[f$params$K],
               PS = mean(L == 0), hit_max_iter = isTRUE(f$result$hit_max_iter))
  }))
  tab$IS <- tab$PEV_sparse * tab$PEV_pca * tab$PS
  tab <- tab[order(tab$K, tab$para), ]
  tab$selected <- stats::ave(tab$IS, tab$K, FUN = function(x) seq_along(x) == which.max(x)) == 1
  tab$at_scan_edge <- NA
  if (!is.null(p)) {
    scan <- tab[tab$K == p$K_scan & round(tab$para, 4) %in% round(p$para_scan, 4), , drop = FALSE]
    if (nrow(scan) > 0) {
      best <- scan$para[which.max(scan$IS)]
      tab$at_scan_edge <- best <= min(p$para_scan) || best >= max(p$para_scan)
    }
  }
  rownames(tab) <- NULL
  tab
}

#' sPCA grid rows for (K, para) pairs
#'
#' Shared by both stages so their rows have the same columns.
#'
#' @param K,para Equal-length vectors of components and sparsity.
#' @param p Resolved sPCA parameters.
#' @return Data frame, one row per (K, para).
spca_grid_rows <- function(K, para, p) {
  data.frame(K = K, para = para, sparse = p$sparse, use.corr = p$use.corr, lambda = p$lambda,
             max.iter = p$max.iter, eps.conv = p$eps.conv, max_hours = p$max_hours,
             stringsAsFactors = FALSE)
}

#' sPCA's first stage: the scan at `K_scan` over `para_scan`
#'
#' @param p Resolved sPCA parameters.
#' @return Data frame, one row per fit.
spca_build_grid <- function(p) {
  if (length(p$lambda) != 1) {
    stop("spca's `lambda` must be a single value, not swept (it has no practical effect under ",
         "sparse = \"varnum\"); got ", length(p$lambda), " values.")
  }
  if (any(p$max.iter %% 10 != 0)) {
    stop("spca's `max.iter` must be a multiple of 10: spca()'s trace prints every 10 iterations, ",
         "and convergence is read from it; got ", paste(p$max.iter, collapse = ", "), ".")
  }
  if (length(p$K_scan) != 1 || !p$K_scan %in% p$K) {
    stop("spca's `K_scan` must be one value among `K` (", paste(p$K, collapse = ", "), "); got ",
         paste(p$K_scan, collapse = ", "), ". Set it in the dataset's methods: spca: block.")
  }
  spca_grid_rows(p$K_scan, sort(unique(p$para_scan)), p)
}

#' sPCA's second stage: every K around the scan's Index-of-Sparseness winner
#'
#' The winner p* among the scan's fits, then each K at p* - para_step, p*
#' and p* + para_step (rounded to 4 digits, kept within [0.02, 1]), minus
#' the (K, para) pairs the scan already fit. A winner at the edge of
#' `para_scan` is refined around all the same, with a warning (the
#' selection table records `at_scan_edge`).
#'
#' @param fits The scan's fits.
#' @param pca_fit The dataset's PCA result.
#' @param p Resolved sPCA parameters.
#' @return Grid rows shaped like spca_build_grid()'s.
spca_refine_grid <- function(fits, pca_fit, p) {
  tab <- spca_select_fits(fits, pca_fit, p)
  if (is.null(tab)) {
    stop("spca: no scan fit finished (K_scan = ", p$K_scan, ", para_scan = ",
         paste(p$para_scan, collapse = ", "), ", max_hours = ", p$max_hours,
         "); every one was stopped or failed, so there is no sparsity to refine around.")
  }
  best <- tab$para[which.max(tab$IS)]
  if (isTRUE(tab$at_scan_edge[1])) {
    warning("spca: the scan's best para (", best, ") is at the edge of para_scan (",
            paste(p$para_scan, collapse = ", "), "); consider widening para_scan for this dataset.")
  }
  para <- round(best + c(-1, 0, 1) * p$para_step, 4)
  para <- unique(pmin(1, pmax(0.02, para)))
  grid <- expand.grid(K = p$K, para = para)
  done <- paste(tab$K, round(tab$para, 4))
  grid <- grid[!paste(grid$K, round(grid$para, 4)) %in% done, , drop = FALSE]
  spca_grid_rows(grid$K, grid$para, p)
}

spca_registry <- list(
  needs_nonneg = FALSE,
  global_object = "mat",
  fn = run_spca_job,
  pkgs = "elasticnet",
  # K_scan, para_scan, para_step: the scan-then-refine sparsity search in
  # this file's header (para is the fraction of genes non-zero per
  # component under sparse = "varnum").
  defaults = list(K = 2:12, K_scan = 5, para_scan = c(0.2, 0.3, 0.4, 0.5), para_step = 0.05,
                   sparse = "varnum", lambda = 1e-6,
                   use.corr = FALSE, max.iter = 200, eps.conv = 1e-3,
                   max_hours = 8),
  build_grid = spca_build_grid,
  refine_grid = spca_refine_grid,
  select_fits = spca_select_fits
)

# Ingest contract (R/targets/ingest.R): kept apart from spca_registry, which
# the fit targets depend on, so editing it never refits. `mse` stays NULL:
# sPCA's natural summary is variance explained (pev_sparse in fit_metrics),
# not a reconstruction error. The IS selection metrics come from
# selection_spca_<dataset> (spca_select_fits()).
spca_ingest <- list(
  family = "param_grid", sign_ambiguous = TRUE, has_loadings = TRUE,
  resource_params = character(0),
  extract = function(result, params) {
    L <- result$loadings
    nonzero <- if (!is.null(L)) colSums(as.matrix(L) != 0) else NULL
    fit_record_from(
      rank = params$K, loadings = L, scores = result$scores,
      diag = c(list(pev = result$pev, var_all = result$var.all, n_nonzero = nonzero,
                    center = result$center, trace = result$trace), result$spca),
      metrics = c(pev_sparse = if (!is.null(result$pev)) sum(result$pev) else NA_real_,
                  iterations = result$iterations %||% NA_real_,
                  converged = if (is.null(result$hit_max_iter)) NA_real_ else as.numeric(!result$hit_max_iter)),
      factor_n_genes = nonzero)
  }
)

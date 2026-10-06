# ICA (fastICA::fastICA()). The fixed-point iteration starts from a random
# unmixing matrix, so ICA is fit across seeds (set.seed() before the call;
# fastICA's `method` is fixed to "C"). `alpha` must be a single value.
#
# Prewhitening: fastICA's own whitening takes the SVD of a genes x genes
# covariance, far too slow at real gene counts, so the data is first reduced
# to n.comp + 5 principal components and ICA runs on those. Loadings in gene
# space are `pc_rotation %*% t(A)`; scores are the sources `S`.
#
# ICASSO's second randomization axis: `bootstrap = TRUE` resamples sample
# columns with replacement before fitting (its scores are dropped, since
# resampled columns repeat sample ids). Kurtosis of each source is the
# non-Gaussianity fastICA maximizes. Fits the matrix `mat`.
#
# The fit targets depend on this file's registry and fit function: editing
# their code refits the method (comments don't count). Docs: Methods.

#' Fit ICA for one (n.comp, seed, bootstrap)
#'
#' @param n.comp Number of components.
#' @param alpha,fun,alg.typ,maxit,tol,... Passed to fastICA::fastICA().
#' @param seed Random seed (set.seed() before fitting and resampling).
#' @param bootstrap Resample sample columns with replacement first.
#' @return `list(rank, seed, bootstrap, mse, loadings, scores, kurtosis, W, K,
#'   prewhiten_sdev, A, S, pc_rotation, pc_center, boot_idx, elapsed)`;
#'   `scores` is `NULL` for bootstrap fits.
run_ica_seed_sweep_job <- function(n.comp, alpha = 1, seed, fun = "logcosh",
                                    alg.typ = "parallel", maxit = 200, tol = 1e-4,
                                    bootstrap = FALSE, ...) capture_fit(
    base = list(rank = n.comp, seed = seed, bootstrap = as.integer(isTRUE(bootstrap)), mse = NA_real_), {
  library(fastICA)
  set.seed(seed)

  mat_use <- mat
  boot_idx <- NULL
  if (isTRUE(bootstrap)) {
    boot_idx <- sample(ncol(mat), ncol(mat), replace = TRUE)
    mat_use <- mat[, boot_idx, drop = FALSE]
  }

  n_pcs <- min(n.comp + 5, ncol(mat_use) - 1, nrow(mat_use) - 1)
  pcfit <- prcomp(t(mat_use), rank. = n_pcs, center = TRUE, scale. = FALSE)
  fit <- fastICA(pcfit$x, n.comp = n.comp, alpha = alpha, fun = fun,
                  alg.typ = alg.typ, method = "C", maxit = maxit, tol = tol, ...)

  loadings <- pcfit$rotation %*% t(fit$A)
  rownames(loadings) <- rownames(mat)
  colnames(loadings) <- paste0("Component_", seq_len(ncol(loadings)))

  # Non-Gaussianity of fastICA's own raw source estimate -- see this
  # file's header. Computed on fit$S directly (before any PCA-rotation
  # mapping), 
  mu <- colMeans(fit$S)
  centered_S <- sweep(fit$S, 2, mu, "-")
  sigma <- sqrt(colMeans(centered_S^2))
  ok <- sigma > 0
  kurtosis <- rep(NA_real_, ncol(fit$S))
  kurtosis[ok] <- colMeans(centered_S[, ok, drop = FALSE]^4) / sigma[ok]^4

  if (isTRUE(bootstrap)) {
    scores <- NULL
  } else {
    scores <- fit$S
    rownames(scores) <- colnames(mat)
    colnames(scores) <- colnames(loadings)
  }

  recon_pc <- fit$S %*% fit$A
  recon <- sweep(recon_pc %*% t(pcfit$rotation), 2, pcfit$center, "+")   # samples x genes
  mse <- mean((t(mat_use) - recon)^2)

  list(rank = n.comp, seed = seed, bootstrap = as.integer(isTRUE(bootstrap)),
       mse = mse, loadings = loadings, scores = scores, kurtosis = kurtosis,
       # W (unmixing matrix): fastICA's own output has no self-reported
       # convergence/quality diagnostic -- W's orthonormality is the cheap
       # real check (a genuinely converged W should be close to
       # orthonormal). K (whitening matrix): the other half of fastICA's
       # own return value, kept for the same reason. pcfit$sdev: the
       # prewhitening PCA step's full spectrum -- the only way to verify
       # how much variance the n_pcs-dimensional reduction (this file's
       # header comment calls it "a standard, well-precedented
       # approximation") actually retained for a given dataset/rank,
       # currently unverifiable after the fact.
       W = fit$W, K = fit$K, prewhiten_sdev = pcfit$sdev,
       # The rest of the model, so ingestion never needs to refit: A (mixing
       # matrix) and S (sources; rows = the fitted samples, i.e. resampled
       # ones for bootstrap fits, which `scores` leaves NULL), the
       # prewhitening PCA's rotation/center (loadings = pc_rotation %*% t(A)),
       # and boot_idx (columns of mat resampled; NULL unless bootstrap).
       A = fit$A, S = fit$S, pc_rotation = pcfit$rotation, pc_center = pcfit$center,
       boot_idx = boot_idx)
})

ica_registry <- list(
  needs_nonneg = FALSE,
  global_object = "mat",
  jobname = "ica_grid",
  fn = run_ica_seed_sweep_job,
  pkgs = "fastICA",
  defaults = list(n.comp = c(5, 7, 10, 15, 20), alpha = 1,
                   seed = c(42, 123, 456, 7, 99, 2024, 8675309, 271828, 31415, 90210),
                   bootstrap = c(FALSE, TRUE),
                   fun = "logcosh", alg.typ = "parallel", maxit = 200, tol = 1e-4),
  build_grid = function(p) {
    if (length(p$alpha) != 1) {
      stop("ica's `alpha` must stay a single static value, not swept -- see this file's header ",
           "for why (breaks the app's same_rank cross-seed comparison). Got length ", length(p$alpha))
    }
    expand.grid(n.comp = p$n.comp, alpha = p$alpha, seed = p$seed, bootstrap = p$bootstrap,
                fun = p$fun, alg.typ = p$alg.typ, maxit = p$maxit, tol = p$tol,
                stringsAsFactors = FALSE)
  }
)

# Ingest contract (R/targets/ingest.R): kept apart from ica_registry, which
# the fit targets depend on, so editing it never refits.
ica_ingest <- list(
  family = "seed_sweep", sign_ambiguous = TRUE, has_loadings = TRUE,
  resource_params = character(0),
  extract = function(result, params) {
    W <- result$W
    # fastICA reports no convergence diagnostic; a converged unmixing
    # matrix is near-orthonormal, so ||W W' - I|| is a cheap check.
    ortho <- if (!is.null(W)) norm(W %*% t(W) - diag(nrow(W)), "F") else NA_real_
    sdev <- result$prewhiten_sdev
    n_pcs <- if (!is.null(result$pc_rotation)) ncol(result$pc_rotation) else NA_integer_
    retained <- if (!is.null(sdev) && !is.na(n_pcs)) sum(sdev[seq_len(n_pcs)]^2) / sum(sdev^2) else NA_real_
    kurt <- result$kurtosis
    fit_record_from(
      rank = params$n.comp, seed = params$seed, bootstrap = isTRUE(as.logical(params$bootstrap)),
      mse = result$mse, loadings = result$loadings, scores = result$scores,
      diag = result[intersect(c("W", "K", "A", "S", "prewhiten_sdev", "pc_rotation", "pc_center",
                                "boot_idx", "kurtosis"), names(result))],
      metrics = c(orthonormality_residual = ortho, prewhiten_var_retained = retained),
      factor_kurtosis = kurt)
  }
)

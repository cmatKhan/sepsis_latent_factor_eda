# Dataset- and fit-level diagnostics, computed as targets (on workers) and
# written by write_dataset_db():
#   matrix_diagnostics_row()    -> matrix_diagnostics
#   wgcna_dataset_diagnostics() -> wgcna_sft, gene_significance
#   wgcna_fit_diagnostics()     -> gene_modules, module_kme_summary,
#                                  module_kme_histograms, the fit's kme artifact

#' Distribution diagnostics of a dataset's matrix
#'
#' Value range and moments, per-gene skewness and kurtosis quantiles, the
#' singular-value spectrum of the sample-centered matrix, and the largest
#' off-diagonal sample-sample correlation (a duplicate check).
#'
#' @param mat Feature x sample matrix.
#' @return One-row data frame: `n_genes`, `n_samples` (for `datasets`) and the
#'   `matrix_diagnostics` columns.
matrix_diagnostics_row <- function(mat) {
  mat <- as.matrix(mat)
  mu <- rowMeans(mat)
  centered <- mat - mu   # per-gene centering (vector recycles down columns)
  sigma <- sqrt(rowMeans(centered^2))
  ok <- sigma > 0
  skew <- rep(NA_real_, nrow(mat)); kurt <- rep(NA_real_, nrow(mat))
  skew[ok] <- rowMeans(centered[ok, , drop = FALSE]^3) / sigma[ok]^3
  kurt[ok] <- rowMeans(centered[ok, , drop = FALSE]^4) / sigma[ok]^4
  sv <- svd(scale(t(mat), center = TRUE, scale = FALSE), nu = 0, nv = 0)$d
  sv <- sv[sv > max(sv) * 1e-8]
  sample_cor <- cor(mat)
  diag(sample_cor) <- NA_real_
  data.frame(
    n_genes = nrow(mat), n_samples = ncol(mat),
    val_min = min(mat), val_max = max(mat), val_mean = mean(mat), val_sd = sd(mat),
    skew_p50 = unname(quantile(skew, 0.5, na.rm = TRUE)),
    skew_p95 = unname(quantile(skew, 0.95, na.rm = TRUE)),
    kurtosis_p50 = unname(quantile(kurt, 0.5, na.rm = TRUE)),
    kurtosis_p95 = unname(quantile(kurt, 0.95, na.rm = TRUE)),
    min_singular_value = min(sv), max_singular_value = max(sv),
    condition_number = max(sv) / min(sv),
    max_sample_cor = max(sample_cor, na.rm = TRUE))
}

#' A dataset's sample metadata (the `sample_meta_<dataset>` target)
#'
#' @param meta The resolved dataset config.
#' @return Data frame with every column a plain vector, or `NULL` if the
#'   config has none.
read_sample_metadata <- function(meta) {
  path <- meta$dataset$sample_metadata_path
  if (is.null(path) || !file.exists(path)) return(NULL)
  df <- as.data.frame(arrow::read_parquet(path))
  df[] <- lapply(df, function(x) x[seq_along(x)])
  df
}

#' Cache sample and feature metadata next to the DB
#'
#' The `metadata_files_<dataset>` file target; the file names are the
#' `dataset_artifacts` kinds.
#'
#' @param meta The resolved dataset config.
#' @param sample_metadata read_sample_metadata() output.
#' @param db_path Path to the DB.
#' @return Paths of the written `sample_metadata.rds` / `feature_metadata.rds`.
write_metadata_files <- function(meta, sample_metadata, db_path) {
  ds <- meta$dataset
  art_dir <- artifacts_dir(db_path, ds$id)
  dir.create(art_dir, recursive = TRUE, showWarnings = FALSE)
  out <- character(0)
  if (!is.null(sample_metadata)) {
    out["sample_metadata"] <- file.path(art_dir, "sample_metadata.rds")
    saveRDS(sample_metadata, out[["sample_metadata"]])
  }
  # (file names are the dataset_artifacts kinds; see write_dataset_db())
  path <- ds$feature_metadata_path
  if (!is.null(path) && file.exists(path)) {
    fm <- as.data.frame(arrow::read_parquet(path))
    fm[] <- lapply(fm, function(x) x[seq_along(x)])
    out["feature_metadata"] <- file.path(art_dir, "feature_metadata.rds")
    saveRDS(fm, out[["feature_metadata"]])
  }
  out
}

#' The samples x genes matrix WGCNA was fit on
#'
#' Applies the same `goodSamplesGenes()` filter as run_wgcna_param_job().
#'
#' @param mat Feature x sample matrix.
#' @return Samples x genes matrix.
wgcna_expr <- function(mat) {
  datExpr <- t(as.matrix(mat))
  gsg <- WGCNA::goodSamplesGenes(datExpr, verbose = 0)
  if (gsg$allOK) datExpr else datExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
}

#' Dataset-level WGCNA diagnostics (the `wgcna_dataset_diag_<dataset>` target)
#'
#' @param mat Feature x sample matrix.
#' @param meta The dataset config (WGCNA's powers and network type).
#' @param sample_metadata read_sample_metadata() output.
#' @return `list(sft, significance)`: `pickSoftThreshold()`'s fit per power
#'   (`wgcna_sft`) and each gene's association with each sample-metadata field
#'   (`gene_significance`).
wgcna_dataset_diagnostics <- function(mat, meta, sample_metadata) {
  wg <- meta$methods$network$wgcna
  datExpr <- wgcna_expr(mat)
  fi <- WGCNA::pickSoftThreshold(datExpr, powerVector = as.integer(wg$power),
                                 networkType = wg$networkType %||% "signed", verbose = 0)$fitIndices
  sft <- data.frame(power = fi$Power, sft_r_sq = fi$SFT.R.sq, slope = fi$slope,
                    truncated_r_sq = fi$truncated.R.sq, mean_k = fi$mean.k.,
                    median_k = fi$median.k., max_k = fi$max.k.)
  id_col <- meta$dataset$sample_id_col %||% "sample_id"
  significance <- NULL
  if (!is.null(sample_metadata) && id_col %in% names(sample_metadata)) {
    gs <- generic_association_scan(datExpr, sample_metadata, id_col = id_col)
    if (nrow(gs)) {
      significance <- data.frame(field = gs$field, gene = gs$component, test = gs$test,
                                 statistic = gs$statistic, p_value = gs$p_value)
    }
  }
  list(sft = sft, significance = significance)
}

N_KME_BINS <- 40L

#' Summary statistics of kME values
#'
#' @param v kME values.
#' @return One-row data frame (n, mean, sd, median, p95, p99) of the finite values.
kme_summary <- function(v) {
  v <- v[is.finite(v)]
  if (length(v) == 0) return(data.frame(n = 0L, mean = NA_real_, sd = NA_real_, median = NA_real_,
                                        p95 = NA_real_, p99 = NA_real_))
  q <- unname(stats::quantile(v, c(0.5, 0.95, 0.99)))
  data.frame(n = length(v), mean = mean(v), sd = if (length(v) > 1) stats::sd(v) else NA_real_,
             median = q[1], p95 = q[2], p99 = q[3])
}

#' Fixed-bin histogram of kME values
#'
#' @param v kME values.
#' @return Data frame (bin, bin_lo, bin_hi, count), `N_KME_BINS` bins over \[-1, 1\].
kme_histogram <- function(v) {
  edges <- seq(-1, 1, length.out = N_KME_BINS + 1)
  bins <- pmin(pmax(findInterval(v[is.finite(v)], edges, rightmost.closed = TRUE), 1L), N_KME_BINS)
  data.frame(bin = seq_len(N_KME_BINS), bin_lo = edges[-length(edges)], bin_hi = edges[-1],
             count = tabulate(bins, N_KME_BINS))
}

#' kME diagnostics for one WGCNA fit (one `wgcna_fit_diag_<dataset>` branch)
#'
#' kME is each gene's correlation with each module eigengene. Per gene: its
#' own-module kME, best competing module and margin, an analytic p-value
#' (Fisher z, n samples) and an empirical p-value against the module's
#' permutation null (the eigengenes' samples shuffled `n_perm` times, seeded
#' from the fit_key). Per module: summaries and histograms for member genes,
#' non-member genes and the permutation null. Saves the genes x modules kME
#' matrix as the fit's `kme` artifact.
#'
#' @param key fit_key.
#' @param ingest The WGCNA ingest_method_fits() result.
#' @param mat Feature x sample matrix.
#' @param db_path Path to the DB.
#' @param n_perm Number of permutations.
#' @return `list(gene_modules, summary, histograms, artifacts)` keyed by
#'   fit_key, or `NULL` if the fit has no eigengenes or modules.
wgcna_fit_diagnostics <- function(key, ingest, mat, db_path, n_perm = 100L) {
  me_file <- fit_artifact(ingest, key, "scores", db_path)
  mods <- ingest$modules[ingest$modules$fit_key == key, c("gene", "module")]
  if (is.na(me_file) || nrow(mods) == 0) return(NULL)
  MEs <- as.matrix(readRDS(me_file))
  ids <- suppressWarnings(as.integer(sub("^ME", "", colnames(MEs))))
  keep <- !is.na(ids) & ids != 0L
  MEs <- MEs[, keep, drop = FALSE]; ids <- ids[keep]
  if (length(ids) == 0) return(NULL)
  datExpr <- t(as.matrix(mat))
  samples <- intersect(rownames(MEs), rownames(datExpr))
  genes <- intersect(mods$gene, colnames(datExpr))
  X <- datExpr[samples, genes, drop = FALSE]
  E <- MEs[samples, , drop = FALSE]
  n <- length(samples)

  kme <- stats::cor(X, E)
  colnames(kme) <- paste0("ME", ids)
  fname <- sprintf("%s_kme.rds", key)
  saveRDS(kme, file.path(artifacts_dir(db_path, ingest$dataset_id), fname))

  module <- mods$module[match(genes, mods$gene)]
  col_of <- match(module, ids)                      # NA for module 0
  kme_own <- ifelse(is.na(col_of), NA_real_, kme[cbind(seq_along(genes), ifelse(is.na(col_of), 1L, col_of))])
  others <- kme
  others[cbind(which(!is.na(col_of)), col_of[!is.na(col_of)])] <- -Inf
  best <- max.col(others, ties.method = "first")
  kme_next <- others[cbind(seq_along(genes), best)]
  kme_next[!is.finite(kme_next)] <- NA_real_
  next_module <- ids[best]
  next_module[is.na(kme_next)] <- NA_integer_
  t_stat <- kme_own * sqrt((n - 2) / pmax(1 - kme_own^2, 1e-12))
  p_analytic <- 2 * stats::pt(-abs(t_stat), df = n - 2)

  # Permutation null: shuffle the samples of all eigengenes together, then
  # correlate every gene with each shuffled eigengene.
  set.seed(strtoi(substr(rlang::hash(key), 1, 7), 16L))
  perms <- lapply(seq_len(n_perm), function(b) stats::cor(X, E[sample(n), , drop = FALSE]))
  p_perm <- rep(NA_real_, length(genes))
  summary <- list(); hist <- list()
  for (m in seq_along(ids)) {
    nul <- sort(unlist(lapply(perms, function(p) p[, m]), use.names = FALSE))
    members <- !is.na(col_of) & col_of == m
    at <- which(members)
    p_perm[at] <- (1 + length(nul) - findInterval(kme_own[at], nul, left.open = TRUE)) / (1 + length(nul))
    groups <- list(member = kme[members, m], nonmember = kme[!members, m], permutation = nul)
    for (g in names(groups)) {
      summary[[length(summary) + 1]] <- cbind(data.frame(fit_key = key, module = ids[m], grp = g),
                                              kme_summary(groups[[g]]))
      hist[[length(hist) + 1]] <- cbind(data.frame(fit_key = key, module = ids[m], grp = g),
                                        kme_histogram(groups[[g]]))
    }
  }
  list(
    gene_modules = data.frame(fit_key = key, gene = genes, module = module, kme_own = kme_own,
                              next_module = next_module, kme_next = kme_next,
                              margin = kme_own - kme_next, p_analytic = p_analytic, p_perm = p_perm),
    summary = do.call(rbind, summary), histograms = do.call(rbind, hist),
    artifacts = data.frame(fit_key = key, kind = "kme", path = artifact_rel(ingest$dataset_id, fname)))
}

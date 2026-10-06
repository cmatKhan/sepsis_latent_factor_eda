# Pairwise stability between the fits of one method within a dataset (all
# seeds and all ranks). Rather than every factor x factor cell, keeps the
# Hungarian-matched pairs plus the context needed to judge them (docs:
# Database, "Pair-stability tables"):
#   fit_pairs   one row per fit pair (WGCNA: + adjusted Rand index)
#   matches     the matched factor pairs (`match` numbers them per call)
#   scores      per match and metric: signed value, runner-up (largest
#               |value| among the unmatched cells in the match's row and
#               column) and margin = |value| - runner-up
#   null        per fit pair and metric: n/mean/sd/median/p95/max of the
#               |unmatched cells|
#   histograms  20-bin counts of |value| over [0, 1], matched vs unmatched,
#               pooled per (level_a <= level_b, metric); level = rank
#               (WGCNA: power)
#   stability   per (fit_key, factor_index): median matched |value| over
#               same-rank pairs, per metric (factors.stability_*)
# Everything is keyed by fit_key; write_dataset_db() resolves keys to ids.
#
# Factor methods compare loadings by cosine, Pearson and Spearman
# (similarity.R) and match on cosine -- on |cosine| for sign-ambiguous
# methods (PCA/sPCA/ICA), whose components can flip sign between fits.
# WGCNA compares modules by gene Jaccard.

N_SIM_BINS <- 20L

#' A histogram accumulator
#'
#' @return An environment collecting matched/unmatched |value| counts per
#'   (level pair, metric); see add_to_histogram().
new_histogram_acc <- function() new.env(parent = emptyenv())

#' Add similarity values to a histogram accumulator
#'
#' @param acc A new_histogram_acc() environment, modified in place.
#' @param level_a,level_b The two fits' levels (rank; power for WGCNA); stored
#'   sorted.
#' @param metric Similarity metric name.
#' @param values Similarity values; binned by |value| into `N_SIM_BINS` bins
#'   over \[0, 1\].
#' @param matched `TRUE` for matched pairs, `FALSE` for unmatched cells.
#' @return `NULL`, invisibly.
add_to_histogram <- function(acc, level_a, level_b, metric, values, matched) {
  if (length(values) == 0) return(invisible(NULL))
  lv <- sort(c(level_a, level_b))
  id <- paste(lv[1], lv[2], metric, sep = "|")
  bins <- pmin(pmax(floor(abs(values) * N_SIM_BINS) + 1L, 1L), N_SIM_BINS)
  cur <- acc[[id]] %||% list(level_a = lv[1], level_b = lv[2], metric = metric,
                             matched = integer(N_SIM_BINS), unmatched = integer(N_SIM_BINS))
  counts <- tabulate(bins, N_SIM_BINS)
  if (matched) cur$matched <- cur$matched + counts else cur$unmatched <- cur$unmatched + counts
  acc[[id]] <- cur
  invisible(NULL)
}

#' Histogram accumulator to `similarity_histograms` rows
#'
#' @param acc A new_histogram_acc() environment.
#' @return Data frame (level_a, level_b, metric, bin, bin_lo, bin_hi,
#'   n_matched, n_unmatched), or `NULL` if empty.
histogram_rows <- function(acc) {
  ids <- ls(acc)
  if (length(ids) == 0) return(NULL)
  edges <- seq(0, 1, length.out = N_SIM_BINS + 1)
  do.call(rbind, lapply(ids, function(id) {
    h <- acc[[id]]
    data.frame(level_a = h$level_a, level_b = h$level_b, metric = h$metric, bin = seq_len(N_SIM_BINS),
               bin_lo = edges[-length(edges)], bin_hi = edges[-1],
               n_matched = h$matched, n_unmatched = h$unmatched)
  }))
}

#' Summary of unmatched similarities
#'
#' @param v Similarity values of the unmatched cells.
#' @return One-row data frame (n, mean, sd, median, p95, max) of |v|.
null_stats <- function(v) {
  v <- abs(v)
  if (length(v) == 0) return(data.frame(n = 0L, mean = NA_real_, sd = NA_real_, median = NA_real_,
                                        p95 = NA_real_, max = NA_real_))
  data.frame(n = length(v), mean = mean(v), sd = if (length(v) > 1) stats::sd(v) else NA_real_,
             median = stats::median(v), p95 = unname(stats::quantile(v, 0.95)), max = max(v))
}

#' Matches and their context for one fit pair
#'
#' @param key_a,key_b The two fits' fit_keys.
#' @param level_a,level_b The two fits' levels (rank, or power for WGCNA).
#' @param sims Named list (by metric) of factors_a x factors_b similarity matrices.
#' @param match_idx Two-column (a, b) matrix of matched indices.
#' @param labels_a,labels_b factor_index values of the rows and columns.
#' @param first_match Number to give the first match (`match` ids run across
#'   the whole method).
#' @param acc Histogram accumulator, updated in place.
#' @return `list(matches, scores, null)`: the matched pairs; per match and
#'   metric the signed value, runner-up (largest |value| among the unmatched
#'   cells in the match's row and column) and margin; per metric the
#'   null_stats() of the unmatched cells.
summarize_pair <- function(key_a, key_b, level_a, level_b, sims, match_idx, labels_a, labels_b,
                           first_match, acc) {
  na <- nrow(sims[[1]]); nb <- ncol(sims[[1]])
  matched <- matrix(FALSE, na, nb)
  matched[match_idx] <- TRUE
  n_match <- nrow(match_idx)
  match_no <- first_match + seq_len(n_match) - 1L
  matches <- data.frame(match = match_no, fit_a = key_a, fit_b = key_b,
                        factor_a = labels_a[match_idx[, 1]], factor_b = labels_b[match_idx[, 2]])
  scores <- list(); nulls <- list()
  for (metric in names(sims)) {
    S <- sims[[metric]]
    A <- abs(S)
    runner <- vapply(seq_len(n_match), function(i) {
      a <- match_idx[i, 1]; b <- match_idx[i, 2]
      others <- c(A[a, -b], A[-a, b])
      if (length(others)) max(others) else NA_real_
    }, numeric(1))
    value <- S[match_idx]
    scores[[metric]] <- data.frame(match = match_no, metric = metric, value = value,
                                   runner_up = runner, margin = abs(value) - runner)
    unmatched <- S[!matched]
    nulls[[metric]] <- cbind(data.frame(fit_a = key_a, fit_b = key_b, metric = metric),
                             null_stats(unmatched))
    add_to_histogram(acc, level_a, level_b, metric, value, TRUE)
    add_to_histogram(acc, level_a, level_b, metric, unmatched, FALSE)
  }
  list(matches = matches, scores = do.call(rbind, scores), null = do.call(rbind, nulls))
}

#' Row-bind one component of several summarize_pair() results
#'
#' @param parts List of summarize_pair() outputs.
#' @param name Component name (`"matches"`, `"scores"` or `"null"`).
#' @return The bound data frame, or `NULL`.
bind_parts <- function(parts, name) {
  rows <- lapply(parts, `[[`, name)
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows)) do.call(rbind, rows) else NULL
}

#' Stability across every pair of a factor method's fits
#'
#' Compares every pair of fits (all ranks, all seeds), Hungarian-matching
#' factors on cosine (|cosine| when `sign_ambiguous`).
#'
#' @param keys,ranks,loadings_files Parallel vectors, one entry per ok fit:
#'   fit_key, rank and absolute path to its genes x factors loadings.
#' @param sign_ambiguous Whether factors are unique only up to sign.
#' @return `list(fit_pairs, matches, scores, null, histograms, stability)`,
#'   keyed by fit_key (`stability` from factor_stability()), or `NULL` for
#'   fewer than 2 fits.
compute_factor_pairs <- function(keys, ranks, loadings_files, sign_ambiguous = FALSE) {
  n <- length(keys)
  if (n < 2) return(NULL)
  prep <- lapply(loadings_files, function(p) prepare_loadings(as.matrix(readRDS(p))))
  acc <- new_histogram_acc()
  parts <- list(); fit_pairs <- list(); next_match <- 1L
  for (i in seq_len(n - 1)) for (j in (i + 1):n) {
    sims <- pair_similarities(prep[[i]], prep[[j]])
    if (is.null(sims)) next
    idx <- if (sign_ambiguous) hungarian_match_abs(sims$cosine) else hungarian_match(sims$cosine)
    p <- summarize_pair(keys[i], keys[j], ranks[i], ranks[j], sims, idx,
                        seq_len(nrow(sims$cosine)), seq_len(ncol(sims$cosine)), next_match, acc)
    next_match <- next_match + nrow(idx)
    parts[[length(parts) + 1]] <- p
    fit_pairs[[length(fit_pairs) + 1]] <- data.frame(
      fit_a = keys[i], fit_b = keys[j], same_rank = as.integer(isTRUE(ranks[i] == ranks[j])),
      n_factors_a = nrow(sims$cosine), n_factors_b = ncol(sims$cosine), ari = NA_real_)
  }
  if (length(parts) == 0) return(NULL)
  out <- list(fit_pairs = do.call(rbind, fit_pairs), matches = bind_parts(parts, "matches"),
              scores = bind_parts(parts, "scores"), null = bind_parts(parts, "null"),
              histograms = histogram_rows(acc))
  out$stability <- factor_stability(out)
  out
}

#' Per-factor stability
#'
#' The median matched |similarity| of each factor across its same-rank fit
#' pairs (`factors.stability_*`).
#'
#' @param pairs A compute_factor_pairs() result.
#' @return Data frame (fit_key, factor_index, cosine, pearson, spearman), or
#'   `NULL` without same-rank pairs.
factor_stability <- function(pairs) {
  same <- pairs$fit_pairs[pairs$fit_pairs$same_rank == 1, c("fit_a", "fit_b")]
  m <- merge(pairs$matches, same, by = c("fit_a", "fit_b"))
  if (nrow(m) == 0) return(NULL)
  s <- merge(m, pairs$scores, by = "match")
  long <- rbind(data.frame(fit_key = s$fit_a, factor_index = s$factor_a, metric = s$metric, v = abs(s$value)),
                data.frame(fit_key = s$fit_b, factor_index = s$factor_b, metric = s$metric, v = abs(s$value)))
  agg <- stats::aggregate(v ~ fit_key + factor_index + metric, data = long, FUN = stats::median)
  wide <- stats::reshape(agg, idvar = c("fit_key", "factor_index"), timevar = "metric", direction = "wide")
  names(wide) <- sub("^v\\.", "", names(wide))
  for (m in c("cosine", "pearson", "spearman")) if (is.null(wide[[m]])) wide[[m]] <- NA_real_
  wide[, c("fit_key", "factor_index", "cosine", "pearson", "spearman")]
}

#' Stability across every pair of WGCNA fits
#'
#' Adjusted Rand index of the two module assignments, and modules
#' Hungarian-matched by gene Jaccard (module 0, unassigned genes, excluded).
#'
#' @param mod_list List (parallel to `keys`) of data frames (gene, module).
#' @param keys The fits' fit_keys.
#' @param levels The fits' soft-thresholding powers (histogram levels).
#' @return Same shape as compute_factor_pairs() with metric `"jaccard"` and
#'   `stability = NULL`, or `NULL` for fewer than 2 fits.
compute_wgcna_pairs <- function(mod_list, keys, levels) {
  n <- length(keys)
  if (n < 2) return(NULL)
  acc <- new_histogram_acc()
  parts <- list(); fit_pairs <- list(); next_match <- 1L
  for (i in seq_len(n - 1)) for (j in (i + 1):n) {
    a <- mod_list[[i]]; b <- mod_list[[j]]
    shared <- intersect(a$gene, b$gene)
    la <- a$module[match(shared, a$gene)]
    lb <- b$module[match(shared, b$gene)]
    mods_a <- setdiff(sort(unique(la)), 0L)
    mods_b <- setdiff(sort(unique(lb)), 0L)
    fit_pairs[[length(fit_pairs) + 1]] <- data.frame(
      fit_a = keys[i], fit_b = keys[j], same_rank = 0L,
      n_factors_a = length(mods_a), n_factors_b = length(mods_b),
      ari = mclust::adjustedRandIndex(la, lb))
    if (length(mods_a) == 0 || length(mods_b) == 0) next
    genes_a <- lapply(mods_a, function(m) shared[la == m])
    genes_b <- lapply(mods_b, function(m) shared[lb == m])
    jac <- outer(seq_along(mods_a), seq_along(mods_b), Vectorize(function(x, y) {
      inter <- length(intersect(genes_a[[x]], genes_b[[y]]))
      uni <- length(genes_a[[x]]) + length(genes_b[[y]]) - inter
      if (uni > 0) inter / uni else 0
    }))
    idx <- hungarian_match(jac)
    parts[[length(parts) + 1]] <- summarize_pair(keys[i], keys[j], levels[i], levels[j],
                                                 list(jaccard = jac), idx, mods_a, mods_b, next_match, acc)
    next_match <- next_match + nrow(idx)
  }
  list(fit_pairs = do.call(rbind, fit_pairs), matches = bind_parts(parts, "matches"),
       scores = bind_parts(parts, "scores"), null = bind_parts(parts, "null"),
       histograms = histogram_rows(acc), stability = NULL)
}

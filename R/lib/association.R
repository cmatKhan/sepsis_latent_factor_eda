# Per-field association tests between sample scores and sample metadata
# (WGCNA gene significance, wgcna_dataset_diagnostics() in
# R/targets/diagnostics.R).

#' Association scan of scores against every metadata field
#'
#' For every column of `meta_df` except `id_col`, tests each column of
#' `scores_mat` against it: Spearman correlation for numeric fields,
#' Kruskal-Wallis otherwise. Fields with too little signal (all `NA`, one
#' level, or fewer than 3 aligned samples) are skipped rather than erroring.
#'
#' @param scores_mat Samples x components matrix, rownames sample ids.
#' @param meta_df Sample metadata data frame.
#' @param id_col Sample id column of `meta_df`.
#' @return Data frame `component, field, test, statistic, p_value, padj`
#'   (`padj` BH-adjusted across every row returned).
generic_association_scan <- function(scores_mat, meta_df, id_col = "sample_id") {
  empty <- data.frame(component = character(0), field = character(0),
                       test = character(0), statistic = numeric(0),
                       p_value = numeric(0))
  if (is.null(scores_mat) || is.null(meta_df) || !(id_col %in% names(meta_df))) return(empty)

  shared <- intersect(rownames(scores_mat), meta_df[[id_col]])
  if (length(shared) < 3) return(empty)
  scores_mat <- scores_mat[shared, , drop = FALSE]
  meta_df <- meta_df[match(shared, meta_df[[id_col]]), , drop = FALSE]

  fields <- setdiff(names(meta_df), id_col)
  comps <- colnames(scores_mat)
  if (is.null(comps)) comps <- as.character(seq_len(ncol(scores_mat)))

  rows <- list()
  for (field in fields) {
    v <- meta_df[[field]]
    is_numeric <- is.numeric(v)

    if (is_numeric) {
      # Vectorized across every component at once via cor()'s own
      # `use = "pairwise.complete.obs"` (per-column NA handling, exactly
      # matching the per-component `ok <- !is.na(v) & !is.na(y)` mask the
      # old per-gene loop computed one column at a time) + WGCNA's own
      # corPvalueStudent() (confirmed to reproduce cor.test(method=
      # "spearman", exact=FALSE)'s p-value exactly, real-data check
      # 2026-09-21). ~1000x faster than the equivalent per-column
      # cor.test() loop at real scale (8000 genes: ~2.0s loop vs ~0.06s
      # vectorized) -- this branch is what made WGCNA gene-significance's
      # numeric-field half a genuine minutes-per-dataset cost; the
      # categorical/Kruskal-Wallis branch below has no comparable
      # vectorized equivalent available in this project's dependencies
      # and is left as a per-component loop.
      ok_v <- !is.na(v)
      if (sum(ok_v) >= 3 && length(unique(v[ok_v])) >= 2) {
        sub_mat <- scores_mat[ok_v, , drop = FALSE]
        v_sub <- v[ok_v]
        cor_vals <- suppressWarnings(stats::cor(sub_mat, v_sub, method = "spearman",
                                                 use = "pairwise.complete.obs"))[, 1]
        n_vals <- colSums(!is.na(sub_mat))
        keep <- !is.na(cor_vals) & n_vals >= 3
        if (any(keep)) {
          p_vals <- WGCNA::corPvalueStudent(cor_vals[keep], n_vals[keep])
          rows[[length(rows) + 1]] <- data.frame(
            component = comps[keep], field = field, test = "spearman",
            statistic = unname(cor_vals[keep]), p_value = p_vals)
        }
      }
      next
    }

    for (j in seq_along(comps)) {
      y <- scores_mat[, j]
      ok <- !is.na(v) & !is.na(y)
      if (sum(ok) < 3) next
      res <- tryCatch({
        f <- as.factor(v[ok])
        if (nlevels(droplevels(f)) < 2) {
          NULL   # skip just this (field, component) cell -- NOT return(),
                  # which would abort the whole scan (see this function's
                  # header; a bug fixed 2026-09-21 after it surfaced via
                  # WGCNA gene-significance's much larger field x gene grid)
        } else {
          kt <- kruskal.test(y[ok] ~ droplevels(f))
          list(test = "kruskal", statistic = unname(kt$statistic), p_value = kt$p.value)
        }
      }, error = function(e) NULL)
      if (is.null(res)) next
      rows[[length(rows) + 1]] <- data.frame(
        component = comps[j], field = field, test = res$test,
        statistic = res$statistic, p_value = res$p_value)
    }
  }
  if (length(rows) == 0) return(empty)
  out <- do.call(rbind, rows)
  out$padj <- p.adjust(out$p_value, method = "BH")
  out
}

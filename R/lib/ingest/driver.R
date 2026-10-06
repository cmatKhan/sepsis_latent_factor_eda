# Pattern drivers: projectR::projectionDriveR() takes ONE factor of a fit and
# two sample groups in the same dataset (e.g. two levels of a metadata
# column) and finds genes whose expression differs between the groups,
# weighted by the factor's loadings -- by confidence intervals ("CI") or
# Welch tests ("PV"). Needs the projectR fork
# (cmatKhan/projectR@feature-add_center_by_loadings).

if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

#' projectionDriveR() for one factor and two sample groups
#'
#' @param loadings Genes x factors matrix.
#' @param factor_index Column of `loadings` to test.
#' @param mat Feature x sample expression matrix.
#' @param ids1,ids2 Sample ids of the two groups.
#' @param mode `"CI"` (confidence intervals) or `"PV"` (Welch tests).
#' @return `list(result, n_considered, n_shared)`, or `NULL` when a group has
#'   fewer than 3 samples, the call fails, or no gene is significant.
pattern_driver <- function(loadings, factor_index, mat, ids1, ids2, mode) {
  pattern_name <- colnames(loadings)[factor_index] %||% paste0("factor_", factor_index)
  colnames(loadings)[factor_index] <- pattern_name   # projectionDriveR looks the pattern up by name
  ids1 <- intersect(ids1, colnames(mat)); ids2 <- intersect(ids2, colnames(mat))
  if (length(ids1) < 3 || length(ids2) < 3) return(NULL)
  result <- tryCatch(
    projectR::projectionDriveR(cellgroup1 = mat[, ids1, drop = FALSE], cellgroup2 = mat[, ids2, drop = FALSE],
                               loadings = loadings, pattern_name = pattern_name, display = FALSE, mode = mode),
    error = function(e) { message("  projectionDriveR failed (", pattern_name, "): ", conditionMessage(e)); NULL })
  if (is.null(result)) return(NULL)
  n_shared <- if (mode == "CI") length(result$sig_genes$significant_shared_genes %||% character(0))
              else length(result$sig_genes$PV_significant_shared_genes %||% character(0))
  if (n_shared == 0) return(NULL)
  list(result = result[setdiff(names(result), "plotted_ci")],   # ggplot objects not worth persisting
       n_considered = if (mode == "CI") nrow(result$mean_ci) else nrow(result$mean_stats),
       n_shared = n_shared)
}

#' Every pattern-driver result for one dataset
#'
#' For each representative fit's first `max_factors` factors, each
#' categorical sample-metadata column with 2-4 levels, each pair of its
#' levels and each mode. Saves each kept result as an artifact.
#'
#' @param ingests Named list (by method) of ingest_method_fits() outputs.
#' @param representatives Named list (by method) of representative fit_keys.
#' @param mat Feature x sample expression matrix.
#' @param sample_metadata Sample-metadata data frame.
#' @param id_col Its sample-id column.
#' @param dataset_id Dataset id.
#' @param db_path Path to the DB (for the artifact directory).
#' @param modes Modes to run.
#' @param max_factors Factors per fit to test.
#' @return `pattern_drivers` rows keyed by (fit_key, factor_index), or `NULL`.
pattern_driver_rows <- function(ingests, representatives, mat, sample_metadata, id_col, dataset_id,
                                db_path, modes = c("CI", "PV"), max_factors = 2) {
  sm <- sample_metadata
  if (is.null(sm) || !(id_col %in% names(sm))) return(NULL)
  cat_cols <- Filter(function(col) {
    col != id_col && (is.character(sm[[col]]) || is.factor(sm[[col]])) &&
      length(unique(na.omit(sm[[col]]))) %in% 2:4
  }, names(sm))
  if (length(cat_cols) == 0) return(NULL)
  mat <- as.matrix(mat)
  art_dir <- artifacts_dir(db_path, dataset_id)
  dir.create(art_dir, recursive = TRUE, showWarnings = FALSE)
  now <- as.character(Sys.time())
  safe <- function(x) gsub("[^A-Za-z0-9.-]+", "_", x)

  rows <- list()
  for (method in names(ingests)) {
    for (key in representatives[[method]]) {
      lf <- fit_artifact(ingests[[method]], key, "loadings", db_path)
      if (is.na(lf)) next
      loadings <- as.matrix(readRDS(lf))
      for (fi in seq_len(min(max_factors, ncol(loadings)))) {
        for (col in cat_cols) {
          for (p in utils::combn(sort(unique(na.omit(sm[[col]]))), 2, simplify = FALSE)) {
            ids1 <- sm[[id_col]][sm[[col]] %in% p[1]]
            ids2 <- sm[[id_col]][sm[[col]] %in% p[2]]
            for (mode in modes) {
              d <- pattern_driver(loadings, fi, mat, ids1, ids2, mode)
              if (is.null(d)) next
              fname <- sprintf("driver_%s_f%d_%s_%s-vs-%s_%s.rds", key, fi, safe(col),
                               safe(p[1]), safe(p[2]), mode)
              saveRDS(d$result, file.path(art_dir, fname))
              rows[[length(rows) + 1]] <- data.frame(
                fit_key = key, factor_index = fi, grouping_col = col,
                group1_level = as.character(p[1]), group2_level = as.character(p[2]), mode = mode,
                n_genes_considered = d$n_considered, n_significant_shared = d$n_shared,
                path = artifact_rel(dataset_id, fname), computed_at = now)
            }
          }
        }
      }
    }
  }
  do.call(rbind, rows)
}

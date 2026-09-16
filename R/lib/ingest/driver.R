# Differential feature identification -- projectR::projectionDriveR()
# (see the projectR vignette's "7 Differential features identification"
# section). A genuinely different analysis from projections.R/projectr_job.R:
# instead of projecting a WHOLE pattern's loadings onto another dataset,
# this takes ONE pattern from an existing fit and asks which genes are
# significantly differentially weighted by that pattern between TWO
# sample groups WITHIN the SAME dataset (e.g. timepoint, disease status,
# responder vs non-responder) -- confidence-interval (default) or
# Welch-test differences in mean expression, weighted by the pattern's
# own loadings.
#
# Requires the `projectR` package (the user's local fork -- see
# R/README.md's container notes) and the dataset's cached matrix
# (cache_dataset_matrix(), login-node only -- see
# R/lib/ingest/ingest_dataset.R) + registered sample metadata
# (dataset_metadata_sources).

if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

#' Run projectionDriveR() for one (fit, factor, grouping column, level
#' pair) combination and store the result. Idempotent -- the UNIQUE
#' constraint on pattern_drivers means calling this twice for the same
#' combination just no-ops the second time (checked explicitly here so we
#' don't rely on a silent constraint-violation swallow).
run_pattern_driver <- function(con, db_path, fit_id, factor_index, dataset_id,
                                grouping_col, group1_level, group2_level, mode = "CI") {
  existing <- DBI::dbGetQuery(con,
    "SELECT driver_id FROM pattern_drivers
     WHERE fit_id = ? AND factor_index = ? AND grouping_col = ? AND group1_level = ? AND group2_level = ? AND mode = ?",
    params = list(fit_id, factor_index, grouping_col, group1_level, group2_level, mode))
  if (nrow(existing) > 0) return(invisible(existing$driver_id[1]))

  f <- DBI::dbGetQuery(con, "SELECT loadings_file FROM fits WHERE fit_id = ?", params = list(fit_id))
  if (nrow(f) == 0 || is.na(f$loadings_file)) return(invisible(NULL))
  loadings <- as.matrix(readRDS(resolve_artifact(f$loadings_file, db_path)))
  if (!(factor_index %in% seq_len(ncol(loadings)))) return(invisible(NULL))
  pattern_name <- colnames(loadings)[factor_index] %||% paste0("factor_", factor_index)
  colnames(loadings)[factor_index] <- pattern_name   # ensure it has a usable name for projectionDriveR's column lookup

  mat_row <- DBI::dbGetQuery(con, "SELECT matrix_file FROM datasets WHERE dataset_id = ?", params = list(dataset_id))
  if (nrow(mat_row) == 0 || is.na(mat_row$matrix_file)) {
    message("  no cached matrix for ", dataset_id, " -- run cache_dataset_matrix() first (login node only)")
    return(invisible(NULL))
  }
  mat <- as.matrix(readRDS(resolve_artifact(mat_row$matrix_file, db_path)))

  meta_row <- DBI::dbGetQuery(con, "SELECT path, id_col FROM dataset_metadata_sources WHERE dataset_id = ? AND kind = 'sample'",
                               params = list(dataset_id))
  if (nrow(meta_row) == 0) return(invisible(NULL))
  sm <- as.data.frame(arrow::read_parquet(meta_row$path[1]))
  if (!(grouping_col %in% names(sm))) return(invisible(NULL))

  ids1 <- sm[[meta_row$id_col[1]]][sm[[grouping_col]] == group1_level]
  ids2 <- sm[[meta_row$id_col[1]]][sm[[grouping_col]] == group2_level]
  ids1 <- intersect(ids1, colnames(mat)); ids2 <- intersect(ids2, colnames(mat))
  if (length(ids1) < 3 || length(ids2) < 3) return(invisible(NULL))   # too few samples per group to be meaningful

  result <- tryCatch(
    projectR::projectionDriveR(
      cellgroup1 = mat[, ids1, drop = FALSE], cellgroup2 = mat[, ids2, drop = FALSE],
      loadings = loadings, pattern_name = pattern_name,
      display = FALSE, mode = mode
    ),
    error = function(e) { message("  projectionDriveR failed (fit ", fit_id, ", ", pattern_name, ", ",
                                   grouping_col, " ", group1_level, " vs ", group2_level, "): ", conditionMessage(e)); NULL }
  )
  if (is.null(result)) return(invisible(NULL))

  n_shared <- if (mode == "CI") length(result$sig_genes$significant_shared_genes %||% character(0))
              else length(result$sig_genes$PV_significant_shared_genes %||% character(0))
  n_considered <- if (mode == "CI") nrow(result$mean_ci) else nrow(result$mean_stats)

  art_dir <- artifacts_dir(db_path, dataset_id)
  dir.create(art_dir, recursive = TRUE, showWarnings = FALSE)
  fname <- sprintf("driver_fit%d_f%d_%s_%s-vs-%s_%s.rds", fit_id, factor_index, grouping_col,
                    make.names(group1_level), make.names(group2_level), mode)
  saveRDS(result[setdiff(names(result), "plotted_ci")], file.path(art_dir, fname))   # ggplot objects not worth persisting

  DBI::dbExecute(con,
    "INSERT INTO pattern_drivers (fit_id, factor_index, grouping_col, group1_level, group2_level, mode,
                                   n_genes_considered, n_significant_shared, result_file, computed_at)
     VALUES (?,?,?,?,?,?,?,?,?,datetime('now'))",
    params = list(fit_id, factor_index, grouping_col, group1_level, group2_level, mode,
                  n_considered, n_shared, file.path("stability_artifacts", dataset_id, fname)))
  invisible(DBI::dbGetQuery(con, "SELECT last_insert_rowid() AS id")$id)
}

#' Batch pass for one dataset: every representative fit of every
#' loadings-bearing method (same scope as fgsea -- pca/nmf/cogaps/spca/ica;
#' cp/tucker's subject/time modes make "sample grouping" ambiguous, wgcna
#' has no continuous pattern to weight by), first 2 factors only (keep it
#' minimal, matching fgsea_job.R's own scope), against every 2-4-level
#' categorical sample-metadata column (skips the id column itself and
#' anything with >4 levels -- avoids combinatorial blowup on high-
#' cardinality columns like patient_id), all pairs of levels.
run_all_pattern_drivers <- function(con, db_path, dataset_id, mode = "CI", max_factors = 2) {
  meta_row <- DBI::dbGetQuery(con, "SELECT path, id_col FROM dataset_metadata_sources WHERE dataset_id = ? AND kind = 'sample'",
                               params = list(dataset_id))
  if (nrow(meta_row) == 0) return(invisible(NULL))
  sm <- as.data.frame(arrow::read_parquet(meta_row$path[1]))

  cat_cols <- Filter(function(col) {
    col != meta_row$id_col[1] && (is.character(sm[[col]]) || is.factor(sm[[col]])) &&
      length(unique(na.omit(sm[[col]]))) %in% 2:4
  }, names(sm))
  if (length(cat_cols) == 0) return(invisible(NULL))

  for (method in c("pca", "nmf", "cogaps", "spca", "ica")) {
    for (fit_id in representative_fit_ids(con, dataset_id, method)) {
      n_factors <- DBI::dbGetQuery(con, "SELECT n_factors FROM fits WHERE fit_id = ?", params = list(fit_id))$n_factors
      if (length(n_factors) == 0 || is.na(n_factors)) next
      for (fi in seq_len(min(max_factors, n_factors))) {
        for (col in cat_cols) {
          levels_ <- sort(unique(na.omit(sm[[col]])))
          pairs <- utils::combn(levels_, 2, simplify = FALSE)
          for (p in pairs) {
            run_pattern_driver(con, db_path, fit_id, fi, dataset_id, col, p[1], p[2], mode = mode)
          }
        }
      }
    }
  }
  invisible(NULL)
}

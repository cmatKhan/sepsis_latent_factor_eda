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
#'
#' @param sample_metadata / id_col: pass a pre-loaded sample-metadata
#'   data.frame + its id column name to use it directly instead of
#'   re-reading `dataset_metadata_sources`. Needed because that table
#'   deliberately stores a LIVE pointer to the config's own
#'   sample_metadata_path (see app/R/metadata_helpers.R's header -- so the
#'   app always reflects the source file's current contents, no re-ingest
#'   needed) -- which is typically an absolute path on whatever machine
#'   holds the raw data, unreachable from inside a cluster container.
#'   run_driver_job() (R/ingest_jobs/driver_job.R) supplies metadata
#'   pre-read on the login node instead of ever touching that raw path
#'   itself from inside the container. Leave both NULL (the default) for
#'   R/ingest_results.R's direct, non-containerized CLI use, where the
#'   raw path DOES resolve fine and the old DB-lookup behavior is kept.
run_pattern_driver <- function(con, db_path, fit_id, factor_index, dataset_id,
                                grouping_col, group1_level, group2_level, mode = "CI",
                                sample_metadata = NULL, id_col = NULL) {
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

  if (is.null(sample_metadata)) {
    meta_row <- DBI::dbGetQuery(con, "SELECT path, id_col FROM dataset_metadata_sources WHERE dataset_id = ? AND kind = 'sample'",
                                 params = list(dataset_id))
    if (nrow(meta_row) == 0) return(invisible(NULL))
    sm <- as.data.frame(arrow::read_parquet(meta_row$path[1]))
    id_col <- meta_row$id_col[1]
  } else {
    sm <- sample_metadata
    if (is.null(id_col)) return(invisible(NULL))
  }
  if (!(grouping_col %in% names(sm))) return(invisible(NULL))

  ids1 <- sm[[id_col]][sm[[grouping_col]] == group1_level]
  ids2 <- sm[[id_col]][sm[[grouping_col]] == group2_level]
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

  # Skip persisting combinations with NO significantly differentially-
  # weighted shared genes ("the length of shared genes are: 0", printed by
  # projectionDriveR() itself) -- these are the vast majority of
  # (fit, factor, grouping column, level pair) combinations
  # run_all_pattern_drivers() tries (most factor/grouping-column pairings
  # just aren't biologically related), and storing them anyway both wastes
  # an artifact file per combo and fills pattern_drivers with rows no
  # downstream consumer (the app's driver browser) can do anything useful
  # with -- n_significant_shared = 0 has no genes to show. NOT gated behind
  # `force`/re-run logic: run_pattern_driver()'s own existing-row check
  # (top of this function) already treats "no row in pattern_drivers" as
  # "not yet computed", so skipping the INSERT here just means a future
  # rerun harmlessly recomputes (and reskips) the same zero-shared-gene
  # combo again, not that it's silently stuck as "done".
  if (n_shared == 0) return(invisible(NULL))

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
#'
#' @param sample_metadata / id_col: see run_pattern_driver()'s matching
#'   doc -- passed straight through to every run_pattern_driver() call
#'   below instead of re-querying dataset_metadata_sources per pair.
run_all_pattern_drivers <- function(con, db_path, dataset_id, mode = "CI", max_factors = 2,
                                     sample_metadata = NULL, id_col = NULL) {
  if (is.null(sample_metadata)) {
    meta_row <- DBI::dbGetQuery(con, "SELECT path, id_col FROM dataset_metadata_sources WHERE dataset_id = ? AND kind = 'sample'",
                                 params = list(dataset_id))
    if (nrow(meta_row) == 0) return(invisible(NULL))
    sm <- as.data.frame(arrow::read_parquet(meta_row$path[1]))
    id_col <- meta_row$id_col[1]
  } else {
    sm <- sample_metadata
    if (is.null(id_col)) return(invisible(NULL))
  }

  cat_cols <- Filter(function(col) {
    col != id_col && (is.character(sm[[col]]) || is.factor(sm[[col]])) &&
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
            run_pattern_driver(con, db_path, fit_id, fi, dataset_id, col, p[1], p[2], mode = mode,
                                sample_metadata = sm, id_col = id_col)
          }
        }
      }
    }
  }
  invisible(NULL)
}

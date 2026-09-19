# Stage-2 SQLite schema + connection/bookkeeping helpers.
#
# The DB path is always a parameter -- nothing here hardcodes a location.
# Opening a nonexistent path bootstraps the schema, so `open_stability_db()`
# on a fresh file is all the "init" there is. Large numeric artifacts
# (loading matrices, wTO edge tables) are NOT stored in SQLite; they live
# under `<db_dir>/stability_artifacts/<dataset_id>/` and are referenced by
# path from `fits.loadings_file`.
#
# Ingest is additive at the granularity of a (dataset_id, jobname) "family"
# (e.g. GSE110487 x nmf_grid): the `ingests` table records what's been
# loaded; a family already present is skipped unless overwrite is
# requested, in which case `delete_family()` removes every row + artifact
# belonging to it before a fresh ingest (plain replacement, no updating).

library(DBI)
library(RSQLite)

open_stability_db <- function(db_path) {
  dir.create(dirname(db_path), recursive = TRUE, showWarnings = FALSE)
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  DBI::dbExecute(con, "PRAGMA foreign_keys = ON")
  ensure_schema(con)
  con
}

artifacts_dir <- function(db_path, dataset_id) {
  file.path(dirname(db_path), "stability_artifacts", dataset_id)
}

#' fits.loadings_file stores paths RELATIVE TO THE DB FILE's directory
#' (artifacts always live next to the DB), so the DB + artifacts move
#' together and resolve correctly no matter where ingest or the app is
#' launched from. This resolves a stored path back to an absolute one.
#'
#' Vectorized over `path` (db_path is a single, shared value) -- most
#' callers already wrap this in vapply()/sapply() for exactly that reason,
#' but a plain vector `path` (e.g. an un-vapply-wrapped `f$loadings_file`
#' with multiple rows) used to hit `is.na(path) || !nzchar(path)`'s
#' length-1-only `||`/`&&`, erroring with "'length = N' in coercion to
#' 'logical(1)'" the moment more than one row came back. Handling vectors
#' directly here fixes every such call site at once, including any future
#' one that forgets to vapply()-wrap it, and is a no-op behavior change for
#' existing scalar/vapply-wrapped callers.
resolve_artifact <- function(path, db_path) {
  out <- rep(NA_character_, length(path))
  ok <- !is.na(path) & nzchar(path)
  is_abs <- ok & startsWith(path, "/")
  out[is_abs] <- path[is_abs]
  rel <- ok & !is_abs
  out[rel] <- file.path(normalizePath(dirname(db_path)), path[rel])
  out
}

ensure_schema <- function(con) {
  statements <- c(
    "CREATE TABLE IF NOT EXISTS datasets (
       dataset_id TEXT PRIMARY KEY,
       description TEXT,
       ingested_at TEXT
     )",
    "CREATE TABLE IF NOT EXISTS ingests (
       ingest_id INTEGER PRIMARY KEY,
       dataset_id TEXT NOT NULL,
       jobname TEXT NOT NULL,
       family TEXT NOT NULL,
       method TEXT NOT NULL,
       ingested_at TEXT,
       n_results INTEGER,
       results_dir TEXT,
       UNIQUE(dataset_id, jobname)
     )",
    "CREATE TABLE IF NOT EXISTS fits (
       fit_id INTEGER PRIMARY KEY,
       dataset_id TEXT NOT NULL,
       method TEXT NOT NULL,
       family TEXT NOT NULL,
       jobname TEXT NOT NULL,
       rank INTEGER, seed INTEGER, alpha REAL,
       power INTEGER, n_boot INTEGER, delta REAL,
       mse REAL, n_factors INTEGER,
       status TEXT NOT NULL DEFAULT 'ok',
       loadings_file TEXT
     )",
    "CREATE INDEX IF NOT EXISTS idx_fits_lookup
       ON fits(dataset_id, method, family, rank, seed)",
    "CREATE TABLE IF NOT EXISTS factors (
       factor_id INTEGER PRIMARY KEY,
       fit_id INTEGER NOT NULL REFERENCES fits(fit_id),
       factor_index INTEGER NOT NULL,
       stability_cosine REAL, stability_pearson REAL, stability_spearman REAL
     )",
    "CREATE INDEX IF NOT EXISTS idx_factors_fit ON factors(fit_id)",
    "CREATE TABLE IF NOT EXISTS factor_pairs (
       fit_a INTEGER NOT NULL REFERENCES fits(fit_id),
       fit_b INTEGER NOT NULL REFERENCES fits(fit_id),
       factor_a INTEGER NOT NULL,
       factor_b INTEGER NOT NULL,
       cosine REAL, pearson REAL, spearman REAL,
       matched INTEGER NOT NULL DEFAULT 0,
       same_rank INTEGER NOT NULL DEFAULT 0
     )",
    "CREATE INDEX IF NOT EXISTS idx_pairs_a ON factor_pairs(fit_a)",
    "CREATE INDEX IF NOT EXISTS idx_pairs_b ON factor_pairs(fit_b)",
    "CREATE TABLE IF NOT EXISTS maskcv_results (
       dataset_id TEXT NOT NULL,
       method TEXT NOT NULL,
       jobname TEXT NOT NULL,
       rank INTEGER, alpha REAL, mse REAL
     )",
    "CREATE TABLE IF NOT EXISTS wgcna_modules (
       fit_id INTEGER NOT NULL REFERENCES fits(fit_id),
       gene TEXT NOT NULL,
       module INTEGER NOT NULL
     )",
    "CREATE INDEX IF NOT EXISTS idx_wgcna_modules_fit ON wgcna_modules(fit_id)",
    "CREATE TABLE IF NOT EXISTS wgcna_fit_pairs (
       fit_a INTEGER NOT NULL REFERENCES fits(fit_id),
       fit_b INTEGER NOT NULL REFERENCES fits(fit_id),
       ari REAL, n_modules_a INTEGER, n_modules_b INTEGER
     )",
    "CREATE TABLE IF NOT EXISTS wgcna_module_pairs (
       fit_a INTEGER NOT NULL, module_a INTEGER NOT NULL,
       fit_b INTEGER NOT NULL, module_b INTEGER NOT NULL,
       jaccard REAL, matched INTEGER NOT NULL DEFAULT 0
     )",
    "CREATE INDEX IF NOT EXISTS idx_wgcna_mp ON wgcna_module_pairs(fit_a, fit_b)",
    "CREATE TABLE IF NOT EXISTS enrichment_cache (
       factor_id INTEGER NOT NULL REFERENCES factors(factor_id),
       query_type TEXT NOT NULL,
       direction TEXT NOT NULL DEFAULT 'pos',
       source TEXT, term_id TEXT, term_name TEXT,
       p_value REAL, intersection_size INTEGER, term_size INTEGER,
       genes TEXT, queried_at TEXT
     )",
    "CREATE INDEX IF NOT EXISTS idx_enrich_factor ON enrichment_cache(factor_id)",
    "CREATE TABLE IF NOT EXISTS enrichment_queried (
       factor_id INTEGER NOT NULL,
       query_type TEXT NOT NULL,
       direction TEXT NOT NULL DEFAULT 'pos',
       queried_at TEXT,
       UNIQUE(factor_id, query_type, direction)
     )",
    "CREATE TABLE IF NOT EXISTS dataset_metadata_sources (
       dataset_id TEXT NOT NULL,
       kind TEXT NOT NULL CHECK (kind IN ('sample','feature')),
       path TEXT NOT NULL,
       id_col TEXT NOT NULL,
       registered_at TEXT,
       UNIQUE(dataset_id, kind)
     )",
    # "Are we discovering too many patterns?" diagnostic (see
    # R/lib/ingest/redundancy.R) -- computed on ONE representative fit per
    # rank (representative_fit_ids()), never across seeds.
    "CREATE TABLE IF NOT EXISTS pattern_markers (
       fit_id INTEGER NOT NULL REFERENCES fits(fit_id),
       factor_index INTEGER NOT NULL,
       gene TEXT NOT NULL,
       score REAL
     )",
    "CREATE INDEX IF NOT EXISTS idx_pattern_markers_fit ON pattern_markers(fit_id)",
    "CREATE TABLE IF NOT EXISTS fit_redundancy (
       fit_id INTEGER PRIMARY KEY REFERENCES fits(fit_id),
       max_offdiag_cosine REAL,
       median_offdiag_cosine REAL,
       n_factors_with_no_markers INTEGER,
       redundancy_file TEXT
     )",
    # projectR results (see R/ingest_jobs/projectr_job.R) -- fits' loadings
    # projected onto OTHER datasets'/timepoints' matrices. projection_type/
    # include_intercept are stored as real columns (not re-derived from the
    # jobname string) even though they're also implied by which of
    # projectr_within_grid/projectr_cross_grid produced the row.
    "CREATE TABLE IF NOT EXISTS projections (
       projection_id INTEGER PRIMARY KEY,
       source_fit_id INTEGER NOT NULL REFERENCES fits(fit_id),
       source_dataset_id TEXT NOT NULL,
       target_dataset_id TEXT NOT NULL,
       method TEXT NOT NULL,
       projection_type TEXT NOT NULL CHECK (projection_type IN ('within_dataset','cross_dataset')),
       include_intercept INTEGER NOT NULL,
       n_genes_matched INTEGER, n_samples INTEGER,
       mean_r_squared REAL, median_r_squared REAL,
       projection_file TEXT
     )",
    "CREATE INDEX IF NOT EXISTS idx_projections_source ON projections(source_fit_id)",
    "CREATE INDEX IF NOT EXISTS idx_projections_type ON projections(projection_type)",
    # Differential feature identification (projectR::projectionDriveR(),
    # vignette section 7 -- see R/lib/ingest/driver.R) -- for a given
    # fit's single pattern, which genes are significantly differentially
    # weighted between two sample groups WITHIN that same dataset (e.g.
    # timepoint, disease status) -- a distinct analysis from projections
    # (which compares across datasets/timepoints via the whole pattern).
    "CREATE TABLE IF NOT EXISTS pattern_drivers (
       driver_id INTEGER PRIMARY KEY,
       fit_id INTEGER NOT NULL REFERENCES fits(fit_id),
       factor_index INTEGER NOT NULL,
       grouping_col TEXT NOT NULL,
       group1_level TEXT NOT NULL,
       group2_level TEXT NOT NULL,
       mode TEXT NOT NULL CHECK (mode IN ('CI','PV')),
       n_genes_considered INTEGER,
       n_significant_shared INTEGER,
       result_file TEXT,
       computed_at TEXT,
       UNIQUE(fit_id, factor_index, grouping_col, group1_level, group2_level, mode)
     )",
    "CREATE INDEX IF NOT EXISTS idx_pattern_drivers_fit ON pattern_drivers(fit_id)"
  )
  for (s in statements) DBI::dbExecute(con, s)
  # additive column migrations for DBs created before this column existed --
  # sample-level scores/eigengenes artifact, same relative-to-DB-dir
  # convention as loadings_file (see resolve_artifact())
  ensure_column(con, "fits", "scores_file", "TEXT")
  # query_size lets the app compute gene-ratio dot plots (intersection_size /
  # query_size) without re-querying g:Profiler
  ensure_column(con, "enrichment_cache", "query_size", "INTEGER")
  # cp/tucker (tensor methods): Tucker's per-mode ranks (CP reuses the
  # existing single `rank` column, like every other rank-1-per-fit
  # method) + the third (time) mode's artifact, same relative-to-DB-dir
  # convention as loadings_file/scores_file
  ensure_column(con, "fits", "rank_genes", "INTEGER")
  ensure_column(con, "fits", "rank_subjects", "INTEGER")
  ensure_column(con, "fits", "rank_time", "INTEGER")
  ensure_column(con, "fits", "time_loadings_file", "TEXT")
  # CoGAPS-only: the FULL raw CogapsResult S4 object, saved as its own
  # artifact at ingest time (from the still-in-scope task result, before
  # extract_result() narrows it down to featureLoadings/sampleFactors) --
  # needed later by R/lib/ingest/redundancy.R's cogaps_pattern_markers()
  # (CoGAPS::patternMarkers() requires the real object, not just loadings).
  ensure_column(con, "fits", "raw_result_file", "TEXT")
  # Cached, dataset-level (not per-family) input matrix artifact -- see
  # R/lib/ingest/ingest_dataset.R::cache_dataset_matrix(). Used by
  # projectr_grid so target-matrix lookups don't re-run preprocessing_script
  # per task.
  ensure_column(con, "datasets", "matrix_file", "TEXT")
  # Cached, dataset-level sample/feature metadata artifacts -- same
  # rationale/convention as matrix_file above (see
  # R/lib/ingest/ingest_dataset.R::cache_dataset_metadata()). Needed so
  # anything reading sample/feature metadata to STAGE cluster jobs (e.g.
  # R/create_ingest_slurm_bundle.R's driver_grid sample_metadata_maps and
  # fgsea_grid/projectr_*_grid's ensembl_maps) can do so
  # from wherever it's actually invoked (typically the cluster login node)
  # without needing raw config paths (sample_metadata_path/
  # feature_metadata_path) that only resolve on whatever machine holds the
  # raw HuggingFace data -- confirmed directly (2026-09-18): with no cache,
  # every dataset's sample_metadata_maps/ensembl_maps entry silently came
  # back NULL when staged from the cluster, since file.exists() on those
  # raw paths is always FALSE there.
  ensure_column(con, "datasets", "sample_metadata_file", "TEXT")
  ensure_column(con, "datasets", "feature_metadata_file", "TEXT")
  # dataset.ensembl_col, for kind='feature' rows -- see
  # register_metadata_source()'s doc above for why this is stored
  # (app-side on-demand enrichment needs it without reading config/*.yml).
  # Datasets ingested before this column existed have NULL here until
  # re-ingested (--stage core / R/ingest_results.R re-registers it fresh
  # every time, no separate backfill step needed).
  ensure_column(con, "dataset_metadata_sources", "ensembl_col", "TEXT")
  # dataset.symbol_col (display-only override) -- see
  # register_metadata_source()'s doc above.
  ensure_column(con, "dataset_metadata_sources", "symbol_col", "TEXT")

  # wTO removed entirely (never had any ingested fits in practice) --
  # drop its table/columns outright rather than leaving dead schema
  # around. Safe/no-op if already dropped or on a DB that never had them.
  drop_table_if_exists(con, "wto_fit_pairs")
  drop_column_if_exists(con, "fits", "n_boot")
  drop_column_if_exists(con, "fits", "delta")

  invisible(con)
}

drop_table_if_exists <- function(con, table) {
  DBI::dbExecute(con, sprintf("DROP TABLE IF EXISTS %s", table))
  invisible(NULL)
}

#' DROP COLUMN requires SQLite >= 3.35 (bundled RSQLite here: 3.53.3) --
#' guarded so it silently no-ops on older SQLite rather than erroring.
drop_column_if_exists <- function(con, table, col) {
  info <- DBI::dbGetQuery(con, sprintf("PRAGMA table_info(%s)", table))
  if (!(col %in% info$name)) return(invisible(NULL))
  ver <- DBI::dbGetQuery(con, "SELECT sqlite_version() AS v")$v
  if (utils::compareVersion(ver, "3.35.0") < 0) {
    warning("SQLite ", ver, " is older than 3.35 -- can't DROP COLUMN ", table, ".", col,
             "; it will remain as a harmless unused column.", call. = FALSE)
    return(invisible(NULL))
  }
  DBI::dbExecute(con, sprintf("ALTER TABLE %s DROP COLUMN %s", table, col))
  invisible(NULL)
}

#' Add a column to an existing table if it isn't already there -- lets DBs
#' created by an earlier schema version upgrade in place, since
#' `CREATE TABLE IF NOT EXISTS` alone never adds columns to a table that
#' already exists.
ensure_column <- function(con, table, col, decl) {
  info <- DBI::dbGetQuery(con, sprintf("PRAGMA table_info(%s)", table))
  if (!(col %in% info$name)) {
    DBI::dbExecute(con, sprintf("ALTER TABLE %s ADD COLUMN %s %s", table, col, decl))
  }
  invisible(NULL)
}

#' Register (or refresh) the sample/feature metadata pointers for a dataset,
#' straight from its config -- pure bookkeeping, no metadata content is
#' copied in. Silently no-ops for a kind whose path/id_col isn't set in the
#' config (e.g. matrix_path-mode datasets with no parquet trio).
#'
#' @param ensembl_col only meaningful for kind = "feature" -- the dataset's
#'   `dataset.ensembl_col` (see R/lib/ingest/symbol_mapping.R), stored here
#'   so the app (which never reads config/*.yml directly -- see
#'   app/R/metadata_helpers.R's header) can build the SAME canonical
#'   Ensembl remap for its own on-demand gprofiler queries as the slurm
#'   pipeline uses, without needing config access. NA for kind = "sample".
#' @param symbol_col only meaningful for kind = "feature" -- the dataset's
#'   `dataset.symbol_col`, REQUIRED to be set explicitly in every config
#'   for the app to show gene symbols (see R/lib/ingest/symbol_mapping.R::
#'   build_symbol_map() and config/dataset_metadata.example.yml).
#'   Display-only (never used computationally -- see that function's
#'   header); deliberately no auto-detection fallback -- an unset or
#'   wrong symbol_col just means the app falls back further down its
#'   display chain (see app/R/metadata_helpers.R::build_display_map()),
#'   never a guess. NA for kind = "sample", or when the config never set
#'   symbol_col.
register_metadata_source <- function(con, dataset_id, kind, path, id_col,
                                      ensembl_col = NA_character_, symbol_col = NA_character_) {
  if (is.null(path) || is.null(id_col) || !nzchar(path) || !nzchar(id_col)) return(invisible(NULL))
  DBI::dbExecute(con,
    "INSERT INTO dataset_metadata_sources (dataset_id, kind, path, id_col, ensembl_col, symbol_col, registered_at)
     VALUES (?, ?, ?, ?, ?, ?, datetime('now'))
     ON CONFLICT(dataset_id, kind) DO UPDATE SET
       path = excluded.path, id_col = excluded.id_col, ensembl_col = excluded.ensembl_col,
       symbol_col = excluded.symbol_col, registered_at = excluded.registered_at",
    params = list(dataset_id, kind, path, id_col, ensembl_col, symbol_col))
  invisible(NULL)
}

ensure_dataset <- function(con, dataset_id, description = NA_character_) {
  DBI::dbExecute(con,
    "INSERT INTO datasets (dataset_id, description, ingested_at)
     VALUES (?, ?, datetime('now'))
     ON CONFLICT(dataset_id) DO UPDATE SET ingested_at = datetime('now')",
    params = list(dataset_id, description))
  invisible(NULL)
}

family_already_ingested <- function(con, dataset_id, jobname) {
  n <- DBI::dbGetQuery(con,
    "SELECT COUNT(*) AS n FROM ingests WHERE dataset_id = ? AND jobname = ?",
    params = list(dataset_id, jobname))$n
  n > 0
}

record_ingest <- function(con, dataset_id, jobname, family, method, n_results, results_dir) {
  DBI::dbExecute(con,
    "INSERT INTO ingests (dataset_id, jobname, family, method, ingested_at, n_results, results_dir)
     VALUES (?, ?, ?, ?, datetime('now'), ?, ?)",
    params = list(dataset_id, jobname, family, method, n_results, results_dir))
  invisible(NULL)
}

#' Remove every row + artifact belonging to one (dataset_id, jobname)
#' family -- the overwrite primitive. Plain replacement, no updating.
delete_family <- function(con, db_path, dataset_id, jobname) {
  fit_ids <- DBI::dbGetQuery(con,
    "SELECT fit_id, loadings_file, scores_file, time_loadings_file, raw_result_file FROM fits
     WHERE dataset_id = ? AND jobname = ?",
    params = list(dataset_id, jobname))

  if (nrow(fit_ids) > 0) {
    for (f in c(fit_ids$loadings_file, fit_ids$scores_file, fit_ids$time_loadings_file, fit_ids$raw_result_file)) {
      fa <- resolve_artifact(f, db_path)
      if (!is.na(fa) && file.exists(fa)) unlink(fa)
    }
    ids_sql <- paste(fit_ids$fit_id, collapse = ",")
    DBI::dbExecute(con, sprintf(
      "DELETE FROM enrichment_cache WHERE factor_id IN
         (SELECT factor_id FROM factors WHERE fit_id IN (%s))", ids_sql))
    DBI::dbExecute(con, sprintf(
      "DELETE FROM enrichment_queried WHERE factor_id IN
         (SELECT factor_id FROM factors WHERE fit_id IN (%s))", ids_sql))
    DBI::dbExecute(con, sprintf(
      "DELETE FROM factor_pairs WHERE fit_a IN (%s) OR fit_b IN (%s)", ids_sql, ids_sql))
    DBI::dbExecute(con, sprintf(
      "DELETE FROM wgcna_module_pairs WHERE fit_a IN (%s) OR fit_b IN (%s)", ids_sql, ids_sql))
    DBI::dbExecute(con, sprintf(
      "DELETE FROM wgcna_fit_pairs WHERE fit_a IN (%s) OR fit_b IN (%s)", ids_sql, ids_sql))
    DBI::dbExecute(con, sprintf("DELETE FROM wgcna_modules WHERE fit_id IN (%s)", ids_sql))

    # redundancy artifacts (see R/lib/ingest/redundancy.R) unlinked + rows
    # dropped -- same "plain replacement" treatment as loadings_file etc.
    red_files <- DBI::dbGetQuery(con, sprintf(
      "SELECT redundancy_file FROM fit_redundancy WHERE fit_id IN (%s)", ids_sql))$redundancy_file
    for (f in red_files) {
      fa <- resolve_artifact(f, db_path)
      if (!is.na(fa) && file.exists(fa)) unlink(fa)
    }
    DBI::dbExecute(con, sprintf("DELETE FROM pattern_markers WHERE fit_id IN (%s)", ids_sql))
    DBI::dbExecute(con, sprintf("DELETE FROM fit_redundancy WHERE fit_id IN (%s)", ids_sql))

    # projections sourced from these fits are stale too -- their source
    # loadings no longer exist
    proj_files <- DBI::dbGetQuery(con, sprintf(
      "SELECT projection_file FROM projections WHERE source_fit_id IN (%s)", ids_sql))$projection_file
    for (f in proj_files) {
      fa <- resolve_artifact(f, db_path)
      if (!is.na(fa) && file.exists(fa)) unlink(fa)
    }
    DBI::dbExecute(con, sprintf("DELETE FROM projections WHERE source_fit_id IN (%s)", ids_sql))

    driver_files <- DBI::dbGetQuery(con, sprintf(
      "SELECT result_file FROM pattern_drivers WHERE fit_id IN (%s)", ids_sql))$result_file
    for (f in driver_files) {
      fa <- resolve_artifact(f, db_path)
      if (!is.na(fa) && file.exists(fa)) unlink(fa)
    }
    DBI::dbExecute(con, sprintf("DELETE FROM pattern_drivers WHERE fit_id IN (%s)", ids_sql))

    DBI::dbExecute(con, sprintf("DELETE FROM factors WHERE fit_id IN (%s)", ids_sql))
    DBI::dbExecute(con, sprintf("DELETE FROM fits WHERE fit_id IN (%s)", ids_sql))
  }

  DBI::dbExecute(con,
    "DELETE FROM maskcv_results WHERE dataset_id = ? AND jobname = ?",
    params = list(dataset_id, jobname))
  DBI::dbExecute(con,
    "DELETE FROM ingests WHERE dataset_id = ? AND jobname = ?",
    params = list(dataset_id, jobname))
  invisible(NULL)
}

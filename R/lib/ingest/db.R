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
resolve_artifact <- function(path, db_path) {
  if (is.na(path) || !nzchar(path)) return(NA_character_)
  if (startsWith(path, "/")) return(path)
  file.path(normalizePath(dirname(db_path)), path)
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
    "CREATE TABLE IF NOT EXISTS wto_fit_pairs (
       fit_a INTEGER NOT NULL REFERENCES fits(fit_id),
       fit_b INTEGER NOT NULL REFERENCES fits(fit_id),
       pearson REAL, spearman REAL,
       jaccard_sig REAL, padj_cutoff REAL, n_edges_common INTEGER
     )",
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
     )"
  )
  for (s in statements) DBI::dbExecute(con, s)
  # additive column migrations for DBs created before this column existed --
  # sample-level scores/eigengenes artifact, same relative-to-DB-dir
  # convention as loadings_file (see resolve_artifact())
  ensure_column(con, "fits", "scores_file", "TEXT")
  # query_size lets the app compute gene-ratio dot plots (intersection_size /
  # query_size) without re-querying g:Profiler
  ensure_column(con, "enrichment_cache", "query_size", "INTEGER")
  invisible(con)
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
register_metadata_source <- function(con, dataset_id, kind, path, id_col) {
  if (is.null(path) || is.null(id_col) || !nzchar(path) || !nzchar(id_col)) return(invisible(NULL))
  DBI::dbExecute(con,
    "INSERT INTO dataset_metadata_sources (dataset_id, kind, path, id_col, registered_at)
     VALUES (?, ?, ?, ?, datetime('now'))
     ON CONFLICT(dataset_id, kind) DO UPDATE SET
       path = excluded.path, id_col = excluded.id_col, registered_at = excluded.registered_at",
    params = list(dataset_id, kind, path, id_col))
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
    "SELECT fit_id, loadings_file FROM fits WHERE dataset_id = ? AND jobname = ?",
    params = list(dataset_id, jobname))

  if (nrow(fit_ids) > 0) {
    for (f in fit_ids$loadings_file) {
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
    DBI::dbExecute(con, sprintf(
      "DELETE FROM wto_fit_pairs WHERE fit_a IN (%s) OR fit_b IN (%s)", ids_sql, ids_sql))
    DBI::dbExecute(con, sprintf("DELETE FROM wgcna_modules WHERE fit_id IN (%s)", ids_sql))
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

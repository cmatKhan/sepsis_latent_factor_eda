# One-off backfill: populate fits.cogaps_diag_file / fits.spca_diag_file
# for fits that were ingested BEFORE those columns existed, using data
# already durably persisted on disk -- no re-fit, no cluster, no re-ingest.
#
#   cogaps: every already-ingested cogaps fit's raw_result_file (the full
#     CogapsResult S4 object, saved at ingest time -- see
#     R/lib/ingest/ingest_dataset.R's header comment on raw_result_file)
#     already contains @loadingStdDev/@factorStdDev. This is now ALSO
#     captured automatically for future ingests by R/lib/ingest/
#     extract.R's cogaps branch -- this script only closes the gap for
#     fits ingested before that change landed.
#   spca: every already-ingested spca fit's loadings_file already lets us
#     compute realized per-component sparsity (nonzero-loadings count).
#     `pev`/`var.all` are NOT recoverable this way (they come from
#     elasticnet::spca()'s own return value, not from the loadings
#     matrix) -- those two fields stay NA in the backfilled diag file
#     until the dataset's spca_grid is re-ingested from its still-present
#     cluster-scratch raw results (Tier 2), which will overwrite this
#     fit's row (and spca_diag_file) with the complete bundle. This
#     backfill is a stopgap for the sparsity part only, safe to run now.
#
# Usage: Rscript R/backfill_diagnostics.R <db_path> [--force]
#   --force: recompute/overwrite even for fits that already have a
#   cogaps_diag_file/spca_diag_file set (default: only fills in NULLs).

library(here)
source(here("R/lib/ingest/db.R"))

args <- commandArgs(trailingOnly = TRUE)
positional <- args[!grepl("^--", args)]
if (length(positional) < 1) stop("Usage: Rscript R/backfill_diagnostics.R <db_path> [--force]")
db_path <- positional[[1]]
force <- "--force" %in% args

con <- open_stability_db(db_path)

backfill_cogaps <- function(con, db_path, force) {
  where <- if (force) "raw_result_file IS NOT NULL" else "raw_result_file IS NOT NULL AND cogaps_diag_file IS NULL"
  rows <- DBI::dbGetQuery(con, sprintf(
    "SELECT fit_id, dataset_id, raw_result_file FROM fits WHERE method = 'cogaps' AND status = 'ok' AND %s", where))
  message("cogaps: backfilling ", nrow(rows), " fits")
  for (i in seq_len(nrow(rows))) {
    fit_id <- rows$fit_id[i]
    raw <- readRDS(resolve_artifact(rows$raw_result_file[i], db_path))
    diag <- list(loading_sd = raw@loadingStdDev, factor_sd = raw@factorStdDev)
    art_dir <- artifacts_dir(db_path, rows$dataset_id[i])
    fname <- sprintf("cogaps_grid_fit%d_diag.rds", fit_id)
    saveRDS(diag, file.path(art_dir, fname))
    rel <- file.path("stability_artifacts", rows$dataset_id[i], fname)
    DBI::dbExecute(con, "UPDATE fits SET cogaps_diag_file = ? WHERE fit_id = ?", params = list(rel, fit_id))
    if (i %% 100 == 0) message("  ", i, "/", nrow(rows))
  }
  invisible(nrow(rows))
}

backfill_spca <- function(con, db_path, force) {
  where <- if (force) "loadings_file IS NOT NULL" else "loadings_file IS NOT NULL AND spca_diag_file IS NULL"
  rows <- DBI::dbGetQuery(con, sprintf(
    "SELECT fit_id, dataset_id, loadings_file FROM fits WHERE method = 'spca' AND status = 'ok' AND %s", where))
  message("spca: backfilling ", nrow(rows), " fits (sparsity only -- pev/var_all need re-ingest, see header)")
  for (i in seq_len(nrow(rows))) {
    fit_id <- rows$fit_id[i]
    loadings <- as.matrix(readRDS(resolve_artifact(rows$loadings_file[i], db_path)))
    diag <- list(pev = NA_real_, var_all = NA_real_, n_nonzero = colSums(loadings != 0))
    art_dir <- artifacts_dir(db_path, rows$dataset_id[i])
    fname <- sprintf("spca_grid_fit%d_diag.rds", fit_id)
    saveRDS(diag, file.path(art_dir, fname))
    rel <- file.path("stability_artifacts", rows$dataset_id[i], fname)
    DBI::dbExecute(con, "UPDATE fits SET spca_diag_file = ? WHERE fit_id = ?", params = list(rel, fit_id))
    if (i %% 200 == 0) message("  ", i, "/", nrow(rows))
  }
  invisible(nrow(rows))
}

DBI::dbExecute(con, "BEGIN")
committed <- FALSE
on.exit(if (!committed) DBI::dbExecute(con, "ROLLBACK"))
n_cogaps <- backfill_cogaps(con, db_path, force)
n_spca   <- backfill_spca(con, db_path, force)
DBI::dbExecute(con, "COMMIT")
committed <- TRUE

message("done: ", n_cogaps, " cogaps + ", n_spca, " spca fits backfilled")
DBI::dbDisconnect(con)

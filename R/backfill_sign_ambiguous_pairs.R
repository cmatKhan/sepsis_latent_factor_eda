# One-off re-derivation: recompute `factor_pairs` (and the `matched` flag
# it stores) for every PCA/sPCA fit using the sign-robust Hungarian match
# introduced 2026-09-29 (see R/lib/ingest/pairs.R's `sign_ambiguous`
# argument) -- PCA and sPCA components are unique only up to sign, and the
# previous signed-cosine Hungarian match either mismatched a sign-flipped
# pair entirely or correctly matched it but recorded a strongly negative
# "stability" for what's actually a rock-solid component. No re-fit
# needed: every fit's `loadings_file` artifact already exists on disk, so
# this only touches `factor_pairs` and the (currently unused by the app,
# but kept consistent) `factors.stability_*` columns.
#
# ICA does NOT need this script -- its fits are being fully regenerated
# (R/methods/ica.R's new bootstrap sweep, see R/README.md) and will pick
# up the same fix automatically on re-ingest.
#
# Usage: Rscript R/backfill_sign_ambiguous_pairs.R <db_path>

library(here)
source(here("R/lib/ingest/db.R"))
source(here("R/lib/ingest/similarity.R"))
source(here("R/lib/ingest/pairs.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript R/backfill_sign_ambiguous_pairs.R <db_path>")
db_path <- args[[1]]

con <- open_stability_db(db_path)

#' Delete and fully recompute one (dataset, method)'s factor_pairs (every
#' fit treated as "new" so every pair is regenerated), then re-derive
#' factors.stability_* from the fresh pairs.
recompute_method <- function(con, db_path, dataset_id, method) {
  fits <- DBI::dbGetQuery(con,
    "SELECT fit_id AS id, rank, loadings_file FROM fits
     WHERE dataset_id = ? AND method = ? AND status = 'ok' AND loadings_file IS NOT NULL",
    params = list(dataset_id, method))
  if (nrow(fits) < 2) return(invisible(0L))
  fits$loadings_file_abs <- vapply(fits$loadings_file, resolve_artifact,
                                    character(1), db_path = db_path)

  # Safe to delete on fit_a alone: compute_factor_pairs() only ever pairs
  # fits WITHIN one dataset+method universe (fit_a is always the
  # lower-fit_id side, per compute_factor_pairs_from_universe()'s i < j
  # loop), so every row here has both sides inside this same fit set.
  DBI::dbExecute(con,
    "DELETE FROM factor_pairs WHERE fit_a IN
       (SELECT fit_id FROM fits WHERE dataset_id = ? AND method = ?)",
    params = list(dataset_id, method))

  rows <- compute_factor_pairs_from_universe(
    fits[, c("id", "rank", "loadings_file_abs")], new_ids = fits$id, sign_ambiguous = TRUE)
  if (!is.null(rows)) {
    for (start in seq(1, nrow(rows), by = 200)) {
      chunk <- rows[start:min(start + 199, nrow(rows)), , drop = FALSE]
      DBI::dbWriteTable(con, "factor_pairs", chunk, append = TRUE)
    }
  }
  update_factor_stability(con, dataset_id, method)
  invisible(if (is.null(rows)) 0L else nrow(rows))
}

datasets <- DBI::dbGetQuery(con,
  "SELECT DISTINCT dataset_id FROM fits WHERE method IN ('pca', 'spca')")$dataset_id

DBI::dbExecute(con, "BEGIN")
committed <- FALSE
on.exit(if (!committed) DBI::dbExecute(con, "ROLLBACK"))

total <- 0L
for (ds in datasets) {
  for (m in c("pca", "spca")) {
    n <- recompute_method(con, db_path, ds, m)
    if (n > 0) message(ds, " / ", m, ": ", n, " factor_pairs rows recomputed (sign-aware)")
    total <- total + n
  }
}

DBI::dbExecute(con, "COMMIT")
committed <- TRUE
message("done: ", total, " total factor_pairs rows recomputed across ", length(datasets), " dataset(s)")
DBI::dbDisconnect(con)

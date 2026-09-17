# Fully purge ONE dataset from the stability DB -- every fit/factor/pair/
# enrichment/redundancy/projection/pattern-driver row belonging to it,
# its `datasets`/`dataset_metadata_sources` rows, and its
# stability_artifacts/<dataset_id>/ directory. Leaves the DB exactly as if
# that dataset had never been ingested -- the next ordinary
# `create_ingest_slurm_bundle.R --stage core` (or R/ingest_results.R) run
# will just re-ingest it fresh, no --overwrite needed.
#
# Built for recovering from a corrupted/partial ingest (e.g. a job that
# died mid-write leaving a dataset half-ingested, or a `database disk image
# is malformed` error where you've decided it's easier to drop one
# dataset's rows and re-ingest than fully rebuild the DB) -- see also the
# sqlite3 `.recover`-based file-level repair playbook in R/README.md.
#
# Usage:
#   Rscript R/purge_dataset.R <dataset_id> [db_path]
# db_path defaults to results/stability.sqlite.
#
# ALWAYS back up the DB file first (cp results/stability.sqlite
# results/stability.sqlite.bak) -- this is a real, irreversible DELETE.

library(here)
source(here("R/lib/ingest/db.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript R/purge_dataset.R <dataset_id> [db_path]")
dataset_id <- args[[1]]
db_path <- if (length(args) >= 2) args[[2]] else "results/stability.sqlite"

con <- open_stability_db(db_path)

jobnames <- DBI::dbGetQuery(con, "SELECT DISTINCT jobname FROM ingests WHERE dataset_id = ?",
                             params = list(dataset_id))$jobname
if (length(jobnames) == 0) {
  message("No `ingests` rows found for dataset_id '", dataset_id, "' -- nothing to delete there ",
          "(still checking `datasets`/dataset_metadata_sources/artifacts below).")
} else {
  message("Purging ", length(jobnames), " job famil", if (length(jobnames) == 1) "y" else "ies",
          " for '", dataset_id, "': ", paste(jobnames, collapse = ", "))
  for (jn in jobnames) delete_family(con, db_path, dataset_id, jn)
}

# projections where this dataset was the TARGET (not the source fit) --
# target_dataset_id is a plain TEXT column, not an FK, so delete_family()
# (fit-scoped) never touches these.
proj_files <- DBI::dbGetQuery(con,
  "SELECT projection_file FROM projections WHERE target_dataset_id = ?", params = list(dataset_id))$projection_file
for (f in proj_files) {
  fa <- resolve_artifact(f, db_path)
  if (!is.na(fa) && file.exists(fa)) unlink(fa)
}
n_target_proj <- DBI::dbExecute(con, "DELETE FROM projections WHERE target_dataset_id = ?", params = list(dataset_id))
if (n_target_proj > 0) message("Deleted ", n_target_proj, " projection row(s) where '", dataset_id, "' was the target dataset")

DBI::dbExecute(con, "DELETE FROM dataset_metadata_sources WHERE dataset_id = ?", params = list(dataset_id))
DBI::dbExecute(con, "DELETE FROM datasets WHERE dataset_id = ?", params = list(dataset_id))

art_dir <- artifacts_dir(db_path, dataset_id)
if (dir.exists(art_dir)) {
  message("Removing ", art_dir)
  unlink(art_dir, recursive = TRUE)
}

DBI::dbDisconnect(con)
message("Done -- '", dataset_id, "' fully purged from ", db_path, ". Re-run the ordinary ingest step to reload it.")

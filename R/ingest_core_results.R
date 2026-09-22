# Phase-1 results -> DB: reads ingest_core's array-job output (one
# compute_ingest_bundle() result per dataset, written by
# R/ingest_jobs/ingest_core_job.R::run_ingest_core_compute_job() to
# rslurm's own `results_*.RDS` files) and writes each bundle into the
# shared stability SQLite DB via R/lib/ingest/ingest_dataset.R::
# write_ingest_bundle(). Same "many parallel compute tasks + one
# lightweight serial merge" shape as R/ingest_enrichment_results.R
# (fgsea_grid/wgcna_ora_grid's merge step) -- this is the ingest_core
# analogue, run AFTER the ingest_core array job completes.
#
# This is the ONLY step that opens a live WRITER connection to the shared
# DB for the core-ingest stage -- one process, one writer, run serially
# over every dataset's already-computed bundle. Since all the expensive
# per-dataset compute (extraction, pair/redundancy math) already happened
# in parallel across the array job's tasks, this step should be
# dramatically faster than the old single non-array `ingest_core` job's
# full serial run.
#
# Usage:
#   Rscript R/ingest_core_results.R --bundle-dir slurm_bundles/ingest --db results/stability.sqlite

library(here); library(optparse); library(DBI)
source(here("R/lib/ingest/db.R"))
source(here("R/lib/ingest/similarity.R"))
source(here("R/lib/ingest/extract.R"))
source(here("R/lib/ingest/pairs.R"))
source(here("R/lib/ingest/redundancy.R"))
source(here("R/lib/ingest/ingest_dataset.R"))

option_list <- list(
  make_option("--bundle-dir", type = "character", help = "e.g. slurm_bundles/ingest"),
  make_option("--db", type = "character", default = "results/stability.sqlite")
)
opt <- parse_args(OptionParser(option_list = option_list))
con <- open_stability_db(opt$db)

fam_dir <- file.path(opt$`bundle-dir`, "_rslurm_ingest_core")
if (!dir.exists(fam_dir)) stop("ingest_core: no results dir found at ", fam_dir)

result_files <- list.files(fam_dir, pattern = "^results_\\d+\\.RDS$", full.names = TRUE)
message("ingest_core: ", length(result_files), " result file(s)")

n_datasets <- 0L
for (f in result_files) {
  x <- readRDS(f)
  # Same either-shape handling as R/ingest_enrichment_results.R --
  # `submit_job_family()`'s max_array_size batches multiple datasets'
  # bundles into one results_<task>.RDS as a plain list when
  # `--core-max-array-size` > 1; `x$dataset_id` is only present on a
  # SINGLE bundle (one dataset per task), never on the outer per-task list.
  bundles <- if (!is.null(x$dataset_id)) list(x) else x

  for (bundle in bundles) {
    if (is.null(bundle) || is.null(bundle$dataset_id)) next   # a failed task (see run_ingest_core_compute_job()'s tryCatch)
    n_datasets <- n_datasets + 1L
    dataset_id <- bundle$dataset_id
    message("dataset: ", dataset_id)
    tryCatch(
      write_ingest_bundle(con, opt$db, bundle),
      error = function(e) message("  FAILED writing bundle for ", dataset_id, ": ", conditionMessage(e))
    )
  }
}
message("\ningest_core: merged ", n_datasets, " dataset bundle(s) across ", length(result_files), " result file(s)")

DBI::dbDisconnect(con)

# Phase-2 slurm job: a SINGLE (non-array) job that loops
# run_all_pattern_drivers() (R/lib/ingest/driver.R, projectR::
# projectionDriveR()-based differential-features pass) over every dataset
# sequentially against the SAME sqlite DB -- one writer, same rationale as
# run_ingest_core_job() (R/ingest_jobs/ingest_core_job.R).
#
# Split out of ingest_core specifically because projectR needs its own
# image (see config/ingest_slurm_config.yml's `driver:` entry, same
# container as projectr_within_grid/projectr_cross_grid), distinct from
# ingest_core's bioconductor-cogaps/fgsea/deseq2 image, which does NOT have
# projectR installed.
#
# Belongs in --stage enrichment, not --stage core: run_all_pattern_drivers()
# calls representative_fit_ids(), which needs fits.mse -- populated only
# after ingest_core has already run (same precondition fgsea_grid/
# gprofiler_grid/projectr_*_grid already have). Also needs each dataset's
# matrix already cached (cache_dataset_matrix(), login-node-only --
# already satisfied by the time --stage enrichment is reached).
#
# Deliberately does NOT `library(here)`/`source()` anything by path -- same
# constraint as every other ingest job (see ingest_core_job.R's header):
# everything it calls must already exist via `global_objects`.
run_driver_job <- function(dataset_ids, db_path, mode = "CI", max_factors = 2) {
  con <- open_stability_db(db_path)
  for (dataset_id in dataset_ids) {
    tryCatch(
      run_all_pattern_drivers(con, db_path, dataset_id, mode = mode, max_factors = max_factors),
      error = function(e) message("  pattern-driver pass failed for ", dataset_id, ": ", conditionMessage(e))
    )
  }
  DBI::dbDisconnect(con)
  invisible(NULL)
}

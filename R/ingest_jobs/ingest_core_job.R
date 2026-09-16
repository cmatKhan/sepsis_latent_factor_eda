# Phase-1 slurm job: a SINGLE (non-array) job that loops every dataset in
# `targets` sequentially against the SAME sqlite DB -- one writer, avoiding
# concurrent-write corruption.
#
# Deliberately does NOT `library(here)`/`source()` anything by path -- the
# container only has slurm_bundles/ingest/ bind-mounted (see
# build_apptainer_rscript_path()), not the project's R/ tree, so
# ingest_one_dataset() and everything it calls must already exist in this
# job's execution environment via `global_objects` (see
# R/create_ingest_slurm_bundle.R, which sources R/lib/ingest/*.R on the
# LOGIN NODE and bundles the resulting function objects alongside
# `targets`/`db_path`/`recompute_redundancy`).
#
# `targets`/`db_path`/`recompute_redundancy` are read as FREE VARIABLES,
# NOT function parameters -- this job is staged via submit_job_family()
# with jobs_df = NULL, which dispatches to rslurm::slurm_call() with no
# `params`. rslurm's generated slurm_run_single_R.txt then calls
# `do.call(f, list())` -- i.e. with ZERO arguments -- relying entirely on
# `add_objects.RData` (loaded into the same global environment this
# function's closure resolves free variables against) to supply them.
# Declaring them as formal parameters instead (as an earlier version of
# this function did) shadows that global lookup and fails with `argument
# "db_path" is missing, with no default` the moment the body references
# it, since nothing ever supplies them positionally. (Contrast with the
# ARRAY jobs in R/ingest_jobs/fgsea_job.R etc., which as jobs_df-driven
# slurm_apply() calls DO get one row's values passed in as real params --
# formal parameters are correct there.)
run_ingest_core_job <- function() {
  con <- open_stability_db(db_path)
  for (i in seq_len(nrow(targets))) {
    force_redundancy <- isTRUE(recompute_redundancy) ||
      (is.character(recompute_redundancy) && targets$dataset_id[i] %in% recompute_redundancy)
    tryCatch(
      # run_pattern_drivers = FALSE -- this container doesn't have projectR
      # installed (a different image -- see R/ingest_jobs/driver_job.R's
      # header and config/ingest_slurm_config.yml's `driver:` entry); that
      # pass is staged as its own driver_grid job family during
      # --stage enrichment instead.
      ingest_one_dataset(con, targets$config_path[i], targets$results_dir[i], db_path,
                          recompute_redundancy = force_redundancy, run_pattern_drivers = FALSE),
      error = function(e) message("  FAILED (", targets$dataset_id[i], "): ", conditionMessage(e))
    )
  }
  DBI::dbDisconnect(con)
  invisible(NULL)
}

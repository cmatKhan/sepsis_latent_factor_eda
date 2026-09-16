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
# `targets`/`db_path`) -- same idiom every R/methods/*.R job function uses
# for `mat`/`mat_nn`/`tnsr`, just for functions instead of data here.
run_ingest_core_job <- function(targets, db_path, recompute_redundancy = FALSE) {
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

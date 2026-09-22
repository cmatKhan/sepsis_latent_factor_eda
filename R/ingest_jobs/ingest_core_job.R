# Phase-1 slurm job: a real ARRAY job (one task per dataset, or a small
# batch of datasets per task via --core-max-array-size -- see
# R/create_ingest_slurm_bundle.R's `--stage core`), each task computing
# ONE dataset's full ingest bundle (R/lib/ingest/ingest_dataset.R::
# compute_ingest_bundle()) with NO shared DB writer at all -- safe to run
# fully in parallel, since compute_ingest_bundle() only ever opens the DB
# for READS (never issues BEGIN/INSERT/UPDATE/DELETE; SQLite's WAL mode
# allows unlimited concurrent readers -- see that function's own header
# for the full audit this design is based on) and every artifact path it
# writes is scoped to that one dataset's own stability_artifacts/<id>/
# directory (no cross-task collisions).
#
# Replaces the old single non-array `run_ingest_core_job()`, which looped
# every dataset serially against ONE shared writer connection -- a real
# seff report showed that approach at only 28.42% CPU efficiency over a
# 1h19m single-core run, and even after fixing the transaction-batching
# portion of that (see ingest_dataset.R's write_ingest_bundle() comment),
# the remaining ~22-25 CPU-minutes summed across ~32 datasets was still a
# genuinely serial, embarrassingly-parallel cost. Each array task here
# just writes its bundle to rslurm's own `results_<task>.RDS` (standard
# `slurm_apply()` output convention -- the actual DB write happens
# afterward, serially, in R/ingest_core_results.R).
#
# Deliberately does NOT `library(here)`/`source()` anything by path -- the
# container only has slurm_bundles/ingest/ bind-mounted, not the
# project's R/ tree, so compute_ingest_bundle() and everything it calls
# must already exist in this job's execution environment via
# `global_objects` (see R/create_ingest_slurm_bundle.R, which sources
# R/lib/ingest/*.R wherever it itself runs -- a plain `Rscript`
# invocation submitted via `srun`, not inside any container -- and bundles
# the resulting function objects alongside `db_path`/`recompute_redundancy`/
# `PROJECT_ROOT`).
#
# `dataset_id`/`config_path`/`results_dir` are real formal parameters
# here (unlike the old single-job version) -- this IS a jobs_df-driven
# `slurm_apply()` array job now, so rslurm supplies one row's values per
# call. `db_path`/`recompute_redundancy`/`PROJECT_ROOT` remain FREE
# VARIABLES (global_objects) -- shared across the whole job family, not
# per-row.
#
# `overwrite` is NOT supported here (always FALSE, same as the old
# run_ingest_core_job() -- it never threaded overwrite through either):
# re-running a family that needs overwriting is a single-dataset,
# deliberate operation better done via the direct CLI
# (R/ingest_results.R), not the bulk array path.
run_ingest_core_compute_job <- function(dataset_id, config_path, results_dir) {
  force_redundancy <- isTRUE(recompute_redundancy) ||
    (is.character(recompute_redundancy) && dataset_id %in% recompute_redundancy)
  tryCatch(
    compute_ingest_bundle(config_path, results_dir, db_path,
                           recompute_redundancy = force_redundancy, project_root = PROJECT_ROOT),
    error = function(e) {
      message("  FAILED (", dataset_id, "): ", conditionMessage(e))
      NULL
    }
  )
}

# TODOs

Deferred changes, grouped by what they would make targets rerun. Costs are for the current
30-dataset build. "Value unchanged" means the rerun target produces the same value, so targets
stops there and nothing downstream reruns.

## Deviations from the docs reorganization (2026-10-06)

These were left as they are because fixing them now would rerun pipeline targets.

1. **Error messages pointing to `R/README.md`.** `R/README.md` is kept as a pointer page,
   rather than deleted, because these code strings still name it:

   | Where | Inside | Fixing it reruns |
   |---|---|---|
   | `R/methods/spca.R` (`spca_build_grid()`) | `spca_registry` | **every sPCA fit (about a day of compute)**, since the registry is a fit-target dependency |
   | `R/lib/matrices.R` (`run_preprocessing_script()`) | the `mat_<dataset>` command | all 30 matrices (preprocessing jobs, minutes each). Downstream stays cached only if each rebuilt matrix is bit-identical -- likely (VST and filtering are deterministic) but unverified; if not, everything refits |
   | `R/lib/metadata.R` (`read_dataset_metadata()`, `validate_dataset_metadata()`) | the `meta_<dataset>` command | the 30 `meta_` targets only (seconds; value unchanged) |
   | `R/lib/method_registry.R` (`discover_method_registries()`) | nothing (runs when `_targets.R` is read) | nothing |

   When fixed, point the messages at `docs/methods.qmd` / `docs/data.qmd` and delete
   `R/README.md`.
2. **`R/db/schema.sql` and `R/db/views.sql` header comments** still say the schema is documented
   in `R/db/SCHEMA.md` (now `docs/database.qmd`). Both files are the `db_schema` file target, so
   any edit -- even a comment -- reruns all 30 `db_<dataset>` writes and `db_cross_dataset`
   (about 1 h 15 m on the main process). Fix with the next real schema change, which also bumps
   `SCHEMA_VERSION` and rebuilds the DB anyway.
3. **App comments.** `app/app.R` and `app/R/db_helpers.R` still reference removed files
   (`R/lib/ingest/db.R`, `extract.R`, `create_ingest_slurm_bundle.R`). No targets impact; fix
   as part of updating the app to the new schema (below).

## Changes that would refit methods

Fit-target dependencies: each `<name>_registry` and fit function, `run_method_fits()` and the
grid helpers, `capture_fit()`, and a controller's `batch_size`/`cpus_per_task`.

- **Remove unused registry fields.** `jobname` (no longer used by ingestion) and the
  `requires_subject_timepoint` check (left from CP/Tucker; no method sets it). Removing
  `jobname` from the registries refits **every method**; the check itself lives in
  `metadata.R` (cheap, see above).
- **sPCA convergence reporting.** elasticnet's `spca()` doesn't return its iteration count, so
  non-convergence is only inferred from the 8-hour stop. Reporting it would need our own
  wrapper of the spca algorithm: refits **every sPCA fit**.
- **Fix the remaining doc-pointer strings in fit-dependency code** (deviation 1 above).
- **Any further batching or thread changes** (`batch_size`, `cpus_per_task` per controller):
  refit that method.

## Changes that would rerun ingestion or the DB (no refits)

Ingestion targets depend on each `<name>_ingest`, `R/targets/*.R` and `R/lib/ingest/*.R`.
`fgsea_*` and `projectr_*` (the expensive ingestion steps, hours) rerun only if their inputs'
values change.

- **Unused argument.** `fit_failed(power = )` in `R/targets/ingest.R` is never used. Removing it
  reruns every `ingest_*` target (minutes); their values are unchanged, so nothing further
  reruns.
- **DB size** (3.6 GB). Options: move enrichment gene lists (~1.1 GB as JSON) to per-dataset
  files, or store only the top N marker genes per factor (`pattern_markers`, 4.6 M rows).
  Either changes the schema: bump `SCHEMA_VERSION`, delete the DB, rerun the DB writes (and
  `redundancy_*` for the markers change).
- **Schema comments** (deviation 2 above).

## Changes that add targets

- **Missing datasets:** `gains` (preprocessing script is a stub; pooling its four accessions
  needs a decision), `ROSE`/`ROSE_DAY0`/`ROSE_DAY2` (data not on the cluster), `MARS` and
  `SHIP-TREND` (preprocessing scripts missing). Adding a dataset builds its own targets, but it
  also changes `matrix_files` and every dataset's projectR job table (a new target dataset), so
  **all projections rerun** (about 3 h) as well as `db_cross_dataset`.

## No targets impact

- **Update the app** to the new DB schema (`docs/database.qmd`), including the stale comments
  (deviation 3).
- **Delete `dev/phase1_parity.R` and `dev/phase2_parity.R`**: they compare against the old
  schema and are obsolete now that the rewrite is verified.
- **`config/pipeline.yml` header** still mentions the deleted rslurm driver and per-dataset
  `slurm:` blocks (the file is read when `_targets.R` loads, not tracked as a target).
- **Commit** the targets rewrite, DB redesign and docs (suggested as separate commits).

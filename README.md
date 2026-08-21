# Sepsis timecourse latent-factor stability

Pipeline for quantifying the stability of latent factors (PCA, NMF, CoGAPS,
WGCNA, wTO) found in the sepsis timecourse expression data, across random
seeds and across rank/parameter choices.

Three stages:

1. **rslurm grid setup + submission** (cluster) -- see [R/README.md](R/README.md)
   for preparing an input matrix, writing a `dataset_metadata.yml`, and
   staging/submitting the batch jobs. Produces `results_<i>.RDS` files per
   job family under `slurm_bundles/<dataset_id>/`, then copied back locally
   into a results directory (e.g. `GSE110487_results/`).
2. **Ingest** -- loads a results directory into a SQLite DB, computing all
   pairwise stability metrics (see below).
3. **Shiny app** -- browses the DB with a drill-down explorer (dataset >
   method > rank/parameter > factor).

Ingest and the app are fully decoupled: the app only ever reads whatever the
DB currently contains, and both take the DB path as a parameter.

## Ingest

```r
Rscript R/ingest_results.R <dataset_config.yml> <results_dir> <db_path> [--overwrite [jobname,...]]
```

Example (the dataset currently in this repo):

```r
Rscript R/ingest_results.R config/GSE110487_config.yml GSE110487_results results/stability.sqlite
```

- Discovers job families as the subdirectories of `<results_dir>`
  containing `results_*.RDS` files.
- **Additive**: re-running the same command only adds families not yet in
  the DB for that dataset (e.g. drop a `wto_grid/` directory in later and
  re-run -- only `wto_grid` gets ingested; everything else is reported as
  skipped).
- **Overwrite**: `--overwrite` (bare) replaces every family found in
  `<results_dir>`; `--overwrite jobname1,jobname2` replaces just those
  families. Overwrite always deletes and re-writes a family's rows and
  artifacts, never updates in place.
- A new dataset (different `dataset.id` in the config) is ingested
  alongside any existing datasets in the same DB file -- point multiple
  configs at the same `<db_path>` to build up a multi-dataset DB.
- Large artifacts (loading matrices, sample/module scores, wTO edge tables)
  are written under `<db_dir>/stability_artifacts/<dataset_id>/` and
  referenced by path from the DB; everything the app filters/aggregates on
  lives in the SQLite tables themselves.
- The dataset config's `sample_metadata_path`/`sample_id_col` and
  `feature_metadata_path`/`feature_id_col` are registered in the DB as
  **pointers** (path + id column), not copied in -- refreshed on every
  ingest run regardless of which job families changed. The app always
  re-reads the pointed-to file live, so editing the metadata file itself
  (e.g. adding a clinical variable) needs no re-ingest. Get `feature_id_col`
  right: it must match the column in `feature_metadata_path` whose values
  equal the fitted matrix's row names (verify with a quick join-coverage
  check by hand -- ingest doesn't validate this for you).

Run interactively instead of via `Rscript` by setting `ingest_config_path`,
`ingest_results_dir`, `ingest_db_path` (and optionally `ingest_overwrite`)
before sourcing `R/ingest_results.R`.

## App

```r
Sys.setenv(STABILITY_DB = "results/stability.sqlite")  # optional; see below
shiny::runApp("app")
```

- `STABILITY_DB` picks the DB to browse; if unset, the app falls back to
  `results/stability.sqlite` (resolved relative to either the repo root or
  `app/`, whichever the app was launched from).
- The dataset list, methods, and job families shown are all queried live
  from the DB, so re-ingesting (e.g. adding `wto_grid`) makes new data
  appear on the next app launch with no app changes needed.

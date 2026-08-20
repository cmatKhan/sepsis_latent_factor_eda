# rslurm stability-grid setup

Metadata-driven setup (not submission) of `rslurm` batch-array jobs that
test stability of latent factors from PCA, NMF, CoGAPS, wTO, and WGCNA
across seeds and across the swept structural parameter (rank / power /
bootstrap-n).

Each config file (`dataset_metadata.yml`) describes **one run against one
matrix**. There's no built-in concept of multiple named subsets ("bases")
within a config -- if you want to analyze several subsets separately (e.g.
Day 0 vs. Day 2 samples), build a separate matrix for each and write a
separate `dataset_metadata.yml` per matrix; each is just another run of
this same framework.

## Quick start

1. **Prepare your input matrix** (see below) and save it as an RDS file.
2. **Write configs**: copy `config/dataset_metadata.example.yml` →
   `config/dataset_metadata.yml`, point `dataset.matrix_path` at your matrix,
   and set each method's parameter grids. Copy
   `config/cluster_config.example.yml` → `config/cluster_config.yml` if your
   cluster's mem/cpus/time/container/module path differs from the example.
3. **Run setup** on the cluster (or source interactively):
   ```r
   Rscript R/run_rslurm_setup.R config/dataset_metadata.yml config/cluster_config.yml
   ```
   `R/run_rslurm_setup.R` reads `dataset_metadata.yml`'s `methods.*.enabled`
   flags and runs whichever of `R/setup_pca_rslurm.R` /
   `setup_nmf_rslurm.R` / `setup_cogaps_rslurm.R` / `setup_network_rslurm.R`
   apply — you don't need to know which individual scripts to call. Each
   underlying setup script builds job grids and calls
   `slurm_apply()`/`slurm_call()` with `submit = FALSE` — it stages the
   sbatch materials without submitting — and writes
   `results/rslurm/<dataset_id>/<method>_sjobs.rds` recording the `sjob_*`
   objects it created. Run a single method directly (e.g.
   `Rscript R/setup_nmf_rslurm.R`) if you only want to redo one.
4. **Submit**: each `_rslurm_<jobname>/` bundle lands under
   `slurm_bundles/<dataset_id>/` (not the directory you ran the script
   from), keyed by `dataset.id` so different datasets/runs never collide or
   overwrite each other's job directories. Either `cd` into one family's
   directory and `sbatch submit.sh` directly (fine on a normal shared
   filesystem), or use the generated `submit_all_<dataset_id>.sh` (see
   "Packaging a bundle as read-only / squashfs" below) to submit every
   family in the bundle at once.

To point at a non-default config location (e.g. running a second subset),
either pass paths as CLI args to `R/run_rslurm_setup.R` (above) or set
these two variables before sourcing any setup script interactively:
```r
dataset_metadata_path <- "config/dataset_metadata_day2.yml"
cluster_config_path   <- "config/cluster_config.yml"
source("R/run_rslurm_setup.R")   # or source("R/setup_nmf_rslurm.R") for just one method
```

## Preparing your input matrix

This framework never bakes in assumptions about normalization or sample
metadata values (e.g. it does not assume a `timepoint` column, "Day 0"/
"Day 2" factor levels, or that VST is the right transform for your data).
`load_input_matrix()` resolves the matrix one of three ways, first match
wins:

| # | Config | Behavior |
|---|---|---|
| 1 | `dataset.matrix_path` set | Reads that RDS file as-is — you built it entirely yourself, however you like. |
| 2 | `dataset.preprocessing_script` set | Runs your script against the raw parquet trio (see contract below). |
| 3 | neither set | Reads `dataset.expression_path` (long format) and pivots it wide — no normalization, no filtering. |

Use (1) if you already have a matrix from some other pipeline. Use (2) if
you need dataset-specific logic (sample/feature filtering, VST/log2/
log-CPM normalization, combining timepoints, anything arbitrary) applied to
the raw HF collection files. Use (3) only if your `value` column is already
in the right numeric form and no filtering is needed — it's a bare
long-to-wide pivot and nothing else.

In every case the resulting matrix must be a plain R matrix with
`rownames` = feature ids and `colnames` = sample ids. There's no separate
post-hoc sample/feature id-list step: if a matrix needs restricting to a
subset, do that inside a preprocessing_script (mode 2) or bake it into the
matrix you point `matrix_path` at (mode 1).

### Preprocessing script contract (mode 2)

Your script is `source()`d with these variables already defined in its
environment:

| Variable | From |
|---|---|
| `expression_path`, `sample_metadata_path`, `feature_metadata_path` | `dataset.*_path` |
| `sample_id_col`, `feature_id_col`, `value_col`, `ensembl_col` | `dataset.*_col` |
| `tmp_dir` | a directory created for this run, guaranteed to persist for the life of the R process running the setup script |

It must write its final feature x sample matrix to
`file.path(tmp_dir, "matrix.rds")`. Any sample/feature filtering the matrix
needs (by timepoint, QC status, candidate gene list, etc.) belongs inside
the script itself. Minimal example:

```r
# R/preprocess_sepsis_day0.R
library(arrow); library(dplyr); library(tidyr); library(DESeq2)

expr_long   <- read_parquet(expression_path)
sample_meta <- read_parquet(sample_metadata_path)

day0_ids <- sample_meta[[sample_id_col]][sample_meta$timepoint == "baseline"]

counts <- expr_long |>
  filter(.data[[sample_id_col]] %in% day0_ids) |>
  pivot_wider(id_cols = all_of(feature_id_col), names_from = all_of(sample_id_col),
              values_from = all_of(value_col)) |>
  tibble::column_to_rownames(feature_id_col) |>
  as.matrix()

dds <- DESeqDataSetFromMatrix(round(counts), data.frame(row.names = colnames(counts)), ~1)
vst_mat <- assay(varianceStabilizingTransformation(dds, blind = TRUE))

saveRDS(vst_mat, file.path(tmp_dir, "matrix.rds"))
```

`tmp_dir` lives under R's per-session temp directory, which isn't cleaned
up until the R process running the setup script exits — long enough for
`submit_job_family()` to serialize the resulting matrix into each job's
input later in that same run. Don't delete it yourself, and don't rely on
it existing across separate `Rscript` invocations.

If you also want to run the script directly (e.g. line-by-line in the
console) to inspect intermediate objects, guard on whether `tmp_dir` is
already defined -- it only is when the framework sourced the script:

```r
framework_run <- exists("tmp_dir", inherits = FALSE)
if (!framework_run) {
  expression_path <- "~/local/path/to/expression.parquet"   # ... etc.
}
# ... your logic ...
if (framework_run) saveRDS(final_mat, file.path(tmp_dir, "matrix.rds"))
```

The only matrix transform this framework performs *for* you regardless of
resolution mode is the non-negativity shift NMF/CoGAPS require
(`shift_nonneg()` in `R/lib/matrices.R`), applied inside
`setup_nmf_rslurm.R`/`setup_cogaps_rslurm.R` — that's an algorithm
requirement, not a normalization choice, so it stays generic.

## Packaging a bundle as read-only / squashfs

Every `_rslurm_<jobname>/` directory rslurm generates is safe to package
read-only, without duplicating any of the resource requests
(mem/cpus/time/array size) rslurm already baked into each family's own
`submit.sh`. Two pieces make that work:

1. **`slurm_run.R` writes results elsewhere.** This project's local
   `config/rslurm_templates/slurm_run_R.txt`/`slurm_run_single_R.txt`
   (used as the default `r_template`) write `results_*.RDS` to
   `Sys.getenv("RSLURM_OUTPUT_DIR")`, falling back to the current
   directory if that's unset.
2. **The container is invoked explicitly, not via a native `--container=`
   SBATCH option.** There's no `#SBATCH --container=` line anymore.
   Instead, `R/lib/submit.R::build_apptainer_rscript_path()` constructs
   the actual `rscript_path` passed to `slurm_apply()`/`slurm_call()` as a
   full command:
   ```
   apptainer run --pwd "/mnt/rslurm_bundle/_rslurm_<jobname>" \
     -B "$SQUASHFS_PATH:/mnt/rslurm_bundle:image-src=/" \
     <container image> Rscript
   ```
   which `config/rslurm_templates/submit_sh.txt`/`submit_single_sh.txt`
   run against `slurm_run.R` (after `eval $(spack load --sh apptainer)`,
   so `apptainer` itself is on `PATH`). Apptainer's native squashfs
   support (`-B archive:dest:image-src=/`, equivalent to
   `--mount type=bind,src=archive,dst=dest,image-src=/` -- see
   [Apptainer's SquashFS docs](https://apptainer.org/docs)) bind-mounts
   the archive into the container itself and sets the container's working
   directory (`--pwd`) to this family's directory inside it, all before
   `Rscript slurm_run.R` runs -- so nothing inside the container needs
   `squashfuse` or loop-mount privilege (which application containers
   essentially never have). `$SQUASHFS_PATH` is a literal, unexpanded
   shell-variable reference in the generated `submit.sh` -- it's resolved
   at job runtime from whatever the job inherited from
   `submit_all_<dataset_id>.sh`'s environment (see below), not baked in at
   bundle-build time (the archive doesn't have a fixed path yet then --
   it hasn't moved to the cluster).

Each `setup_*_rslurm.R` script also writes
`slurm_bundles/<dataset_id>/submit_all_<dataset_id>.sh` -- a small,
dataset-agnostic bash script (not R, since it runs standalone on the
cluster after the bundle has moved) written *inside* the bundle directory
so it travels with it into the squashfs archive. It doesn't know or care
about any family's resources or container: it mounts the archive briefly
on the (uncontainerized) submission host to discover every
`_rslurm_*/submit.sh` inside it, `export`s `SQUASHFS_PATH` and
`RSLURM_OUTPUT_DIR` in its own shell, then submits each family's existing
`submit.sh` as-is via plain `sbatch --output=... submit.sh` -- no
`--export=` flag needed, since Slurm propagates the submitting shell's
environment to the job by default. Workflow:

```bash
# 1. Build the bundle locally (as above), then package it:
mksquashfs slurm_bundles/GSE110487 GSE110487.sqfs

# 2. Move GSE110487.sqfs to the cluster, then copy just the submit
#    script out of it to wherever you want to submit/collect results from:
unsquashfs -d /tmp/x GSE110487.sqfs submit_all_GSE110487.sh
cp /tmp/x/submit_all_GSE110487.sh /scratch/mblab/chasem/latent_factor_eda/

# 3. Run it from there, pointing at the archive:
cd /scratch/mblab/chasem/latent_factor_eda
./submit_all_GSE110487.sh GSE110487.sqfs
# results land in ./GSE110487_results by default, or pass a second
# argument to override: ./submit_all_GSE110487.sh GSE110487.sqfs /some/other/dir

# To (re)submit just one job family instead of all of them, use --only
# with its _rslurm_<jobname> subdirectory name -- a flag rather than a
# positional argument specifically so it can't be silently mistaken for
# results-dir if given in a different order:
./submit_all_GSE110487.sh GSE110487.sqfs --only _rslurm_wgcna_grid
./submit_all_GSE110487.sh GSE110487.sqfs /some/other/dir --only _rslurm_wgcna_grid
```

This assumes `apptainer` is available via `spack load` and `squashfuse`
(or loop-mount privilege) is available on the submission host -- adjust
`config/rslurm_templates/submit_sh.txt`/`submit_single_sh.txt` and
`R/lib/submit_all_script.R` if your cluster's setup differs (e.g. a
different module system than spack, or a non-apptainer container
runtime, which would need a different bind-mount mechanism entirely).

## Container images

Each `cluster_config.yml` method entry's `container` value is passed
straight through to `apptainer run` (see above) -- it can be a local
`.sif` file path or a pullable image reference (`docker://...`,
`oras://...`), whichever your `apptainer` build/version and network
access support. `config/cluster_config.yml` currently points at local
`.sif` files pre-downloaded to `/ref/mblab/containers/chasem/` (faster,
and works without outbound network access from compute nodes), built
from these sources:

| Method | Local `.sif` | Original source |
|---|---|---|
| pca, nmf | `tidyverse_latest.sif` | `docker://rocker/tidyverse:latest` |
| cogaps | `cogaps_sha-3b3e002.sif` | `docker://ghcr.io/fertiglab/cogaps:sha-3b3e002` |
| network.wgcna | `r-wgcna_1.74--5149c638df2976dd.sif` | `oras://community.wave.seqera.io/library/r-wgcna:1.74--5149c638df2976dd` |
| network.wto | `wto_image.sif` | custom-built (installs the `wTO` package on top of a base R image; no public registry source) |

`config/cluster_config.example.yml` still uses the pullable
`docker://`/`oras://` references above, since a new dataset/cluster won't
have these `.sif` files pre-staged -- copy them down with
`apptainer pull <name>.sif <source>` (or reuse an existing `.sif` cache
like this one) and update `container:` accordingly.

### Minimal images: `rscript_path`

Some container images (the WGCNA one here) don't put R on `PATH` --
`rscript_path` in `cluster_config.yml` can be a full path (e.g.
`/opt/conda/bin/Rscript`) instead of just `"Rscript"` if the image's R
lives in an unactivated conda env or somewhere else not on `PATH`. Find it
with `apptainer exec <image> bash -lc 'command -v Rscript'` (a bash
builtin, so it works even if `which` itself is missing) inside an
interactive `srun --container=<image> ...` session.

(A separate issue with this same WGCNA image -- `utils`'s `.onLoad`
crashing because the image has no `which`/coreutils at all -- was worked
around here previously via a bundled `which` shim, but that broke other
methods and has been reverted; it's a known upstream Seqera Wave
container bug being fixed via a rebuilt image instead.)

## Layout

| Path | What |
|---|---|
| `R/lib/` | generic helpers: config validation, matrix loading, job-grid builders, slurm submission wrapper, submit_all script generation |
| `R/methods/*.R` | per-method job functions + which stability designs apply |
| `R/run_rslurm_setup.R` | wrapper — runs whichever `setup_*_rslurm.R` scripts apply, based on `methods.*.enabled` |
| `R/setup_*_rslurm.R` | per-method orchestrators (can be run standalone too) |
| `R/legacy/` | original hardcoded sepsis-only scripts, kept for reference |
| `config/*.example.yml` | copy these to `config/dataset_metadata.yml` / `config/cluster_config.yml` |
| `config/rslurm_templates/` | local copies of rslurm's sh/R templates (R ones patched for `RSLURM_OUTPUT_DIR`; sh ones load apptainer via spack and drop the native `--container=` directive, since the container is invoked explicitly instead -- see squashfs section) |
| `slurm_bundles/<dataset_id>/` | generated `_rslurm_<jobname>/` dirs + `submit_all_<dataset_id>.sh` |

## Stability designs

| Design | Question | Methods |
|---|---|---|
| seed-sweep | Do factors change across random seeds? | PCA (deterministic, degenerate), NMF, CoGAPS |
| masking-CV | Which parameter value best reconstructs held-out entries? | PCA, NMF, CoGAPS |
| param-grid | How do structures change as power/bootstrap-n increases? | wTO, WGCNA |

Not yet generalized (still dataset-specific scripts in `R/legacy/`):
subject-subsample-composition stability, WGCNA cross-basis module
preservation.

## Adding a new method

Add `R/methods/<method>.R` with its `run_*_job()` function(s) (operating on
the global `mat`/`mat_nn` object, no basis argument) + a
`<method>_stability_designs` character vector, then a
`R/setup_<method>_rslurm.R` orchestrator following the existing ones as a
template (read configs → load matrix → build grid via `R/lib/grids.R` →
`submit_job_family()`).

## Stage 2: ingesting results into the stability DB

Once a dataset's job outputs are copied back from the cluster (e.g.
`results/GSE110487_results/`, one subdirectory of `results_<i>.RDS` files
per job family), ingest them into a SQLite DB + external artifacts:

```bash
Rscript R/ingest_results.R <dataset_config.yml> <results_dir> <db_path> [--overwrite [jobname,...]]

# e.g.
Rscript R/ingest_results.R config/GSE110487_config.yml results/GSE110487_results results/stability.sqlite
```

Semantics (fully decoupled from the app; DB path is always a parameter):

- **Additive across datasets**: a second dataset ingests into the same DB
  alongside the first (every table carries `dataset_id`).
- **Additive across job families**: families are discovered as the
  subdirectories of `<results_dir>` containing `results_*.RDS`; families
  already in the DB are skipped with a message. So: ingest now without
  `wto_grid/`, drop it into the results dir when its jobs finish, re-run
  the same command -- only `wto_grid` gets added.
- **Overwrite is explicit and whole-family**: `--overwrite` (bare = every
  family present in the results dir) or `--overwrite nmf_grid,nmf_maskcv`
  deletes and re-ingests those families -- plain replacement, never
  row-level updating.
- Requires the matching `slurm_bundles/<dataset_id>/_rslurm_<jobname>/params.RDS`
  (aligns array task indices to their parameters); missing/failed tasks
  are recorded as `status = 'missing'`/`'failed'` fits, not silently
  dropped.

What gets computed at ingest (see `R/lib/ingest/`):

- **Factorization methods** (PCA/NMF/CoGAPS): loading matrices saved as
  artifacts under `<db_dir>/stability_artifacts/<dataset_id>/`; cosine,
  Pearson, and Spearman similarity for EVERY factor pair across EVERY fit
  pair of a method (all rank pairs x all seed pairs), with Hungarian
  1-to-1 matching (`clue::solve_LSAP` on the cosine matrix); per-factor
  stability summaries (median matched same-rank similarity per metric).
  Incremental: adding a family later computes new x (new + existing)
  pairs only.
- **WGCNA**: gene -> module tables; ARI between module assignments for
  every power pair (`mclust::adjustedRandIndex`); module x module Jaccard
  with Hungarian matching.
- **wTO**: edge tables as parquet artifacts; per-run-pair Pearson/Spearman
  of edge wTO values + Jaccard of significant edges (padj < 0.05).
- **Masking-CV** families: rank (x alpha) -> held-out MSE table.

Artifact paths are stored relative to the DB file's directory, so the DB
and its `stability_artifacts/` folder move together as a unit.

## Stage 3: the Shiny stability explorer

```bash
STABILITY_DB=/path/to/stability.sqlite Rscript -e "shiny::runApp('app')"
# or from R: Sys.setenv(STABILITY_DB = "..."); shiny::runApp("app")
```

`STABILITY_DB` defaults to `results/stability.sqlite` if unset. The app
reads whatever the DB currently contains -- datasets, methods, and job
families all queried live, so newly ingested data appears on the next
launch with zero app changes.

Drill-down levels (breadcrumb at the top navigates back up; a global
metric toggle switches cosine/Pearson/Spearman everywhere):

| Level | Scope | Shows | Drill via |
|---|---|---|---|
| 0 | dataset | per-method headline stability, fit/failure counts, ingested-family inventory | Explore button |
| 1 | method | seed-stability-vs-rank boxplots; masking-CV curves; cross-rank persistence heatmap + factor-tracking trajectories (WGCNA: ARI heatmap, module counts; wTO: run-pair correlation/Jaccard heatmaps) | click a rank / select a power / pick two runs |
| 2 | rank/parameter | seed x seed matched-similarity matrix; per-factor stability strips; factor x factor heatmaps with Hungarian matches outlined (WGCNA: module sizes + cross-power Jaccard; wTO: per-edge scatter) | click a factor |
| 3 | factor | top-loading genes; this factor's Hungarian match in every other fit (all seeds + ranks) with loading scatter; gprofiler2 ORA/GSEA enrichment, run on demand and cached into the DB | -- |

PCA has no cross-seed stability (deterministic, one fit per rank); its
Level-1 seed-stability tab says so and points at masking-CV/cross-rank
instead.

## Deferred

- Subsample-size stability (stage 1 doesn't produce it yet; the schema's
  `family` column accommodates it later).
- Multi-dataset comparison views in the app (schema-ready via `dataset_id`).
- Fetching parquet files directly from the HF Hub (local paths only for now).

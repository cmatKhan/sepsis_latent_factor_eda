# rslurm stability-grid setup

Metadata-driven setup (not submission) of `rslurm` batch-array jobs that
test stability of latent factors from PCA, NMF, CoGAPS, sPCA, CP, Tucker,
and WGCNA across seeds and across the swept structural parameter
(rank / power / per-mode rank).

CP and Tucker (both via the `rTensor` package) decompose a THIRD kind of
input -- a genes x subjects x timepoints array, built from the usual
feature x sample matrix plus two dataset-level sample-metadata columns
(`subject_id_col`/`timepoint_col`) -- rather than the plain feature x
sample matrix every other method uses. They're treated as two entirely
separate methods (own config blocks, own `fits.method`, own ingest
branches) despite sharing `R/lib/tensors.R::build_tensor()` and the
`rTensor` package, since CP (one shared rank across all three modes) and
Tucker (an independent rank per mode) are genuinely different analyses.

Each config file (`dataset_metadata.yml`) describes **one run against one
matrix**. There's no built-in concept of multiple named subsets ("bases")
within a config -- if you want to analyze several subsets separately (e.g.
Day 0 vs. Day 2 samples), build a separate matrix for each and write a
separate `dataset_metadata.yml` per matrix; each is just another run of
this same framework.

## Quick start

1. **Prepare your input matrix** (see below) and save it as an RDS file,
   or point at a `preprocessing_script` (see below).
2. **Write a config**: copy `config/dataset_metadata.example.yml` to
   `config/<dataset>_config.yml` -- this ONE file holds the dataset block,
   SLURM resource settings (`slurm:`, formerly a separate
   `cluster_config.yml`), and every method's parameter grid (`methods:`).
   See "Config reference" below for the full field-by-field explanation.
3. **Run setup** on the cluster (or source interactively):
   ```r
   Rscript R/create_slurm_bundle.R config/<dataset>_config.yml
   # or just some methods:
   Rscript R/create_slurm_bundle.R config/<dataset>_config.yml --only pca,cogaps
   ```
   `R/create_slurm_bundle.R` sets up every method PRESENT in `methods:` --
   presence = enabled, there is no `enabled:` flag; delete/comment out a
   method's block to skip it. It builds each method's job grid and calls
   `slurm_apply()`/`slurm_call()` with `submit = FALSE` -- it stages the
   sbatch materials without submitting -- and writes
   `results/rslurm/<dataset_id>/<method>_sjobs.rds` recording the `sjob_*`
   objects it created.
4. **Submit**: each `_rslurm_<jobname>/` bundle lands under
   `slurm_bundles/<dataset_id>/` (not the directory you ran the script
   from), keyed by `dataset.id` so different datasets/runs never collide or
   overwrite each other's job directories. Either `cd` into one family's
   directory and `sbatch submit.sh` directly (fine on a normal shared
   filesystem), or use the generated `submit_all_<dataset_id>.sh` (see
   "Packaging a bundle for the cluster" below) to submit every
   family in the bundle at once.

To point at a non-default config location (e.g. running a second subset),
either pass the path as a CLI arg to `R/create_slurm_bundle.R` (above) or set
this variable before sourcing it interactively:
```r
dataset_metadata_path <- "config/dataset_metadata_day2.yml"
source("R/create_slurm_bundle.R")
```

## Config reference

Each method is enabled purely by being **present** under `methods:` --
there is no `enabled:` flag anywhere. Each method's block IS its
overrides, directly (e.g. `methods.pca.rank`) -- there's no wrapper key,
since a method now describes exactly one kind of run (a stability-design
second family -- masking-CV, then a subsample-projection design -- and
the `full`/`enabled`/`params` wrapper that came with it, were explored
and deliberately removed; see git history if reviving either).

Every key is the underlying tool's OWN argument name, but a key only
takes effect if that method's `build_grid()` function actually references
it -- there is no automatic generic pass-through for unlisted keys. A
**scalar** value is a static setting; a **list** value is a swept grid
dimension, crossed via `expand.grid()` inside `build_grid()` with every
other swept dimension in the same call. See each `R/methods/<name>.R`
file for exactly which keys it wires through and what its `defaults` are.

| Method | Config keys | Real argument of |
|---|---|---|
| `pca` | `rank` (swept) | *documented exception* -- `prcomp()`'s real argument is `rank.` (trailing dot, a base-R naming quirk not worth reproducing) |
| `nmf` | `k` (swept), `seed` (swept, framework-applied) | `NNLM::nnmf()` -- `seed` has no `nnmf()` equivalent at all; applied via `set.seed()` before the call, per `R/methods/nmf.R` |
| `cogaps.params` | `nPatterns`, `seed`, `nIterations`, `distributed` | `CogapsParams(...)` |
| `cogaps.distributed_params` | `nSets`, `cut`, `minNS`, `maxNS` | `setDistributedParams(params, ...)` -- a sibling key to `params`, not nested inside it |
| `cogaps.run` | `nThreads`, `uncertainty`, ... | `CoGAPS(data, params, ...)` -- a sibling key to `params`, not nested inside it |
| `spca` | `K` (swept), `para` (swept, CROSSED with `K`), `type`, `sparse`, `use.corr`, `lambda`, `max.iter`, `eps.conv` | `elasticnet::spca()` |
| `cp` | `num_components` (swept) | `rTensor::cp()` |
| `tucker` | `rank_genes`/`rank_subjects`/`rank_time` (each independently swept) | *documented exception* -- `rTensor::tucker()`'s real argument is one length-3 vector `ranks = c(r1, r2, r3)`; assembled internally by `run_tucker_job()` from these three independent keys |
| `network.wgcna` | `power` (swept), `minModuleSize`, `mergeCutHeight`, `networkType` | `WGCNA::blockwiseModules()` |

**CoGAPS's `params`/`distributed_params`/`run`** exist as separate
sibling keys directly under `methods.cogaps:` (unlike every other method
here, this nesting is intrinsic to CoGAPS's three call sites, not a
schema wrapper -- there's no flatter form available) because CoGAPS
genuinely spans three call sites with distinct (and in one case,
colliding) argument namespaces -- `checkpointInFile`/`checkpointOutFile`/
`checkpointInterval` are real argument names on BOTH `CogapsParams` and
`CoGAPS()` itself. Put an argument under the sub-block matching which
function it belongs to. `nSets`/`nThreads` default to
`slurm.cogaps.cpus_per_task` unless set explicitly (see
`cogaps_resource_defaults()` in `R/methods/cogaps.R`).

**Watch for YAML boolean-token keys**: R's `yaml` package (2.3.12,
confirmed) parses an UNQUOTED key that happens to be a legacy YAML 1.1
boolean token (`y`/`Y`/`yes`/`Yes`/`YES`/`n`/`N`/`no`/`No`/`NO`/`true`/
`True`/`TRUE`/`false`/`False`/`FALSE`/`on`/`On`/`ON`/`off`/`Off`/`OFF`) as
a literal `TRUE`/`FALSE`, not the string you wrote -- this bit the
now-removed wTO method's real argument name `n` (bootstrap count): it
silently became the key `FALSE`, which downstream became the column name
`FALSE.` (R's `make.names()` escaping the reserved word) and broke the
tool call with `unused argument (FALSE. = 100)`. `read_dataset_metadata()`
checks for a literal `TRUE`/`FALSE` key anywhere in `methods:` and errors
with a pointer to this section, but quote proactively (`"n": [...]`) if
you're ever adding a new single-letter or word-like argument name.

**Tensor methods (`cp`/`tucker`)** additionally need, at the `dataset:`
level (not inside `methods:`): `subject_id_col` and `timepoint_col`, both
columns in `dataset.sample_metadata_path`. `R/lib/tensors.R::build_tensor()`
uses them (plus `sample_id_col`) to reshape the usual feature x sample
matrix into a genes x subjects x timepoints array. `rTensor::cp()`/
`tucker()` require a COMPLETE (dense) array -- there's no native
missing-entry handling -- so any subject missing a timepoint, or with a
duplicate sample at one, is dropped automatically, and this is always
reported (never silent) via a `message()` at setup time listing exactly
which subjects were dropped and why.

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
(`shift_nonneg()` in `R/lib/matrices.R`), applied by `R/create_slurm_bundle.R`
for any method whose registry sets `needs_nonneg = TRUE` -- that's an
algorithm requirement, not a normalization choice, so it stays generic.

## Packaging a bundle for the cluster

Every `_rslurm_<jobname>/` directory rslurm generates travels to the
cluster as a **plain directory copy** (e.g. via `rsync`/`scp`) -- no
squashfs/archive packaging step at all. Two pieces make that work:

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
     -B "$RSLURM_BUNDLE_DIR:/mnt/rslurm_bundle" \
     <container image> Rscript
   ```
   which `config/rslurm_templates/submit_sh.txt`/`submit_single_sh.txt`
   run against `slurm_run.R` (after `eval $(spack load --sh apptainer)`,
   so `apptainer` itself is on `PATH`). This is an ordinary directory
   bind-mount (no squashfs image involved) -- `$RSLURM_BUNDLE_DIR` is a
   literal, unexpanded shell-variable reference in the generated
   `submit.sh` -- it's resolved at job runtime from whatever the job
   inherited from `submit_all_<dataset_id>.sh`'s environment (see below),
   not baked in at bundle-build time (the bundle doesn't have a fixed
   path yet then -- it hasn't moved to the cluster).

Each `setup_*_rslurm.R` script also writes
`slurm_bundles/<dataset_id>/submit_all_<dataset_id>.sh` -- a small,
dataset-agnostic bash script (not R, since it runs standalone on the
cluster after the bundle has moved) written *inside* the bundle directory
so it travels with it wherever it's copied. It doesn't know or care about
any family's resources or container: it discovers every
`_rslurm_*/submit.sh` directly inside the (already-copied) bundle
directory, `export`s `RSLURM_BUNDLE_DIR` and `RSLURM_OUTPUT_DIR` in its
own shell, then submits each family's existing `submit.sh` as-is via
plain `sbatch --output=... submit.sh` -- no `--export=` flag needed,
since Slurm propagates the submitting shell's environment to the job by
default. Workflow:

```bash
# 1. Build the bundle locally (as above), then copy the whole directory
#    to wherever you want to submit/collect results from -- no
#    packaging step:
rsync -av slurm_bundles/GSE110487/ user@cluster:/scratch/mblab/chasem/latent_factor_eda/GSE110487_bundle/

# 2. Run it from there, pointing at the copied bundle directory:
ssh user@cluster
cd /scratch/mblab/chasem/latent_factor_eda
./GSE110487_bundle/submit_all_GSE110487.sh ./GSE110487_bundle
# results land in ./GSE110487_results by default, or pass a second
# argument to override: ./GSE110487_bundle/submit_all_GSE110487.sh ./GSE110487_bundle /some/other/dir

# To (re)submit just one job family instead of all of them, use --only
# with its _rslurm_<jobname> subdirectory name -- a flag rather than a
# positional argument specifically so it can't be silently mistaken for
# results-dir if given in a different order:
./GSE110487_bundle/submit_all_GSE110487.sh ./GSE110487_bundle --only _rslurm_wgcna_grid
./GSE110487_bundle/submit_all_GSE110487.sh ./GSE110487_bundle /some/other/dir --only _rslurm_wgcna_grid
```

This assumes `apptainer` is available via `spack load` on the submission
host -- adjust `config/rslurm_templates/submit_sh.txt`/
`submit_single_sh.txt` and `R/lib/submit_all_script.R` if your cluster's
setup differs (e.g. a different module system than spack, or a
non-apptainer container runtime, which would need a different
bind-mount mechanism entirely).

## Container images

Each `slurm:` method entry's `container` value is passed straight through
to `apptainer run` (see above) -- it can be a local `.sif` file path or a
pullable image reference (`docker://...`, `oras://...`), whichever your
`apptainer` build/version and network access support. This project's
dataset configs point at local `.sif` files pre-downloaded to
`/ref/mblab/containers/chasem/` (faster, and works without outbound
network access from compute nodes), built from these sources:

| Method | Local `.sif` | Original source |
|---|---|---|
| pca, nmf | `tidyverse_latest.sif` | `docker://rocker/tidyverse:latest` |
| cogaps | `cogaps_sha-3b3e002.sif` | `docker://ghcr.io/fertiglab/cogaps:sha-3b3e002` |
| network.wgcna | `r-wgcna_1.74--5149c638df2976dd.sif` | `oras://community.wave.seqera.io/library/r-wgcna:1.74--5149c638df2976dd` |
| spca, cp, tucker | **TODO** -- not built yet | need an image with `elasticnet` (spca) / `rTensor` (cp, tucker) installed; `container: "TODO"` in every dataset config's `slurm.spca`/`slurm.cp`/`slurm.tucker` is a deliberate placeholder, not a bug -- both packages are only confirmed available in local/interactive R for now |

`config/dataset_metadata.example.yml`'s `slurm:` section uses these same
local paths as a starting point -- a new dataset/cluster without them
pre-staged should use the pullable `docker://`/`oras://` references above
instead: copy them down with `apptainer pull <name>.sif <source>` (or
reuse an existing `.sif` cache like this one) and update `container:`
accordingly.

### Minimal images: `rscript_path`

Some container images (the WGCNA one here) don't put R on `PATH` --
`rscript_path` in a dataset config's `slurm:` entry can be a full path
(e.g. `/opt/conda/bin/Rscript`) instead of just `"Rscript"` if the image's R
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
| `R/lib/` | generic helpers: config validation, method-registry discovery, matrix loading, slurm submission wrapper, submit_all script generation. No generic job-grid builder anymore -- each method builds its own (see `R/methods/*.R`) |
| `R/methods/*.R` | per-method job function(s) + a flat registry list (`<method>_registry`) with its own `build_grid()` function -- the generic driver just calls what each method provides. See "Adding a new method" |
| `R/create_slurm_bundle.R` | THE generic driver -- sets up every method present in a config's `methods:` block; `--only a,b` restricts to specific methods |
| `R/legacy/` | retired scripts, kept for reference: the original hardcoded sepsis-only scripts, plus the pre-config-driven-registry `setup_*_rslurm.R`/`run_rslurm_setup.R`/`cluster_config*.yml` |
| `config/dataset_metadata.example.yml` | copy this to `config/<dataset>_config.yml` -- one file per dataset, includes dataset/slurm/methods |
| `config/rslurm_templates/` | local copies of rslurm's sh/R templates (R ones patched for `RSLURM_OUTPUT_DIR`; sh ones load apptainer via spack and drop the native `--container=` directive, since the container is invoked explicitly instead -- see "Packaging a bundle for the cluster") |
| `slurm_bundles/<dataset_id>/` | generated `_rslurm_<jobname>/` dirs + `submit_all_<dataset_id>.sh` |

## Stability designs

Every method's registry describes exactly one kind of run -- but it
serves conceptually different stability questions depending on the
method:

| Design | Question | Methods |
|---|---|---|
| seed-sweep | Do factors change across random seeds? | PCA (deterministic, degenerate), NMF, CoGAPS |
| param-grid | How do structures change as the swept parameter(s) increase? | WGCNA (power), sPCA (K/para), CP (num_components), Tucker (rank_genes/rank_subjects/rank_time) -- all deterministic given their parameters, so (like PCA) there's no genuine cross-seed stability question |

A held-out-entry reconstruction design (masking-CV) and a held-out-subject
projection design (fitting on a subsample and projecting the rest onto the
learned basis) were both explored and deliberately removed for now -- see
git history if reviving either. Not yet implemented at all: WGCNA
cross-basis module preservation.

## Adding a new method

`R/create_slurm_bundle.R` discovers methods automatically, via
`discover_method_registries()` (see `R/lib/method_registry.R`): it
`source()`s every `R/methods/*.R` file and builds its manifest from
whatever registry object each one defines. **Adding a method means adding
one file -- nothing in `R/create_slurm_bundle.R` or `R/lib/metadata.R`
needs to change.**

Write ONE `R/methods/<name>.R` file:
1. Job function(s) operating on the global `mat`/`mat_nn`/`tnsr` object (no
   basis argument) -- formals are whatever your `build_grid()` puts in its
   output columns (see below). If your tool has multiple call sites like
   CoGAPS, give your job function one formal per named sub-block instead
   (each arrives as a proper named list -- see `R/methods/cogaps.R`).
2. A `<name>_registry` list -- **the file name and this variable name must
   match** (`R/methods/foo.R` -> `foo_registry`; discovery errors loudly at
   bundle-creation time if they don't). A flat list (no per-family
   nesting -- each method describes exactly one kind of run) with
   required fields:
   - `global_object` -- the matrix/tensor variable name your job function
     expects.
   - `jobname` -- passed to `submit_job_family()`/rslurm (e.g. `"pca_grid"`).
   - `fn` -- the job function.
   - `build_grid` -- `function(params)` returning a data.frame ready for
     `submit_job_family()`'s `jobs_df`: YOUR OWN sweep/cross logic
     (typically just `expand.grid(...)` over the specific arguments your
     job function needs -- see `R/methods/pca.R` for the simple case,
     `R/methods/cogaps.R` for one with a nested sub-block). `params` is
     this method's `defaults` merged with the dataset config's overrides
     (`merge_named_list()` in `R/lib/grids.R` recurses one level into any
     key that's a named list in both, e.g. CoGAPS's `params` sub-block, so
     overriding one nested key doesn't wipe out its siblings' defaults).

   Optional fields:
   - `defaults` (default `list()`) -- named list of default parameter
     values, overridable per-key by the dataset config.
   - `pkgs` -- packages the job function needs.
   - `resource_defaults` -- `function(params, slurm_cfg)` for defaults
     derived from cluster config (e.g. CoGAPS's nSets/nThreads).

   Two further optional fields (default `FALSE` if omitted):
   - `needs_nonneg` / `needs_tensor` -- which cached input
     (`mat`/`mat_nn`/`tnsr`) `global_object` refers to.
   - `network` -- `TRUE` if this method's config/slurm entries nest under
     `methods.network.<name>`/`slurm.network.<name>` instead of the flat
     `methods.<name>`/`slurm.<name>` every other method uses (WGCNA today).
   - `requires_subject_timepoint` -- `TRUE` if this method needs
     `dataset.subject_id_col`/`timepoint_col`/`sample_metadata_path`
     (CP/Tucker today); checked generically by `R/lib/metadata.R` off this
     flag rather than a hardcoded method-name list.

   Any existing `R/methods/*.R` is a template -- `R/methods/pca.R` for the
   simplest case, `R/methods/cogaps.R` for nested sub-blocks,
   `R/methods/wgcna.R` for `network`,
   `R/methods/cp.R`/`R/methods/tucker.R` for `requires_subject_timepoint`.
3. That's it -- run `R/create_slurm_bundle.R` (or restart your R session
   first, so the new file gets sourced) and the method shows up. If the
   registry variable is missing or misnamed, or a required field is
   missing, `discover_method_registries()` stops with a message naming the
   exact problem before any slurm job is built.

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
  `cp_grid/`, drop it into the results dir when its jobs finish, re-run
  the same command -- only `cp_grid` gets added.
- **Overwrite is explicit and whole-family**: `--overwrite` (bare = every
  family present in the results dir) or `--overwrite nmf_grid,nmf_maskcv`
  deletes and re-ingests those families -- plain replacement, never
  row-level updating.
- Requires the matching `slurm_bundles/<dataset_id>/_rslurm_<jobname>/params.RDS`
  (aligns array task indices to their parameters); missing/failed tasks
  are recorded as `status = 'missing'`/`'failed'` fits, not silently
  dropped.

What gets computed at ingest (see `R/lib/ingest/`):

- **Factorization methods** (PCA/NMF/CoGAPS/sPCA/CP/Tucker): loading
  matrices saved as artifacts under
  `<db_dir>/stability_artifacts/<dataset_id>/`; cosine, Pearson, and
  Spearman similarity for EVERY factor pair across EVERY fit pair of a
  method (all rank pairs x all seed pairs, where applicable), with
  Hungarian 1-to-1 matching (`clue::solve_LSAP` on the cosine matrix);
  per-factor stability summaries (median matched same-rank similarity per
  metric). Incremental: adding a family later computes new x (new +
  existing) pairs only. **CP/Tucker caveat**: `factor_pairs.same_rank` is
  keyed off the single `fits.rank` column, which Tucker never sets (it has
  `rank_genes`/`rank_subjects`/`rank_time` instead) -- pairs are still
  computed (no crash), just never flagged `same_rank` against each other;
  this framework doesn't yet have a "same 3-tuple of ranks" concept.
- **WGCNA**: gene -> module tables; ARI between module assignments for
  every power pair (`mclust::adjustedRandIndex`); module x module Jaccard
  with Hungarian matching.
- **Masking-CV** families: rank (x alpha) -> held-out MSE table.
- **CP/Tucker only**: a third, time-mode loading matrix (rows = timepoint
  levels) saved alongside the usual feature-loadings/sample-scores
  artifacts -- `fits.time_loadings_file`. Their "scores" artifact is
  SUBJECT-mode, not sample-mode (one row per subject, not per sample) --
  the app's sample-metadata-association machinery doesn't yet reconcile
  this distinction (deferred, see below).

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
| 1 | method | seed-stability-vs-rank boxplots; masking-CV curves; cross-rank persistence heatmap + factor-tracking trajectories (WGCNA: ARI heatmap, module counts) | click a rank / select a power |
| 2 | rank/parameter | seed x seed matched-similarity matrix; per-factor stability strips; factor x factor heatmaps with Hungarian matches outlined (WGCNA: module sizes + cross-power Jaccard) | click a factor |
| 3 | factor | top-loading genes; this factor's Hungarian match in every other fit (all seeds + ranks) with loading scatter; gprofiler2 ORA/GSEA enrichment, run on demand and cached into the DB | -- |

PCA has no cross-seed stability (deterministic, one fit per rank) -- its
Level 1 is reduced to just the masking-CV rank-selection curve, and it
serves instead as a comparison BASELINE from NMF/CoGAPS's (and WGCNA's)
Level 2 "Compare to PCA" tab: held-out MSE side by side, signed factor/
eigengene similarity (Hungarian-matched by absolute value, so a strong
match to a PCA component's negative side shows up correctly), and a
side-by-side listing of whatever enrichment has already been queried for
each side.

## Deferred

- Subsample-size stability (stage 1 doesn't produce it yet; the schema's
  `family` column accommodates it later).
- Multi-dataset comparison views in the app (schema-ready via `dataset_id`).
- Fetching parquet files directly from the HF Hub (local paths only for now).
- **App (Stage 3) UI for sPCA/CP/Tucker**: stage 1 (setup) and stage 2
  (ingest) fully support all three (verified end to end against real
  GSE110487 data), but `app/app.R` doesn't have Level 1/2/3 views wired up
  for them yet -- a method card would appear at Level 0 once ingested
  (method/fit counts are queried live from the DB), but "Explore" won't do
  anything (no `open_<method>` observer registered). CP/Tucker in
  particular need real design work before that's straightforward: their
  factor view has a THIRD (time-mode) matrix the current Level 3
  "Loadings"/"Sample scores" pair doesn't have a slot for, and their
  `scores` artifact is subject-mode rather than sample-mode (see the
  ingest section above).
- CP/Tucker container images (`slurm.{spca,cp,tucker}.container` is
  `"TODO"` in every dataset config) -- `elasticnet`/`rTensor` are only
  confirmed available in local/interactive R for now.

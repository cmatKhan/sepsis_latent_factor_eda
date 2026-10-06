# Sepsis latent-factor stability

Fits six latent-factor methods (PCA, NMF, CoGAPS, sPCA, ICA, WGCNA) to sepsis gene-expression
datasets across seeds and rank/sparsity choices, measures how stable the factors are, enriches
and projects them across datasets, and writes it all to a SQLite database that a Shiny app
browses. One [targets](https://books.ropensci.org/targets/) pipeline runs everything, with the
heavy steps as SLURM jobs.

**Full documentation:** the Quarto site in [`docs/`](docs/) -- `quarto preview docs` to browse it.

## Prerequisites

- A SLURM cluster (HTCF here). The main process submits worker jobs and must stay running.
- R 4.6 with system libraries: on HTCF, the spack environment
  `/ref/mblab/software/chasem/spack_envs/rstudio-4.6.1` (the one RStudio Server runs from).
  Elsewhere, see [docs/setup.qmd](docs/setup.qmd).
- The raw data (HuggingFace sepsis collection). On HTCF:
  `/scratch/mblab/chasem/hf_sepsis_collection`.

## Set up (once)

```sh
git clone <this repo> && cd sepsis_latent_factor_eda
# R from the spack environment (RStudio Server on HTCF already has it):
SPACK_ENV=/ref/mblab/software/chasem/spack_envs/rstudio-4.6.1
eval "$(spack env activate --sh "$SPACK_ENV")"
export LD_LIBRARY_PATH="$SPACK_ENV/view/lib:$SPACK_ENV/view/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
unset R_LIBS R_LIBS_USER R_LIBS_SITE
Rscript -e 'options(renv.install.timeout = 6 * 3600); renv::restore()'   # hours: builds from source
```

CoGAPS needs a patched build for OpenMP, and projectR comes from a fork; both are covered in
[docs/setup.qmd](docs/setup.qmd). Point `data_root` in `config/pipeline.yml` (or the
`SEPSIS_DATA_ROOT` environment variable) at the data.

## Add a dataset

1. Write `config/<dataset>_config.yml` (copy `config/dataset_metadata.example.yml`): data paths,
   how to build the matrix, each method's parameter grid. See [docs/data.qmd](docs/data.qmd).
2. Add its `dataset.id` to `datasets:` in `config/pipeline.yml`.

## Run

```sh
sbatch run_pipeline.sbatch                                     # build everything outdated
sbatch run_pipeline.sbatch 'tidyselect::starts_with("fit_")'   # fits only, no ingestion
tail -f logs/targets_main_<jobid>.out
```

Or from an RStudio Server session on the cluster: `targets::tar_make()` (and
`targets::tar_visnetwork()` to see what's outdated). Only what a change invalidates is rebuilt.
See [docs/workflow.qmd](docs/workflow.qmd).

## Results

The database is at `db_path` in `config/pipeline.yml` (default
`results/targets/stability.sqlite`), with large artifacts under
`results/targets/stability_artifacts/`. Every table is described in
[docs/database.qmd](docs/database.qmd).

## Browse

```r
shiny::runApp("app")
```

The app has not been updated to the current database schema yet; see
[docs/app.qmd](docs/app.qmd).

## Documentation

```sh
quarto render docs     # builds docs/_site/ (Quarto ships with RStudio Server)
quarto preview docs    # live preview
```

Pages: [setup](docs/setup.qmd), [data](docs/data.qmd), [methods](docs/methods.qmd),
[database](docs/database.qmd), [running the workflow](docs/workflow.qmd),
[launching the app](docs/app.qmd), and a function reference generated from the code's roxygen
comments.

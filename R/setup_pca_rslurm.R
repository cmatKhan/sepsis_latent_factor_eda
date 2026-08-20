# PCA rslurm job setup -- metadata-driven, operates on ONE matrix (see
# config/dataset_metadata.example.yml and R/README.md's "Preparing your
# input matrix"). To run PCA on a different sample/feature subset, build a
# separate matrix and dataset_metadata.yml and rerun this script against it.
#
# Builds two job families:
#   1. seed-sweep (degenerate single fit per rank -- PCA is deterministic,
#      kept for shape-parity with the other methods' rank-vs-MSE curves)
#   2. masking-CV rank selection (held-out-entry reconstruction MSE)
#
# Run out-of-band: `Rscript R/setup_pca_rslurm.R`.

library(here)
source(here("R/lib/metadata.R"))
source(here("R/lib/matrices.R"))
source(here("R/lib/grids.R"))
source(here("R/lib/submit.R"))
source(here("R/methods/pca.R"))

if (!exists("dataset_metadata_path")) dataset_metadata_path <- here("config/dataset_metadata.yml")
if (!exists("cluster_config_path")) cluster_config_path <- here("config/cluster_config.yml")

dataset_meta <- read_dataset_metadata(dataset_metadata_path)
cluster_cfg  <- read_cluster_config(cluster_config_path)

mat <- load_input_matrix(dataset_meta)

pca_meta   <- dataset_meta$methods$pca
rank_range <- rank_seq(pca_meta$rank_range)

# All _rslurm_<jobname> bundles land under slurm_bundles/<dataset_id>/ (not
# the working directory this script happens to be run from), so different
# datasets/runs never collide or overwrite each other's job directories.
job_output_dir <- here("slurm_bundles", dataset_meta$dataset$id)

## ---- Family 1: seed-sweep (single fit per rank) ---------------------------

seed_sweep_grid <- build_grid_seed_sweep(rank_range, seeds = NA_integer_)
seed_sweep_grid$seed <- NULL   # PCA is deterministic; drop the degenerate seed column

sjob_pca_grid <- submit_job_family(
  f              = run_pca_seed_sweep_job,
  jobs_df        = seed_sweep_grid,
  jobname        = "pca_grid",
  global_objects = c("mat"),
  pkgs           = character(0),
  cluster_cfg    = cluster_cfg$pca,
  output_dir     = job_output_dir
)

## ---- Family 2: masking-CV rank selection ----------------------------------

sjob_pca_maskcv <- NULL
if (isTRUE(pca_meta$masking_cv$enabled)) {
  mask_idx <- make_mask_index(mat, mask_frac = pca_meta$masking_cv$mask_frac,
                               mask_seed = pca_meta$masking_cv$mask_seed)

  maskcv_grid <- build_grid_masking_cv(rank_range)

  sjob_pca_maskcv <- submit_job_family(
    f              = run_pca_masking_cv_job,
    jobs_df        = maskcv_grid,
    jobname        = "pca_maskcv",
    global_objects = c("mat", "mask_idx"),
    pkgs           = character(0),
    cluster_cfg    = cluster_cfg$pca,
    output_dir     = job_output_dir
  )
}

save_sjobs(list(grid = sjob_pca_grid, maskcv = sjob_pca_maskcv),
           dataset_id = dataset_meta$dataset$id, method = "pca")

# Regenerates slurm_bundles/<dataset_id>/submit_all_<dataset_id>.sh from
# every job family recorded so far (this method's plus any others already
# run for this dataset) -- safe/expected to call after every setup script.
write_submit_all_script(dataset_meta$dataset$id, job_output_dir)

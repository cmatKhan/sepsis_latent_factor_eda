# CoGAPS rslurm job setup -- metadata-driven, operates on ONE matrix (see
# config/dataset_metadata.example.yml and R/README.md's "Preparing your
# input matrix"). To run CoGAPS on a different sample/feature subset, build
# a separate matrix and dataset_metadata.yml and rerun this script against it.
#
# Builds two job families:
#   1. seed-sweep (rank x seed grid, full in-sample reconstruction error)
#   2. masking-CV rank/alpha selection (held-out-entry reconstruction MSE
#      via CoGAPS's uncertainty-matrix masking mechanism)
#
# Run out-of-band: `Rscript R/setup_cogaps_rslurm.R`.

library(here)
source(here("R/lib/metadata.R"))
source(here("R/lib/matrices.R"))
source(here("R/lib/grids.R"))
source(here("R/lib/submit.R"))
source(here("R/methods/cogaps.R"))

if (!exists("dataset_metadata_path")) dataset_metadata_path <- here("config/dataset_metadata.yml")
if (!exists("cluster_config_path")) cluster_config_path <- here("config/cluster_config.yml")

dataset_meta <- read_dataset_metadata(dataset_metadata_path)
cluster_cfg  <- read_cluster_config(cluster_config_path)

mat_nn <- shift_nonneg(load_input_matrix(dataset_meta))   # CoGAPS requires non-negative input

cogaps_meta  <- dataset_meta$methods$cogaps
rank_range   <- rank_seq(cogaps_meta$rank_range)
n_iterations <- cogaps_meta$n_iterations
nsets        <- cluster_cfg$cogaps$cpus_per_task

# All _rslurm_<jobname> bundles land under slurm_bundles/<dataset_id>/ (not
# the working directory this script happens to be run from), so different
# datasets/runs never collide or overwrite each other's job directories.
job_output_dir <- here("slurm_bundles", dataset_meta$dataset$id)

## ---- Family 1: seed-sweep --------------------------------------------------

seed_sweep_grid <- build_grid_seed_sweep(rank_range, cogaps_meta$seeds)
seed_sweep_grid$n_iterations <- n_iterations
seed_sweep_grid$nsets        <- nsets

sjob_cogaps_grid <- submit_job_family(
  f              = run_cogaps_seed_sweep_job,
  jobs_df        = seed_sweep_grid,
  jobname        = "cogaps_grid",
  global_objects = c("mat_nn"),
  pkgs           = c("CoGAPS"),
  cluster_cfg    = cluster_cfg$cogaps,
  output_dir     = job_output_dir
)

## ---- Family 2: masking-CV rank/alpha selection ----------------------------

sjob_cogaps_maskcv <- NULL
if (isTRUE(cogaps_meta$masking_cv$enabled)) {
  mask_idx    <- make_mask_index(mat_nn, mask_frac = cogaps_meta$masking_cv$mask_frac,
                                  mask_seed = cogaps_meta$masking_cv$mask_seed)
  uncertainty <- build_uncertainty(mat_nn)

  maskcv_grid <- build_grid_masking_cv(
    rank_range,
    extra = list(alpha = cogaps_meta$masking_cv$alpha_range)
  )
  maskcv_grid$maskcv_seed  <- cogaps_meta$masking_cv$maskcv_seed
  maskcv_grid$n_iterations <- n_iterations
  maskcv_grid$nsets        <- nsets

  sjob_cogaps_maskcv <- submit_job_family(
    f              = run_cogaps_masking_cv_job,
    jobs_df        = maskcv_grid,
    jobname        = "cogaps_maskcv",
    global_objects = c("mat_nn", "mask_idx", "uncertainty"),
    pkgs           = c("CoGAPS"),
    cluster_cfg    = cluster_cfg$cogaps,
    output_dir     = job_output_dir
  )
}

save_sjobs(list(grid = sjob_cogaps_grid, maskcv = sjob_cogaps_maskcv),
           dataset_id = dataset_meta$dataset$id, method = "cogaps")

# Regenerates slurm_bundles/<dataset_id>/submit_all_<dataset_id>.sh from
# every job family recorded so far (this method's plus any others already
# run for this dataset) -- safe/expected to call after every setup script.
write_submit_all_script(dataset_meta$dataset$id, job_output_dir)

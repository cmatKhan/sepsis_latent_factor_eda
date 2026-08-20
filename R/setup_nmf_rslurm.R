# NMF rslurm job setup -- metadata-driven, operates on ONE matrix (see
# config/dataset_metadata.example.yml and R/README.md's "Preparing your
# input matrix"). To run NMF on a different sample/feature subset, build a
# separate matrix and dataset_metadata.yml and rerun this script against it.
#
# Builds two job families:
#   1. seed-sweep (rank x seed grid, full in-sample reconstruction error)
#   2. masking-CV rank selection (held-out-entry reconstruction MSE, via
#      NNLM::nnmf's native NA handling)
#
# Run out-of-band: `Rscript R/setup_nmf_rslurm.R`.

library(here)
source(here("R/lib/metadata.R"))
source(here("R/lib/matrices.R"))
source(here("R/lib/grids.R"))
source(here("R/lib/submit.R"))
source(here("R/methods/nmf.R"))

if (!exists("dataset_metadata_path")) dataset_metadata_path <- here("config/dataset_metadata.yml")
if (!exists("cluster_config_path")) cluster_config_path <- here("config/cluster_config.yml")

dataset_meta <- read_dataset_metadata(dataset_metadata_path)
cluster_cfg  <- read_cluster_config(cluster_config_path)

mat_nn <- shift_nonneg(load_input_matrix(dataset_meta))   # NMF requires non-negative input

nmf_meta   <- dataset_meta$methods$nmf
rank_range <- rank_seq(nmf_meta$rank_range)

# All _rslurm_<jobname> bundles land under slurm_bundles/<dataset_id>/ (not
# the working directory this script happens to be run from), so different
# datasets/runs never collide or overwrite each other's job directories.
job_output_dir <- here("slurm_bundles", dataset_meta$dataset$id)

## ---- Family 1: seed-sweep --------------------------------------------------

seed_sweep_grid <- build_grid_seed_sweep(rank_range, nmf_meta$seeds)

sjob_nmf_grid <- submit_job_family(
  f              = run_nmf_seed_sweep_job,
  jobs_df        = seed_sweep_grid,
  jobname        = "nmf_grid",
  global_objects = c("mat_nn"),
  pkgs           = c("NNLM"),
  cluster_cfg    = cluster_cfg$nmf,
  output_dir     = job_output_dir
)

## ---- Family 2: masking-CV rank selection ----------------------------------

sjob_nmf_maskcv <- NULL
if (isTRUE(nmf_meta$masking_cv$enabled)) {
  mask_idx <- make_mask_index(mat_nn, mask_frac = nmf_meta$masking_cv$mask_frac,
                               mask_seed = nmf_meta$masking_cv$mask_seed)

  maskcv_grid <- build_grid_masking_cv(rank_range)

  sjob_nmf_maskcv <- submit_job_family(
    f              = run_nmf_masking_cv_job,
    jobs_df        = maskcv_grid,
    jobname        = "nmf_maskcv",
    global_objects = c("mat_nn", "mask_idx"),
    pkgs           = c("NNLM"),
    cluster_cfg    = cluster_cfg$nmf,
    output_dir     = job_output_dir
  )
}

save_sjobs(list(grid = sjob_nmf_grid, maskcv = sjob_nmf_maskcv),
           dataset_id = dataset_meta$dataset$id, method = "nmf")

# Regenerates slurm_bundles/<dataset_id>/submit_all_<dataset_id>.sh from
# every job family recorded so far (this method's plus any others already
# run for this dataset) -- safe/expected to call after every setup script.
write_submit_all_script(dataset_meta$dataset$id, job_output_dir)

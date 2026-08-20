# Network-method rslurm job setup -- metadata-driven, operates on ONE matrix
# (see config/dataset_metadata.example.yml and R/README.md's "Preparing your
# input matrix"), covering both backends declared in dataset_metadata.yml's
# `methods$network` block:
#   - wTO:   n (bootstrap replicate count) x [seed] param grid
#   - WGCNA: power (soft-threshold) param grid
#
# To run on a different sample/feature subset, build a separate matrix and
# dataset_metadata.yml and rerun this script against it.
#
# Neither backend has a masking-CV analogue (see R/methods/wto.R and
# R/methods/wgcna.R header comments for why), so each gets a single
# `param_grid` family.
#
# Run out-of-band: `Rscript R/setup_network_rslurm.R`.

library(here)
source(here("R/lib/metadata.R"))
source(here("R/lib/matrices.R"))
source(here("R/lib/grids.R"))
source(here("R/lib/submit.R"))
source(here("R/methods/wto.R"))
source(here("R/methods/wgcna.R"))

if (!exists("dataset_metadata_path")) dataset_metadata_path <- here("config/dataset_metadata.yml")
if (!exists("cluster_config_path")) cluster_config_path <- here("config/cluster_config.yml")

dataset_meta <- read_dataset_metadata(dataset_metadata_path)
cluster_cfg  <- read_cluster_config(cluster_config_path)

mat <- load_input_matrix(dataset_meta)

network_meta <- dataset_meta$methods$network

# All _rslurm_<jobname> bundles land under slurm_bundles/<dataset_id>/ (not
# the working directory this script happens to be run from), so different
# datasets/runs never collide or overwrite each other's job directories.
job_output_dir <- here("slurm_bundles", dataset_meta$dataset$id)

## ---- wTO -------------------------------------------------------------------

sjob_wto <- NULL
if (isTRUE(network_meta$wto$enabled)) {
  wto_k <- cluster_cfg$network$wto$cpus_per_task

  seeds <- network_meta$wto$seeds
  if (is.null(seeds)) seeds <- NA_integer_

  wto_grid <- build_grid_param(list(
    n     = network_meta$wto$n_range,
    delta = network_meta$wto$delta,
    seed  = seeds
  ))
  wto_grid$wto_k <- wto_k

  sjob_wto <- submit_job_family(
    f              = run_wto_param_job,
    jobs_df        = wto_grid,
    jobname        = "wto_grid",
    global_objects = c("mat"),
    pkgs           = c("wTO"),
    cluster_cfg    = cluster_cfg$network$wto,
    output_dir     = job_output_dir
  )
}

## ---- WGCNA -------------------------------------------------------------------

sjob_wgcna <- NULL
if (isTRUE(network_meta$wgcna$enabled)) {
  wgcna_cpu_per_task <- cluster_cfg$network$wgcna$cpus_per_task

  wgcna_grid <- build_grid_param(list(power = network_meta$wgcna$power_grid))
  wgcna_grid$min_module_size    <- network_meta$wgcna$min_module_size
  wgcna_grid$merge_cut_height   <- network_meta$wgcna$merge_cut_height
  wgcna_grid$network_type       <- network_meta$wgcna$network_type
  wgcna_grid$wgcna_cpu_per_task <- wgcna_cpu_per_task

  sjob_wgcna <- submit_job_family(
    f              = run_wgcna_param_job,
    jobs_df        = wgcna_grid,
    jobname        = "wgcna_grid",
    global_objects = c("mat"),
    pkgs           = c("WGCNA"),
    cluster_cfg    = cluster_cfg$network$wgcna,
    output_dir     = job_output_dir
  )
}

save_sjobs(list(wto = sjob_wto, wgcna = sjob_wgcna),
           dataset_id = dataset_meta$dataset$id, method = "network")

# Regenerates slurm_bundles/<dataset_id>/submit_all_<dataset_id>.sh from
# every job family recorded so far (this method's plus any others already
# run for this dataset) -- safe/expected to call after every setup script.
write_submit_all_script(dataset_meta$dataset$id, job_output_dir)

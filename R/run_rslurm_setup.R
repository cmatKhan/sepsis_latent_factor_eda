# Wrapper: runs whichever R/setup_*_rslurm.R scripts are relevant for a
# given dataset_metadata.yml, based on each method's `enabled` flag --
# so you don't need to know (or separately invoke) setup_pca_rslurm.R /
# setup_nmf_rslurm.R / setup_cogaps_rslurm.R / setup_network_rslurm.R
# yourself. Every underlying setup script still runs standalone too, if you
# ever want to rerun just one method.
#
# Usage (CLI):
#   Rscript R/run_rslurm_setup.R config/GSE110487_config.yml config/cluster_config.yml
#
# Usage (interactive, matching the individual setup scripts' convention):
#   dataset_metadata_path <- "config/GSE110487_config.yml"
#   cluster_config_path   <- "config/cluster_config.yml"
#   source("R/run_rslurm_setup.R")

library(here)
source(here("R/lib/metadata.R"))

args <- commandArgs(trailingOnly = TRUE)
if (!exists("dataset_metadata_path")) {
  dataset_metadata_path <- if (length(args) >= 1) args[[1]] else here("config/dataset_metadata.yml")
}
if (!exists("cluster_config_path")) {
  cluster_config_path <- if (length(args) >= 2) args[[2]] else here("config/cluster_config.yml")
}

dataset_meta <- read_dataset_metadata(dataset_metadata_path)
methods_meta <- dataset_meta$methods

run_if_enabled <- function(enabled, script, label) {
  if (isTRUE(enabled)) {
    message("== ", label, " (", script, ") ==")
    source(here(script))
  } else {
    message("-- skipping ", label, " (not enabled in ", basename(dataset_metadata_path), ") --")
  }
}

run_if_enabled(methods_meta$pca$enabled,    "R/setup_pca_rslurm.R",    "PCA")
run_if_enabled(methods_meta$nmf$enabled,    "R/setup_nmf_rslurm.R",    "NMF")
run_if_enabled(methods_meta$cogaps$enabled, "R/setup_cogaps_rslurm.R", "CoGAPS")

# setup_network_rslurm.R covers both backends and already gates each
# internally on its own enabled flag, so it only needs one combined check.
network_enabled <- isTRUE(methods_meta$network$wto$enabled) ||
  isTRUE(methods_meta$network$wgcna$enabled)
run_if_enabled(network_enabled, "R/setup_network_rslurm.R", "network (wTO/WGCNA)")

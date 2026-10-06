# targets pipeline: latent factor method fitting.
#
# Run from the project root (RStudio Server on the cluster, project open):
#   targets::tar_visnetwork()                     # see the graph / what's outdated
#   targets::tar_make()                           # build everything outdated
#   targets::tar_make(names = starts_with("fit_pca_"))  # just some targets
#   targets::tar_make(as_job = TRUE)              # as an RStudio background job
#   targets::tar_read(fit_pca_EARLI)              # results (list of branches)
# For runs longer than the RStudio session: sbatch run_pipeline.sbatch
#
# Each fit_<method>_<dataset> target runs on SLURM workers launched by
# crew.cluster (one controller per method, config/pipeline.yml). The main
# process (this one) must stay alive for the run: workers connect back to it.
# See README.md's "Running the pipeline (targets)" section.

library(targets)
library(tarchetypes)
library(crew)
library(crew.cluster)

# Not a bare tar_source(): R/ also holds side-effecting driver scripts.
tar_source(c(
  "R/lib/grids.R",
  "R/lib/data_paths.R",
  "R/lib/matrices.R",
  "R/lib/metadata.R",
  "R/lib/method_registry.R",
  "R/db/connect.R",
  "R/lib/ingest",
  "R/ingest_jobs/fgsea_job.R",
  "R/ingest_jobs/wgcna_ora_job.R",
  "R/ingest_jobs/projectr_job.R",
  "app/R/metadata_helpers.R",   # generic_association_scan() (WGCNA gene significance)
  "R/targets"
))

# Sources every R/methods/*.R and defines its <name>_registry -- the fit
# targets reference these registries, so editing R/methods/pca.R re-runs
# only the pca fits.
METHOD_MANIFEST <- discover_method_registries()
# `methods` table rows (R/db/schema.sql), from each <name>_ingest contract.
METHOD_TABLE <- method_table(METHOD_MANIFEST)

pipeline <- read_pipeline_config("config/pipeline.yml")
data_root <- default_data_root()

tar_option_set(
  controller = build_controller_group(pipeline),
  # Shared filesystem: workers read their inputs from and write their
  # results to _targets/ directly instead of sending them through the main
  # process.
  storage = "worker",
  retrieval = "worker",
  # One failed fit shouldn't stop the other branches.
  error = "continue"
)

list(
  # MSigDB gene sets for enrichment, fetched once on the main process.
  tar_target(msigdb_pathways, fetch_msigdb_pathways(), deployment = "main"),
  lapply(pipeline$datasets, tar_dataset_fits,
         manifest = METHOD_MANIFEST, pipeline = pipeline, data_root = data_root),
  # Projections and gene-space agreement across datasets.
  tar_cross_dataset(pipeline, METHOD_MANIFEST, data_root)
)

# Standalone matrix-caching step, split out of R/create_ingest_slurm_bundle.R's
# --stage core (see that script's line ~118-127). Exists because
# cache_dataset_matrix() sources each dataset's preprocessing_script and
# reads the raw parquet trio via config$dataset$expression_path/
# sample_metadata_path/feature_metadata_path -- paths that, in this
# project's configs, are absolute paths on the machine that holds the raw
# HuggingFace data (typically your laptop), NOT on the cluster where the
# actual `*_results/` directories (and the slurm/apptainer tooling) live.
#
# Notably, this step needs NONE of that: cache_dataset_matrix(con,
# dataset_id, dataset_yaml, db_path) only reads dataset_yaml and writes to
# the DB + stability_artifacts/ -- it never looks at a results directory.
# So run this HERE (wherever the parquet resolves), then rsync/scp the two
# outputs (`results/stability_artifacts/<dataset_id>/matrix.rds` and the
# updated `results/stability.sqlite`) to the cluster before running
# R/create_ingest_slurm_bundle.R --stage core there. On the cluster,
# cache_dataset_matrix()'s staleness check (see R/lib/ingest/
# ingest_dataset.R) filters out any source path that doesn't
# file.exists() -- so with the parquet paths absent there, it can never
# detect staleness and will just short-circuit on your synced cache,
# never touching parquet again.
#
# Usage:
#   Rscript R/cache_dataset_matrices.R --datasets datasets.txt \
#     [--db results/stability.sqlite] [--recache-matrix [id1,id2,...]]
#
# `datasets.txt` is the SAME file used by --stage core: one results-dir
# path per line (e.g. "/scratch/.../GSE110487_T2_results") -- only its
# basename (with a trailing "_results" stripped) is used here to derive
# dataset_id and look up config/<dataset_id>_config.yml; the directory
# itself does not need to exist on this machine.

library(here); library(optparse); library(yaml); library(DBI)
source(here("R/lib/matrices.R"))
source(here("R/lib/ingest/db.R"))
source(here("R/lib/ingest/ingest_dataset.R"))   # cache_dataset_matrix()

option_list <- list(
  make_option("--datasets", type = "character", help = "file listing results dir paths, one per line"),
  make_option("--db", type = "character", default = "results/stability.sqlite"),
  make_option("--recache-matrix", type = "character", default = NULL,
              help = "bare flag = recache every dataset's matrix; or a comma-separated list of dataset ids")
)
opt <- parse_args(OptionParser(option_list = option_list))

result_dirs <- readLines(opt$datasets) |> trimws()
result_dirs <- result_dirs[nzchar(result_dirs)]
targets <- data.frame(
  results_dir = result_dirs,
  dataset_id  = sub("_results$", "", basename(result_dirs)),
  stringsAsFactors = FALSE
)
# Same case-insensitive config lookup as R/create_ingest_slurm_bundle.R --
# see its comment for why this can't just be file.exists(paste0(...)).
available_cfg <- list.files("config", pattern = "_config\\.yml$")
cfg_match <- match(tolower(paste0(targets$dataset_id, "_config.yml")), tolower(available_cfg))
targets$config_path <- ifelse(is.na(cfg_match), NA_character_, file.path("config", available_cfg[cfg_match]))
missing_cfg <- is.na(targets$config_path)
if (any(missing_cfg)) stop("No config found for: ", paste(targets$dataset_id[missing_cfg], collapse = ", "))

parse_flag_list <- function(x) {
  if (is.null(x)) return(FALSE)
  if (nzchar(x)) strsplit(x, ",")[[1]] else TRUE
}
recache_matrix <- parse_flag_list(opt$`recache-matrix`)

con <- open_stability_db(opt$db)
for (i in seq_len(nrow(targets))) {
  force_i <- isTRUE(recache_matrix) || (is.character(recache_matrix) && targets$dataset_id[i] %in% recache_matrix)
  ds_yaml <- yaml::read_yaml(targets$config_path[i])
  cache_dataset_matrix(con, targets$dataset_id[i], ds_yaml, opt$db, force = force_i)
  # Also cache sample/feature metadata -- needed so
  # create_ingest_slurm_bundle.R's --stage enrichment can build
  # sample_metadata_maps/ensembl_maps from the cache instead of these same
  # raw (laptop-only) paths once it's run on the cluster -- see
  # cache_dataset_metadata()'s header for the full story.
  cache_dataset_metadata(con, targets$dataset_id[i], ds_yaml, opt$db, force = force_i)
}
DBI::dbDisconnect(con)

message("\nDone -- cached matrices + sample/feature metadata for ", nrow(targets), " dataset(s) into ",
        dirname(opt$db), "/stability_artifacts/ and recorded in ", opt$db)
message("Sync both of those to the cluster before running create_ingest_slurm_bundle.R there.")

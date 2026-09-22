# Post-ingest_core step: compute WGCNA kME (intramodular connectivity /
# module membership) and Gene Significance (per-gene x sample-trait
# correlation) for every dataset -- see R/lib/ingest/ingest_dataset.R's
# compute_wgcna_kme()/compute_wgcna_gene_significance() headers for the
# full rationale.
#
# Split out as its own script, run AFTER ingest_core (unlike
# R/cache_dataset_matrices.R, which runs BEFORE it), because:
#   - compute_wgcna_kme() needs each WGCNA fit's `scores_file` (module
#     eigengenes), which only exists once ingest_core has ingested that
#     dataset's wgcna_grid results.
#   - Neither function is safe to call from inside
#     run_ingest_core_compute_job() (R/ingest_jobs/ingest_core_job.R):
#     that job runs inside ingest_core's own container, which -- like
#     compute_wgcna_sft()'s identical constraint -- does not have WGCNA
#     installed. This script, like
#     cache_dataset_matrices.R and R/ingest_results.R, is a plain
#     `Rscript` invocation meant to be submitted via `srun` (an
#     interactive scheduled allocation, with WGCNA + the raw
#     sample-metadata path both reachable), never folded into a
#     container-bound `sbatch`/array job.
#
# Usage:
#   Rscript R/compute_wgcna_diagnostics.R --datasets datasets.txt \
#     [--db results/stability.sqlite] [--recompute-kme [id1,id2,...]] \
#     [--recompute-gs [id1,id2,...]]
#
# `datasets.txt` is the SAME file used by cache_dataset_matrices.R/
# --stage core: one results-dir path per line -- only its basename (with
# a trailing "_results" stripped) is used here to derive dataset_id and
# look up config/<dataset_id>_config.yml.
#
# `--recompute-kme`/`--recompute-gs`: bare flag = recompute for every
# dataset even if already computed; or a comma-separated list of dataset
# ids. Otherwise additive/no-op per compute_wgcna_kme()'s (per-fit,
# additive) and compute_wgcna_gene_significance()'s (per-dataset,
# no-op-if-already-present) own semantics -- see their headers in
# R/lib/ingest/ingest_dataset.R.

# Deliberately does NOT `library(here)`/`library(optparse)` -- confirmed
# directly (2026-09-22) that the dedicated WGCNA container this script
# needs (r-wgcna, a minimal biocontainers-style image -- see this file's
# header) has WGCNA/DBI/RSQLite/yaml but NOT `here` or `optparse`
# (`library(here)` failing immediately was the actual, reproduced error).
# Hand-rolled arg parsing instead, same idiom as R/ingest_results.R's
# `--overwrite` handling. Also no longer sources R/lib/matrices.R --
# compute_wgcna_kme()/compute_wgcna_gene_significance() never call
# anything from it (that file exists for cache_dataset_matrix(), a
# different step entirely -- see this file's header), so it was dead
# weight that would ALSO have needed `arrow` (missing from this container
# too) the moment any of its functions were actually invoked.
library(yaml); library(DBI)
source("R/lib/ingest/db.R")
source("R/lib/ingest/ingest_dataset.R")   # compute_wgcna_kme()/compute_wgcna_gene_significance()

args <- commandArgs(trailingOnly = TRUE)
get_opt <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 1 && i < length(args)) args[[i + 1]] else default
}
get_bare_or_list <- function(flag) {
  i <- which(args == flag)
  if (length(i) != 1) return(FALSE)
  nxt <- if (i < length(args)) args[[i + 1]] else ""
  if (nzchar(nxt) && !startsWith(nxt, "--")) strsplit(nxt, ",")[[1]] else TRUE
}

datasets_file <- get_opt("--datasets")
if (is.null(datasets_file)) {
  stop("Usage: Rscript R/compute_wgcna_diagnostics.R --datasets datasets.txt ",
       "[--db results/stability.sqlite] [--recompute-kme [id1,id2,...]] [--recompute-gs [id1,id2,...]]")
}
db_path <- get_opt("--db", "results/stability.sqlite")
recompute_kme <- get_bare_or_list("--recompute-kme")
recompute_gs  <- get_bare_or_list("--recompute-gs")

result_dirs <- readLines(datasets_file) |> trimws()
result_dirs <- result_dirs[nzchar(result_dirs)]
targets <- data.frame(
  results_dir = result_dirs,
  dataset_id  = sub("_results$", "", basename(result_dirs)),
  stringsAsFactors = FALSE
)
# Same case-insensitive config lookup as R/create_ingest_slurm_bundle.R /
# R/cache_dataset_matrices.R.
available_cfg <- list.files("config", pattern = "_config\\.yml$")
cfg_match <- match(tolower(paste0(targets$dataset_id, "_config.yml")), tolower(available_cfg))
targets$config_path <- ifelse(is.na(cfg_match), NA_character_, file.path("config", available_cfg[cfg_match]))
missing_cfg <- is.na(targets$config_path)
if (any(missing_cfg)) stop("No config found for: ", paste(targets$dataset_id[missing_cfg], collapse = ", "))

con <- open_stability_db(db_path)
for (i in seq_len(nrow(targets))) {
  dataset_id <- targets$dataset_id[i]
  kme_force_i <- isTRUE(recompute_kme) || (is.character(recompute_kme) && dataset_id %in% recompute_kme)
  gs_force_i  <- isTRUE(recompute_gs)  || (is.character(recompute_gs)  && dataset_id %in% recompute_gs)
  ds_yaml <- yaml::read_yaml(targets$config_path[i])
  message("dataset: ", dataset_id)
  compute_wgcna_kme(con, dataset_id, ds_yaml, db_path, force = kme_force_i)
  compute_wgcna_gene_significance(con, dataset_id, ds_yaml, db_path, force = gs_force_i)
}
DBI::dbDisconnect(con)

message("\nDone -- computed WGCNA kME + gene significance for ", nrow(targets),
        " dataset(s) into ", db_path)

# Stage-2 ingest CLI: load ONE dataset's rslurm job outputs into the
# stability SQLite DB + external artifacts.
#
# Thin wrapper around R/lib/ingest/ingest_dataset.R::ingest_one_dataset()
# -- see that file for the actual per-family ingest logic (shared with
# R/ingest_jobs/ingest_core_job.R, which loops this over many datasets in
# a single slurm job; see R/create_ingest_slurm_bundle.R).
#
# Usage (CLI):
#   Rscript R/ingest_results.R <dataset_config.yml> <results_dir> <db_path> \
#     [--overwrite [jobname,...]] [--recache-matrix] [--recompute-redundancy] [--run-enrichment]
#
# `--overwrite` (bare = every family present in <results_dir>; or
# `--overwrite jobname1,jobname2` = just those) deletes and re-ingests --
# plain replacement, never row-level updating.
# `--recache-matrix` forces the cached input-matrix artifact to rebuild
# even if not auto-detected as stale (see cache_dataset_matrix()).
# `--recompute-redundancy` forces the pattern-redundancy diagnostic to
# recompute even if not auto-detected as stale (see run_all_redundancy()).
# `--run-enrichment` (legacy, on-demand path; the slurm pipeline's
# fgsea_grid/gprofiler_grid jobs are the preferred route now -- see
# R/create_ingest_slurm_bundle.R) runs gprofiler2 ORA/GSEA for every ok
# fit's every factor across the WHOLE db, skipping anything already cached.
#
# Usage (interactive): set `ingest_config_path`, `ingest_results_dir`,
# `ingest_db_path` (and optionally `ingest_overwrite`, `ingest_recache_matrix`,
# `ingest_recompute_redundancy`, `ingest_run_enrichment`) then source this file.

library(here)
library(yaml)
source(here("R/lib/matrices.R"))
source(here("R/lib/ingest/db.R"))
source(here("R/lib/ingest/similarity.R"))
source(here("R/lib/ingest/extract.R"))
source(here("R/lib/ingest/pairs.R"))
source(here("R/lib/ingest/redundancy.R"))
source(here("R/lib/ingest/driver.R"))
source(here("R/lib/ingest/ingest_dataset.R"))

## ---- argument handling ------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
positional <- args[!grepl("^--", args)]
if (!exists("ingest_config_path")) {
  if (length(positional) < 3) {
    stop("Usage: Rscript R/ingest_results.R <dataset_config.yml> <results_dir> <db_path> [--overwrite [jobname,...]]")
  }
  ingest_config_path <- positional[[1]]
  ingest_results_dir <- positional[[2]]
  ingest_db_path     <- positional[[3]]
}
if (!exists("ingest_overwrite")) {
  ow_flag <- which(args == "--overwrite")
  ingest_overwrite <- if (length(ow_flag) == 1) {
    nxt <- if (ow_flag < length(args)) args[[ow_flag + 1]] else ""
    if (nzchar(nxt) && !grepl("^--", nxt) && !(nxt %in% positional[1:3])) {
      strsplit(nxt, ",")[[1]]
    } else {
      TRUE
    }
  } else {
    FALSE
  }
}
if (!exists("ingest_recache_matrix")) {
  ingest_recache_matrix <- "--recache-matrix" %in% args
}
if (!exists("ingest_recompute_redundancy")) {
  ingest_recompute_redundancy <- "--recompute-redundancy" %in% args
}
if (!exists("ingest_run_enrichment")) {
  ingest_run_enrichment <- "--run-enrichment" %in% args
}

con <- open_stability_db(ingest_db_path)
# cache_dataset_matrix() sources preprocessing_script by path and reads
# the raw parquet trio -- only ever call it from an ordinary (non-
# containerized) R session like this CLI, never from inside a slurm job's
# container (see R/ingest_jobs/ingest_core_job.R's header for why).
dataset_yaml_for_cache <- yaml::read_yaml(ingest_config_path)
cache_dataset_matrix(con, dataset_yaml_for_cache$dataset$id, dataset_yaml_for_cache,
                      ingest_db_path, force = ingest_recache_matrix)
ingest_one_dataset(con, ingest_config_path, ingest_results_dir, ingest_db_path,
                    overwrite = ingest_overwrite, recompute_redundancy = ingest_recompute_redundancy)

dataset_id <- yaml::read_yaml(ingest_config_path)$dataset$id
message("\n===== ingest report: ", dataset_id, " -> ", ingest_db_path, " =====")
counts <- DBI::dbGetQuery(con, "
  SELECT 'fits' AS tbl, COUNT(*) AS n FROM fits WHERE dataset_id = :d
  UNION ALL SELECT 'factors', COUNT(*) FROM factors f JOIN fits ft ON ft.fit_id = f.fit_id WHERE ft.dataset_id = :d
  UNION ALL SELECT 'factor_pairs', COUNT(*) FROM factor_pairs fp JOIN fits ft ON ft.fit_id = fp.fit_a WHERE ft.dataset_id = :d
  UNION ALL SELECT 'maskcv_results', COUNT(*) FROM maskcv_results WHERE dataset_id = :d
  UNION ALL SELECT 'wgcna_fit_pairs', COUNT(*) FROM wgcna_fit_pairs wp JOIN fits ft ON ft.fit_id = wp.fit_a WHERE ft.dataset_id = :d
  UNION ALL SELECT 'fit_redundancy', COUNT(*) FROM fit_redundancy fr JOIN fits ft ON ft.fit_id = fr.fit_id WHERE ft.dataset_id = :d",
  params = list(d = dataset_id))
for (r in seq_len(nrow(counts))) message(sprintf("  %-22s %d rows", counts$tbl[r], counts$n[r]))

if (isTRUE(ingest_run_enrichment)) {
  message("\n===== --run-enrichment: running gprofiler2 ORA/GSEA for the WHOLE db =====")
  source(here("R/lib/ingest/enrichment.R"))
  run_all_enrichment(con, ingest_db_path)
}

DBI::dbDisconnect(con)

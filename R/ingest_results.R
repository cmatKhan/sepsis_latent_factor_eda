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
#     [--overwrite [jobname,...]] [--recache-matrix] [--recompute-redundancy] \
#     [--recompute-sft] [--recompute-kme] [--recompute-gs]
#
# `--overwrite` (bare = every family present in <results_dir>; or
# `--overwrite jobname1,jobname2` = just those) deletes and re-ingests --
# plain replacement, never row-level updating.
# `--recache-matrix` forces the cached input-matrix artifact to rebuild
# even if not auto-detected as stale (see cache_dataset_matrix()).
# `--recompute-redundancy` forces the pattern-redundancy diagnostic to
# recompute even if not auto-detected as stale (see run_all_redundancy()).
# `--recompute-sft` forces WGCNA::pickSoftThreshold()'s scale-free-topology
# fit diagnostic to recompute even if the dataset's configured power grid
# hasn't changed (see compute_wgcna_sft()). No-ops for datasets with no
# methods.network.wgcna config.
# `--recompute-kme` forces WGCNA::signedKME() module-membership recompute
# for every wgcna fit of this dataset, even for fits that already have
# wgcna_kme rows (see compute_wgcna_kme()) -- otherwise only NEW fits (zero
# existing rows) get computed, additive across fits.
# `--recompute-gs` forces WGCNA Gene Significance (per-gene x sample-trait
# correlation) to recompute even if this dataset already has
# wgcna_gene_significance rows (see compute_wgcna_gene_significance()) --
# otherwise a no-op once any rows exist, replace-on-recompute like SFT.
#
# (`--run-enrichment`/R/lib/ingest/enrichment.R's run_all_enrichment(), a
# legacy whole-DB gprofiler2 ORA/GSEA pass, was retired 2026-09-19 alongside
# the slurm pipeline's gprofiler_grid job family -- both had the same
# gprofiler2-API-rate-limit scaling problem; see R/ingest_jobs/fgsea_job.R's
# header for the local fora()/fgsea()-based replacement, now folded into
# fgsea_grid. The app's own per-factor, on-demand queries in app/app.R were
# later retired ENTIRELY (not just switched off gprofiler2) -- the app has
# no live-compute enrichment path at all now, it only ever reads what this
# ingest pipeline has already computed, so gprofiler2 is no longer called
# anywhere in this project.)
#
# Usage (interactive): set `ingest_config_path`, `ingest_results_dir`,
# `ingest_db_path` (and optionally `ingest_overwrite`, `ingest_recache_matrix`,
# `ingest_recompute_redundancy`) then source this file.

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
if (!exists("ingest_recompute_sft")) {
  ingest_recompute_sft <- "--recompute-sft" %in% args
}
if (!exists("ingest_recompute_kme")) {
  ingest_recompute_kme <- "--recompute-kme" %in% args
}
if (!exists("ingest_recompute_gs")) {
  ingest_recompute_gs <- "--recompute-gs" %in% args
}

con <- open_stability_db(ingest_db_path)
# cache_dataset_matrix() sources preprocessing_script by path and reads
# the raw parquet trio -- only ever call it from an ordinary (non-
# containerized) R session like this CLI, never from inside a slurm job's
# container (see R/ingest_jobs/ingest_core_job.R's header for why).
dataset_yaml_for_cache <- yaml::read_yaml(ingest_config_path)
cache_dataset_matrix(con, dataset_yaml_for_cache$dataset$id, dataset_yaml_for_cache,
                      ingest_db_path, force = ingest_recache_matrix)
# compute_wgcna_sft() only needs the matrix just cached above + WGCNA
# (installed here, same non-containerized-session constraint as
# cache_dataset_matrix() -- see that function's header) -- never callable
# from inside ingest_core's container.
compute_wgcna_sft(con, dataset_yaml_for_cache$dataset$id, dataset_yaml_for_cache,
                   ingest_db_path, force = ingest_recompute_sft)
# kME needs each wgcna fit's scores_file (module eigengenes) to already be
# ingested, so it runs AFTER ingest_one_dataset() below, not alongside SFT.
# Gene significance only needs the cached matrix + registered metadata
# (both already in place by this point), so its ordering doesn't matter --
# kept next to kME for readability.
ingest_one_dataset(con, ingest_config_path, ingest_results_dir, ingest_db_path,
                    overwrite = ingest_overwrite, recompute_redundancy = ingest_recompute_redundancy)
compute_wgcna_kme(con, dataset_yaml_for_cache$dataset$id, dataset_yaml_for_cache,
                   ingest_db_path, force = ingest_recompute_kme)
compute_wgcna_gene_significance(con, dataset_yaml_for_cache$dataset$id, dataset_yaml_for_cache,
                                 ingest_db_path, force = ingest_recompute_gs)

dataset_id <- yaml::read_yaml(ingest_config_path)$dataset$id
message("\n===== ingest report: ", dataset_id, " -> ", ingest_db_path, " =====")
counts <- DBI::dbGetQuery(con, "
  SELECT 'fits' AS tbl, COUNT(*) AS n FROM fits WHERE dataset_id = :d
  UNION ALL SELECT 'factors', COUNT(*) FROM factors f JOIN fits ft ON ft.fit_id = f.fit_id WHERE ft.dataset_id = :d
  UNION ALL SELECT 'factor_pairs', COUNT(*) FROM factor_pairs fp JOIN fits ft ON ft.fit_id = fp.fit_a WHERE ft.dataset_id = :d
  UNION ALL SELECT 'wgcna_sft', COUNT(*) FROM wgcna_sft WHERE dataset_id = :d
  UNION ALL SELECT 'wgcna_kme', COUNT(*) FROM wgcna_kme wk JOIN fits ft ON ft.fit_id = wk.fit_id WHERE ft.dataset_id = :d
  UNION ALL SELECT 'wgcna_gene_significance', COUNT(*) FROM wgcna_gene_significance WHERE dataset_id = :d
  UNION ALL SELECT 'wgcna_fit_pairs', COUNT(*) FROM wgcna_fit_pairs wp JOIN fits ft ON ft.fit_id = wp.fit_a WHERE ft.dataset_id = :d
  UNION ALL SELECT 'fit_redundancy', COUNT(*) FROM fit_redundancy fr JOIN fits ft ON ft.fit_id = fr.fit_id WHERE ft.dataset_id = :d",
  params = list(d = dataset_id))
for (r in seq_len(nrow(counts))) message(sprintf("  %-22s %d rows", counts$tbl[r], counts$n[r]))

DBI::dbDisconnect(con)

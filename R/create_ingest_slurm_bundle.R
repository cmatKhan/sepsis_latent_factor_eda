# Two-phase slurm ingest pipeline driver.
#
# Phase 1 (--stage core): a SINGLE (non-array) "ingest_core" job that
# loops every dataset's job-family results (from a user-supplied list of
# results directories, one per dataset -- see --datasets) into the shared
# stability DB, one writer at a time (avoiding concurrent-write hazards).
#
# Phase 2 (--stage enrichment): run ONLY after Phase 1's job has finished
# and the DB reflects its results (representative-fit selection needs
# fits.mse, which doesn't exist until ingest_core has run). Queries the DB
# for representative fits per (dataset, method) and stages three more job
# families: fgsea_grid (array, parallel), gprofiler_grid (single, serial
# -- API rate-limited), and projectr_within_grid/projectr_cross_grid
# (array, "reduced all-pairs" -- see R/lib/ingest/projectr_pairs.R).
#
# Usage:
#   Rscript R/create_ingest_slurm_bundle.R --datasets datasets.txt --stage core
#   # ... wait for ingest_core to finish ...
#   Rscript R/create_ingest_slurm_bundle.R --datasets datasets.txt --stage enrichment
#   # ... wait for fgsea_grid/gprofiler_grid/projectr_*_grid ...
#   Rscript R/ingest_enrichment_results.R --bundle-dir slurm_bundles/ingest --db results/stability.sqlite
#   Rscript R/ingest_projectr_results.R   --bundle-dir slurm_bundles/ingest --db results/stability.sqlite
#
# `datasets.txt` lists results directory paths, one per line, e.g.:
#   results/EARLI_results
#   results/GSE110487_results
# dataset_id is derived as basename(path) with a trailing "_results"
# stripped; config/<dataset_id>_config.yml must exist.

library(here); library(optparse); library(yaml); library(DBI)
source(here("R/lib/matrices.R"))
source(here("R/lib/ingest/db.R"))
source(here("R/lib/ingest/similarity.R"))
source(here("R/lib/ingest/extract.R"))
source(here("R/lib/ingest/pairs.R"))
source(here("R/lib/ingest/redundancy.R"))
source(here("R/lib/ingest/driver.R"))
source(here("R/lib/ingest/ingest_dataset.R"))
source(here("R/lib/ingest/symbol_mapping.R"))
source(here("R/lib/ingest/projectr_pairs.R"))
source(here("R/lib/submit.R"))
source(here("R/ingest_jobs/ingest_core_job.R"))
source(here("R/ingest_jobs/fgsea_job.R"))
source(here("R/ingest_jobs/gprofiler_job.R"))
source(here("R/ingest_jobs/projectr_job.R"))

`%||%` <- function(a, b) if (is.null(a)) b else a

# Every function defined by the sourced lib/job files above -- baked into
# EVERY submit_job_family() call below as `global_objects`, since none of
# these compute-node jobs can `source()` project files by path (only
# slurm_bundles/ingest/ is bind-mounted inside the container -- see
# R/ingest_jobs/ingest_core_job.R's header). Over-including per job is
# harmless (plain function objects are cheap to serialize); hand-
# maintaining a separate minimal list per job family is not worth the
# fragility. Captured here, before `opt`/`targets`/etc. exist, so it's
# exactly "everything sourced above" and nothing else.
FRAMEWORK_FUNCS <- Filter(function(n) is.function(get(n, envir = .GlobalEnv)), ls(envir = .GlobalEnv))

# Bind-mounted (read-write, same absolute path in and out of the
# container) + set as --pwd for every job family below, so
# config/<id>_config.yml, results/<id>_results/, results/stability.sqlite
# (+ stability_artifacts/), and slurm_bundles/<id>/_rslurm_<jobname>/
# params.RDS all resolve via their ordinary project-relative paths exactly
# as they do outside the container -- no path translation, no baking
# config/params content into global_objects. Mounting the whole project
# root (rather than just `results/`) is the simplest way to cover all
# three trees ingest_one_dataset() touches in one bind; narrow this to
# just `here("results")` (+ read-only `here("config")`/`here("slurm_bundles")`
# binds) later if you'd rather not expose the rest of the project tree.
PROJECT_ROOT <- here::here()

option_list <- list(
  make_option("--datasets", type = "character", help = "file listing results dir paths, one per line"),
  make_option("--db", type = "character", default = "results/stability.sqlite"),
  make_option("--slurm-config", type = "character", default = "config/ingest_slurm_config.yml"),
  make_option("--dataset-families", type = "character", default = "config/dataset_families.yml"),
  make_option("--stage", type = "character", default = "core", help = "'core' or 'enrichment'"),
  make_option("--output", type = "character", default = here("slurm_bundles/ingest")),
  make_option("--recache-matrix", type = "character", default = NULL,
              help = "bare flag = recache every dataset's matrix; or a comma-separated list of dataset ids"),
  make_option("--recompute-redundancy", type = "character", default = NULL,
              help = "bare flag = recompute for every dataset; or a comma-separated list of dataset ids")
)
opt <- parse_args(OptionParser(option_list = option_list))
slurm_cfg <- yaml::read_yaml(opt$`slurm-config`)

result_dirs <- readLines(opt$datasets) |> trimws()
result_dirs <- result_dirs[nzchar(result_dirs)]
targets <- data.frame(
  results_dir = result_dirs,
  dataset_id  = sub("_results$", "", basename(result_dirs)),
  stringsAsFactors = FALSE
)
# Case-insensitive match against actual config/ filenames: config file
# names lowercase their day/timepoint suffix (e.g. "GSE110487_t2_config.yml",
# "CORTICUS_pre_config.yml") while dataset_id here -- derived from the
# results directory basename -- keeps whatever casing that directory (and
# the config's own internal `dataset: id:` field, which DOES match
# dataset_id exactly) uses, e.g. "GSE110487_T2", "CORTICUS_PRE". Matching
# file.exists() directly against paste0(dataset_id, "_config.yml") is
# case-SENSITIVE on Linux and fails for every such dataset.
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
recompute_redundancy <- parse_flag_list(opt$`recompute-redundancy`)

if (opt$stage == "core") {

  # cache_dataset_matrix() sources each dataset's preprocessing_script and
  # reads the raw parquet trio -- ONLY safe to run here, on the login node,
  # where project-relative paths actually resolve (see
  # R/ingest_jobs/ingest_core_job.R's header). Never inside the slurm job.
  con_cache <- open_stability_db(opt$db)
  for (i in seq_len(nrow(targets))) {
    force_i <- isTRUE(recache_matrix) || (is.character(recache_matrix) && targets$dataset_id[i] %in% recache_matrix)
    cache_dataset_matrix(con_cache, targets$dataset_id[i], yaml::read_yaml(targets$config_path[i]), opt$db, force = force_i)
  }
  DBI::dbDisconnect(con_cache)

  # ingest_one_dataset() reads targets$config_path/results_dir (project-
  # relative, e.g. "config/EARLI_config.yml", "results/EARLI_results/")
  # and slurm_bundles/<dataset_id>/_rslurm_<jobname>/params.RDS for every
  # dataset being ingested -- resolved via the PROJECT_ROOT bind +
  # pwd_override above, not via global_objects.
  assign("targets", targets, envir = .GlobalEnv)
  assign("db_path", opt$db, envir = .GlobalEnv)
  assign("recompute_redundancy", recompute_redundancy, envir = .GlobalEnv)
  submit_job_family(
    f = run_ingest_core_job, jobs_df = NULL, jobname = "ingest_core",
    global_objects = c(FRAMEWORK_FUNCS, "targets", "db_path", "recompute_redundancy"),
    # projectR added here (not just the projectr container) -- run_all_pattern_drivers()
    # (R/lib/ingest/driver.R) now runs inside ingest_core itself, calling
    # projectR::projectionDriveR() for the differential-features pass.
    pkgs = c("DBI", "RSQLite", "arrow", "yaml", "CoGAPS", "clue", "matrixStats", "projectR"),
    cluster_cfg = slurm_cfg$ingest_core, output_dir = opt$output,
    extra_binds = PROJECT_ROOT, pwd_override = PROJECT_ROOT
  )
  message("Staged ingest_core -- run this FIRST, wait for it to finish, then re-run with --stage enrichment")

} else if (opt$stage == "enrichment") {

  con <- open_stability_db(opt$db)
  all_dataset_ids <- DBI::dbGetQuery(con, "SELECT DISTINCT dataset_id FROM fits")$dataset_id
  if (length(all_dataset_ids) == 0) stop("No fits in the DB yet -- run --stage core first")

  # symbol maps + dataset matrix file lookups, built once on the login node
  symbol_maps <- setNames(lapply(all_dataset_ids, function(ds) {
    cfg_path <- file.path("config", paste0(ds, "_config.yml"))
    if (!file.exists(cfg_path)) return(NULL)
    build_symbol_map(yaml::read_yaml(cfg_path))
  }), all_dataset_ids)

  target_matrix_path_for <- function(dataset_id) {
    f <- DBI::dbGetQuery(con, "SELECT matrix_file FROM datasets WHERE dataset_id = ?",
                          params = list(dataset_id))$matrix_file
    if (length(f) != 1 || is.na(f)) stop("No cached matrix for dataset '", dataset_id,
                                          "' -- run --stage core first (cache_dataset_matrix() populates this).")
    resolve_artifact(f, opt$db)
  }

  # ---- fgsea + gprofiler: representative fits only ----
  rep_rows <- do.call(rbind, lapply(all_dataset_ids, function(ds) {
    do.call(rbind, lapply(c("pca", "nmf", "cogaps", "spca", "ica"), function(m) {
      fids <- representative_fit_ids(con, ds, m)
      if (length(fids) == 0) return(NULL)
      f <- DBI::dbGetQuery(con, sprintf("SELECT fit_id, loadings_file FROM fits WHERE fit_id IN (%s)",
                                         paste(fids, collapse = ",")))
      data.frame(dataset_id = ds, method = m, fit_id = f$fit_id,
                 loadings_file = resolve_artifact(f$loadings_file, opt$db), stringsAsFactors = FALSE)
    }))
  }))

  if (!is.null(rep_rows) && nrow(rep_rows) > 0) {
    if (!requireNamespace("msigdbr", quietly = TRUE)) stop("Package 'msigdbr' is required")
    msig <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
    pathways <- split(msig$gene_symbol, msig$gs_name)
    assign("pathways", pathways, envir = .GlobalEnv)

    # per-row cogaps marker genes (already computed by run_all_redundancy() during --stage core)
    rep_rows$cogaps_marker_genes <- vector("list", nrow(rep_rows))
    cogaps_idx <- which(rep_rows$method == "cogaps")
    for (i in cogaps_idx) {
      mk <- DBI::dbGetQuery(con, "SELECT factor_index, gene FROM pattern_markers WHERE fit_id = ?",
                             params = list(rep_rows$fit_id[i]))
      rep_rows$cogaps_marker_genes[[i]] <- if (nrow(mk) > 0) split(mk$gene, mk$factor_index) else NULL
    }
    rep_rows$symbol_map <- lapply(rep_rows$dataset_id, function(ds) symbol_maps[[ds]])

    submit_job_family(
      f = run_fgsea_job,
      jobs_df = rep_rows[, c("dataset_id", "method", "fit_id", "loadings_file", "symbol_map", "cogaps_marker_genes")],
      jobname = "fgsea_grid", global_objects = c(FRAMEWORK_FUNCS, "pathways"),
      pkgs = c("CoGAPS", "BiocParallel", "arrow"),
      cluster_cfg = slurm_cfg$fgsea, output_dir = opt$output,
      extra_binds = PROJECT_ROOT, pwd_override = PROJECT_ROOT
    )

    gprofiler_targets <- do.call(rbind, lapply(seq_len(nrow(rep_rows)), function(i) {
      L <- as.matrix(readRDS(rep_rows$loadings_file[i]))
      dirs <- if (rep_rows$method[i] %in% c("pca", "ica", "spca")) c("pos", "neg") else "pos"
      expand.grid(dataset_id = rep_rows$dataset_id[i], method = rep_rows$method[i],
                  fit_id = rep_rows$fit_id[i], loadings_file = rep_rows$loadings_file[i],
                  factor_index = seq_len(ncol(L)), direction = dirs, stringsAsFactors = FALSE)
    }))
    assign("gprofiler_targets", gprofiler_targets, envir = .GlobalEnv)
    assign("symbol_maps", symbol_maps, envir = .GlobalEnv)
    submit_job_family(
      f = run_gprofiler_job, jobs_df = NULL, jobname = "gprofiler_grid",
      global_objects = c(FRAMEWORK_FUNCS, "gprofiler_targets", "symbol_maps"), pkgs = c("gprofiler2", "arrow"),
      cluster_cfg = slurm_cfg$gprofiler, output_dir = opt$output,
      extra_binds = PROJECT_ROOT, pwd_override = PROJECT_ROOT
    )
  } else {
    message("No representative fits found for fgsea/gprofiler -- skipping both")
  }

  # ---- projectr: reduced all-pairs, two explicitly-named families ----
  families <- yaml::read_yaml(opt$`dataset-families`)
  pairs <- build_projectr_pairs(con, all_dataset_ids, families)
  message("projectr: ", nrow(pairs), " total (source fit x target dataset) pairs -- ",
          sum(pairs$projection_type == "within_dataset"), " within-dataset, ",
          sum(pairs$projection_type == "cross_dataset"), " cross-dataset")

  for (ptype in c("within_dataset", "cross_dataset")) {
    sub <- pairs[pairs$projection_type == ptype, ]
    if (nrow(sub) == 0) next
    sub$loadings_file <- vapply(sub$source_fit_id, function(fid) {
      resolve_artifact(DBI::dbGetQuery(con, "SELECT loadings_file FROM fits WHERE fit_id = ?",
                                        params = list(fid))$loadings_file, opt$db)
    }, character(1))
    sub$target_matrix_file <- vapply(sub$target_dataset_id, target_matrix_path_for, character(1))
    sub$source_mat_file <- vapply(sub$source_dataset_id, function(ds) {
      f <- DBI::dbGetQuery(con, "SELECT matrix_file FROM datasets WHERE dataset_id = ?",
                            params = list(ds))$matrix_file
      if (length(f) == 1 && !is.na(f)) resolve_artifact(f, opt$db) else NA_character_
    }, character(1))
    sub$source_scores_file <- vapply(sub$source_fit_id, function(fid) {
      f <- DBI::dbGetQuery(con, "SELECT scores_file FROM fits WHERE fit_id = ?", params = list(fid))$scores_file
      if (length(f) == 1 && !is.na(f)) resolve_artifact(f, opt$db) else NA_character_
    }, character(1))

    assign("symbol_maps", symbol_maps, envir = .GlobalEnv)
    submit_job_family(
      f = run_projectr_job, jobs_df = sub,
      jobname = paste0("projectr_", sub("_dataset$", "", ptype), "_grid"),
      global_objects = c(FRAMEWORK_FUNCS, "symbol_maps"), pkgs = c("projectR", "arrow"),
      cluster_cfg = slurm_cfg$projectr, output_dir = opt$output,
      extra_binds = PROJECT_ROOT, pwd_override = PROJECT_ROOT
    )
  }
  DBI::dbDisconnect(con)
}

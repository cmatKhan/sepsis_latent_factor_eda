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
source(here("R/ingest_jobs/driver_job.R"))

`%||%` <- function(a, b) if (is.null(a)) b else a

# EVERY object (functions AND plain data constants, e.g. extract.R's
# PARAM_GRID_METHODS) defined by the sourced lib/job files above -- baked
# into EVERY submit_job_family() call below as `global_objects`, since none
# of these compute-node jobs can `source()` project files by path (only
# slurm_bundles/ingest/ is bind-mounted inside the container -- see
# R/ingest_jobs/ingest_core_job.R's header). Over-including per job is
# harmless (these are all cheap to serialize); hand-maintaining a separate
# minimal list per job family is not worth the fragility. Captured here,
# before `opt`/`targets`/etc. exist, so it's exactly "everything sourced
# above" and nothing else -- deliberately NOT filtered down to
# is.function() only (an earlier version did that and silently dropped
# PARAM_GRID_METHODS, breaking classify_jobname() -> ingest_one_dataset()
# inside run_ingest_core_job() with "object 'PARAM_GRID_METHODS' not
# found" -- any future non-function constant added to these lib files
# would hit the same gap if this were re-narrowed).
FRAMEWORK_FUNCS <- ls(envir = .GlobalEnv)

# Bind-mounted (read-write, same absolute path in and out of the
# container) via extra_binds on every submit_job_family() call below, so
# config/<id>_config.yml, results/<id>_results/, results/stability.sqlite
# (+ stability_artifacts/), and slurm_bundles/<id>/_rslurm_<jobname>/
# params.RDS are all REACHABLE inside the container at this same absolute
# path. NOT used as --pwd (never pass this as submit_job_family()'s
# pwd_override) -- that collides with rslurm's own slurm_run.R, which loads
# f.RDS/params.RDS/add_objects.RData via relative paths from ITS bundle
# directory before calling the job function (see run_ingest_core_job()'s
# header for the full story). Consequently, every path a containerized job
# actually touches under this tree must be built as an ABSOLUTE path baked
# in via global_objects (e.g. ingest_one_dataset()'s project_root
# argument, or normalizePath()-ed config_path/db_path below) -- cwd inside
# the container can't be relied on to resolve anything project-relative.
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
# Absolute from here on: run_ingest_core_job() (inside the ingest_core
# container) reads this back -- see the --stage core block below for why
# it must NOT be project-relative there.
targets$config_path <- normalizePath(targets$config_path, mustWork = TRUE)

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

  # ingest_one_dataset() reads targets$config_path/results_dir and
  # slurm_bundles/<dataset_id>/_rslurm_<jobname>/params.RDS for every
  # dataset being ingested. These must be ABSOLUTE, not project-relative:
  # rslurm's own slurm_run.R (see config/rslurm_templates/
  # slurm_run_single_R.txt) loads f.RDS/params.RDS/add_objects.RData via
  # relative paths from its own _rslurm_ingest_core/ directory BEFORE
  # calling run_ingest_core_job() -- so the container's cwd at that point
  # must stay the default `$RSLURM_BUNDLE_DIR/_rslurm_ingest_core` (never
  # override it via pwd_override, which breaks that load with "cannot open
  # file 'slurm_run.R'"). extra_binds (below) is still what makes
  # PROJECT_ROOT reachable/writable inside the container at its ordinary
  # absolute path; run_ingest_core_job() just needs paths that already
  # point there regardless of its own cwd, which is what
  # normalizePath()-ing config_path (above) and db_path (below) achieves
  # -- results_dir is already absolute, since datasets.txt itself lists
  # absolute cluster paths.
  assign("targets", targets, envir = .GlobalEnv)
  assign("db_path", normalizePath(opt$db, mustWork = FALSE), envir = .GlobalEnv)
  assign("recompute_redundancy", recompute_redundancy, envir = .GlobalEnv)
  # PROJECT_ROOT itself was defined AFTER FRAMEWORK_FUNCS was captured
  # (above), so it isn't already among those baked-in globals -- needed by
  # run_ingest_core_job() to pass ingest_one_dataset()'s project_root arg.
  assign("PROJECT_ROOT", PROJECT_ROOT, envir = .GlobalEnv)
  # ingest_one_dataset() also list.dirs()/reads results_*.RDS straight out
  # of each targets$results_dir -- these commonly live OUTSIDE PROJECT_ROOT
  # entirely (e.g. as sibling directories of the project checkout, not
  # nested under it: /scratch/.../GSE110487_T2_results next to /scratch/
  # .../sepsis_latent_factor_eda/), so binding PROJECT_ROOT alone leaves
  # them invisible inside the container -- list.dirs(results_dir) silently
  # returns character(0) there even though the same path is populated on
  # the login node, surfacing as "no job family subdirectories ... --
  # skipping" for every single dataset. Bind each results_dir's PARENT
  # (deduplicated -- typically just one shared parent covering everything)
  # in addition to PROJECT_ROOT.
  results_parents <- unique(dirname(targets$results_dir))
  submit_job_family(
    f = run_ingest_core_job, jobs_df = NULL, jobname = "ingest_core",
    global_objects = c(FRAMEWORK_FUNCS, "targets", "db_path", "recompute_redundancy", "PROJECT_ROOT"),
    # NOTE: no projectR here -- ingest_core's image doesn't have it (a
    # different image than driver_grid needs, see run_ingest_core_job()'s
    # header); run_pattern_drivers = FALSE there, staged as its own
    # driver_grid job family during --stage enrichment instead.
    pkgs = c("DBI", "RSQLite", "arrow", "yaml", "CoGAPS", "clue", "matrixStats", "mclust"),
    cluster_cfg = slurm_cfg$ingest_core, output_dir = opt$output,
    extra_binds = unique(c(PROJECT_ROOT, results_parents))
  )
  message("Staged ingest_core -- run this FIRST, wait for it to finish, then re-run with --stage enrichment")

} else if (opt$stage == "enrichment") {

  con <- open_stability_db(opt$db)
  all_dataset_ids <- DBI::dbGetQuery(con, "SELECT DISTINCT dataset_id FROM fits")$dataset_id
  if (length(all_dataset_ids) == 0) stop("No fits in the DB yet -- run --stage core first")

  # ---- driver_grid: projectR::projectionDriveR() pattern-driver pass ----
  # Split out of ingest_core (see R/ingest_jobs/driver_job.R's header) --
  # needs the projectr image, not ingest_core's. One single (non-array) job,
  # looping every dataset sequentially against the SAME db (one writer),
  # same rationale as ingest_core itself.
  assign("all_dataset_ids", all_dataset_ids, envir = .GlobalEnv)
  assign("db_path", normalizePath(opt$db, mustWork = FALSE), envir = .GlobalEnv)
  submit_job_family(
    f = run_driver_job, jobs_df = NULL, jobname = "driver_grid",
    global_objects = c(FRAMEWORK_FUNCS, "all_dataset_ids", "db_path"),
    pkgs = c("DBI", "RSQLite", "arrow", "projectR"),
    cluster_cfg = slurm_cfg$driver, output_dir = opt$output,
    extra_binds = PROJECT_ROOT
  )

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
      extra_binds = PROJECT_ROOT
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
      extra_binds = PROJECT_ROOT
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
      extra_binds = PROJECT_ROOT
    )
  }
  DBI::dbDisconnect(con)
}

# Every _rslurm_<jobname>/submit.sh generated above expects
# $RSLURM_BUNDLE_DIR (and, per slurm_run.R, $RSLURM_OUTPUT_DIR) to already
# be exported in the submitting shell -- see build_apptainer_rscript_path()'s
# header and config/rslurm_templates/submit_sh.txt's `set -euo pipefail`.
# Running a job family's submit.sh directly (e.g. a bare `sbatch submit.sh`)
# fails with "RSLURM_BUNDLE_DIR: unbound variable" for exactly that reason.
# write_submit_all_script() (safe to call repeatedly -- see its own header)
# writes/refreshes slurm_bundles/ingest/submit_all_ingest.sh, which exports
# both variables and discovers every _rslurm_* directory dynamically; use
# IT to submit, not `sbatch` directly against a job family's own submit.sh.
write_submit_all_script("ingest", opt$output)
message("Wrote/refreshed ", file.path(opt$output, "submit_all_ingest.sh"),
        " -- use it to actually submit the job(s) staged above, e.g.:\n  ",
        "cd ", opt$output, " && ./submit_all_ingest.sh .")

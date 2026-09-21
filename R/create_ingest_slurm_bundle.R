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
# for representative fits per (dataset, method) and stages two more job
# families: fgsea_grid (array, parallel -- ALSO covers the local
# fora()/fgsea()-based ORA/GSEA pass across GO/KEGG/Reactome/WikiPathways/
# Hallmark that used to be a separate, API-rate-limited gprofiler_grid job
# family; see R/ingest_jobs/fgsea_job.R's header and this script's
# `pathways_by_source` comment for why gprofiler_grid was retired
# 2026-09-19), and projectr_within_grid/projectr_cross_grid (array,
# "reduced all-pairs" -- see R/lib/ingest/projectr_pairs.R).
#
# Usage:
#   Rscript R/create_ingest_slurm_bundle.R --datasets datasets.txt --stage core
#   # ... wait for ingest_core to finish ...
#   Rscript R/create_ingest_slurm_bundle.R --datasets datasets.txt -nd -stage enrichment
#   # ... wait for fgsea_grid/projectr_*_grid ...
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
source(here("R/ingest_jobs/wgcna_ora_job.R"))
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
# would hit the same gap if this were re-narrowed). ALSO deliberately
# `all.names = TRUE` -- plain ls() silently excludes dot-prefixed names,
# which is a SEPARATE exclusion from the is.function() one above and bit
# us the same way: R/lib/ingest/symbol_mapping.R's private `.remap_ids()`
# helper (called internally by remap_to_ensembl()/remap_to_symbol(), both
# ordinary non-dot names that WERE captured) was silently missing from
# every job's add_objects.RData, failing at runtime -- not submission
# time -- with "could not find function '.remap_ids'" the moment
# fgsea_grid/projectr_*_grid actually called it. Any future
# dot-prefixed ("private") helper added to these lib files needs this too.
FRAMEWORK_FUNCS <- ls(envir = .GlobalEnv, all.names = TRUE)

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
              help = "bare flag = recompute for every dataset; or a comma-separated list of dataset ids"),
  make_option("--max-array-size", type = "integer", default = 1000,
              help = paste("max rows any ONE array task processes sequentially for",
                            "projectr_within_grid/projectr_cross_grid (fgsea_grid has its own",
                            "--fgsea-max-array-size below -- its per-row cost is minutes of real",
                            "GSEA/ORA compute, nothing like projectr's cheap per-row cost, so the two",
                            "were decoupled) -- a single small array job is",
                            "still submitted regardless of total grid size (see submit_job_family()'s",
                            "doc): a 3000-row grid with the default 1000 becomes one 3-task array job,",
                            "not 3000 array tasks (or the old behavior of three separate 1000-task",
                            "array-job submissions). Lower this to shorten each task's runtime, or",
                            "raise it to shrink the array further; it is NOT a Slurm MaxArraySize limit",
                            "to stay under (the resulting array is always small) -- pass 1 to fall back",
                            "to one row per array task, closest to the pre-2026-09-19 default.")),
  make_option("--fgsea-max-array-size", type = "integer", default = 20,
              help = paste("same idea as --max-array-size, but for fgsea_grid specifically. Default",
                            "20 (not 1000): confirmed against this project's real DB that fgsea_grid's",
                            "row count (~1000-2200, dominated by sPCA before representative_fit_ids()",
                            "collapsed it per-K) at the old shared default of 1000 became just 2-3",
                            "array tasks, each sequentially grinding through 700+ fits inside one",
                            "long-running R process -- the actual cause of observed fgsea_grid",
                            "timeouts/OOMs, not underprovisioned mem/time. A single 20-factor PCA fit",
                            "measured ~350s locally; 20 such fits/task stays comfortably inside the",
                            "slurm config's time/mem budget with real margin, and any one task failing",
                            "only costs ~20 fits of re-work instead of 700+."))
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
    ds_yaml_i <- yaml::read_yaml(targets$config_path[i])
    cache_dataset_matrix(con_cache, targets$dataset_id[i], ds_yaml_i, opt$db, force = force_i)
    # Also cache sample/feature metadata -- see cache_dataset_metadata()'s
    # header. Harmless/no-op here if this script itself is run somewhere
    # that can't reach the raw paths (e.g. the cluster) and a cache from
    # R/cache_dataset_matrices.R already exists; genuinely populates it
    # when run wherever the raw HuggingFace data resolves.
    cache_dataset_metadata(con_cache, targets$dataset_id[i], ds_yaml_i, opt$db, force = force_i)
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
  #
  # sample_metadata_maps: read from the CACHED artifact (datasets.
  # sample_metadata_file, populated by cache_dataset_metadata() -- see its
  # header) rather than the raw sample_metadata_path in each dataset's
  # config, which is typically an absolute path on whatever machine holds
  # the raw HuggingFace data -- NOT wherever this script actually gets
  # invoked (commonly the cluster login node, confirmed directly
  # (2026-09-18) to leave every dataset's map silently NULL otherwise,
  # since file.exists() on that raw path is always FALSE there). A
  # dataset with no cached artifact yet (e.g. never run through
  # cache_dataset_matrices.R/--stage core's own caching loop since this
  # fix landed) just gets NULL here; run_driver_job() skips it rather
  # than ever falling back to a raw-path read itself.
  sample_metadata_maps <- setNames(lapply(all_dataset_ids, function(ds) {
    f <- DBI::dbGetQuery(con, "SELECT sample_metadata_file FROM datasets WHERE dataset_id = ?",
                          params = list(ds))$sample_metadata_file
    if (length(f) != 1 || is.na(f)) return(NULL)
    path <- resolve_artifact(f, opt$db)
    if (!file.exists(path)) return(NULL)
    readRDS(path)
  }), all_dataset_ids)
  sample_id_cols <- setNames(lapply(all_dataset_ids, function(ds) {
    cp <- targets$config_path[match(ds, targets$dataset_id)]
    if (is.na(cp) || !file.exists(cp)) return(NULL)
    yaml::read_yaml(cp)$dataset$sample_id_col %||% "sample_id"
  }), all_dataset_ids)

  assign("all_dataset_ids", all_dataset_ids, envir = .GlobalEnv)
  assign("db_path", normalizePath(opt$db, mustWork = FALSE), envir = .GlobalEnv)
  assign("sample_metadata_maps", sample_metadata_maps, envir = .GlobalEnv)
  assign("sample_id_cols", sample_id_cols, envir = .GlobalEnv)
  submit_job_family(
    f = run_driver_job, jobs_df = NULL, jobname = "driver_grid",
    global_objects = c(FRAMEWORK_FUNCS, "all_dataset_ids", "db_path", "sample_metadata_maps", "sample_id_cols"),
    pkgs = c("DBI", "RSQLite", "arrow", "projectR"),
    cluster_cfg = slurm_cfg$driver, output_dir = opt$output,
    extra_binds = PROJECT_ROOT
  )

  # Ensembl maps (THE canonical cross-dataset identifier -- see
  # R/lib/ingest/symbol_mapping.R's header) + dataset matrix file lookups,
  # built once on the login node. Uses targets$config_path (already
  # case-corrected/absolutized -- see its own construction above) rather
  # than re-deriving "config/<id>_config.yml" from scratch, which is
  # case-SENSITIVE and silently misses every dataset whose config
  # filename lowercases a day/timepoint suffix (GSE110487_T2 ->
  # GSE110487_t2_config.yml, etc.). Config files themselves ARE reachable
  # here (checked into git) -- only feature_metadata_path (the raw parquet
  # build_ensembl_map() would otherwise read directly) commonly isn't, so
  # the CACHED artifact (datasets.feature_metadata_file, populated by
  # cache_dataset_metadata()) is read and passed in via `fm` instead. Same
  # rationale/confirmed failure mode as sample_metadata_maps above.
  ensembl_maps <- setNames(lapply(all_dataset_ids, function(ds) {
    cp <- targets$config_path[match(ds, targets$dataset_id)]
    if (is.na(cp) || !file.exists(cp)) return(NULL)
    f <- DBI::dbGetQuery(con, "SELECT feature_metadata_file FROM datasets WHERE dataset_id = ?",
                          params = list(ds))$feature_metadata_file
    fm <- NULL
    if (length(f) == 1 && !is.na(f)) {
      path <- resolve_artifact(f, opt$db)
      if (file.exists(path)) fm <- readRDS(path)
    }
    build_ensembl_map(yaml::read_yaml(cp), fm = fm)
  }), all_dataset_ids)

  target_matrix_path_for <- function(dataset_id) {
    f <- DBI::dbGetQuery(con, "SELECT matrix_file FROM datasets WHERE dataset_id = ?",
                          params = list(dataset_id))$matrix_file
    if (length(f) != 1 || is.na(f)) stop("No cached matrix for dataset '", dataset_id,
                                          "' -- run --stage core first (cache_dataset_matrix() populates this).")
    resolve_artifact(f, opt$db)
  }

  # ---- fgsea (+ local ORA/GSEA replacing gprofiler_grid): representative
  # fits only. Includes cp/tucker (added when the app's live-compute
  # enrichment path was retired -- they have the same loadings-based
  # shape as every other method here and representative_fit_ids() already
  # has a working fallback for them: "every ok fit is representative",
  # same idea as sPCA before it got its own per-K collapse, appropriate
  # since cp/tucker grids are small). WGCNA is NOT here -- no loadings to
  # rank, see R/ingest_jobs/wgcna_ora_job.R's separate job family below. ----
  rep_rows <- do.call(rbind, lapply(all_dataset_ids, function(ds) {
    do.call(rbind, lapply(c("pca", "nmf", "cogaps", "spca", "ica", "cp", "tucker"), function(m) {
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
    # collection = (not the deprecated category = -- msigdbr >= 10.0.0),
    # ensembl_gene = (not gene_symbol =) to match this pipeline's
    # canonical cross-dataset identifier space (see R/lib/ingest/
    # symbol_mapping.R's header) -- fgsea/fora both need `pathways` in the
    # SAME id space as the Ensembl-remapped loadings they're compared
    # against.
    fetch_msig <- function(collection, subcollection = NULL) {
      msig <- msigdbr::msigdbr(species = "Homo sapiens", collection = collection, subcollection = subcollection)
      msig <- msig[!is.na(msig$ensembl_gene) & nzchar(msig$ensembl_gene), ]
      split(msig$ensembl_gene, msig$gs_name)
    }
    pathways <- fetch_msig("H")
    assign("pathways", pathways, envir = .GlobalEnv)

    # NEW (2026-09-19): local, parallel replacement for gprofiler_grid's
    # per-request g:Profiler API calls (retired -- see git history/project
    # notes; g:Profiler itself confirmed via their own FAQ that they don't
    # offer a self-hosted/local instance for high query volumes, so a local
    # mirror of their DB isn't possible). `fora()`/`fgsea()` against these
    # same msigdbr collections reproduce gprofiler's ORA/GSEA modes
    # entirely locally, run as part of fgsea_grid (R/ingest_jobs/
    # fgsea_job.R) instead of a separate rate-limited job family -- see
    # that file's header for exactly how each collection is used.
    #
    # Collection/subcollection choices mirror gprofiler's default `sources
    # = c("GO:BP", "GO:MF", "REAC", "KEGG", "WP")`: GO:BP/GO:MF come from
    # msigdbr's C5 collection (subcollection strings match verbatim);
    # REACTOME/WIKIPATHWAYS from C2's CP:REACTOME/CP:WIKIPATHWAYS.  KEGG
    # is NOT a single subcollection in the installed msigdbr version
    # (confirmed via `msigdbr::msigdbr_collections()`: C2 splits it into
    # CP:KEGG_LEGACY and CP:KEGG_MEDICUS) -- CP:KEGG_LEGACY is used here as
    # the closer analog to gprofiler's classic KEGG pathway source.
    # HALLMARK reuses the same `pathways` object built above (already
    # Hallmark-only) rather than re-querying msigdbr for it.
    pathways_by_source <- list(
      HALLMARK = pathways,
      `GO:BP` = fetch_msig("C5", "GO:BP"),
      `GO:MF` = fetch_msig("C5", "GO:MF"),
      KEGG    = fetch_msig("C2", "CP:KEGG_LEGACY"),
      REAC    = fetch_msig("C2", "CP:REACTOME"),
      WP      = fetch_msig("C2", "CP:WIKIPATHWAYS")
    )
    assign("pathways_by_source", pathways_by_source, envir = .GlobalEnv)

    # per-row cogaps marker genes (already computed by run_all_redundancy() during --stage core)
    rep_rows$cogaps_marker_genes <- vector("list", nrow(rep_rows))
    cogaps_idx <- which(rep_rows$method == "cogaps")
    for (i in cogaps_idx) {
      mk <- DBI::dbGetQuery(con, "SELECT factor_index, gene FROM pattern_markers WHERE fit_id = ?",
                             params = list(rep_rows$fit_id[i]))
      rep_rows$cogaps_marker_genes[[i]] <- if (nrow(mk) > 0) split(mk$gene, mk$factor_index) else NULL
    }
    rep_rows$ensembl_map <- lapply(rep_rows$dataset_id, function(ds) ensembl_maps[[ds]])

    # gprofiler_grid (retired 2026-09-19) used to be staged here as a
    # separate job family covering the SAME (fit, factor, direction) grid
    # via live gprofiler2::gost() calls -- see fgsea_job.R's header and this
    # block's `pathways_by_source` comment above for why that's now folded
    # into fgsea_grid instead (local fora()/fgsea() against the same
    # msigdbr collections, no API rate limit, no separate container image).
    # The app's own per-factor, on-demand enrichment path (app/app.R) was
    # retired ENTIRELY, not just switched off gprofiler2 -- the app is now
    # a pure read-only viewer of whatever this pipeline has already
    # computed; see app/app.R's FGSEA_GRID_METHODS comment.
    submit_job_family(
      f = run_fgsea_job,
      jobs_df = rep_rows[, c("dataset_id", "method", "fit_id", "loadings_file", "ensembl_map", "cogaps_marker_genes")],
      jobname = "fgsea_grid", global_objects = c(FRAMEWORK_FUNCS, "pathways", "pathways_by_source"),
      pkgs = c("CoGAPS", "BiocParallel", "arrow"),
      cluster_cfg = slurm_cfg$fgsea, output_dir = opt$output,
      extra_binds = PROJECT_ROOT, max_array_size = opt$`fgsea-max-array-size`
    )
  } else {
    message("No representative fits found for fgsea -- skipping")
  }

  # ---- wgcna_ora: every WGCNA fit's every module, ORA only (no loadings
  # to rank, so no GSEA equivalent -- see R/ingest_jobs/wgcna_ora_job.R's
  # header). Reuses the SAME pathways_by_source built above for fgsea_grid.
  wgcna_fits_df <- DBI::dbGetQuery(con,
    "SELECT fit_id, dataset_id FROM fits WHERE method = 'wgcna' AND status = 'ok'")
  if (nrow(wgcna_fits_df) > 0 && exists("pathways_by_source", inherits = FALSE)) {
    wgcna_rows <- wgcna_fits_df
    # module 0 = WGCNA's "unassigned" -- excluded as a queryable gene set
    # (never gets a `factors` row either, see ingest_dataset.R), but its
    # genes still belong in the universe (real tested network genes).
    wgcna_rows$module_genes <- lapply(wgcna_rows$fit_id, function(fid) {
      mods <- DBI::dbGetQuery(con, "SELECT gene, module FROM wgcna_modules WHERE fit_id = ? AND module != 0",
                               params = list(fid))
      if (nrow(mods) == 0) return(NULL)
      split(mods$gene, mods$module)
    })
    wgcna_rows$universe_genes <- lapply(wgcna_rows$fit_id, function(fid) {
      DBI::dbGetQuery(con, "SELECT DISTINCT gene FROM wgcna_modules WHERE fit_id = ?", params = list(fid))$gene
    })
    wgcna_rows <- wgcna_rows[lengths(wgcna_rows$module_genes) > 0, ]
    wgcna_rows$ensembl_map <- lapply(wgcna_rows$dataset_id, function(ds) ensembl_maps[[ds]])

    if (nrow(wgcna_rows) > 0) {
      submit_job_family(
        f = run_wgcna_ora_job,
        jobs_df = wgcna_rows[, c("dataset_id", "fit_id", "module_genes", "universe_genes", "ensembl_map")],
        jobname = "wgcna_ora_grid", global_objects = c(FRAMEWORK_FUNCS, "pathways_by_source"),
        pkgs = c("BiocParallel"),
        cluster_cfg = slurm_cfg$wgcna_ora, output_dir = opt$output,
        extra_binds = PROJECT_ROOT, max_array_size = opt$`fgsea-max-array-size`
      )
    } else {
      message("No non-empty WGCNA module sets found -- skipping wgcna_ora")
    }
  } else {
    message("No WGCNA fits found (or fgsea_grid was skipped, so pathways_by_source was never built) -- skipping wgcna_ora")
  }

  # ---- projectr: reduced all-pairs, two explicitly-named families ----
  families <- yaml::read_yaml(opt$`dataset-families`)
  pairs <- build_projectr_pairs(con, all_dataset_ids, families)
  message("projectr: ", nrow(pairs), " total (source fit x target dataset) pairs -- ",
          sum(pairs$projection_type == "within_dataset"), " within-dataset, ",
          sum(pairs$projection_type == "cross_dataset"), " cross-dataset")

  # Looks up one column (via `query_fn`, called ONCE per DISTINCT key, not
  # once per row) and broadcasts it back across every row of `keys`. The
  # "reduced all-pairs" cross-dataset grid repeats the same handful of
  # source_fit_id/source_dataset_id/target_dataset_id values across tens of
  # thousands of rows (e.g. 64321 cross-dataset rows here, but only a few
  # hundred distinct source fits and ~30 distinct datasets) -- querying
  # per-row instead of per-distinct-key means the same tiny result gets
  # re-fetched thousands of times over, and on a networked/scratch
  # filesystem (where each SQLite round-trip can cost far more than on
  # local disk) that turns a sub-second lookup into a run that looks hung
  # for many minutes.
  lookup_by_distinct_key <- function(keys, query_fn) {
    uniq <- unique(keys)
    vals <- vapply(uniq, query_fn, character(1))
    unname(vals[match(keys, uniq)])
  }

  for (ptype in c("within_dataset", "cross_dataset")) {
    sub <- pairs[pairs$projection_type == ptype, ]
    if (nrow(sub) == 0) next
    sub$loadings_file <- lookup_by_distinct_key(sub$source_fit_id, function(fid) {
      resolve_artifact(DBI::dbGetQuery(con, "SELECT loadings_file FROM fits WHERE fit_id = ?",
                                        params = list(fid))$loadings_file, opt$db)
    })
    sub$target_matrix_file <- lookup_by_distinct_key(sub$target_dataset_id, target_matrix_path_for)
    sub$source_mat_file <- lookup_by_distinct_key(sub$source_dataset_id, function(ds) {
      f <- DBI::dbGetQuery(con, "SELECT matrix_file FROM datasets WHERE dataset_id = ?",
                            params = list(ds))$matrix_file
      if (length(f) == 1 && !is.na(f)) resolve_artifact(f, opt$db) else NA_character_
    })
    sub$source_scores_file <- lookup_by_distinct_key(sub$source_fit_id, function(fid) {
      f <- DBI::dbGetQuery(con, "SELECT scores_file FROM fits WHERE fit_id = ?", params = list(fid))$scores_file
      if (length(f) == 1 && !is.na(f)) resolve_artifact(f, opt$db) else NA_character_
    })

    assign("ensembl_maps", ensembl_maps, envir = .GlobalEnv)
    submit_job_family(
      f = run_projectr_job, jobs_df = sub,
      jobname = paste0("projectr_", sub("_dataset$", "", ptype), "_grid"),
      global_objects = c(FRAMEWORK_FUNCS, "ensembl_maps"), pkgs = c("projectR", "arrow"),
      cluster_cfg = slurm_cfg$projectr, output_dir = opt$output,
      extra_binds = PROJECT_ROOT, max_array_size = opt$`max-array-size`
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

# Core per-dataset ingest logic, extracted from R/ingest_results.R so both
# the single-dataset CLI (R/ingest_results.R) and the multi-dataset slurm
# job (R/ingest_jobs/ingest_core_job.R) share ONE implementation.
#
# Depends on R/lib/ingest/db.R, similarity.R, extract.R, pairs.R,
# redundancy.R, and R/lib/matrices.R (for cache_dataset_matrix()) already
# being sourced.

if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

#' Build (if not already cached, or if stale) and record this dataset's
#' feature x sample matrix as a stability_artifacts artifact -- the SAME
#' matrix every factorization job for this dataset was built from
#' (deterministic given its preprocessing_script/expression_path).
#'
#' Staleness is auto-detected: if ANY of preprocessing_script/
#' expression_path/sample_metadata_path/feature_metadata_path has an mtime
#' newer than the cached matrix.rds, it's rebuilt automatically (loudly,
#' never silently). `force = TRUE` (--recache-matrix) rebuilds regardless
#' -- the escape hatch for changes mtime can't see, e.g. editing a helper
#' script your preprocessing_script source()s.
cache_dataset_matrix <- function(con, dataset_id, dataset_yaml, db_path, force = FALSE) {
  ds <- dataset_yaml$dataset
  existing <- DBI::dbGetQuery(con, "SELECT matrix_file FROM datasets WHERE dataset_id = ?",
                               params = list(dataset_id))$matrix_file
  cached_path <- if (length(existing) == 1 && !is.na(existing)) resolve_artifact(existing, db_path) else NA_character_
  have_cache <- !is.na(cached_path) && file.exists(cached_path)

  stale <- FALSE
  if (!force && have_cache) {
    cached_mtime <- file.mtime(cached_path)
    source_files <- c(ds$preprocessing_script, ds$expression_path, ds$sample_metadata_path, ds$feature_metadata_path)
    source_files <- source_files[!is.null(source_files)]
    source_files <- source_files[file.exists(source_files)]
    if (length(source_files) > 0) {
      newer <- source_files[file.mtime(source_files) > cached_mtime]
      if (length(newer) > 0) {
        message("  cached matrix for ", dataset_id, " is STALE -- newer than cache: ",
                paste(basename(newer), collapse = ", "), " -- rebuilding")
        stale <- TRUE
      }
    }
  }

  if (!force && !stale && have_cache) return(invisible(existing))   # up to date, nothing to do

  message("  ", if (force) "force-recaching" else if (stale) "recaching (stale)" else "caching",
          " input matrix for ", dataset_id,
          if (!is.null(ds$preprocessing_script)) " (running preprocessing_script)" else "", "...")
  mat <- load_input_matrix(dataset_yaml)
  art_dir <- artifacts_dir(db_path, dataset_id)
  dir.create(art_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(mat, file.path(art_dir, "matrix.rds"))
  rel_path <- file.path("stability_artifacts", dataset_id, "matrix.rds")
  # ensure_dataset() first: a plain UPDATE silently affects 0 rows (and
  # never errors) for a dataset_id with no `datasets` row yet -- true for
  # any dataset that hasn't gone through ingest_one_dataset() before, which
  # is exactly the case the FIRST time this function ever runs for it. Without
  # this, matrix.rds gets written to disk but matrix_file is never actually
  # recorded, so every subsequent call sees "not cached" and rebuilds from
  # raw parquet again -- including on a machine (e.g. the cluster) that
  # doesn't have that parquet at all.
  ensure_dataset(con, dataset_id, ds$description %||% NA_character_)
  DBI::dbExecute(con, "UPDATE datasets SET matrix_file = ? WHERE dataset_id = ?",
                 params = list(rel_path, dataset_id))
  invisible(rel_path)
}

#' Ingest ONE dataset's job-family results (+ redundancy) into `con`.
#' Identical behavior to what R/ingest_results.R's body used to do inline;
#' extracted so R/ingest_jobs/ingest_core_job.R can loop this over many
#' datasets in a single slurm job (one SQLite writer, avoiding concurrent
#' writes -- see R/README.md's ingest-slurm design).
#'
#' NOTE: deliberately does NOT call cache_dataset_matrix() -- that function
#' runs `preprocessing_script`/`load_input_matrix()`, which `source()`s a
#' project-relative file path and reads the raw parquet trio. Neither is
#' reachable from inside the ingest_core slurm job's container (only
#' slurm_bundles/ingest/ is bind-mounted). cache_dataset_matrix() must be
#' called on the LOGIN NODE instead, before staging ingest_core -- see
#' R/create_ingest_slurm_bundle.R's --stage core.
#'
#' @param overwrite FALSE, TRUE (every family), or a character vector of
#'   jobnames to replace.
#' @param recompute_redundancy passed straight through to
#'   run_all_redundancy()'s `force` argument.
#' @param run_pattern_drivers whether to also run the
#'   projectR::projectionDriveR()-based pattern-driver pass (see
#'   R/lib/ingest/driver.R) inline. Default TRUE for R/ingest_results.R's
#'   direct (non-slurm) CLI use, where projectR is just whatever's in that
#'   R session. FALSE for run_ingest_core_job() specifically -- ingest_core
#'   runs inside a container that does NOT have projectR installed (a
#'   different image than the one it needs -- see config/
#'   ingest_slurm_config.yml's `driver:` entry); the pattern-driver pass is
#'   staged as its own `driver_grid` job family instead, during
#'   --stage enrichment (see R/ingest_jobs/driver_job.R).
#' @param project_root absolute project-root path, used (a) as the
#'   fallback location for this dataset's rslurm bundle when it isn't
#'   colocated with results_dir -- see bundle_dir below -- and (b) to
#'   locate redundancy.R's own source file for staleness detection (see
#'   run_all_redundancy()). Defaults to getwd(), correct for
#'   R/ingest_results.R's direct CLI use (run from the project root).
#'   run_ingest_core_job() passes the container's bind-mounted PROJECT_ROOT
#'   explicitly instead, since its own cwd is rslurm's bundle directory,
#'   NOT the project root (see that function's header for why cwd can't
#'   just be overridden to fix this the other way around).
ingest_one_dataset <- function(con, config_path, results_dir, db_path,
                                overwrite = FALSE, recompute_redundancy = FALSE,
                                run_pattern_drivers = TRUE, project_root = getwd()) {
  dataset_yaml <- yaml::read_yaml(config_path)
  dataset_id <- dataset_yaml$dataset$id
  stopifnot(!is.null(dataset_id))
  message("dataset: ", dataset_id)

  family_dirs <- Filter(
    function(d) length(list.files(d, pattern = "^results_\\d+\\.RDS$")) > 0,
    list.dirs(results_dir, recursive = FALSE)
  )
  if (length(family_dirs) == 0) {
    message("  no job family subdirectories with results_*.RDS in ", results_dir, " -- skipping")
    return(invisible(dataset_id))
  }
  jobnames <- basename(family_dirs)
  message("  families found: ", paste(jobnames, collapse = ", "))

  # Each family's ORIGINAL rslurm bundle (_rslurm_<jobname>/params.RDS) --
  # NOT the same tree as family_dirs/results_dir above (that's where
  # RESULTS landed; the bundle is where params.RDS/f.RDS/etc. that
  # PRODUCED those results still live). Two conventions in the wild here:
  # (1) results_dir's own sibling, dropping its "_results" suffix -- e.g.
  # /scratch/.../GSE110487_T2_results (results) next to /scratch/.../
  # GSE110487_T2 (bundle) -- what R/lib/submit_all_script.R's rsync
  # workflow actually produces when deployed to a cluster, one dataset's
  # bundle/results copied in as two independent top-level directories; (2)
  # project_root/slurm_bundles/<dataset_id> -- what R/create_slurm_bundle.R
  # writes locally when bundle-building and ingest happen from the SAME
  # project checkout, never separately relocated. Try (1) first since it's
  # colocated with results_dir (the input we actually have in hand), (2) as
  # a fallback for co-located/local use.
  bundle_dir <- sub("_results$", "", results_dir)
  if (!dir.exists(bundle_dir)) bundle_dir <- file.path(project_root, "slurm_bundles", dataset_id)

  ensure_dataset(con, dataset_id, dataset_yaml$dataset$description %||% NA_character_)
  art_dir <- artifacts_dir(db_path, dataset_id)
  dir.create(art_dir, recursive = TRUE, showWarnings = FALSE)

  ds_cfg <- dataset_yaml$dataset
  register_metadata_source(con, dataset_id, "sample", ds_cfg$sample_metadata_path, ds_cfg$sample_id_col)
  register_metadata_source(con, dataset_id, "feature", ds_cfg$feature_metadata_path, ds_cfg$feature_id_col,
                            ensembl_col = ds_cfg$ensembl_col %||% "ensembl",
                            symbol_col = ds_cfg$symbol_col %||% NA_character_)

  overwrite_requested <- function(jobname) {
    isTRUE(overwrite) || (is.character(overwrite) && jobname %in% overwrite)
  }

  new_fits_by_method <- list()
  report <- list()

  for (k in seq_along(jobnames)) {
    jobname <- jobnames[k]
    fam_dir <- family_dirs[k]
    cls <- classify_jobname(jobname)

    if (family_already_ingested(con, dataset_id, jobname)) {
      if (overwrite_requested(jobname)) {
        message("  [", jobname, "] already ingested -- OVERWRITING")
        delete_family(con, db_path, dataset_id, jobname)
      } else {
        message("  [", jobname, "] already ingested -- skipping (use overwrite = '", jobname, "' to replace)")
        report[[jobname]] <- "skipped"
        next
      }
    }

    params_path <- file.path(bundle_dir, paste0("_rslurm_", jobname), "params.RDS")
    if (!file.exists(params_path)) {
      message("  [", jobname, "] bundle params not found at ", params_path, " -- skipping this family")
      report[[jobname]] <- "skipped (no params.RDS)"
      next
    }
    params <- readRDS(params_path)

    result_files <- list.files(fam_dir, pattern = "^results_\\d+\\.RDS$", full.names = TRUE)
    task_ids <- as.integer(sub("^results_(\\d+)\\.RDS$", "\\1", basename(result_files)))

    n_ok <- 0L; n_failed <- 0L; n_missing <- 0L; n_scores <- 0L
    new_ids <- integer(0)

    DBI::dbExecute(con, "BEGIN")
    for (task in seq_len(nrow(params)) - 1L) {
      params_row <- params[task + 1L, , drop = FALSE]
      f <- result_files[match(task, task_ids)]
      result <- if (!is.na(f)) {
        x <- readRDS(f)
        if (is.list(x) && length(x) == 1 && is.null(names(x))) x[[1]] else x
      } else {
        NULL
      }

      ext <- extract_result(jobname, result, params_row)
      fit <- ext$fit

      loadings_file <- NA_character_
      if (!is.null(ext$loadings)) {
        fname <- sprintf("%s_task%03d_loadings.rds", jobname, task)
        saveRDS(ext$loadings, file.path(art_dir, fname))
        loadings_file <- file.path("stability_artifacts", dataset_id, fname)
      }

      scores_file <- NA_character_
      if (!is.null(ext$scores)) {
        fname <- sprintf("%s_task%03d_scores.rds", jobname, task)
        saveRDS(ext$scores, file.path(art_dir, fname))
        scores_file <- file.path("stability_artifacts", dataset_id, fname)
      }

      time_loadings_file <- NA_character_
      if (!is.null(ext$time_loadings)) {
        fname <- sprintf("%s_task%03d_time_loadings.rds", jobname, task)
        saveRDS(ext$time_loadings, file.path(art_dir, fname))
        time_loadings_file <- file.path("stability_artifacts", dataset_id, fname)
      }

      # CoGAPS only: stash the FULL raw CogapsResult while it's still in
      # scope (before it's discarded) -- see R/lib/ingest/redundancy.R's
      # header for why patternMarkers() needs the real object.
      raw_result_file <- NA_character_
      if (ext$method == "cogaps" && fit$status == "ok" && !is.null(result) && !is.null(result$result)) {
        fname <- sprintf("%s_task%03d_cograw.rds", jobname, task)
        saveRDS(result$result, file.path(art_dir, fname))
        raw_result_file <- file.path("stability_artifacts", dataset_id, fname)
      }

      DBI::dbExecute(con,
        "INSERT INTO fits (dataset_id, method, family, jobname, rank, seed, alpha, power,
                           rank_genes, rank_subjects, rank_time,
                           mse, n_factors, status, loadings_file, scores_file, time_loadings_file, raw_result_file)
         VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        params = list(dataset_id, ext$method, ext$family, jobname,
                      fit$rank, fit$seed, fit$alpha, fit$power,
                      fit$rank_genes, fit$rank_subjects, fit$rank_time,
                      fit$mse, fit$n_factors, fit$status, loadings_file, scores_file, time_loadings_file,
                      raw_result_file))
      fit_id <- DBI::dbGetQuery(con, "SELECT last_insert_rowid() AS id")$id

      if (fit$status == "ok") {
        n_ok <- n_ok + 1L
        if (!is.na(scores_file)) n_scores <- n_scores + 1L
        new_ids <- c(new_ids, fit_id)
        if (!is.null(ext$loadings)) {
          DBI::dbWriteTable(con, "factors", data.frame(
            fit_id = fit_id, factor_index = seq_len(ncol(ext$loadings)),
            stability_cosine = NA_real_, stability_pearson = NA_real_,
            stability_spearman = NA_real_
          ), append = TRUE)
        }
        if (!is.null(ext$modules)) {
          DBI::dbWriteTable(con, "wgcna_modules",
                            cbind(fit_id = fit_id, ext$modules), append = TRUE)
          mod_ids <- setdiff(sort(unique(ext$modules$module)), 0L)
          if (length(mod_ids) > 0) {
            DBI::dbWriteTable(con, "factors", data.frame(
              fit_id = fit_id, factor_index = mod_ids,
              stability_cosine = NA_real_, stability_pearson = NA_real_,
              stability_spearman = NA_real_
            ), append = TRUE)
          }
        }
      } else if (fit$status == "failed") {
        n_failed <- n_failed + 1L
      } else {
        n_missing <- n_missing + 1L
      }

      if (ext$family == "maskcv") {
        DBI::dbExecute(con,
          "INSERT INTO maskcv_results (dataset_id, method, jobname, rank, alpha, mse)
           VALUES (?, ?, ?, ?, ?, ?)",
          params = list(dataset_id, ext$method, jobname, fit$rank, fit$alpha, fit$mse))
      }
    }
    record_ingest(con, dataset_id, jobname, cls$family, cls$method,
                  n_results = length(result_files), results_dir = normalizePath(results_dir))
    DBI::dbExecute(con, "COMMIT")

    if (cls$family != "maskcv" && length(new_ids) > 0) {
      key <- cls$method
      new_fits_by_method[[key]] <- c(new_fits_by_method[[key]], new_ids)
    }
    report[[jobname]] <- sprintf("ingested (%d ok, %d failed, %d missing of %d tasks)",
                                  n_ok, n_failed, n_missing, nrow(params))
    message("  [", jobname, "] ", report[[jobname]])
  }

  for (method in names(new_fits_by_method)) {
    ids <- new_fits_by_method[[method]]
    message("  [pairs:", method, "] computing similarities for ", length(ids), " new fits ...")
    if (method == "wgcna") {
      compute_wgcna_pairs(con, dataset_id, ids)
    } else {
      n <- compute_factor_pairs(con, db_path, dataset_id, method, ids)
      message("  [pairs:", method, "] wrote ", n, " factor-pair rows")
      update_factor_stability(con, dataset_id, method)
    }
  }

  run_all_redundancy(con, db_path, dataset_id, force = recompute_redundancy, project_root = project_root)

  # Differential feature identification (projectR::projectionDriveR(), see
  # R/lib/ingest/driver.R) -- batch pass over representative fits x
  # 2-4-level categorical sample-metadata columns. Silently no-ops if the
  # dataset's matrix isn't cached yet (cache_dataset_matrix() is
  # login-node-only) or has no small categorical columns -- never blocks
  # the rest of ingest. Skipped entirely when run_pattern_drivers = FALSE
  # (see this function's @param doc) -- staged as its own driver_grid job
  # family instead in that case.
  if (run_pattern_drivers) {
    tryCatch(run_all_pattern_drivers(con, db_path, dataset_id),
             error = function(e) message("  pattern-driver pass failed for ", dataset_id, ": ", conditionMessage(e)))
  }

  invisible(dataset_id)
}

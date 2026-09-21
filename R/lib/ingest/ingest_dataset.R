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

#' Build (if not already cached, or if stale) and record this dataset's
#' sample AND feature metadata as stability_artifacts artifacts -- same
#' login-node-only rationale as cache_dataset_matrix() (reads raw
#' sample_metadata_path/feature_metadata_path, which typically only
#' resolve on whatever machine holds the raw HuggingFace data, NOT
#' wherever create_ingest_slurm_bundle.R actually gets invoked -- commonly
#' the cluster login node).
#'
#' Needed because R/create_ingest_slurm_bundle.R's --stage enrichment
#' builds `sample_metadata_maps` (driver_grid) and `ensembl_maps`
#' (fgsea_grid/projectr_*_grid) by reading each dataset's
#' metadata -- confirmed directly (2026-09-18): without a cached artifact
#' to fall back on, EVERY dataset's map silently came back NULL when
#' staged from the cluster, since file.exists() on the raw config paths is
#' always FALSE there (identical failure mode to the original
#' cache_dataset_matrix() bug this mirrors). Both `build_ensembl_map()`/
#' `build_symbol_map()` (R/lib/ingest/symbol_mapping.R) accept an optional
#' pre-loaded `fm` to read the CACHED artifact instead of re-reading raw
#' feature_metadata_path -- see their docs.
#'
#' Gracefully keeps whatever's already cached (rather than erroring) if
#' the raw path is unreachable from wherever this happens to run --
#' running this on the machine that actually HAS the raw data (same as
#' cache_dataset_matrix()) is what actually populates the cache in the
#' first place; running it elsewhere afterward is a safe no-op.
cache_dataset_metadata <- function(con, dataset_id, dataset_yaml, db_path, force = FALSE) {
  ds <- dataset_yaml$dataset
  art_dir <- artifacts_dir(db_path, dataset_id)
  dir.create(art_dir, recursive = TRUE, showWarnings = FALSE)

  cache_one <- function(kind, raw_path, col) {
    if (is.null(raw_path)) return(invisible(NULL))
    existing <- DBI::dbGetQuery(con, sprintf("SELECT %s FROM datasets WHERE dataset_id = ?", col),
                                 params = list(dataset_id))[[col]]
    cached_path <- if (length(existing) == 1 && !is.na(existing)) resolve_artifact(existing, db_path) else NA_character_
    have_cache <- !is.na(cached_path) && file.exists(cached_path)
    raw_reachable_here <- file.exists(raw_path)
    stale <- have_cache && !force && raw_reachable_here && file.mtime(raw_path) > file.mtime(cached_path)

    if (!force && !stale && have_cache) return(invisible(existing))
    if (!raw_reachable_here) {
      if (have_cache) return(invisible(existing))   # keep existing cache -- can't rebuild from here
      message("  no ", kind, "_metadata cached for ", dataset_id, " and raw path unreachable here (",
              raw_path, ") -- skipping")
      return(invisible(NULL))
    }

    message("  ", if (force) "force-recaching" else if (stale) "recaching (stale)" else "caching",
            " ", kind, " metadata for ", dataset_id, "...")
    df <- as.data.frame(arrow::read_parquet(raw_path))
    fname <- paste0(kind, "_metadata.rds")
    saveRDS(df, file.path(art_dir, fname))
    rel_path <- file.path("stability_artifacts", dataset_id, fname)
    ensure_dataset(con, dataset_id, ds$description %||% NA_character_)
    DBI::dbExecute(con, sprintf("UPDATE datasets SET %s = ? WHERE dataset_id = ?", col),
                   params = list(rel_path, dataset_id))
    invisible(rel_path)
  }

  cache_one("sample", ds$sample_metadata_path, "sample_metadata_file")
  cache_one("feature", ds$feature_metadata_path, "feature_metadata_file")
  invisible(NULL)
}

#' Compute (or refresh) WGCNA::pickSoftThreshold()'s scale-free-topology fit
#' diagnostic for one dataset, over its configured power grid
#' (methods.network.wgcna.power) -- the tutorial-standard justification for a
#' soft-threshold power choice, which this pipeline otherwise skips (it
#' instead sweeps the full grid and judges power only by cross-power module
#' stability -- see R/methods/wgcna.R's header). A dataset-level diagnostic,
#' not tied to any one wgcna_grid fit: it needs only the cached input matrix
#' (cache_dataset_matrix() must have already run) + the power grid, so it's
#' stored in its own `wgcna_sft` table keyed by dataset_id, not by fit_id.
#'
#' Must be called from an ordinary (non-containerized) R session with the
#' WGCNA package installed -- same constraint as cache_dataset_matrix(),
#' and for the same reason blockwiseModules() itself needs its own
#' container (methods.network.wgcna.container): the `ingest_core` slurm
#' container does not have WGCNA installed. Call this from
#' R/ingest_results.R or R/cache_dataset_matrices.R (both already
#' non-containerized, matrix-caching entry points), never from
#' run_ingest_core_job().
#'
#' Silently no-ops if `methods.network.wgcna` isn't configured for this
#' dataset, or if the matrix hasn't been cached yet. `force = TRUE` always
#' recomputes; otherwise recomputes only if the configured power grid has
#' changed since the last computation (replace-on-recompute, like
#' matrix_file -- not additive).
compute_wgcna_sft <- function(con, dataset_id, dataset_yaml, db_path, force = FALSE) {
  wg <- dataset_yaml$methods$network$wgcna
  if (is.null(wg) || is.null(wg$power)) return(invisible(NULL))

  existing <- DBI::dbGetQuery(con, "SELECT power FROM wgcna_sft WHERE dataset_id = ?",
                               params = list(dataset_id))$power
  if (!force && length(existing) > 0 && setequal(existing, as.integer(wg$power))) {
    return(invisible(NULL))   # already computed for this exact grid
  }

  mat_row <- DBI::dbGetQuery(con, "SELECT matrix_file FROM datasets WHERE dataset_id = ?",
                              params = list(dataset_id))
  if (nrow(mat_row) == 0 || is.na(mat_row$matrix_file)) {
    message("  no cached matrix for ", dataset_id, " -- run cache_dataset_matrix() first -- skipping SFT fit")
    return(invisible(NULL))
  }
  mat <- as.matrix(readRDS(resolve_artifact(mat_row$matrix_file, db_path)))

  message("  computing scale-free-topology fit for ", dataset_id, " (", length(wg$power), " powers)...")
  datExpr <- t(mat)
  gsg <- WGCNA::goodSamplesGenes(datExpr, verbose = 0)   # mirrors run_wgcna_param_job()'s own filter
  if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]

  sft <- WGCNA::pickSoftThreshold(datExpr, powerVector = as.integer(wg$power),
                                   networkType = wg$networkType %||% "signed", verbose = 0)
  fi <- sft$fitIndices

  DBI::dbExecute(con, "DELETE FROM wgcna_sft WHERE dataset_id = ?", params = list(dataset_id))
  DBI::dbWriteTable(con, "wgcna_sft", data.frame(
    dataset_id = dataset_id, power = fi$Power,
    sft_r_sq = fi$SFT.R.sq, slope = fi$slope, truncated_r_sq = fi$truncated.R.sq,
    mean_k = fi$mean.k., median_k = fi$median.k., max_k = fi$max.k.,
    computed_at = as.character(Sys.time()),
    stringsAsFactors = FALSE
  ), append = TRUE)
  invisible(NULL)
}

#' Compute WGCNA intramodular connectivity / module membership (kME) for
#' every already-ingested WGCNA fit of this dataset that has module
#' eigengenes (fits.scores_file -- see extract.R's wgcna branch; older
#' fits ingested before `samples` was captured have none and are silently
#' skipped). `WGCNA::signedKME(datExpr, MEs)` correlates every gene's
#' expression against every module eigengene -- this is the standard
#' mechanism for identifying hub genes (highest |kME| within their
#' assigned module) that `blockwiseModules()` computes INTERNALLY during
#' module trimming/merging but never returns (confirmed via
#' ?blockwiseModules's Details section) -- so it must be recomputed
#' separately, exactly like compute_wgcna_sft() recomputes the SFT fit
#' blockwiseModules() also never returns.
#'
#' Per-FIT (not per-dataset, unlike SFT/gene-significance below): each
#' power's module structure/eigengenes differ, so kME differs per fit.
#' Same non-containerized login-node constraint as compute_wgcna_sft() --
#' needs WGCNA installed and the cached matrix already on disk.
#'
#' `force = TRUE` recomputes for every fit even if wgcna_kme already has
#' rows for it; otherwise only fits with zero existing rows are computed
#' (additive across fits, unlike matrix_file/wgcna_sft's replace-on-
#' recompute -- a NEW wgcna fit ingested later just adds its own rows).
compute_wgcna_kme <- function(con, dataset_id, dataset_yaml, db_path, force = FALSE) {
  wg <- dataset_yaml$methods$network$wgcna
  if (is.null(wg)) return(invisible(NULL))

  mat_row <- DBI::dbGetQuery(con, "SELECT matrix_file FROM datasets WHERE dataset_id = ?",
                              params = list(dataset_id))
  if (nrow(mat_row) == 0 || is.na(mat_row$matrix_file)) {
    message("  no cached matrix for ", dataset_id, " -- skipping WGCNA kME")
    return(invisible(NULL))
  }
  datExpr_full <- t(as.matrix(readRDS(resolve_artifact(mat_row$matrix_file, db_path))))

  fits <- DBI::dbGetQuery(con,
    "SELECT fit_id, scores_file FROM fits
     WHERE dataset_id = ? AND method = 'wgcna' AND status = 'ok' AND scores_file IS NOT NULL",
    params = list(dataset_id))
  if (nrow(fits) == 0) return(invisible(NULL))

  for (i in seq_len(nrow(fits))) {
    fit_id <- fits$fit_id[i]
    n_existing <- DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM wgcna_kme WHERE fit_id = ?",
                                   params = list(fit_id))$n
    if (!force && n_existing > 0) next

    MEs <- as.matrix(readRDS(resolve_artifact(fits$scores_file[i], db_path)))
    common_samples <- intersect(rownames(datExpr_full), rownames(MEs))
    if (length(common_samples) < 3) next
    datExpr <- datExpr_full[common_samples, , drop = FALSE]
    MEs <- MEs[common_samples, , drop = FALSE]

    kme <- WGCNA::signedKME(datExpr, as.data.frame(MEs), outputColumnName = "kME")
    module_ids <- suppressWarnings(as.integer(sub("^ME", "", colnames(MEs))))

    message("  [wgcna_kme] fit_id ", fit_id, ": ", ncol(datExpr), " genes x ", ncol(MEs), " modules")
    DBI::dbExecute(con, "DELETE FROM wgcna_kme WHERE fit_id = ?", params = list(fit_id))
    DBI::dbWriteTable(con, "wgcna_kme", data.frame(
      fit_id = fit_id,
      gene = rep(colnames(datExpr), times = ncol(kme)),
      module = rep(module_ids, each = nrow(kme)),
      kme = as.vector(as.matrix(kme)),
      stringsAsFactors = FALSE
    ), append = TRUE)
  }
  invisible(NULL)
}

#' Compute WGCNA Gene Significance (GS) -- per-gene correlation of
#' expression against each sample-metadata trait -- for this dataset.
#' Dataset-level, not per-fit: GS depends only on the cached expression
#' matrix + registered sample metadata, neither of which varies by
#' module/power choice, unlike kME above. Combined with kME, reconstructs
#' the classic GS-vs-MM hub-gene scatter the WGCNA workflow is built
#' around (a plot neither table alone supports).
#'
#' Reuses app/R/metadata_helpers.R's generic_association_scan() (Spearman
#' for numeric fields, Kruskal-Wallis otherwise, no field names hardcoded)
#' verbatim -- passing the full expression matrix in place of a
#' factor/module scores matrix works unchanged, since that function only
#' ever treats its first argument as "samples x things to correlate
#' against metadata," genes being no different from factors/eigengenes
#' for that purpose. Sourced lazily (once) since the ingest pipeline
#' doesn't otherwise depend on any app/R file.
#'
#' Same non-containerized login-node constraint as compute_wgcna_sft()/
#' compute_wgcna_kme() (needs the cached matrix + registered metadata
#' pointer already on disk here). `force = TRUE` always recomputes;
#' otherwise no-ops if this dataset already has any rows (replace-on-
#' recompute like matrix_file/wgcna_sft, not additive).
compute_wgcna_gene_significance <- function(con, dataset_id, dataset_yaml, db_path, force = FALSE) {
  wg <- dataset_yaml$methods$network$wgcna
  if (is.null(wg)) return(invisible(NULL))
  if (!exists("generic_association_scan")) {
    source(here::here("app/R/metadata_helpers.R"))
  }

  n_existing <- DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM wgcna_gene_significance WHERE dataset_id = ?",
                                 params = list(dataset_id))$n
  if (!force && n_existing > 0) return(invisible(NULL))

  mat_row <- DBI::dbGetQuery(con, "SELECT matrix_file FROM datasets WHERE dataset_id = ?",
                              params = list(dataset_id))
  if (nrow(mat_row) == 0 || is.na(mat_row$matrix_file)) {
    message("  no cached matrix for ", dataset_id, " -- skipping WGCNA gene significance")
    return(invisible(NULL))
  }
  datExpr <- t(as.matrix(readRDS(resolve_artifact(mat_row$matrix_file, db_path))))

  meta <- dataset_metadata(con, dataset_id, "sample")
  if (is.null(meta)) {
    message("  no sample metadata registered for ", dataset_id, " -- skipping WGCNA gene significance")
    return(invisible(NULL))
  }

  message("  computing WGCNA gene significance for ", dataset_id, " (", ncol(datExpr), " genes x ",
          length(setdiff(names(meta), "sample_id")), " fields)...")
  gs <- generic_association_scan(datExpr, meta, id_col = "sample_id")
  if (nrow(gs) == 0) return(invisible(NULL))

  DBI::dbExecute(con, "DELETE FROM wgcna_gene_significance WHERE dataset_id = ?", params = list(dataset_id))
  DBI::dbWriteTable(con, "wgcna_gene_significance", data.frame(
    dataset_id = dataset_id, gene = gs$component, field = gs$field,
    statistic = gs$statistic, p_value = gs$p_value, stringsAsFactors = FALSE
  ), append = TRUE)
  invisible(NULL)
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

      # method-specific diagnostics bundle (see db.R's *_diag_file
      # columns / extract_result()'s `diag`) -- one small named-list
      # artifact per method, saved into the ONE column that method uses;
      # every other method's *_diag_file column stays NULL for this row.
      diag_cols <- c(cogaps = "cogaps_diag_file", pca = "pca_diag_file",
                      spca = "spca_diag_file", ica = "ica_diag_file",
                      nmf = "nmf_diag_file", cp = "cp_diag_file", tucker = "tucker_diag_file")
      diag_file_val <- list(cogaps_diag_file = NA_character_, pca_diag_file = NA_character_,
                             spca_diag_file = NA_character_, ica_diag_file = NA_character_,
                             nmf_diag_file = NA_character_, cp_diag_file = NA_character_,
                             tucker_diag_file = NA_character_)
      if (!is.null(ext$diag) && ext$method %in% names(diag_cols)) {
        fname <- sprintf("%s_task%03d_diag.rds", jobname, task)
        saveRDS(ext$diag, file.path(art_dir, fname))
        diag_file_val[[diag_cols[[ext$method]]]] <- file.path("stability_artifacts", dataset_id, fname)
      }

      DBI::dbExecute(con,
        "INSERT INTO fits (dataset_id, method, family, jobname, rank, seed, alpha, power,
                           rank_genes, rank_subjects, rank_time,
                           mse, n_factors, status, converged, loadings_file, scores_file, time_loadings_file, raw_result_file,
                           cogaps_diag_file, pca_diag_file, spca_diag_file, ica_diag_file, nmf_diag_file, cp_diag_file, tucker_diag_file)
         VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        params = list(dataset_id, ext$method, ext$family, jobname,
                      fit$rank, fit$seed, fit$alpha, fit$power,
                      fit$rank_genes, fit$rank_subjects, fit$rank_time,
                      fit$mse, fit$n_factors, fit$status, fit$converged, loadings_file, scores_file, time_loadings_file,
                      raw_result_file,
                      diag_file_val$cogaps_diag_file, diag_file_val$pca_diag_file, diag_file_val$spca_diag_file,
                      diag_file_val$ica_diag_file, diag_file_val$nmf_diag_file, diag_file_val$cp_diag_file,
                      diag_file_val$tucker_diag_file))
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

    }
    record_ingest(con, dataset_id, jobname, cls$family, cls$method,
                  n_results = length(result_files), results_dir = normalizePath(results_dir))
    DBI::dbExecute(con, "COMMIT")

    if (length(new_ids) > 0) {
      key <- cls$method
      new_fits_by_method[[key]] <- c(new_fits_by_method[[key]], new_ids)
    }
    report[[jobname]] <- sprintf("ingested (%d ok, %d failed, %d missing of %d tasks)",
                                  n_ok, n_failed, n_missing, nrow(params))
    message("  [", jobname, "] ", report[[jobname]])
  }

  # Wrapped in ONE transaction: compute_wgcna_pairs()'s full pairwise
  # double loop (R/lib/ingest/pairs.R) and run_all_redundancy()'s
  # per-representative-fit loop (R/lib/ingest/redundancy.R) both issue
  # many small autocommit dbExecute()/dbWriteTable() calls with no
  # batching of their own -- confirmed the real cause of a genuinely
  # measured ingest_core seff report showing only 28.42% CPU efficiency
  # over a 1h19m single-core run (i.e. ~56 minutes of the wall-clock spent
  # blocked on per-statement fsyncs, not computing) -- same "death by a
  # thousand fsyncs" pattern already diagnosed and fixed for driver_grid
  # earlier (see run_all_pattern_drivers()'s own transaction-wrapping
  # comment in driver.R for the original real-seff-backed diagnosis this
  # mirrors). on.exit()'s rollback-if-not-committed guard ensures an
  # error partway through never leaves an open transaction for the NEXT
  # dataset's writes (run_ingest_core_job() loops many datasets against
  # the same connection) to land inside.
  DBI::dbExecute(con, "BEGIN")
  pairs_committed <- FALSE
  on.exit(if (!pairs_committed) DBI::dbExecute(con, "ROLLBACK"), add = TRUE)

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

  DBI::dbExecute(con, "COMMIT")
  pairs_committed <- TRUE

  # Differential feature identification (projectR::projectionDriveR(), see
  # R/lib/ingest/driver.R) -- batch pass over representative fits x
  # 2-4-level categorical sample-metadata columns. Silently no-ops if the
  # dataset's matrix isn't cached yet (cache_dataset_matrix() is
  # login-node-only) or has no small categorical columns -- never blocks
  # the rest of ingest. Skipped entirely when run_pattern_drivers = FALSE
  # (see this function's @param doc) -- staged as its own driver_grid job
  # family instead in that case.
  if (run_pattern_drivers) {
    # CI (default) and PV are the projectR fork's own paired standard
    # modes (its vignette runs both) -- see R/ingest_jobs/driver_job.R's
    # matching comment for why PV is a second, independently useful pass
    # rather than a redundant re-run of CI.
    for (mode in c("CI", "PV")) {
      tryCatch(run_all_pattern_drivers(con, db_path, dataset_id, mode = mode),
               error = function(e) message("  pattern-driver pass (", mode, ") failed for ", dataset_id, ": ", conditionMessage(e)))
    }
  }

  invisible(dataset_id)
}

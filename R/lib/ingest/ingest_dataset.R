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
#' run_ingest_core_compute_job() (its container lacks WGCNA).
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

#' Compute ONE dataset's full ingest bundle -- everything ingest_one_dataset()
#' used to do LIVE against a writer connection, EXCEPT the actual DB writes,
#' so this can run inside a parallel slurm ARRAY task with no shared writer
#' at all (see R/ingest_jobs/ingest_core_job.R::run_ingest_core_compute_job()
#' and R/create_ingest_slurm_bundle.R's `--stage core`). Returns a `bundle`
#' list consumed by write_ingest_bundle() below.
#'
#' Every fit-referencing id in the bundle uses this convention: POSITIVE =
#' a real, already-merged `fits.fit_id` (read from the DB); NEGATIVE =
#' `-local_id`, a placeholder for one of THIS bundle's own new fits, which
#' doesn't have a real fit_id yet (only assigned once write_ingest_bundle()
#' actually INSERTs it). `bundle$fits$local_id` itself is always a small
#' positive integer (1, 2, 3, ... in insertion order within this bundle) --
#' `write_ingest_bundle()` inserts fits in that exact order and remaps
#' `-local_id` references to the real fit_id each INSERT returns.
#'
#' Opens its own DB connection to `db_path`, but ONLY for reads (existing
#' fits/pairs/redundancy-cache lookups needed for "new x existing"
#' comparisons and incremental staleness checks -- see
#' R/lib/ingest/pairs.R/redundancy.R's `_from_universe` functions) --
#' never issues BEGIN/INSERT/UPDATE/DELETE. This is safe to run from many
#' concurrent array tasks against the SAME db_path at once: SQLite's WAL
#' mode allows unlimited concurrent readers, and no writer lock is ever
#' requested here (see R/lib/ingest/db.R's WAL/busy_timeout pragmas and
#' this project's own audit of SQLite's documented WAL concurrency model).
#'
#' Cross-dataset coupling is a non-issue -- every artifact path, SQL
#' filter, and pair/redundancy comparison here is scoped to THIS ONE
#' dataset_id, so one array task per dataset never needs to see another
#' task's in-flight work (confirmed directly: compute_factor_pairs()/
#' compute_wgcna_pairs() only ever compare fits WITHIN the same dataset).
#'
#' @param overwrite FALSE, TRUE (every family), or a character vector of
#'   jobnames to replace -- same semantics as the old ingest_one_dataset().
#'   The actual delete_family() call (a write) is deferred to
#'   write_ingest_bundle(); this function only records WHICH families need
#'   it in `bundle$families[[jobname]]$overwrite`.
#' @param recompute_redundancy passed through to the same staleness/force
#'   logic run_all_redundancy() used to apply, now inlined here since it
#'   needs the read-only `con` this function already holds.
#' @param project_root see run_all_redundancy()'s matching doc -- used
#'   identically here for redundancy.R's own staleness mtime check.
compute_ingest_bundle <- function(config_path, results_dir, db_path,
                                    overwrite = FALSE, recompute_redundancy = FALSE,
                                    project_root = getwd()) {
  dataset_yaml <- yaml::read_yaml(config_path)
  dataset_id <- dataset_yaml$dataset$id
  stopifnot(!is.null(dataset_id))
  message("dataset: ", dataset_id)

  bundle <- list(dataset_id = dataset_id,
                 description = dataset_yaml$dataset$description %||% NA_character_,
                 sample_metadata = dataset_yaml$dataset[c("sample_metadata_path", "sample_id_col")],
                 feature_metadata = list(path = dataset_yaml$dataset$feature_metadata_path,
                                         id_col = dataset_yaml$dataset$feature_id_col,
                                         ensembl_col = dataset_yaml$dataset$ensembl_col %||% "ensembl",
                                         symbol_col = dataset_yaml$dataset$symbol_col %||% NA_character_),
                 families = list(), report = list(),
                 fits = NULL, factors = NULL, wgcna_modules = NULL,
                 factor_pairs = NULL, wgcna_fit_pairs = NULL, wgcna_module_pairs = NULL,
                 factor_stability_updates = NULL,
                 fit_redundancy = NULL, pattern_markers = NULL, redundancy_deletes = integer(0))

  family_dirs <- Filter(
    function(d) length(list.files(d, pattern = "^results_\\d+\\.RDS$")) > 0,
    list.dirs(results_dir, recursive = FALSE)
  )
  if (length(family_dirs) == 0) {
    message("  no job family subdirectories with results_*.RDS in ", results_dir, " -- skipping")
    return(invisible(bundle))
  }
  jobnames <- basename(family_dirs)
  message("  families found: ", paste(jobnames, collapse = ", "))

  # See ingest_one_dataset()'s original comment (now here) for the two
  # bundle_dir conventions this tries.
  bundle_dir <- sub("_results$", "", results_dir)
  if (!dir.exists(bundle_dir)) bundle_dir <- file.path(project_root, "slurm_bundles", dataset_id)

  art_dir <- artifacts_dir(db_path, dataset_id)
  dir.create(art_dir, recursive = TRUE, showWarnings = FALSE)

  con <- open_stability_db(db_path)   # READ-ONLY IN PRACTICE -- see this function's header
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  overwrite_requested <- function(jobname) {
    isTRUE(overwrite) || (is.character(overwrite) && jobname %in% overwrite)
  }

  local_counter <- 0L
  fits_rows <- list(); factors_rows <- list(); wgcna_modules_rows <- list()
  new_local_ids_by_method <- list()

  for (k in seq_along(jobnames)) {
    jobname <- jobnames[k]
    fam_dir <- family_dirs[k]
    cls <- classify_jobname(jobname)

    already <- family_already_ingested(con, dataset_id, jobname)
    if (already) {
      if (!overwrite_requested(jobname)) {
        message("  [", jobname, "] already ingested -- skipping (use overwrite = '", jobname, "' to replace)")
        bundle$report[[jobname]] <- "skipped"
        next
      }
      message("  [", jobname, "] already ingested -- will OVERWRITE at merge time")
    }

    params_path <- file.path(bundle_dir, paste0("_rslurm_", jobname), "params.RDS")
    if (!file.exists(params_path)) {
      message("  [", jobname, "] bundle params not found at ", params_path, " -- skipping this family")
      bundle$report[[jobname]] <- "skipped (no params.RDS)"
      next
    }
    params <- readRDS(params_path)

    result_files <- list.files(fam_dir, pattern = "^results_\\d+\\.RDS$", full.names = TRUE)
    task_ids <- as.integer(sub("^results_(\\d+)\\.RDS$", "\\1", basename(result_files)))

    bundle$families[[jobname]] <- list(family = cls$family, method = cls$method,
                                        overwrite = already, n_results = length(result_files),
                                        results_dir = normalizePath(results_dir))

    n_ok <- 0L; n_failed <- 0L; n_missing <- 0L

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

      local_counter <- local_counter + 1L
      local_id <- local_counter
      fits_rows[[length(fits_rows) + 1]] <- data.frame(
        local_id = local_id, dataset_id = dataset_id, method = ext$method, family = ext$family,
        jobname = jobname, rank = fit$rank, seed = fit$seed, alpha = fit$alpha, power = fit$power,
        rank_genes = fit$rank_genes, rank_subjects = fit$rank_subjects, rank_time = fit$rank_time,
        mse = fit$mse, n_factors = fit$n_factors, status = fit$status, converged = fit$converged,
        loadings_file = loadings_file, scores_file = scores_file, time_loadings_file = time_loadings_file,
        raw_result_file = raw_result_file,
        cogaps_diag_file = diag_file_val$cogaps_diag_file, pca_diag_file = diag_file_val$pca_diag_file,
        spca_diag_file = diag_file_val$spca_diag_file, ica_diag_file = diag_file_val$ica_diag_file,
        nmf_diag_file = diag_file_val$nmf_diag_file, cp_diag_file = diag_file_val$cp_diag_file,
        tucker_diag_file = diag_file_val$tucker_diag_file,
        stringsAsFactors = FALSE)

      if (fit$status == "ok") {
        n_ok <- n_ok + 1L
        new_local_ids_by_method[[ext$method]] <- c(new_local_ids_by_method[[ext$method]], -local_id)
        if (!is.null(ext$loadings)) {
          factors_rows[[length(factors_rows) + 1]] <- data.frame(
            fit_ref = -local_id, factor_index = seq_len(ncol(ext$loadings)),
            stability_cosine = NA_real_, stability_pearson = NA_real_, stability_spearman = NA_real_)
        }
        if (!is.null(ext$modules)) {
          wgcna_modules_rows[[length(wgcna_modules_rows) + 1]] <- cbind(fit_ref = -local_id, ext$modules)
          mod_ids <- setdiff(sort(unique(ext$modules$module)), 0L)
          if (length(mod_ids) > 0) {
            factors_rows[[length(factors_rows) + 1]] <- data.frame(
              fit_ref = -local_id, factor_index = mod_ids,
              stability_cosine = NA_real_, stability_pearson = NA_real_, stability_spearman = NA_real_)
          }
        }
      } else if (fit$status == "failed") {
        n_failed <- n_failed + 1L
      } else {
        n_missing <- n_missing + 1L
      }
    }
    bundle$report[[jobname]] <- sprintf("ingested (%d ok, %d failed, %d missing of %d tasks)",
                                         n_ok, n_failed, n_missing, nrow(params))
    message("  [", jobname, "] ", bundle$report[[jobname]])
  }

  bundle$fits <- if (length(fits_rows)) do.call(rbind, fits_rows) else NULL
  bundle$factors <- if (length(factors_rows)) do.call(rbind, factors_rows) else NULL
  bundle$wgcna_modules <- if (length(wgcna_modules_rows)) do.call(rbind, wgcna_modules_rows) else NULL

  resolve_local <- function(rel_paths) vapply(rel_paths, resolve_artifact, character(1), db_path = db_path)

  # ---- pairs (new x existing, per method) ----
  for (method in names(new_local_ids_by_method)) {
    new_ids <- new_local_ids_by_method[[method]]
    if (method == "wgcna") {
      existing <- DBI::dbGetQuery(con, "SELECT fit_id FROM fits WHERE dataset_id = ? AND method = 'wgcna' AND status = 'ok'",
                                   params = list(dataset_id))
      existing_ids <- existing$fit_id
      existing_mods <- if (length(existing_ids)) DBI::dbGetQuery(con, sprintf(
        "SELECT fit_id, gene, module FROM wgcna_modules WHERE fit_id IN (%s)", paste(existing_ids, collapse = ","))
      ) else data.frame(fit_id = integer(0), gene = character(0), module = integer(0))
      mod_list <- split(existing_mods[, c("gene", "module")], existing_mods$fit_id)
      if (!is.null(bundle$wgcna_modules)) {
        new_mod_list <- split(bundle$wgcna_modules[, c("gene", "module")], bundle$wgcna_modules$fit_ref)
        mod_list <- c(mod_list, new_mod_list)
      }
      all_ids <- c(existing_ids, if (!is.null(bundle$wgcna_modules)) unique(bundle$wgcna_modules$fit_ref) else integer(0))
      out <- compute_wgcna_pairs_from_universe(mod_list, all_ids, new_ids)
      if (!is.null(out)) {
        bundle$wgcna_fit_pairs <- rbind(bundle$wgcna_fit_pairs, out$fit_pairs)
        bundle$wgcna_module_pairs <- rbind(bundle$wgcna_module_pairs, out$module_pairs)
      }
    } else {
      existing <- DBI::dbGetQuery(con,
        "SELECT fit_id AS id, rank, loadings_file FROM fits
         WHERE dataset_id = ? AND method = ? AND status = 'ok' AND loadings_file IS NOT NULL",
        params = list(dataset_id, method))
      existing_u <- if (nrow(existing)) {
        data.frame(id = existing$id, rank = existing$rank, loadings_file_abs = resolve_local(existing$loadings_file))
      } else NULL

      new_sub <- bundle$fits[bundle$fits$method == method & bundle$fits$status == "ok" &
                                !is.na(bundle$fits$loadings_file), , drop = FALSE]
      new_u <- if (nrow(new_sub)) {
        data.frame(id = -new_sub$local_id, rank = new_sub$rank, loadings_file_abs = resolve_local(new_sub$loadings_file))
      } else NULL

      universe <- rbind(existing_u, new_u)
      rows <- compute_factor_pairs_from_universe(universe, new_ids)
      if (!is.null(rows)) {
        message("  [pairs:", method, "] computed ", nrow(rows), " factor-pair rows for ", length(new_ids), " new fits")
        bundle$factor_pairs <- rbind(bundle$factor_pairs, rows)
      }

      existing_pairs <- DBI::dbGetQuery(con,
        "SELECT fp.fit_a, fp.fit_b, fp.factor_a, fp.factor_b, fp.cosine, fp.pearson, fp.spearman, fp.matched, fp.same_rank
         FROM factor_pairs fp JOIN fits fa ON fa.fit_id = fp.fit_a
         WHERE fa.dataset_id = ? AND fa.method = ?", params = list(dataset_id, method))
      combined_pairs <- rbind(existing_pairs, rows)
      agg <- if (!is.null(combined_pairs)) update_factor_stability_from_pairs(combined_pairs) else NULL
      if (!is.null(agg)) bundle$factor_stability_updates <- rbind(bundle$factor_stability_updates, agg)
    }
  }

  # ---- redundancy (representative fits, new x existing universe) ----
  this_file <- file.path(project_root, "R/lib/ingest/redundancy.R")
  redundancy_src_mtime <- if (file.exists(this_file)) file.mtime(this_file) else Sys.time()

  for (method in c("nmf", "cogaps", "spca", "ica")) {
    existing_all <- DBI::dbGetQuery(con,
      "SELECT fit_id AS id, method, family, rank, alpha, mse, status FROM fits WHERE dataset_id = ? AND method = ?",
      params = list(dataset_id, method))
    new_all <- if (!is.null(bundle$fits)) bundle$fits[bundle$fits$method == method,
      c("local_id", "method", "family", "rank", "alpha", "mse", "status"), drop = FALSE] else NULL
    if (!is.null(new_all) && nrow(new_all)) {
      new_all$id <- -new_all$local_id
      new_all$local_id <- NULL
    } else new_all <- NULL
    universe_fits <- rbind(existing_all, new_all)
    if (is.null(universe_fits) || nrow(universe_fits) == 0) next
    rep_ids <- select_representative_ids_from_universe(universe_fits, method)

    for (rid in rep_ids) {
      if (rid < 0) {
        row <- bundle$fits[bundle$fits$local_id == -rid, , drop = FALSE]
        if (nrow(row) == 0 || is.na(row$loadings_file)) next
        loadings_path <- resolve_artifact(row$loadings_file, db_path)
        raw_path <- if (is.na(row$raw_result_file)) NA_character_ else resolve_artifact(row$raw_result_file, db_path)
        out <- run_redundancy_for_fit(method, loadings_path, raw_path)   # new fit -- always compute, no cache possible
        if (is.null(out)) next
        fname <- sprintf("%s_local%d_redundancy.rds", method, -rid)
        saveRDS(out$summary$matrix, file.path(art_dir, fname))
      } else {
        existing_cache <- DBI::dbGetQuery(con, "SELECT redundancy_file FROM fit_redundancy WHERE fit_id = ?", params = list(rid))
        have_cache <- nrow(existing_cache) > 0 && !is.na(existing_cache$redundancy_file[1])
        stale <- FALSE
        if (!recompute_redundancy && have_cache) {
          art_path <- resolve_artifact(existing_cache$redundancy_file[1], db_path)
          if (file.exists(art_path) && file.mtime(art_path) < redundancy_src_mtime) {
            message("  fit_redundancy for fit_id ", rid, " is STALE (redundancy.R edited since) -- recomputing")
            stale <- TRUE
          }
        }
        if (!recompute_redundancy && !stale && have_cache) next

        f <- DBI::dbGetQuery(con, "SELECT * FROM fits WHERE fit_id = ?", params = list(rid))
        if (nrow(f) == 0 || is.na(f$loadings_file)) next
        out <- run_redundancy_for_fit(method, resolve_artifact(f$loadings_file, db_path),
                                       if (is.na(f$raw_result_file)) NA_character_ else resolve_artifact(f$raw_result_file, db_path))
        if (is.null(out)) next
        fname <- sprintf("%s_fit%d_redundancy.rds", method, rid)
        saveRDS(out$summary$matrix, file.path(art_dir, fname))
        if (stale || have_cache) bundle$redundancy_deletes <- c(bundle$redundancy_deletes, rid)
      }
      bundle$pattern_markers <- rbind(bundle$pattern_markers, cbind(fit_ref = rid, out$markers[, c("factor_index", "gene", "score")]))
      bundle$fit_redundancy <- rbind(bundle$fit_redundancy, data.frame(
        fit_ref = rid, max_offdiag_cosine = out$summary$max_offdiag, median_offdiag_cosine = out$summary$median_offdiag,
        n_factors_with_no_markers = out$summary$n_factors_with_no_markers,
        redundancy_file = file.path("stability_artifacts", dataset_id, fname)))
    }
  }

  invisible(bundle)
}

#' Apply one dataset's compute_ingest_bundle() output to the shared DB.
#' Wrapped in ONE transaction -- confirmed the fix for a real ingest_core
#' seff report showing only 28.42% CPU efficiency over a 1h19m single-core
#' run (i.e. ~56 minutes of wall-clock spent blocked on per-statement
#' fsyncs, not computing) -- same "death by a thousand fsyncs" pattern
#' already diagnosed and fixed for driver_grid. on.exit()'s
#' rollback-if-not-committed guard ensures an error partway through never
#' leaves an open transaction for the NEXT dataset's writes (the merge
#' script loops many datasets' bundles against the same connection) to
#' land inside.
write_ingest_bundle <- function(con, db_path, bundle) {
  dataset_id <- bundle$dataset_id
  if (length(bundle$families) == 0 && is.null(bundle$fits)) return(invisible(dataset_id))

  DBI::dbExecute(con, "BEGIN")
  committed <- FALSE
  on.exit(if (!committed) DBI::dbExecute(con, "ROLLBACK"), add = TRUE)

  ensure_dataset(con, dataset_id, bundle$description)
  register_metadata_source(con, dataset_id, "sample", bundle$sample_metadata$sample_metadata_path,
                            bundle$sample_metadata$sample_id_col)
  register_metadata_source(con, dataset_id, "feature", bundle$feature_metadata$path, bundle$feature_metadata$id_col,
                            ensembl_col = bundle$feature_metadata$ensembl_col, symbol_col = bundle$feature_metadata$symbol_col)

  for (jobname in names(bundle$families)) {
    if (isTRUE(bundle$families[[jobname]]$overwrite)) {
      message("  [", jobname, "] OVERWRITING")
      delete_family(con, db_path, dataset_id, jobname)
    }
  }

  local_to_fit_id <- integer(0)   # names = local_id (character), values = real fit_id
  if (!is.null(bundle$fits)) {
    fits <- bundle$fits[order(bundle$fits$local_id), , drop = FALSE]
    for (i in seq_len(nrow(fits))) {
      r <- fits[i, ]
      DBI::dbExecute(con,
        "INSERT INTO fits (dataset_id, method, family, jobname, rank, seed, alpha, power,
                           rank_genes, rank_subjects, rank_time,
                           mse, n_factors, status, converged, loadings_file, scores_file, time_loadings_file, raw_result_file,
                           cogaps_diag_file, pca_diag_file, spca_diag_file, ica_diag_file, nmf_diag_file, cp_diag_file, tucker_diag_file)
         VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        params = list(r$dataset_id, r$method, r$family, r$jobname, r$rank, r$seed, r$alpha, r$power,
                      r$rank_genes, r$rank_subjects, r$rank_time, r$mse, r$n_factors, r$status, r$converged,
                      r$loadings_file, r$scores_file, r$time_loadings_file, r$raw_result_file,
                      r$cogaps_diag_file, r$pca_diag_file, r$spca_diag_file, r$ica_diag_file,
                      r$nmf_diag_file, r$cp_diag_file, r$tucker_diag_file))
      fit_id <- DBI::dbGetQuery(con, "SELECT last_insert_rowid() AS id")$id
      local_to_fit_id[as.character(r$local_id)] <- fit_id
    }
  }
  # Remap: positive ids (already-merged, from an earlier run) pass through
  # unchanged; negative ids (-local_id, this bundle's own new fits) map to
  # the real fit_id just assigned above.
  remap <- function(ids) ifelse(ids < 0, local_to_fit_id[as.character(-ids)], ids)

  if (!is.null(bundle$factors)) {
    df <- bundle$factors
    df$fit_id <- remap(df$fit_ref); df$fit_ref <- NULL
    DBI::dbWriteTable(con, "factors", df, append = TRUE)
  }
  if (!is.null(bundle$wgcna_modules)) {
    df <- bundle$wgcna_modules
    df$fit_id <- remap(df$fit_ref); df$fit_ref <- NULL
    DBI::dbWriteTable(con, "wgcna_modules", df, append = TRUE)
  }
  if (!is.null(bundle$factor_pairs)) {
    df <- bundle$factor_pairs
    df$fit_a <- remap(df$fit_a); df$fit_b <- remap(df$fit_b)
    for (start in seq(1, nrow(df), by = 200)) {
      DBI::dbWriteTable(con, "factor_pairs", df[start:min(start + 199, nrow(df)), , drop = FALSE], append = TRUE)
    }
  }
  if (!is.null(bundle$wgcna_fit_pairs)) {
    df <- bundle$wgcna_fit_pairs
    df$fit_a <- remap(df$fit_a); df$fit_b <- remap(df$fit_b)
    DBI::dbWriteTable(con, "wgcna_fit_pairs", df, append = TRUE)
  }
  if (!is.null(bundle$wgcna_module_pairs)) {
    df <- bundle$wgcna_module_pairs
    df$fit_a <- remap(df$fit_a); df$fit_b <- remap(df$fit_b)
    DBI::dbWriteTable(con, "wgcna_module_pairs", df, append = TRUE)
  }
  if (!is.null(bundle$factor_stability_updates)) {
    agg <- bundle$factor_stability_updates
    agg$fit_id <- remap(agg$fit_id)
    for (r in seq_len(nrow(agg))) {
      DBI::dbExecute(con,
        "UPDATE factors SET stability_cosine = ?, stability_pearson = ?, stability_spearman = ?
         WHERE fit_id = ? AND factor_index = ?",
        params = list(agg$cosine[r], agg$pearson[r], agg$spearman[r], agg$fit_id[r], agg$factor_index[r]))
    }
  }
  for (rid in bundle$redundancy_deletes) {
    DBI::dbExecute(con, "DELETE FROM pattern_markers WHERE fit_id = ?", params = list(rid))
    DBI::dbExecute(con, "DELETE FROM fit_redundancy WHERE fit_id = ?", params = list(rid))
  }
  if (!is.null(bundle$pattern_markers)) {
    df <- bundle$pattern_markers
    df$fit_id <- remap(df$fit_ref); df$fit_ref <- NULL
    DBI::dbWriteTable(con, "pattern_markers", df, append = TRUE)
  }
  if (!is.null(bundle$fit_redundancy)) {
    df <- bundle$fit_redundancy
    df$fit_id <- remap(df$fit_ref); df$fit_ref <- NULL
    DBI::dbWriteTable(con, "fit_redundancy", df, append = TRUE)
  }

  for (jobname in names(bundle$families)) {
    fam <- bundle$families[[jobname]]
    record_ingest(con, dataset_id, jobname, fam$family, fam$method, n_results = fam$n_results, results_dir = fam$results_dir)
  }

  DBI::dbExecute(con, "COMMIT")
  committed <- TRUE
  invisible(dataset_id)
}

#' Ingest ONE dataset's job-family results (+ redundancy) into `con`.
#' Thin wrapper -- compute_ingest_bundle() + write_ingest_bundle() --
#' kept for the single-dataset direct CLI (R/ingest_results.R) and any
#' other caller that just wants "ingest this one dataset now," without
#' the array-job/merge-script split R/ingest_jobs/ingest_core_job.R uses
#' for the multi-dataset production path.
#'
#' NOTE: deliberately does NOT call cache_dataset_matrix() -- that function
#' runs `preprocessing_script`/`load_input_matrix()`, which `source()`s a
#' project-relative file path and reads the raw parquet trio. Neither is
#' reachable from inside the ingest_core slurm job's container (only
#' slurm_bundles/ingest/ is bind-mounted). cache_dataset_matrix() must be
#' called separately first (a plain `Rscript` invocation, submitted via
#' `srun`, never inside ingest_core's container) -- see
#' R/create_ingest_slurm_bundle.R's --stage core.
#'
#' @param overwrite FALSE, TRUE (every family), or a character vector of
#'   jobnames to replace.
#' @param recompute_redundancy passed straight through to
#'   compute_ingest_bundle()'s redundancy staleness/force logic.
#' @param run_pattern_drivers whether to also run the
#'   projectR::projectionDriveR()-based pattern-driver pass (see
#'   R/lib/ingest/driver.R) inline. Default TRUE for R/ingest_results.R's
#'   direct (non-slurm) CLI use, where projectR is just whatever's in that
#'   R session. FALSE for the multi-dataset compute-job path -- ingest_core
#'   runs inside a container that does NOT have projectR installed (a
#'   different image than the one it needs -- see config/
#'   ingest_slurm_config.yml's `driver:` entry); the pattern-driver pass is
#'   staged as its own `driver_grid` job family instead, during
#'   --stage enrichment (see R/ingest_jobs/driver_job.R).
#' @param project_root absolute project-root path, used (a) as the
#'   fallback location for this dataset's rslurm bundle when it isn't
#'   colocated with results_dir -- see bundle_dir below -- and (b) to
#'   locate redundancy.R's own source file for staleness detection.
#'   Defaults to getwd(), correct for R/ingest_results.R's direct CLI use
#'   (run from the project root).
ingest_one_dataset <- function(con, config_path, results_dir, db_path,
                                overwrite = FALSE, recompute_redundancy = FALSE,
                                run_pattern_drivers = TRUE, project_root = getwd()) {
  bundle <- compute_ingest_bundle(config_path, results_dir, db_path,
                                   overwrite = overwrite, recompute_redundancy = recompute_redundancy,
                                   project_root = project_root)
  write_ingest_bundle(con, db_path, bundle)
  dataset_id <- bundle$dataset_id

  # Differential feature identification (projectR::projectionDriveR(), see
  # R/lib/ingest/driver.R) -- batch pass over representative fits x
  # 2-4-level categorical sample-metadata columns. Silently no-ops if the
  # dataset's matrix isn't cached yet or has no small categorical columns
  # -- never blocks the rest of ingest. Skipped entirely when
  # run_pattern_drivers = FALSE (see this function's @param doc) --
  # staged as its own driver_grid job family instead in that case.
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

# Stage-2 ingest: load a dataset's rslurm job outputs into the stability
# SQLite DB + external artifacts, computing all pairwise stability metrics.
#
# Fully decoupled from the Shiny app (which just reads whatever the DB
# currently contains) and ADDITIVE:
#   - a new dataset ingests alongside existing ones in the same DB;
#   - job families (subdirectories of <results_dir> containing
#     results_*.RDS) not yet in the DB for this dataset are added; families
#     already present are skipped with a message;
#   - `--overwrite` (bare = every family present in <results_dir>;
#     or `--overwrite jobname1,jobname2` = just those) deletes and
#     re-ingests -- plain replacement, never row-level updating.
#
# Usage (CLI):
#   Rscript R/ingest_results.R <dataset_config.yml> <results_dir> <db_path> [--overwrite [jobname,...]]
#
# Usage (interactive): set `ingest_config_path`, `ingest_results_dir`,
# `ingest_db_path` (and optionally `ingest_overwrite`, a character vector
# of jobnames or TRUE for all) then source this file.
#
# Example workflow (the acceptance case this was built around):
#   Rscript R/ingest_results.R config/GSE110487_config.yml results/GSE110487_results results/stability.sqlite
#   # ... later, wto_grid/ lands in results/GSE110487_results/ ...
#   Rscript R/ingest_results.R config/GSE110487_config.yml results/GSE110487_results results/stability.sqlite
#   # -> only wto_grid is added; everything else reported as skipped
#   # ... nmf re-run on the cluster ...
#   Rscript R/ingest_results.R ... --overwrite nmf_grid,nmf_maskcv

library(here)
library(yaml)
source(here("R/lib/ingest/db.R"))
source(here("R/lib/ingest/similarity.R"))
source(here("R/lib/ingest/extract.R"))
source(here("R/lib/ingest/pairs.R"))

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
      TRUE   # bare --overwrite = everything present in results_dir
    }
  } else {
    FALSE
  }
}

if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

dataset_yaml <- yaml::read_yaml(ingest_config_path)
dataset_id   <- dataset_yaml$dataset$id
stopifnot(!is.null(dataset_id))
description  <- dataset_yaml$dataset$description %||% NA_character_

## ---- discovery --------------------------------------------------------------

family_dirs <- Filter(
  function(d) length(list.files(d, pattern = "^results_\\d+\\.RDS$")) > 0,
  list.dirs(ingest_results_dir, recursive = FALSE)
)
if (length(family_dirs) == 0) stop("no job family subdirectories with results_*.RDS found in ", ingest_results_dir)
jobnames <- basename(family_dirs)

message("dataset: ", dataset_id)
message("families found in ", ingest_results_dir, ": ", paste(jobnames, collapse = ", "))

# NOTE: no on.exit() here -- at the top level of a source()d script it
# would fire as soon as its own expression completes, closing the
# connection immediately. Disconnect happens explicitly at the end.
con <- open_stability_db(ingest_db_path)
ensure_dataset(con, dataset_id, description)
art_dir <- artifacts_dir(ingest_db_path, dataset_id)
dir.create(art_dir, recursive = TRUE, showWarnings = FALSE)

overwrite_requested <- function(jobname) {
  isTRUE(ingest_overwrite) ||
    (is.character(ingest_overwrite) && jobname %in% ingest_overwrite)
}

## ---- per-family ingest ------------------------------------------------------

# track fits added this run, per method, for incremental pair computation
new_fits_by_method <- list()
report <- list()

for (k in seq_along(jobnames)) {
  jobname <- jobnames[k]
  fam_dir <- family_dirs[k]
  cls <- classify_jobname(jobname)

  if (family_already_ingested(con, dataset_id, jobname)) {
    if (overwrite_requested(jobname)) {
      message("[", jobname, "] already ingested -- OVERWRITING")
      delete_family(con, ingest_db_path, dataset_id, jobname)
    } else {
      message("[", jobname, "] already ingested -- skipping (use --overwrite ", jobname, " to replace)")
      report[[jobname]] <- "skipped"
      next
    }
  }

  params_path <- here("slurm_bundles", dataset_id, paste0("_rslurm_", jobname), "params.RDS")
  if (!file.exists(params_path)) {
    stop("[", jobname, "] bundle params not found at ", params_path,
         " -- the slurm bundle used to generate these results must be present ",
         "to align array indices to parameters")
  }
  params <- readRDS(params_path)

  result_files <- list.files(fam_dir, pattern = "^results_\\d+\\.RDS$", full.names = TRUE)
  task_ids <- as.integer(sub("^results_(\\d+)\\.RDS$", "\\1", basename(result_files)))

  n_ok <- 0L; n_failed <- 0L; n_missing <- 0L
  new_ids <- integer(0)

  DBI::dbExecute(con, "BEGIN")
  for (task in seq_len(nrow(params)) - 1L) {   # array task ids are 0-based; row = task + 1
    params_row <- params[task + 1L, , drop = FALSE]
    f <- result_files[match(task, task_ids)]
    result <- if (!is.na(f)) {
      x <- readRDS(f)
      # rslurm wraps each task's chunk in a list (nchunk = 1 here)
      if (is.list(x) && length(x) == 1 && is.null(names(x))) x[[1]] else x
    } else {
      NULL
    }

    ext <- extract_result(jobname, result, params_row)
    fit <- ext$fit

    # loadings_file is stored RELATIVE TO THE DB's directory (see
    # resolve_artifact() in R/lib/ingest/db.R), so DB + artifacts stay
    # portable as a unit regardless of where ingest/app are launched from
    loadings_file <- NA_character_
    if (!is.null(ext$loadings)) {
      fname <- sprintf("%s_task%03d_loadings.rds", jobname, task)
      saveRDS(ext$loadings, file.path(art_dir, fname))
      loadings_file <- file.path("stability_artifacts", dataset_id, fname)
    } else if (!is.null(ext$edges)) {
      fname <- sprintf("%s_task%03d_edges.parquet", jobname, task)
      arrow::write_parquet(ext$edges, file.path(art_dir, fname))
      loadings_file <- file.path("stability_artifacts", dataset_id, fname)
    }

    DBI::dbExecute(con,
      "INSERT INTO fits (dataset_id, method, family, jobname, rank, seed, alpha,
                         power, n_boot, delta, mse, n_factors, status, loadings_file)
       VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
      params = list(dataset_id, ext$method, ext$family, jobname,
                    fit$rank, fit$seed, fit$alpha, fit$power, fit$n_boot, fit$delta,
                    fit$mse, fit$n_factors, fit$status, loadings_file))
    fit_id <- DBI::dbGetQuery(con, "SELECT last_insert_rowid() AS id")$id

    if (fit$status == "ok") {
      n_ok <- n_ok + 1L
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
                n_results = length(result_files), results_dir = normalizePath(ingest_results_dir))
  DBI::dbExecute(con, "COMMIT")

  if (cls$family != "maskcv" && length(new_ids) > 0) {
    key <- cls$method
    new_fits_by_method[[key]] <- c(new_fits_by_method[[key]], new_ids)
  }
  report[[jobname]] <- sprintf("ingested (%d ok, %d failed, %d missing of %d tasks)",
                                n_ok, n_failed, n_missing, nrow(params))
  message("[", jobname, "] ", report[[jobname]])
}

## ---- pairwise stability (incremental: new fits x all same-method fits) ------

for (method in names(new_fits_by_method)) {
  ids <- new_fits_by_method[[method]]
  message("[pairs:", method, "] computing similarities for ", length(ids), " new fits ...")
  if (method == "wgcna") {
    compute_wgcna_pairs(con, dataset_id, ids)
  } else if (method == "wto") {
    compute_wto_pairs(con, ingest_db_path, dataset_id, ids)
  } else {
    n <- compute_factor_pairs(con, ingest_db_path, dataset_id, method, ids)
    message("[pairs:", method, "] wrote ", n, " factor-pair rows")
    update_factor_stability(con, dataset_id, method)
  }
}

## ---- report ------------------------------------------------------------------

message("\n===== ingest report: ", dataset_id, " -> ", ingest_db_path, " =====")
for (jn in names(report)) message(sprintf("  %-22s %s", jn, report[[jn]]))
counts <- DBI::dbGetQuery(con, "
  SELECT 'fits' AS tbl, COUNT(*) AS n FROM fits WHERE dataset_id = :d
  UNION ALL SELECT 'factors', COUNT(*) FROM factors f JOIN fits ft ON ft.fit_id = f.fit_id WHERE ft.dataset_id = :d
  UNION ALL SELECT 'factor_pairs', COUNT(*) FROM factor_pairs fp JOIN fits ft ON ft.fit_id = fp.fit_a WHERE ft.dataset_id = :d
  UNION ALL SELECT 'maskcv_results', COUNT(*) FROM maskcv_results WHERE dataset_id = :d
  UNION ALL SELECT 'wgcna_fit_pairs', COUNT(*) FROM wgcna_fit_pairs wp JOIN fits ft ON ft.fit_id = wp.fit_a WHERE ft.dataset_id = :d
  UNION ALL SELECT 'wto_fit_pairs', COUNT(*) FROM wto_fit_pairs wp JOIN fits ft ON ft.fit_id = wp.fit_a WHERE ft.dataset_id = :d",
  params = list(d = dataset_id))
for (r in seq_len(nrow(counts))) message(sprintf("  %-22s %d rows", counts$tbl[r], counts$n[r]))

DBI::dbDisconnect(con)

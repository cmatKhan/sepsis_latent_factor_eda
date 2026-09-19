# Phase-2 results -> DB: reads fgsea_grid/gprofiler_grid's plain
# results_*.RDS files (written by R/ingest_jobs/fgsea_job.R /
# gprofiler_job.R, staged by R/create_ingest_slurm_bundle.R --stage
# enrichment) and writes them into enrichment_cache/enrichment_queried --
# the SAME tables/shape the app and R/lib/ingest/enrichment.R's
# gprofiler2-only --run-enrichment path already use, so nothing downstream
# (app, "Compare methods" screen) needs to know which route produced a
# given row. Never touches a compute node -- this runs on the login node
# after both job families have finished.
#
# Usage:
#   Rscript R/ingest_enrichment_results.R --bundle-dir slurm_bundles/ingest --db results/stability.sqlite

library(here); library(optparse); library(DBI)
source(here("R/lib/ingest/db.R"))
`%||%` <- function(a, b) if (is.null(a)) b else a

option_list <- list(
  make_option("--bundle-dir", type = "character", help = "e.g. slurm_bundles/ingest"),
  make_option("--db", type = "character", default = "results/stability.sqlite")
)
opt <- parse_args(OptionParser(option_list = option_list))
con <- open_stability_db(opt$db)

get_factor_id <- function(fit_id, factor_index) {
  DBI::dbGetQuery(con, "SELECT factor_id FROM factors WHERE fit_id = ? AND factor_index = ?",
                   params = list(fit_id, factor_index))$factor_id
}

store_enrichment <- function(factor_id, query_type, direction, rows) {
  DBI::dbExecute(con,
    "INSERT OR REPLACE INTO enrichment_queried (factor_id, query_type, direction, queried_at)
     VALUES (?, ?, ?, datetime('now'))",
    params = list(factor_id, query_type, direction))
  if (!is.null(rows) && nrow(rows) > 0) {
    DBI::dbWriteTable(con, "enrichment_cache", rows, append = TRUE)
  }
}

## ---- fgsea_grid --------------------------------------------------------------

fam_dir <- file.path(opt$`bundle-dir`, "_rslurm_fgsea_grid")
if (dir.exists(fam_dir)) {
  result_files <- list.files(fam_dir, pattern = "^results_\\d+\\.RDS$", full.names = TRUE)
  message("fgsea_grid: ", length(result_files), " result file(s)")
  n_entries <- 0
  for (f in result_files) {
    x <- readRDS(f)
    # Each results_<task>.RDS holds either ONE run_fgsea_job() return value
    # (task processed a single row -- `submit_job_family()`'s max_array_size
    # = NULL / 1 case) or a plain list of several (batched -- see
    # submit_job_family()'s `max_array_size` doc). `x$fit_id` is only
    # present on the former (a single row's own named-list result never
    # has a top-level "fit_id" name if it's actually the outer per-task
    # list instead).
    entries <- if (!is.null(x$fit_id)) list(x) else x
    n_entries <- n_entries + length(entries)

    for (x in entries) {
      fit_id <- x$fit_id
      if (is.null(fit_id)) next

      for (g in x$gsea %||% list()) {
        res <- g$result
        factor_id <- get_factor_id(fit_id, g$factor_index)
        if (length(factor_id) != 1) next
        rows <- if (!is.null(res) && nrow(res) > 0) {
          sig <- res[res$padj < 0.05, ]
          if (nrow(sig) == 0) NULL else data.frame(
            factor_id = factor_id, query_type = "fgsea",
            direction = ifelse(sig$NES > 0, "pos", "neg"),
            source = "MSigDB", term_id = sig$pathway, term_name = sig$pathway,
            p_value = sig$padj, intersection_size = lengths(sig$leadingEdge), term_size = sig$size,
            query_size = NA_integer_,
            genes = vapply(sig$leadingEdge, paste, character(1), collapse = ","),
            queried_at = as.character(Sys.time())
          )
        } else NULL
        store_enrichment(factor_id, "fgsea", "pos", rows)
      }

      if (!is.null(x$fora)) {
        # x$fora is named by factor_index (a character), one fora() result per CoGAPS pattern
        for (fac_idx_chr in names(x$fora)) {
          res <- x$fora[[fac_idx_chr]]
          factor_id <- get_factor_id(fit_id, as.integer(fac_idx_chr))
          if (length(factor_id) != 1) next
          rows <- if (!is.null(res) && nrow(res) > 0) {
            sig <- res[res$padj < 0.05, ]
            if (nrow(sig) == 0) NULL else data.frame(
              factor_id = factor_id, query_type = "cogaps_fora", direction = "pos",
              source = "MSigDB", term_id = sig$pathway, term_name = sig$pathway,
              p_value = sig$padj, intersection_size = sig$overlap, term_size = sig$size,
              query_size = NA_integer_, genes = NA_character_, queried_at = as.character(Sys.time())
            )
          } else NULL
          store_enrichment(factor_id, "cogaps_fora", "pos", rows)
        }
      }
    }
  }
  message("fgsea_grid: ", n_entries, " fit-level result(s) across those file(s)")
} else {
  message("fgsea_grid: no results dir found, skipping")
}

## ---- gprofiler_grid ------------------------------------------------------------

fam_dir <- file.path(opt$`bundle-dir`, "_rslurm_gprofiler_grid")
if (dir.exists(fam_dir)) {
  result_files <- list.files(fam_dir, pattern = "^results_\\d+\\.RDS$", full.names = TRUE)
  message("gprofiler_grid: ", length(result_files), " result file(s)")
  for (f in result_files) {
    x <- readRDS(f)
    if (is.list(x) && length(x) == 1 && is.null(names(x))) x <- x[[1]]
    # run_gprofiler_job() returns a LIST of per-row results (one job, many rows)
    entries <- if (!is.null(x$fit_id)) list(x) else x
    for (e in entries) {
      factor_id <- get_factor_id(e$fit_id, e$factor_index)
      if (length(factor_id) != 1) next
      for (qt in c("ora", "gsea")) {
        gost_result <- e[[qt]]
        rows <- if (!is.null(gost_result) && !is.null(gost_result$result) && nrow(gost_result$result) > 0) {
          r <- gost_result$result
          data.frame(
            factor_id = factor_id, query_type = qt, direction = e$direction,
            source = r$source, term_id = r$term_id, term_name = r$term_name,
            p_value = r$p_value, intersection_size = r$intersection_size, term_size = r$term_size,
            query_size = if (!is.null(r$query_size)) r$query_size else NA_integer_,
            genes = if (!is.null(r$intersection)) as.character(r$intersection) else NA_character_,
            queried_at = as.character(Sys.time())
          )
        } else NULL
        store_enrichment(factor_id, qt, e$direction, rows)
      }
    }
  }
} else {
  message("gprofiler_grid: no results dir found, skipping")
}

DBI::dbDisconnect(con)

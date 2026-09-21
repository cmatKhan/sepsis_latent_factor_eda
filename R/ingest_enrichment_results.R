# Phase-2 results -> DB: reads fgsea_grid's plain results_*.RDS files
# (written by R/ingest_jobs/fgsea_job.R, staged by
# R/create_ingest_slurm_bundle.R --stage enrichment) and writes them into
# enrichment_cache/enrichment_queried -- the SAME tables/shape the app's
# own live, on-demand gprofiler2 queries (app/app.R) already use, so
# nothing downstream (app, "Compare methods" screen) needs to know which
# route produced a given row. Never touches a compute node -- this runs on
# the login node after the job family has finished.
#
# gprofiler_grid (a separate job family that made live, rate-limited
# gprofiler2::gost() API calls) and R/lib/ingest/enrichment.R's legacy
# whole-DB --run-enrichment CLI path (same API, same scaling problem) were
# both retired 2026-09-19: fgsea_job.R now runs the equivalent ORA/GSEA
# passes locally via fora()/fgsea() against msigdbr collections instead
# (see that file's header) -- `x$local_gsea`/`x$local_ora` below are that
# replacement's output, written under the SAME `query_type` values
# ("gsea"/"ora") the app's live queries already use, distinguished by the
# `source` column (GO:BP/GO:MF/KEGG/REAC/WP/HALLMARK) exactly like
# gprofiler2's own `result$source` field did.
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
          # `res$padj < 0.05` is NA (not FALSE) for any pathway fgsea
          # couldn't compute a p-value for -- fgsea's own "unbalanced
          # (positive and negative) gene-level statistic values" case,
          # confirmed to occur on real data. Indexing a data.frame with a
          # logical vector containing NA does NOT drop those rows, it
          # inserts a literal all-NA row per NA (every column, including
          # factor_id below via recycling into that row), which crashes the
          # INSERT on enrichment_cache's NOT NULL factor_id -- must exclude
          # NA explicitly. Same fix applies to every `sig <- res[...]` line
          # below (fora/local_gsea/local_ora).
          sig <- res[!is.na(res$padj) & res$padj < 0.05, ]
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
            sig <- res[!is.na(res$padj) & res$padj < 0.05, ]   # see x$gsea's comment above
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

      # ---- local GSEA/ORA replacing gprofiler_grid (see this file's header
      # and R/ingest_jobs/fgsea_job.R's header) -- written under the same
      # query_type values ("gsea"/"ora") the app's live gprofiler2 queries
      # already use, one `source` per msigdbr collection.
      for (g in x$local_gsea %||% list()) {
        res <- g$result
        factor_id <- get_factor_id(fit_id, g$factor_index)
        if (length(factor_id) != 1) next
        rows <- if (!is.null(res) && nrow(res) > 0) {
          sig <- res[!is.na(res$padj) & res$padj < 0.05, ]   # see x$gsea's comment above
          if (nrow(sig) == 0) NULL else data.frame(
            factor_id = factor_id, query_type = "gsea",
            direction = ifelse(sig$NES > 0, "pos", "neg"),
            source = g$source, term_id = sig$pathway, term_name = sig$pathway,
            p_value = sig$padj, intersection_size = lengths(sig$leadingEdge), term_size = sig$size,
            query_size = NA_integer_,
            genes = vapply(sig$leadingEdge, paste, character(1), collapse = ","),
            queried_at = as.character(Sys.time())
          )
        } else NULL
        # GSEA rows span both directions (NES sign) in one fgsea() call --
        # mark the query as "queried" for both, same as gprofiler_grid's
        # single ordered_query = TRUE call used to (it also produced
        # mixed-sign NES results under one nominal `direction`).
        store_enrichment(factor_id, "gsea", "pos", if (!is.null(rows)) rows[rows$direction == "pos", ] else NULL)
        store_enrichment(factor_id, "gsea", "neg", if (!is.null(rows)) rows[rows$direction == "neg", ] else NULL)
      }

      for (o in x$local_ora %||% list()) {
        res <- o$result
        factor_id <- get_factor_id(fit_id, o$factor_index)
        if (length(factor_id) != 1) next
        rows <- if (!is.null(res) && nrow(res) > 0) {
          sig <- res[!is.na(res$padj) & res$padj < 0.05, ]   # see x$gsea's comment above
          if (nrow(sig) == 0) NULL else data.frame(
            factor_id = factor_id, query_type = "ora", direction = o$direction,
            source = o$source, term_id = sig$pathway, term_name = sig$pathway,
            p_value = sig$padj, intersection_size = sig$overlap, term_size = sig$size,
            query_size = o$n_genes %||% NA_integer_,
            # fora() output has no `leadingEdge` column (that's GSEA-only,
            # used above for x$gsea/x$local_gsea) -- fora()'s genes column
            # is `overlapGenes`. Using leadingEdge here crashed the whole
            # script on the first ORA hit (character(0) from an always-NULL
            # column against every other length-N column), silently
            # aborting every ingest attempt before most rows were written.
            genes = vapply(sig$overlapGenes, paste, character(1), collapse = ","),
            queried_at = as.character(Sys.time())
          )
        } else NULL
        store_enrichment(factor_id, "ora", o$direction, rows)
      }
    }
  }
  message("fgsea_grid: ", n_entries, " fit-level result(s) across those file(s)")
} else {
  message("fgsea_grid: no results dir found, skipping")
}

DBI::dbDisconnect(con)

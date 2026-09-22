# Phase-2 results -> DB: reads projectr_within_grid/projectr_cross_grid's
# plain results_*.RDS files (written by R/ingest_jobs/projectr_job.R) and
# writes them into the `projections` table + stability_artifacts/.
#
# Usage:
#   Rscript R/ingest_projectr_results.R --bundle-dir slurm_bundles/ingest --db results/stability.sqlite

library(here); library(optparse); library(DBI)
source(here("R/lib/ingest/db.R"))

option_list <- list(
  make_option("--bundle-dir", type = "character", help = "e.g. slurm_bundles/ingest"),
  make_option("--db", type = "character", default = "results/stability.sqlite")
)
opt <- parse_args(OptionParser(option_list = option_list))
con <- open_stability_db(opt$db)

# Accepts EITHER rslurm's own raw `_rslurm_<jobname>` naming (pass the
# directory that CONTAINS it) OR a directory that already IS the results
# folder itself, OR a shared parent with a plain `<jobname>/` subfolder
# (results_*.RDS directly inside, no `_rslurm_` prefix) -- confirmed in
# practice this project's actual cluster output can land any of these
# ways depending on how it's collected/organized after a job completes.
# Same fallback pattern as R/ingest_core_results.R.
resolve_family_dir <- function(bundle_dir, jobname) {
  candidates <- c(file.path(bundle_dir, paste0("_rslurm_", jobname)), file.path(bundle_dir, jobname))
  for (d in candidates) if (length(list.files(d, pattern = "^results_\\d+\\.RDS$")) > 0) return(d)
  if (length(list.files(bundle_dir, pattern = "^results_\\d+\\.RDS$")) > 0) return(bundle_dir)
  NA_character_
}

for (jobname in c("projectr_within_grid", "projectr_cross_grid")) {
  fam_dir <- resolve_family_dir(opt$`bundle-dir`, jobname)
  if (is.na(fam_dir)) { message(jobname, ": no results dir found, skipping"); next }
  result_files <- list.files(fam_dir, pattern = "^results_\\d+\\.RDS$", full.names = TRUE)
  message(jobname, ": ", length(result_files), " result file(s)")

  DBI::dbExecute(con, "BEGIN")
  n_entries <- 0
  for (f in result_files) {
    x <- readRDS(f)
    # Same either-shape handling as R/ingest_enrichment_results.R's
    # fgsea_grid block -- see submit_job_family()'s `max_array_size` doc:
    # a results_<task>.RDS is either one run_projectr_job() return value
    # or a plain list of several (batched). `x$source_fit_id` is only
    # present on the former.
    entries <- if (!is.null(x$source_fit_id)) list(x) else x
    n_entries <- n_entries + length(entries)

    for (x in entries) {
      if (is.null(x$result) || is.null(x$result$projection)) {
        message("  ", basename(f), ": no projection (failed task) -- skipping")
        next
      }
      proj <- x$result$projection
      r2 <- x$result$r_squared
      # pval (generic/non-pca mode: a components x samples significance
      # matrix) and pvar (pca mode only: a per-component % variance-
      # explained-in-target vector) are mutually exclusive -- projectR's
      # dispatch returns exactly one of the two depending on whether
      # `loadings` was a plain matrix or a `prcomp`-classed object (see
      # R/ingest_jobs/projectr_job.R's pca branch). Both already sat in
      # the saved projection_file artifact, just never summarized into a
      # queryable column the way r_squared already was.
      pval_mat <- x$result$pval
      pvar_vec <- x$result$pvar

      art_dir <- artifacts_dir(opt$db, x$source_dataset_id)
      dir.create(art_dir, recursive = TRUE, showWarnings = FALSE)
      fname <- sprintf("%s_src%d_tgt%s_projection.rds", jobname, x$source_fit_id, x$target_dataset_id)
      saveRDS(x$result, file.path(art_dir, fname))

      DBI::dbExecute(con,
        "INSERT INTO projections (source_fit_id, source_dataset_id, target_dataset_id, method,
                                   projection_type, include_intercept, n_genes_matched, n_samples,
                                   mean_r_squared, median_r_squared,
                                   mean_pval, median_pval, mean_pvar, median_pvar, projection_file)
         VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        params = list(x$source_fit_id, x$source_dataset_id, x$target_dataset_id, x$method,
                      x$projection_type, isTRUE(x$include_intercept), x$n_genes_matched, ncol(proj),
                      mean(r2, na.rm = TRUE), median(r2, na.rm = TRUE),
                      if (!is.null(pval_mat)) mean(pval_mat, na.rm = TRUE) else NA_real_,
                      if (!is.null(pval_mat)) stats::median(pval_mat, na.rm = TRUE) else NA_real_,
                      if (!is.null(pvar_vec)) mean(pvar_vec, na.rm = TRUE) else NA_real_,
                      if (!is.null(pvar_vec)) stats::median(pvar_vec, na.rm = TRUE) else NA_real_,
                      file.path("stability_artifacts", x$source_dataset_id, fname)))
    }
  }
  DBI::dbExecute(con, "COMMIT")
  message(jobname, ": ", n_entries, " pair-level result(s) across those file(s)")
}
DBI::dbDisconnect(con)

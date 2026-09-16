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

for (jobname in c("projectr_within_grid", "projectr_cross_grid")) {
  fam_dir <- file.path(opt$`bundle-dir`, paste0("_rslurm_", jobname))
  if (!dir.exists(fam_dir)) { message(jobname, ": no results dir found, skipping"); next }
  result_files <- list.files(fam_dir, pattern = "^results_\\d+\\.RDS$", full.names = TRUE)
  message(jobname, ": ", length(result_files), " result file(s)")

  DBI::dbExecute(con, "BEGIN")
  for (f in result_files) {
    x <- readRDS(f)
    if (is.list(x) && length(x) == 1 && is.null(names(x))) x <- x[[1]]
    if (is.null(x$result) || is.null(x$result$projection)) {
      message("  ", basename(f), ": no projection (failed task) -- skipping")
      next
    }
    proj <- x$result$projection
    r2 <- x$result$r_squared

    art_dir <- artifacts_dir(opt$db, x$source_dataset_id)
    dir.create(art_dir, recursive = TRUE, showWarnings = FALSE)
    fname <- sprintf("%s_src%d_tgt%s_projection.rds", jobname, x$source_fit_id, x$target_dataset_id)
    saveRDS(x$result, file.path(art_dir, fname))

    DBI::dbExecute(con,
      "INSERT INTO projections (source_fit_id, source_dataset_id, target_dataset_id, method,
                                 projection_type, include_intercept, n_genes_matched, n_samples,
                                 mean_r_squared, median_r_squared, projection_file)
       VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
      params = list(x$source_fit_id, x$source_dataset_id, x$target_dataset_id, x$method,
                    x$projection_type, isTRUE(x$include_intercept), x$n_genes_matched, ncol(proj),
                    mean(r2, na.rm = TRUE), median(r2, na.rm = TRUE),
                    file.path("stability_artifacts", x$source_dataset_id, fname)))
  }
  DBI::dbExecute(con, "COMMIT")
}
DBI::dbDisconnect(con)

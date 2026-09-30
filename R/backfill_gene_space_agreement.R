# Offline backfill: gene-space agreement between the SAME method's basis
# learned on two different datasets -- Overview Tab 4's gene-space matrix
# (app/app.R, nav$mode == "overview"). Aggregates the exact same live
# computation the app's "Compare bases (projectR)" screen's Gene-space
# alignment tab already does for one pair at a time
# (project_basis_onto_basis() + Hungarian-matched reordering, both in
# app/R/comparison_helpers.R) into a single scalar per
# (source_dataset_id, target_dataset_id, method) -- cheap enough to show
# as a full dataset x dataset heatmap, unlike running it live for
# hundreds of pairs inside the app.
#
# Each row costs a live projectR() call + Ensembl remapping (unlike
# `projections`, which is a free aggregate read of data already computed
# offline by R/lib/ingest/projectr_pairs.R) -- deliberately scoped down by
# default (see --scope below), not run for all 992 x 5-method pairs at
# once in one sitting. Additive: re-running the same command only fills
# in rows not already present, unless --force.
#
# Usage:
#   Rscript R/backfill_gene_space_agreement.R <db_path> [--scope=within|all] [--methods=pca,spca,...] [--force]
#     --scope=within (default): only within-family pairs (same dataset
#       family, e.g. ANEMONES -> ANEMONES_DAY1, matched by the
#       "<family>_<suffix>" dataset_id naming convention this project's
#       configs already use) -- cheapest and smallest in count, per the
#       plan's own sequencing (build this first, `--scope=all` is a real
#       cluster job, not an interactive one).
#     --scope=all: every (source,target) pair that already has
#       sample-space projections for this method (i.e. every row already
#       in `projections`) -- much larger; intended to be run via sbatch,
#       not interactively.
#     --methods: comma-separated subset of pca,spca,nmf,cogaps,ica
#       (default: all five -- wgcna is never included, it has no loadings).
#     --force: recompute even pairs already in gene_space_agreement
#       (default: additive, skips existing rows).

library(here)
source(here("R/lib/ingest/db.R"))
suppressMessages(library(DBI))

args <- commandArgs(trailingOnly = TRUE)
positional <- args[!grepl("^--", args)]
if (length(positional) < 1) {
  stop("Usage: Rscript R/backfill_gene_space_agreement.R <db_path> [--scope=within|all] [--methods=...] [--force]")
}
db_path <- positional[[1]]
force <- "--force" %in% args
scope_arg <- args[grepl("^--scope=", args)]
scope <- if (length(scope_arg) > 0) sub("^--scope=", "", scope_arg[1]) else "within"
methods_arg <- args[grepl("^--methods=", args)]
methods <- if (length(methods_arg) > 0) strsplit(sub("^--methods=", "", methods_arg[1]), ",")[[1]] else
  c("pca", "spca", "nmf", "cogaps", "ica")

con <- open_stability_db(db_path)
options(stability.db_dir = normalizePath(dirname(db_path)))
# Reuses the app's own helpers directly -- optimal_fit_for_method()
# (db_helpers.R), project_basis_onto_basis()/hungarian_match_abs()
# (comparison_helpers.R) -- all plain functions with no Shiny-reactive
# dependencies, safe to call from a bare Rscript. app/R/db_helpers.R's
# single-arg resolve_artifact() (relying on the stability.db_dir option
# just set above) intentionally shadows R/lib/ingest/db.R's two-arg
# version sourced above it -- same idiom used throughout this session's
# verification Rscript calls.
source(here("app/R/db_helpers.R"))
source(here("app/R/metadata_helpers.R"))
source(here("app/R/comparison_helpers.R"))
source(here("R/lib/ingest/similarity.R"))
source(here("R/lib/ingest/symbol_mapping.R"))

is_within_family <- function(a, b) {
  startsWith(a, paste0(b, "_")) || startsWith(b, paste0(a, "_"))
}

pairs <- DBI::dbGetQuery(con, "SELECT DISTINCT source_dataset_id, target_dataset_id, method FROM projections")
pairs <- pairs[pairs$method %in% methods, ]
if (scope == "within") {
  keep <- mapply(is_within_family, pairs$source_dataset_id, pairs$target_dataset_id)
  pairs <- pairs[keep, ]
} else if (scope != "all") {
  stop("--scope must be 'within' or 'all'")
}
if (!force) {
  done <- DBI::dbGetQuery(con, "SELECT source_dataset_id, target_dataset_id, method FROM gene_space_agreement")
  done_key <- paste(done$source_dataset_id, done$target_dataset_id, done$method)
  pairs <- pairs[!paste(pairs$source_dataset_id, pairs$target_dataset_id, pairs$method) %in% done_key, ]
}
message(sprintf("gene_space_agreement: %d pairs to compute (scope=%s, methods=%s, force=%s)",
                 nrow(pairs), scope, paste(methods, collapse = ","), force))

for (i in seq_len(nrow(pairs))) {
  src_ds <- pairs$source_dataset_id[i]; tgt_ds <- pairs$target_dataset_id[i]; m <- pairs$method[i]
  fit_a <- optimal_fit_for_method(con, src_ds, m)
  fit_b <- optimal_fit_for_method(con, tgt_ds, m)
  if (is.null(fit_a) || is.null(fit_b)) {
    message(sprintf("  [%d/%d] %s -> %s (%s): SKIP (no optimal fit on one side)", i, nrow(pairs), src_ds, tgt_ds, m))
    next
  }
  proj <- tryCatch(project_basis_onto_basis(con, fit_a$fit_id, src_ds, fit_b$fit_id, tgt_ds), error = function(e) NULL)
  if (is.null(proj)) {
    message(sprintf("  [%d/%d] %s -> %s (%s): SKIP (no usable loadings/shared genes)", i, nrow(pairs), src_ds, tgt_ds, m))
    next
  }
  P <- proj$projection
  match_idx <- hungarian_match_abs(P)
  mean_abs_diag <- mean(abs(P[match_idx]), na.rm = TRUE)
  DBI::dbExecute(con,
    "INSERT INTO gene_space_agreement
       (source_dataset_id, target_dataset_id, method, source_fit_id, target_fit_id, n_genes_matched, mean_abs_diagonal, computed_at)
     VALUES (?, ?, ?, ?, ?, ?, ?, ?)
     ON CONFLICT(source_dataset_id, target_dataset_id, method) DO UPDATE SET
       source_fit_id = excluded.source_fit_id, target_fit_id = excluded.target_fit_id,
       n_genes_matched = excluded.n_genes_matched, mean_abs_diagonal = excluded.mean_abs_diagonal,
       computed_at = excluded.computed_at",
    params = list(src_ds, tgt_ds, m, fit_a$fit_id, fit_b$fit_id, proj$n_genes_matched, mean_abs_diag, as.character(Sys.time())))
  message(sprintf("  [%d/%d] %s -> %s (%s): mean_abs_diagonal=%.3f (%d genes)",
                   i, nrow(pairs), src_ds, tgt_ds, m, mean_abs_diag, proj$n_genes_matched))
}

DBI::dbDisconnect(con)
message("Done.")

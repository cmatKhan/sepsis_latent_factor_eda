# Parameterized queries for the stability explorer app. Every function
# takes the DBI connection; the app reads whatever the DB currently
# contains (see app.R for the STABILITY_DB path parameterization).

METRICS <- c("cosine", "pearson", "spearman")

check_metric <- function(metric) {
  match.arg(metric, METRICS)
}

list_datasets <- function(con) {
  DBI::dbGetQuery(con, "SELECT dataset_id, description FROM datasets ORDER BY dataset_id")
}

list_ingests <- function(con, dataset_id) {
  DBI::dbGetQuery(con,
    "SELECT jobname, family, method, n_results, ingested_at
     FROM ingests WHERE dataset_id = ? ORDER BY method, family",
    params = list(dataset_id))
}

#' Per-method overview: fit counts by status + headline stability
#' (mean matched same-rank similarity) per metric.
method_overview <- function(con, dataset_id) {
  counts <- DBI::dbGetQuery(con,
    "SELECT method,
            SUM(CASE WHEN status = 'ok' THEN 1 ELSE 0 END) AS n_ok,
            SUM(CASE WHEN status = 'failed' THEN 1 ELSE 0 END) AS n_failed,
            SUM(CASE WHEN status = 'missing' THEN 1 ELSE 0 END) AS n_missing
     FROM fits WHERE dataset_id = ? GROUP BY method",
    params = list(dataset_id))
  stab <- DBI::dbGetQuery(con,
    "SELECT ft.method,
            AVG(fp.cosine) AS cosine, AVG(fp.pearson) AS pearson, AVG(fp.spearman) AS spearman
     FROM factor_pairs fp JOIN fits ft ON ft.fit_id = fp.fit_a
     WHERE fp.matched = 1 AND fp.same_rank = 1 AND ft.dataset_id = ?
     GROUP BY ft.method",
    params = list(dataset_id))
  wgcna <- DBI::dbGetQuery(con,
    "SELECT 'wgcna' AS method, AVG(wp.ari) AS ari
     FROM wgcna_fit_pairs wp JOIN fits ft ON ft.fit_id = wp.fit_a
     WHERE ft.dataset_id = ?",
    params = list(dataset_id))
  wto <- DBI::dbGetQuery(con,
    "SELECT 'wto' AS method, AVG(wp.pearson) AS pearson, AVG(wp.spearman) AS spearman
     FROM wto_fit_pairs wp JOIN fits ft ON ft.fit_id = wp.fit_a
     WHERE ft.dataset_id = ?",
    params = list(dataset_id))
  list(counts = counts, stability = stab, wgcna = wgcna, wto = wto)
}

#' Matched same-rank factor similarities with rank/seed info -- the Level-1
#' "seed stability vs rank" data.
seed_stability_by_rank <- function(con, dataset_id, method) {
  DBI::dbGetQuery(con,
    "SELECT fa.rank, fa.seed AS seed_a, fb.seed AS seed_b,
            fp.factor_a, fp.factor_b, fp.cosine, fp.pearson, fp.spearman
     FROM factor_pairs fp
     JOIN fits fa ON fa.fit_id = fp.fit_a
     JOIN fits fb ON fb.fit_id = fp.fit_b
     WHERE fp.matched = 1 AND fp.same_rank = 1
       AND fa.dataset_id = ? AND fa.method = ?",
    params = list(dataset_id, method))
}

maskcv_curve <- function(con, dataset_id, method) {
  DBI::dbGetQuery(con,
    "SELECT rank, alpha, mse FROM maskcv_results
     WHERE dataset_id = ? AND method = ? ORDER BY rank, alpha",
    params = list(dataset_id, method))
}

#' Cross-rank persistence: mean matched similarity between every rank pair.
crossrank_matrix <- function(con, dataset_id, method) {
  DBI::dbGetQuery(con,
    "SELECT fa.rank AS rank_a, fb.rank AS rank_b,
            AVG(fp.cosine) AS cosine, AVG(fp.pearson) AS pearson, AVG(fp.spearman) AS spearman,
            COUNT(*) AS n
     FROM factor_pairs fp
     JOIN fits fa ON fa.fit_id = fp.fit_a
     JOIN fits fb ON fb.fit_id = fp.fit_b
     WHERE fp.matched = 1 AND fp.same_rank = 0
       AND fa.dataset_id = ? AND fa.method = ?
     GROUP BY fa.rank, fb.rank",
    params = list(dataset_id, method))
}

#' Seed x seed mean matched similarity at one rank (Level 2).
seedpair_matrix <- function(con, dataset_id, method, rank) {
  DBI::dbGetQuery(con,
    "SELECT fa.seed AS seed_a, fb.seed AS seed_b,
            AVG(fp.cosine) AS cosine, AVG(fp.pearson) AS pearson, AVG(fp.spearman) AS spearman
     FROM factor_pairs fp
     JOIN fits fa ON fa.fit_id = fp.fit_a
     JOIN fits fb ON fb.fit_id = fp.fit_b
     WHERE fp.matched = 1 AND fp.same_rank = 1
       AND fa.dataset_id = ? AND fa.method = ? AND fa.rank = ?
     GROUP BY fa.seed, fb.seed",
    params = list(dataset_id, method, rank))
}

fits_at_rank <- function(con, dataset_id, method, rank) {
  DBI::dbGetQuery(con,
    "SELECT fit_id, seed, mse, n_factors, loadings_file FROM fits
     WHERE dataset_id = ? AND method = ? AND rank = ? AND status = 'ok'
       AND family = 'seed_sweep'
     ORDER BY seed",
    params = list(dataset_id, method, rank))
}

#' Per-factor matched similarities for one fit vs all other same-rank fits
#' (handles both storage directions of the pair).
factor_stability_for_fit <- function(con, fit_id) {
  DBI::dbGetQuery(con,
    "SELECT fp.factor_a AS factor_index, fb.seed AS other_seed,
            fp.cosine, fp.pearson, fp.spearman
     FROM factor_pairs fp
     JOIN fits fb ON fb.fit_id = fp.fit_b
     WHERE fp.fit_a = ?1 AND fp.matched = 1 AND fp.same_rank = 1
     UNION ALL
     SELECT fp.factor_b AS factor_index, fa.seed AS other_seed,
            fp.cosine, fp.pearson, fp.spearman
     FROM factor_pairs fp
     JOIN fits fa ON fa.fit_id = fp.fit_a
     WHERE fp.fit_b = ?1 AND fp.matched = 1 AND fp.same_rank = 1",
    params = list(fit_id))
}

#' All factor x factor similarities for one specific fit pair.
factor_pair_heatmap_data <- function(con, fit_a, fit_b) {
  d <- DBI::dbGetQuery(con,
    "SELECT factor_a, factor_b, cosine, pearson, spearman, matched
     FROM factor_pairs WHERE fit_a = ?1 AND fit_b = ?2
     UNION ALL
     SELECT factor_b AS factor_a, factor_a AS factor_b, cosine, pearson, spearman, matched
     FROM factor_pairs WHERE fit_a = ?2 AND fit_b = ?1",
    params = list(fit_a, fit_b))
  d
}

#' This factor's Hungarian match in every other fit (all seeds + ranks).
factor_matches_everywhere <- function(con, fit_id, factor_index) {
  DBI::dbGetQuery(con,
    "SELECT fb.fit_id AS other_fit, fb.rank AS other_rank, fb.seed AS other_seed,
            fp.factor_b AS other_factor, fp.cosine, fp.pearson, fp.spearman
     FROM factor_pairs fp JOIN fits fb ON fb.fit_id = fp.fit_b
     WHERE fp.fit_a = ?1 AND fp.factor_a = ?2 AND fp.matched = 1
     UNION ALL
     SELECT fa.fit_id, fa.rank, fa.seed, fp.factor_a, fp.cosine, fp.pearson, fp.spearman
     FROM factor_pairs fp JOIN fits fa ON fa.fit_id = fp.fit_a
     WHERE fp.fit_b = ?1 AND fp.factor_b = ?2 AND fp.matched = 1",
    params = list(fit_id, factor_index))
}

get_fit <- function(con, fit_id) {
  DBI::dbGetQuery(con, "SELECT * FROM fits WHERE fit_id = ?", params = list(fit_id))
}

get_factor_id <- function(con, fit_id, factor_index) {
  DBI::dbGetQuery(con,
    "SELECT factor_id FROM factors WHERE fit_id = ? AND factor_index = ?",
    params = list(fit_id, factor_index))$factor_id
}

#' loadings_file paths are stored relative to the DB file's directory
#' (see R/lib/ingest/db.R) -- resolve via the option set at app startup.
resolve_artifact <- function(path) {
  if (is.na(path) || !nzchar(path)) return(NA_character_)
  if (startsWith(path, "/")) return(path)
  db_dir <- getOption("stability.db_dir")
  stopifnot(!is.null(db_dir))
  file.path(db_dir, path)
}

load_loadings <- function(con, fit_id) {
  f <- get_fit(con, fit_id)
  if (nrow(f) == 0 || is.na(f$loadings_file)) return(NULL)
  path <- resolve_artifact(f$loadings_file)
  if (!file.exists(path)) return(NULL)
  readRDS(path)
}

## ---- network methods ---------------------------------------------------------

wgcna_fits <- function(con, dataset_id) {
  DBI::dbGetQuery(con,
    "SELECT fit_id, power, n_factors FROM fits
     WHERE dataset_id = ? AND method = 'wgcna' AND status = 'ok' ORDER BY power",
    params = list(dataset_id))
}

wgcna_ari_matrix <- function(con, dataset_id) {
  DBI::dbGetQuery(con,
    "SELECT fa.power AS power_a, fb.power AS power_b, wp.ari
     FROM wgcna_fit_pairs wp
     JOIN fits fa ON fa.fit_id = wp.fit_a
     JOIN fits fb ON fb.fit_id = wp.fit_b
     WHERE fa.dataset_id = ?",
    params = list(dataset_id))
}

wgcna_module_sizes <- function(con, fit_id) {
  DBI::dbGetQuery(con,
    "SELECT module, COUNT(*) AS n_genes FROM wgcna_modules
     WHERE fit_id = ? GROUP BY module ORDER BY module",
    params = list(fit_id))
}

wgcna_module_jaccard <- function(con, fit_a, fit_b) {
  DBI::dbGetQuery(con,
    "SELECT module_a, module_b, jaccard, matched FROM wgcna_module_pairs
     WHERE fit_a = ?1 AND fit_b = ?2
     UNION ALL
     SELECT module_b AS module_a, module_a AS module_b, jaccard, matched
     FROM wgcna_module_pairs WHERE fit_a = ?2 AND fit_b = ?1",
    params = list(fit_a, fit_b))
}

wgcna_module_best_matches <- function(con, dataset_id, fit_id) {
  DBI::dbGetQuery(con,
    "SELECT p.module_a AS module, fb.power AS other_power, p.module_b AS other_module, p.jaccard
     FROM wgcna_module_pairs p JOIN fits fb ON fb.fit_id = p.fit_b
     WHERE p.fit_a = ?1 AND p.matched = 1
     UNION ALL
     SELECT p.module_b, fa.power, p.module_a, p.jaccard
     FROM wgcna_module_pairs p JOIN fits fa ON fa.fit_id = p.fit_a
     WHERE p.fit_b = ?1 AND p.matched = 1",
    params = list(fit_id))
}

wto_fits <- function(con, dataset_id) {
  DBI::dbGetQuery(con,
    "SELECT fit_id, n_boot, seed, delta, loadings_file FROM fits
     WHERE dataset_id = ? AND method = 'wto' AND status = 'ok' ORDER BY n_boot, seed",
    params = list(dataset_id))
}

wto_pair_matrix <- function(con, dataset_id) {
  DBI::dbGetQuery(con,
    "SELECT fa.fit_id AS fit_a, fb.fit_id AS fit_b,
            fa.n_boot AS n_a, fa.seed AS seed_a, fb.n_boot AS n_b, fb.seed AS seed_b,
            wp.pearson, wp.spearman, wp.jaccard_sig
     FROM wto_fit_pairs wp
     JOIN fits fa ON fa.fit_id = wp.fit_a
     JOIN fits fb ON fb.fit_id = wp.fit_b
     WHERE fa.dataset_id = ?",
    params = list(dataset_id))
}

## ---- enrichment cache ---------------------------------------------------------

enrichment_cached <- function(con, factor_id, query_type, direction) {
  hit <- DBI::dbGetQuery(con,
    "SELECT COUNT(*) AS n FROM enrichment_queried
     WHERE factor_id = ? AND query_type = ? AND direction = ?",
    params = list(factor_id, query_type, direction))$n > 0
  if (!hit) return(NULL)
  DBI::dbGetQuery(con,
    "SELECT source, term_id, term_name, p_value, intersection_size, term_size
     FROM enrichment_cache
     WHERE factor_id = ? AND query_type = ? AND direction = ?
     ORDER BY p_value",
    params = list(factor_id, query_type, direction))
}

enrichment_store <- function(con, factor_id, query_type, direction, gost_result) {
  DBI::dbExecute(con,
    "INSERT OR REPLACE INTO enrichment_queried (factor_id, query_type, direction, queried_at)
     VALUES (?, ?, ?, datetime('now'))",
    params = list(factor_id, query_type, direction))
  if (!is.null(gost_result) && !is.null(gost_result$result) && nrow(gost_result$result) > 0) {
    r <- gost_result$result
    DBI::dbWriteTable(con, "enrichment_cache", data.frame(
      factor_id = factor_id, query_type = query_type, direction = direction,
      source = r$source, term_id = r$term_id, term_name = r$term_name,
      p_value = r$p_value, intersection_size = r$intersection_size,
      term_size = r$term_size,
      genes = if (!is.null(r$intersection)) as.character(r$intersection) else NA_character_,
      queried_at = as.character(Sys.time())
    ), append = TRUE)
  }
  invisible(NULL)
}

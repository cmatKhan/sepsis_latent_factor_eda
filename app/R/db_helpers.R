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
  list(counts = counts, stability = stab, wgcna = wgcna)
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

#' Mean +/- SD of the seed-sweep family's own (in-sample) reconstruction
#' `mse` per rank -- distinct from scree_mse_by_rank() (PCA/sPCA have no
#' seed dimension at all, so there's nothing to average over). Answers
#' "how much does reconstruction quality vary across random seeds at this
#' rank?", a genuine stability question a plain scree plot doesn't
#' address.
seed_sweep_mse_by_rank <- function(con, dataset_id, method) {
  d <- DBI::dbGetQuery(con,
    "SELECT rank, mse FROM fits
     WHERE dataset_id = ? AND method = ? AND family = 'seed_sweep' AND status = 'ok' AND mse IS NOT NULL",
    params = list(dataset_id, method))
  if (nrow(d) == 0) return(d)
  agg <- aggregate(mse ~ rank, data = d, FUN = function(x) c(mean = mean(x), sd = sd(x), n = length(x)))
  out <- data.frame(rank = agg$rank, mean_mse = agg$mse[, "mean"], sd_mse = agg$mse[, "sd"], n = agg$mse[, "n"])
  out[order(out$rank), ]
}

#' Scree-plot data for PCA/sPCA (SCREE_RANK_METHODS): each is deterministic
#' given its rank (+ sPCA's para), one `fits` row per rank already carrying
#' an in-sample reconstruction `mse` -- no separate masking-CV family
#' needed (that one was never actually wired up -- see app.R's
#' SCREE_RANK_METHODS comment). `alpha` is sPCA's `para` sparsity penalty,
#' crossed with rank; always NA for PCA.
scree_mse_by_rank <- function(con, dataset_id, method) {
  DBI::dbGetQuery(con,
    "SELECT rank, alpha, mse FROM fits
     WHERE dataset_id = ? AND method = ? AND status = 'ok' AND mse IS NOT NULL
     ORDER BY rank, alpha",
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

#' Distinct ranks with at least one ok fit for a method -- used by PCA's
#' reduced Level 1 rank-select (falls back to this when a method has no
#' masking-CV family at all) and by the Level 2+ breadcrumb's rank
#' dropdown for every FACTORIZATION_METHODS method. No `family` filter
#' (see fits_at_rank()'s matching comment): sPCA's fits live under
#' 'param_grid', not 'seed_sweep', so filtering on the latter always
#' returned zero rows for sPCA and crashed the breadcrumb (setNames() on
#' a length-0 vector against paste0()'s length-1 "rank " -- paste0()
#' treats a zero-length argument as "" for recycling, not as "propagate
#' zero length", so the mismatch only surfaces here, not as an empty
#' dropdown).
distinct_ranks <- function(con, dataset_id, method) {
  DBI::dbGetQuery(con,
    "SELECT DISTINCT rank FROM fits
     WHERE dataset_id = ? AND method = ? AND status = 'ok'
     ORDER BY rank",
    params = list(dataset_id, method))$rank
}

#' No `family` filter needed here (every ok fit at this rank is a real
#' fit -- both seed_sweep methods (pca/nmf/cogaps/ica) and sPCA's
#' 'param_grid' family (see PARAM_GRID_METHODS in R/lib/ingest/extract.R)
#' since K/para have no genuine seed dimension) -- `alpha` is included
#' because sPCA's `para` value is stored there (see extract_result()'s
#' spca branch) and is what the app labels sPCA's per-rank fit selector
#' with (there being no real `seed` to show).
fits_at_rank <- function(con, dataset_id, method, rank) {
  DBI::dbGetQuery(con,
    "SELECT fit_id, seed, alpha, mse, n_factors, loadings_file FROM fits
     WHERE dataset_id = ? AND method = ? AND rank = ? AND status = 'ok'
     ORDER BY seed, alpha",
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

#' Every ok fit for a method, in whatever shape that method's Level 1
#' already uses to build its own fit choices (seed_sweep across all ranks
#' for pca/nmf/cogaps; direct_fits() for sPCA/CP/Tucker; wgcna_fits() for
#' WGCNA) -- used by the standalone "Compare methods" screen's fit
#' pickers, which need to offer EVERY fit up front rather than drilling
#' down one rank at a time.
all_fits_for_compare <- function(con, dataset_id, method) {
  if (method == "wgcna") return(wgcna_fits(con, dataset_id)$fit_id)
  if (method %in% c("spca", "cp", "tucker")) return(direct_fits(con, dataset_id, method)$fit_id)
  DBI::dbGetQuery(con,
    "SELECT fit_id FROM fits
     WHERE dataset_id = ? AND method = ? AND family = 'seed_sweep' AND status = 'ok'
     ORDER BY rank, seed",
    params = list(dataset_id, method))$fit_id
}

get_fit <- function(con, fit_id) {
  DBI::dbGetQuery(con, "SELECT * FROM fits WHERE fit_id = ?", params = list(fit_id))
}

#' Human-readable descriptor for one fit, method-aware -- used in
#' breadcrumbs, Level 1's fit-select labels, and Level 3 plot titles so
#' every method (including the "direct fit" ones with no seed dimension:
#' sPCA/CP/Tucker) gets a sensible label instead of a literal "rank NA
#' seed NA".
fit_descriptor <- function(con, method, fit_id) {
  f <- get_fit(con, fit_id)
  if (nrow(f) == 0) return("")
  if (method == "wgcna") {
    sprintf("power %s", f$power)
  } else if (method == "spca") {
    sprintf("K=%s", f$rank)
  } else if (method == "cp") {
    sprintf("num_components=%s", f$rank)
  } else if (method == "tucker") {
    sprintf("rank_genes=%s, rank_subjects=%s, rank_time=%s", f$rank_genes, f$rank_subjects, f$rank_time)
  } else if (!is.na(f$seed)) {
    sprintf("rank %s seed %s", f$rank, f$seed)
  } else {
    sprintf("rank %s", f$rank)
  }
}

#' All ok fits for a "direct fit" method (sPCA/CP/Tucker) -- no seed
#' dimension to group by, so (unlike fits_at_rank()) this returns every
#' fit for the method directly.
direct_fits <- function(con, dataset_id, method) {
  DBI::dbGetQuery(con,
    "SELECT fit_id, rank, rank_genes, rank_subjects, rank_time, mse, n_factors
     FROM fits WHERE dataset_id = ? AND method = ? AND status = 'ok'
     ORDER BY rank, rank_genes, rank_subjects, rank_time",
    params = list(dataset_id, method))
}

#' Lowest-MSE fit for a direct-fit method -- Level 0's headline, standing
#' in for the masking-CV-based headlines other methods use (sPCA/CP/Tucker
#' have no masking_cv family at all).
best_direct_fit <- function(con, dataset_id, method) {
  d <- direct_fits(con, dataset_id, method)
  d <- d[!is.na(d$mse), ]
  if (nrow(d) == 0) return(NULL)
  d[which.min(d$mse), ]
}

get_factor_id <- function(con, fit_id, factor_index) {
  DBI::dbGetQuery(con,
    "SELECT factor_id FROM factors WHERE fit_id = ? AND factor_index = ?",
    params = list(fit_id, factor_index))$factor_id
}

#' Resolve which factor_id's cached enrichment should represent
#' (fit_id, factor_index) -- either that fit's own factor_id (if
#' anything has been queried for it already), or an exact-loadings-match
#' sibling fit's factor_id at some other rank/seed within the SAME
#' method + dataset.
#'
#' This generalizes what used to be a PCA-only special case in app.R's
#' find_or_reuse_enrichment(). PCA's rank truncation is a genuine
#' mathematical guarantee (prcomp() with a smaller `rank.` never changes
#' earlier components -- a rank-10 fit's factor 3 IS, bit for bit,
#' rank-20's factor 3), so representative-fit selection
#' (R/lib/ingest/redundancy.R::representative_fit_ids() collapsing PCA to
#' a single max-rank fit before running fgsea_grid) never actually omits
#' any other rank's enrichment -- it's identical, just never redundantly
#' recomputed. Confirmed directly (2026-09-22): the pre-existing PCA-only
#' code already searched siblings across ALL ranks (no rank filter), so
#' this was already correct for PCA -- the gap was that (a) it was
#' hardcoded to PCA only, and (b) cached_enrichment_summary() (used by
#' the Compare tab) never called it at all, for any method.
#'
#' Other seed-sweep methods (nmf/cogaps/ica/spca) have NO such guarantee
#' across DIFFERENT seeds/para values -- each seed's factorization is
#' genuinely its own answer, not a truncation of anything else -- but
#' this still checks for an exact match there (fastICA/NNLM/CoGAPS
#' occasionally converge to bit-identical solutions from different
#' seeds), a real if less common win, at negligible cost since it's just
#' comparing already-loaded loading vectors.
#'
#' Returns NULL if the fit/factor doesn't exist; otherwise always returns
#' SOME factor_id (falling back to the fit's own, even if nothing is
#' cached for it, so callers can distinguish "genuinely nothing
#' computed anywhere" from "this factor doesn't exist").
resolve_enrichment_factor_id <- function(con, fit_id, factor_index) {
  factor_id <- get_factor_id(con, fit_id, factor_index)
  if (length(factor_id) != 1) return(NULL)

  was_queried <- function(fid) {
    DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM enrichment_queried WHERE factor_id = ?",
                     params = list(fid))$n > 0
  }
  if (was_queried(factor_id)) return(factor_id)

  fit <- get_fit(con, fit_id)
  if (nrow(fit) == 0 || !(fit$method %in% c("pca", "nmf", "cogaps", "ica", "spca"))) return(factor_id)

  L <- load_loadings(con, fit_id)
  if (is.null(L) || !(factor_index %in% seq_len(ncol(L)))) return(factor_id)
  v <- L[, factor_index]

  siblings <- DBI::dbGetQuery(con,
    "SELECT fit_id FROM fits WHERE dataset_id = ? AND method = ? AND fit_id != ? AND status = 'ok'",
    params = list(fit$dataset_id, fit$method, fit_id))$fit_id
  for (other_fit in siblings) {
    L2 <- load_loadings(con, other_fit)
    if (is.null(L2)) next
    for (fi2 in seq_len(ncol(L2))) {
      v2 <- L2[, fi2]
      if (length(v) != length(v2) || !identical(names(v), names(v2))) next
      if (!isTRUE(all.equal(as.numeric(v), as.numeric(v2), tolerance = 1e-8))) next
      other_factor_id <- get_factor_id(con, other_fit, fi2)
      if (length(other_factor_id) == 1 && was_queried(other_factor_id)) return(other_factor_id)
    }
  }
  factor_id
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

#' CP/Tucker's third (time-mode) loading matrix -- rows = timepoint
#' levels, columns = component. NULL for any other method (no
#' time_loadings_file).
load_time_loadings <- function(con, fit_id) {
  f <- get_fit(con, fit_id)
  if (nrow(f) == 0 || is.na(f$time_loadings_file)) return(NULL)
  path <- resolve_artifact(f$time_loadings_file)
  if (!file.exists(path)) return(NULL)
  readRDS(path)
}

#' sPCA's diagnostics bundle -- list(pev, var_all, n_nonzero), see
#' R/lib/ingest/db.R's spca_diag_file column comment. NULL for any other
#' method (no spca_diag_file) or a fit predating this column.
load_spca_diag <- function(con, fit_id) {
  f <- get_fit(con, fit_id)
  if (nrow(f) == 0 || is.na(f$spca_diag_file)) return(NULL)
  path <- resolve_artifact(f$spca_diag_file)
  if (!file.exists(path)) return(NULL)
  readRDS(path)
}

#' This fit's gene-clustering tree (merge/height/order/labels -- see
#' R/lib/ingest/ingest_dataset.R::compute_wgcna_dendro()), reconstructed
#' as a real `hclust`-classed object ready for WGCNA::plotDendroAndColors().
#' NULL for any fit compute_wgcna_dendro() hasn't been run for yet
#' (real, expected -- it's an explicit opt-in step, not part of normal
#' ingest; see --recompute-dendro).
load_wgcna_dendro <- function(con, fit_id) {
  f <- get_fit(con, fit_id)
  if (nrow(f) == 0 || is.na(f$wgcna_dendro_file)) return(NULL)
  path <- resolve_artifact(f$wgcna_dendro_file)
  if (!file.exists(path)) return(NULL)
  tree <- readRDS(path)
  hc <- list(merge = tree$merge, height = tree$height, order = tree$order,
             labels = tree$labels, method = "average", dist.method = "1 - TOM")
  class(hc) <- "hclust"
  hc
}

#' Classic PCA/sPCA biplot data: ALL sample scores + a SUBSET of gene
#' loadings (the `top_n` genes by combined magnitude on the two chosen
#' components -- showing every gene as an arrow is unreadable at
#' genomics scale, thousands of features vs. the handful typical in a
#' textbook biplot). Loadings are pre-scaled so arrow tips land within
#' ~80% of the score cloud's radius -- the standard ad hoc biplot
#' convention (e.g. factoextra::fviz_pca_biplot, ggbiplot) since raw
#' loadings and raw scores live on very different numeric scales and a
#' biplot's real content is DIRECTION/relative magnitude, not a shared
#' numeric axis.
#'
#' `comp_x`/`comp_y` are 1-based component indices (matching
#' colnames(loadings)/colnames(scores)). Returns NULL if the fit has no
#' loadings/scores, or either component index is out of range.
biplot_data <- function(con, fit_id, comp_x, comp_y, top_n = 15) {
  L <- load_loadings(con, fit_id)
  S <- load_scores(con, fit_id)
  if (is.null(L) || is.null(S)) return(NULL)
  if (!(comp_x %in% seq_len(ncol(L))) || !(comp_y %in% seq_len(ncol(L)))) return(NULL)

  scores_df <- data.frame(sample_id = rownames(S), x = S[, comp_x], y = S[, comp_y])

  lx <- L[, comp_x]; ly <- L[, comp_y]
  mag <- sqrt(lx^2 + ly^2)
  keep <- order(-mag)[seq_len(min(top_n, sum(mag > 0)))]

  score_radius <- suppressWarnings(max(sqrt(scores_df$x^2 + scores_df$y^2), na.rm = TRUE))
  loading_radius <- suppressWarnings(max(mag[keep], na.rm = TRUE))
  scale_factor <- if (is.finite(loading_radius) && loading_radius > 0) 0.8 * score_radius / loading_radius else 1

  arrows_df <- data.frame(
    gene = rownames(L)[keep],
    x = lx[keep] * scale_factor,
    y = ly[keep] * scale_factor
  )
  list(scores = scores_df, arrows = arrows_df)
}

## ---- pattern drivers (differential features, projectR::projectionDriveR()) ---
##
## Read-only helpers for R/lib/ingest/driver.R's ingest-time batch pass
## (pattern_drivers table) plus a self-contained on-demand runner for
## combinations that pass wasn't scoped to cover (arbitrary grouping
## column / factor / mode) -- deliberately NOT sourcing R/lib/ingest/
## driver.R itself. Unlike enrichment (all computed by the cluster ingest
## pipeline now, R/ingest_jobs/fgsea_job.R -- this app has no live-compute
## enrichment path at all), pattern-driver combinations genuinely can't
## all be pre-enumerated at ingest time (arbitrary grouping column x level
## pair x factor), so this one on-demand runner stays.

#' Already-cached projectionDriveR() results for one fit -- populates a
#' selector of (factor, grouping column, level pair, mode) combinations
#' already computed at ingest time.
pattern_drivers_for_fit <- function(con, fit_id) {
  DBI::dbGetQuery(con,
    "SELECT driver_id, factor_index, grouping_col, group1_level, group2_level, mode,
            n_genes_considered, n_significant_shared, computed_at
     FROM pattern_drivers WHERE fit_id = ? ORDER BY factor_index, grouping_col, group1_level, group2_level",
    params = list(fit_id))
}

load_pattern_driver_result <- function(con, driver_id) {
  f <- DBI::dbGetQuery(con, "SELECT result_file FROM pattern_drivers WHERE driver_id = ?", params = list(driver_id))
  if (nrow(f) == 0 || is.na(f$result_file)) return(NULL)
  path <- resolve_artifact(f$result_file)
  if (!file.exists(path)) return(NULL)
  readRDS(path)
}

#' 2-6 level categorical columns of a dataset's registered sample
#' metadata -- same eligibility rule as driver.R's batch pass (excludes
#' the id column and anything higher-cardinality), but a slightly wider
#' 2-6 range since this is a user-driven on-demand selector, not an
#' exhaustive ingest sweep.
available_grouping_columns <- function(con, dataset_id) {
  sm <- dataset_metadata(con, dataset_id, "sample")
  if (is.null(sm)) return(character(0))
  Filter(function(col) {
    col != "sample_id" && (is.character(sm[[col]]) || is.factor(sm[[col]])) &&
      length(unique(stats::na.omit(sm[[col]]))) %in% 2:6
  }, names(sm))
}

#' On-demand projectionDriveR() run for a combination not already cached
#' -- checks pattern_drivers first (cache hit -> just loads it), otherwise
#' computes live and stores the result the SAME way the ingest-time batch
#' pass does, so it shows up in pattern_drivers_for_fit() on next visit.
run_pattern_driver_on_demand <- function(con, db_path, fit_id, factor_index, dataset_id,
                                          grouping_col, group1_level, group2_level, mode = "CI") {
  cached <- DBI::dbGetQuery(con,
    "SELECT driver_id FROM pattern_drivers
     WHERE fit_id=? AND factor_index=? AND grouping_col=? AND group1_level=? AND group2_level=? AND mode=?",
    params = list(fit_id, factor_index, grouping_col, group1_level, group2_level, mode))
  if (nrow(cached) > 0) return(load_pattern_driver_result(con, cached$driver_id[1]))

  loadings <- load_loadings(con, fit_id)
  if (is.null(loadings) || !(factor_index %in% seq_len(ncol(loadings)))) return(NULL)
  pattern_name <- colnames(loadings)[factor_index] %||% paste0("factor_", factor_index)
  colnames(loadings)[factor_index] <- pattern_name

  mat_file <- DBI::dbGetQuery(con, "SELECT matrix_file FROM datasets WHERE dataset_id = ?", params = list(dataset_id))
  if (nrow(mat_file) == 0 || is.na(mat_file$matrix_file)) return(NULL)
  mat <- as.matrix(readRDS(resolve_artifact(mat_file$matrix_file)))

  sm <- dataset_metadata(con, dataset_id, "sample")
  if (is.null(sm) || !(grouping_col %in% names(sm))) return(NULL)
  ids1 <- intersect(sm$sample_id[sm[[grouping_col]] == group1_level], colnames(mat))
  ids2 <- intersect(sm$sample_id[sm[[grouping_col]] == group2_level], colnames(mat))
  if (length(ids1) < 3 || length(ids2) < 3) return(NULL)

  result <- tryCatch(
    projectR::projectionDriveR(cellgroup1 = mat[, ids1, drop = FALSE], cellgroup2 = mat[, ids2, drop = FALSE],
                                loadings = loadings, pattern_name = pattern_name, display = FALSE, mode = mode),
    error = function(e) NULL)
  if (is.null(result)) return(NULL)

  n_shared <- if (mode == "CI") length(result$sig_genes$significant_shared_genes %||% character(0))
              else length(result$sig_genes$PV_significant_shared_genes %||% character(0))
  n_considered <- if (mode == "CI") nrow(result$mean_ci) else nrow(result$mean_stats)
  art_dir <- file.path(getOption("stability.db_dir"), "stability_artifacts", dataset_id)
  dir.create(art_dir, recursive = TRUE, showWarnings = FALSE)
  fname <- sprintf("driver_fit%d_f%d_%s_%s-vs-%s_%s_ondemand.rds", fit_id, factor_index, grouping_col,
                    make.names(group1_level), make.names(group2_level), mode)
  saveRDS(result[setdiff(names(result), "plotted_ci")], file.path(art_dir, fname))
  DBI::dbExecute(con,
    "INSERT OR IGNORE INTO pattern_drivers (fit_id, factor_index, grouping_col, group1_level, group2_level, mode,
                                             n_genes_considered, n_significant_shared, result_file, computed_at)
     VALUES (?,?,?,?,?,?,?,?,?,datetime('now'))",
    params = list(fit_id, factor_index, grouping_col, group1_level, group2_level, mode,
                  n_considered, n_shared, file.path("stability_artifacts", dataset_id, fname)))
  result
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

#' One module's member genes ranked by kME (intramodular connectivity /
#' module membership, WGCNA::signedKME()) -- the standard hub-gene
#' ranking blockwiseModules() computes internally but never returns (see
#' compute_wgcna_kme()'s header in R/lib/ingest/ingest_dataset.R). Genes
#' with no wgcna_kme row (older fits ingested before kME support, or a
#' sample-overlap edge case at compute time) still appear with kme = NA,
#' sorted last -- never silently dropped from the membership list.
wgcna_module_kme <- function(con, fit_id, module) {
  DBI::dbGetQuery(con,
    "SELECT m.gene, k.kme FROM wgcna_modules m
     LEFT JOIN wgcna_kme k ON k.fit_id = m.fit_id AND k.gene = m.gene AND k.module = m.module
     WHERE m.fit_id = ? AND m.module = ?
     ORDER BY k.kme DESC",
    params = list(fit_id, module))
}

#' Every gene's module assignment for this fit, module 0 (unassigned)
#' included -- for coloring a full gene dendrogram by module, where every
#' leaf needs a color regardless of whether it landed in a real module.
wgcna_all_modules <- function(con, fit_id) {
  DBI::dbGetQuery(con, "SELECT gene, module FROM wgcna_modules WHERE fit_id = ?",
                   params = list(fit_id))
}

#' Every gene's kME to its OWN assigned module, for every module in this
#' fit -- the "how well-defined is each module" diagnostic: a module
#' where members' kME clusters tightly near 1 is a real, coherent
#' co-expression unit; one with a broad/low kME spread is diffuse. Module
#' 0 (WGCNA's "unassigned genes" bucket) has no real eigengene to belong
#' to, so it's excluded here, matching the app's existing convention
#' (e.g. gene_sets_from_wgcna()) of treating module 0 as not a real module.
wgcna_kme_all <- function(con, fit_id) {
  DBI::dbGetQuery(con,
    "SELECT k.module, k.gene, k.kme FROM wgcna_kme k
     WHERE k.fit_id = ? AND k.module != 0",
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

#' WGCNA::pickSoftThreshold()'s scale-free-topology fit per power -- see
#' R/lib/ingest/ingest_dataset.R::compute_wgcna_sft(). Dataset-level (not
#' per-fit), populated as a side effect of the normal ingest matrix-caching
#' step; NULL/empty for datasets ingested before this was added, until
#' re-ingested.
wgcna_sft_fit <- function(con, dataset_id) {
  DBI::dbGetQuery(con,
    "SELECT power, sft_r_sq, slope, mean_k FROM wgcna_sft WHERE dataset_id = ? ORDER BY power",
    params = list(dataset_id))
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

#' Numeric-metadata field names available for this dataset's Gene
#' Significance (see compute_wgcna_gene_significance()) -- restricted to
#' test='spearman' rows, since only those carry a SIGNED statistic
#' usable in a GS-vs-kME scatter (Kruskal-Wallis/categorical fields have
#' no sign). Empty for datasets whose sample metadata is all-categorical
#' (real, expected -- e.g. ANEMONES), or that predate this table
#' (compute_wgcna_gene_significance() hasn't been run for them yet).
wgcna_gene_significance_fields <- function(con, dataset_id) {
  DBI::dbGetQuery(con,
    "SELECT DISTINCT field FROM wgcna_gene_significance
     WHERE dataset_id = ? AND test = 'spearman' ORDER BY field",
    params = list(dataset_id))$field
}

#' The classic WGCNA hub-gene validation join: for one module, each
#' member gene's kME (module membership) against its Gene Significance
#' for one chosen trait. A real kME-vs-GS relationship within a module
#' means that module is phenotype-relevant, not just a network-structure
#' artifact -- see R/lib/ingest/ingest_dataset.R::compute_wgcna_gene_significance()'s
#' header. `dataset_id` is needed because wgcna_gene_significance is
#' dataset-level (not per-fit, unlike wgcna_kme).
#'
#' Joins through wgcna_modules (real membership), NOT just `wgcna_kme
#' WHERE module = ?` -- signedKME() correlates every gene against EVERY
#' module's eigengene, not just its own, so wgcna_kme alone has a
#' (module) row for every gene in the dataset regardless of actual
#' assignment; filtering on it directly silently returns all genes
#' instead of this module's real members (confirmed directly: module 1
#' of one real ROSE fit has 2,236 assigned genes, but `wgcna_kme WHERE
#' module = 1` alone returns all 7,080).
wgcna_module_gs_kme <- function(con, fit_id, module, dataset_id, field) {
  DBI::dbGetQuery(con,
    "SELECT m.gene, k.kme, g.statistic AS gs
     FROM wgcna_modules m
     JOIN wgcna_kme k ON k.fit_id = m.fit_id AND k.gene = m.gene AND k.module = m.module
     JOIN wgcna_gene_significance g ON g.dataset_id = ?3 AND g.gene = m.gene AND g.field = ?4
     WHERE m.fit_id = ?1 AND m.module = ?2 AND g.test = 'spearman'",
    params = list(fit_id, module, dataset_id, field))
}

## ---- enrichment cache ---------------------------------------------------------
## Read-only from this app's side -- all enrichment computation happens in
## the cluster ingest pipeline (R/ingest_jobs/fgsea_job.R,
## R/ingest_jobs/wgcna_ora_job.R + R/ingest_enrichment_results.R), which
## writes enrichment_cache/enrichment_queried directly.

enrichment_cached <- function(con, factor_id, query_type, direction) {
  hit <- DBI::dbGetQuery(con,
    "SELECT COUNT(*) AS n FROM enrichment_queried
     WHERE factor_id = ? AND query_type = ? AND direction = ?",
    params = list(factor_id, query_type, direction))$n > 0
  if (!hit) return(NULL)
  DBI::dbGetQuery(con,
    "SELECT source, term_id, term_name, p_value, intersection_size, term_size, query_size, is_main_pathway
     FROM enrichment_cache
     WHERE factor_id = ? AND query_type = ? AND direction = ?
     ORDER BY p_value",
    params = list(factor_id, query_type, direction))
}

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
#' (mean matched same-rank similarity) per metric. ABS() on cosine/
#' pearson/spearman: a magnitude claim ("how similar"), not a sign claim
#' -- needed for pca/spca/ica (unique only up to sign; see
#' R/lib/ingest/pairs.R's `sign_ambiguous`), a no-op for nmf/cogaps
#' (already non-negative).
method_overview <- function(con, dataset_id) {
  counts <- DBI::dbGetQuery(con,
    "SELECT method,
            SUM(CASE WHEN status = 'ok' THEN 1 ELSE 0 END) AS n_ok,
            SUM(CASE WHEN status = 'failed' THEN 1 ELSE 0 END) AS n_failed,
            SUM(CASE WHEN status = 'missing' THEN 1 ELSE 0 END) AS n_missing
     FROM fits WHERE dataset_id = ? GROUP BY method",
    params = list(dataset_id))
  stab <- DBI::dbGetQuery(con,
    "SELECT fa.method,
            AVG(ABS(fp.cosine)) AS cosine, AVG(ABS(fp.pearson)) AS pearson, AVG(ABS(fp.spearman)) AS spearman
     FROM factor_pairs fp
     JOIN fits fa ON fa.fit_id = fp.fit_a
     JOIN fits fb ON fb.fit_id = fp.fit_b
     WHERE fp.matched = 1 AND fp.same_rank = 1 AND fa.dataset_id = ?
       AND (fa.bootstrap IS NULL OR fa.bootstrap = 0)
       AND (fb.bootstrap IS NULL OR fb.bootstrap = 0)
     GROUP BY fa.method",
    params = list(dataset_id))
  wgcna <- DBI::dbGetQuery(con,
    "SELECT 'wgcna' AS method, AVG(wp.ari) AS ari
     FROM wgcna_fit_pairs wp JOIN fits ft ON ft.fit_id = wp.fit_a
     WHERE ft.dataset_id = ?",
    params = list(dataset_id))
  list(counts = counts, stability = stab, wgcna = wgcna)
}

#' Matched same-rank factor similarities with rank/seed info -- the Level-1
#' "seed stability vs rank" data. ABS() on cosine/pearson/spearman: see
#' method_overview()'s comment -- a no-op for nmf/cogaps, a real fix for
#' ica (its only STABILITY_METHODS sibling that's sign-ambiguous). Excludes
#' ICA's bootstrap fits on both sides of the pair -- this panel's meaning
#' ("how much do plain reseeded reruns agree") stays exactly what it
#' always was; the new ICASSO clustering (R/lib/ingest/icasso.R,
#' app.R's "ICASSO cluster quality" tab) is where bootstrap+reinit runs
#' are deliberately pooled together instead.
seed_stability_by_rank <- function(con, dataset_id, method) {
  DBI::dbGetQuery(con,
    "SELECT fa.rank, fa.seed AS seed_a, fb.seed AS seed_b,
            fp.factor_a, fp.factor_b, ABS(fp.cosine) AS cosine, ABS(fp.pearson) AS pearson, ABS(fp.spearman) AS spearman
     FROM factor_pairs fp
     JOIN fits fa ON fa.fit_id = fp.fit_a
     JOIN fits fb ON fb.fit_id = fp.fit_b
     WHERE fp.matched = 1 AND fp.same_rank = 1
       AND fa.dataset_id = ? AND fa.method = ?
       AND (fa.bootstrap IS NULL OR fa.bootstrap = 0)
       AND (fb.bootstrap IS NULL OR fb.bootstrap = 0)",
    params = list(dataset_id, method))
}

#' Mean +/- SD of the seed-sweep family's own (in-sample) reconstruction
#' `mse` per rank -- distinct from scree_mse_by_rank() (PCA/sPCA have no
#' seed dimension at all, so there's nothing to average over). Answers
#' "how much does reconstruction quality vary across random seeds at this
#' rank?", a genuine stability question a plain scree plot doesn't
#' address. Excludes ICA's bootstrap fits (`bootstrap = 1`): their `mse`
#' is reconstruction error against a RESAMPLED dataset, not a comparable
#' claim about the real one (see R/lib/ingest/db.R's `fits.bootstrap`
#' comment) -- a no-op filter for every other method (bootstrap is always
#' NULL there).
seed_sweep_mse_by_rank <- function(con, dataset_id, method) {
  d <- DBI::dbGetQuery(con,
    "SELECT rank, mse FROM fits
     WHERE dataset_id = ? AND method = ? AND family = 'seed_sweep' AND status = 'ok' AND mse IS NOT NULL
       AND (bootstrap IS NULL OR bootstrap = 0)",
    params = list(dataset_id, method))
  if (nrow(d) == 0) return(d)
  agg <- aggregate(mse ~ rank, data = d, FUN = function(x) c(mean = mean(x), sd = sd(x), n = length(x)))
  out <- data.frame(rank = agg$rank, mean_mse = agg$mse[, "mean"], sd_mse = agg$mse[, "sd"], n = agg$mse[, "n"])
  out[order(out$rank), ]
}

#' Scree-plot data for PCA/sPCA (SCREE_RANK_METHODS): in-sample
#' reconstruction `mse` per rank. PCA has one fit per dataset at the max
#' rank, so its curve is that fit's mse_by_rank (pca_rank_curve()); sPCA
#' has one `fits` row per (rank, para), with para in `alpha`.
scree_mse_by_rank <- function(con, dataset_id, method) {
  if (method == "pca") {
    # One PCA fit per dataset at the max rank (components are nested): the
    # per-rank curve is its stored mse_by_rank (see pca_rank_curve()).
    curve <- pca_rank_curve(con, dataset_id)
    return(if (is.null(curve)) data.frame(rank = integer(0), alpha = numeric(0), mse = numeric(0))
           else data.frame(rank = curve$rank, alpha = NA_real_, mse = curve$mse))
  }
  DBI::dbGetQuery(con,
    "SELECT rank, alpha, mse FROM fits
     WHERE dataset_id = ? AND method = ? AND status = 'ok' AND mse IS NOT NULL
     ORDER BY rank, alpha",
    params = list(dataset_id, method))
}

#' Cross-rank persistence: mean matched similarity between every rank pair.
#' ABS(): see method_overview()'s comment -- used only by pca/spca here
#' (both sign-ambiguous), so this is a real fix, not a no-op.
crossrank_matrix <- function(con, dataset_id, method) {
  DBI::dbGetQuery(con,
    "SELECT fa.rank AS rank_a, fb.rank AS rank_b,
            AVG(ABS(fp.cosine)) AS cosine, AVG(ABS(fp.pearson)) AS pearson, AVG(ABS(fp.spearman)) AS spearman,
            COUNT(*) AS n
     FROM factor_pairs fp
     JOIN fits fa ON fa.fit_id = fp.fit_a
     JOIN fits fb ON fb.fit_id = fp.fit_b
     WHERE fp.matched = 1 AND fp.same_rank = 0
       AND fa.dataset_id = ? AND fa.method = ?
     GROUP BY fa.rank, fb.rank",
    params = list(dataset_id, method))
}

#' Seed x seed mean matched similarity at one rank (Level 2). ABS(): see
#' method_overview()'s comment. Excludes ICA's bootstrap fits -- they
#' share the same `seed` values as the plain reinit-only fits (same 10-
#' seed list, see R/methods/ica.R), so without this exclusion a bootstrap
#' fit would silently merge into the same seed x seed cell as its
#' non-bootstrap namesake via GROUP BY fa.seed, fb.seed.
seedpair_matrix <- function(con, dataset_id, method, rank) {
  DBI::dbGetQuery(con,
    "SELECT fa.seed AS seed_a, fb.seed AS seed_b,
            AVG(ABS(fp.cosine)) AS cosine, AVG(ABS(fp.pearson)) AS pearson, AVG(ABS(fp.spearman)) AS spearman
     FROM factor_pairs fp
     JOIN fits fa ON fa.fit_id = fp.fit_a
     JOIN fits fb ON fb.fit_id = fp.fit_b
     WHERE fp.matched = 1 AND fp.same_rank = 1
       AND fa.dataset_id = ? AND fa.method = ? AND fa.rank = ?
       AND (fa.bootstrap IS NULL OR fa.bootstrap = 0)
       AND (fb.bootstrap IS NULL OR fb.bootstrap = 0)
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
  if (method == "pca") {
    f <- pca_fit_row(con, dataset_id)
    return(if (nrow(f) == 0) integer(0) else seq_len(f$rank))
  }
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
#' with (there being no real `seed` to show). Excludes ICA's bootstrap
#' fits: they'd otherwise show up as a confusing duplicate "seed" entry
#' (same seed value as their non-bootstrap namesake, NULL scores_file --
#' no real sample identity) in this per-rank fit picker; they're pooled
#' into ICASSO's clustering instead (R/lib/ingest/icasso.R), not
#' individually selectable here.
fits_at_rank <- function(con, dataset_id, method, rank) {
  if (method == "pca") {
    # The single PCA fit, viewed at `rank` components -- callers load it
    # with load_loadings()/load_scores(n = rank).
    f <- pca_fit_row(con, dataset_id)
    f <- f[!is.na(f$rank) & f$rank >= rank, c("fit_id", "seed", "alpha", "mse", "n_factors", "loadings_file"),
           drop = FALSE]
    if (nrow(f) == 1) {
      curve <- pca_rank_curve(con, dataset_id)
      f$n_factors <- rank
      if (!is.null(curve)) f$mse <- curve$mse[match(rank, curve$rank)]
    }
    return(f)
  }
  DBI::dbGetQuery(con,
    "SELECT fit_id, seed, alpha, mse, n_factors, loadings_file FROM fits
     WHERE dataset_id = ? AND method = ? AND rank = ? AND status = 'ok'
       AND (bootstrap IS NULL OR bootstrap = 0)
     ORDER BY seed, alpha",
    params = list(dataset_id, method, rank))
}

#' Per-factor matched similarities for one fit vs all other same-rank fits
#' (handles both storage directions of the pair). ABS(): see
#' method_overview()'s comment. `fit_id` itself is always a non-bootstrap
#' fit in practice (fits_at_rank()'s selector already excludes bootstrap
#' fits), but the OTHER side of a stored pair could still be one --
#' excluded on both sides for the same "keep this panel's meaning to
#' plain reseeded reruns" reasoning as seed_stability_by_rank().
factor_stability_for_fit <- function(con, fit_id) {
  DBI::dbGetQuery(con,
    "SELECT fp.factor_a AS factor_index, fb.seed AS other_seed,
            ABS(fp.cosine) AS cosine, ABS(fp.pearson) AS pearson, ABS(fp.spearman) AS spearman
     FROM factor_pairs fp
     JOIN fits fb ON fb.fit_id = fp.fit_b
     WHERE fp.fit_a = ?1 AND fp.matched = 1 AND fp.same_rank = 1
       AND (fb.bootstrap IS NULL OR fb.bootstrap = 0)
     UNION ALL
     SELECT fp.factor_b AS factor_index, fa.seed AS other_seed,
            ABS(fp.cosine) AS cosine, ABS(fp.pearson) AS pearson, ABS(fp.spearman) AS spearman
     FROM factor_pairs fp
     JOIN fits fa ON fa.fit_id = fp.fit_a
     WHERE fp.fit_b = ?1 AND fp.matched = 1 AND fp.same_rank = 1
       AND (fa.bootstrap IS NULL OR fa.bootstrap = 0)",
    params = list(fit_id))
}

#' All factor x factor similarities for one specific fit pair.
#' Deliberately SIGNED, unlike the magnitude-style aggregates above (this
#' feeds the Level-2 raw heatmap, scale c(-1, 1), where seeing an actual
#' anti-correlation is the point) -- the `matched` column already reflects
#' the sign-robust Hungarian assignment from R/lib/ingest/pairs.R, so the
#' correct pair is outlined even though the displayed value stays signed.
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
#' Deliberately SIGNED, same reasoning as factor_pair_heatmap_data() --
#' feeds Level 3's loading-scatter view, where the actual sign relationship
#' is informative, not a "how similar" magnitude claim.
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
#' for pca/nmf/cogaps; direct_fits() for sPCA; wgcna_fits() for WGCNA) --
#' used by the standalone "Compare methods" screen's fit pickers, which
#' need to offer EVERY fit up front rather than drilling down one rank at
#' a time. Excludes ICA's bootstrap fits -- they have no `scores_file`
#' (no real per-sample identity to compare against another fit) and
#' aren't meant to be individually inspected, only pooled into ICASSO's
#' clustering (see R/lib/ingest/icasso.R); a no-op filter for every other
#' method.
all_fits_for_compare <- function(con, dataset_id, method) {
  if (method == "wgcna") return(wgcna_fits(con, dataset_id)$fit_id)
  if (method == "spca") return(direct_fits(con, dataset_id, method)$fit_id)
  DBI::dbGetQuery(con,
    "SELECT fit_id FROM fits
     WHERE dataset_id = ? AND method = ? AND family = 'seed_sweep' AND status = 'ok'
       AND (bootstrap IS NULL OR bootstrap = 0)
     ORDER BY rank, seed",
    params = list(dataset_id, method))$fit_id
}

get_fit <- function(con, fit_id) {
  DBI::dbGetQuery(con, "SELECT * FROM fits WHERE fit_id = ?", params = list(fit_id))
}

#' Human-readable descriptor for one fit, method-aware -- used in
#' breadcrumbs, Level 1's fit-select labels, and Level 3 plot titles so
#' every method (including the "direct fit" one with no seed dimension,
#' sPCA) gets a sensible label instead of a literal "rank NA seed NA".
fit_descriptor <- function(con, method, fit_id) {
  f <- get_fit(con, fit_id)
  if (nrow(f) == 0) return("")
  if (method == "wgcna") {
    sprintf("power %s", f$power)
  } else if (method == "spca") {
    sprintf("K=%s", f$rank)
  } else if (!is.na(f$seed)) {
    sprintf("rank %s seed %s", f$rank, f$seed)
  } else {
    sprintf("rank %s", f$rank)
  }
}

#' All ok fits for a "direct fit" method (sPCA) -- no seed dimension to
#' group by, so (unlike fits_at_rank()) this returns every fit for the
#' method directly. `rank_genes`/`rank_subjects`/`rank_time` are still
#' selected for backward compatibility with historical CP/Tucker rows
#' (removed 2026-09-29, see FACTORIZATION_METHODS's header in app.R) that
#' remain in the DB but are no longer surfaced anywhere in the app.
direct_fits <- function(con, dataset_id, method) {
  DBI::dbGetQuery(con,
    "SELECT fit_id, rank, rank_genes, rank_subjects, rank_time, mse, n_factors
     FROM fits WHERE dataset_id = ? AND method = ? AND status = 'ok'
     ORDER BY rank, rank_genes, rank_subjects, rank_time",
    params = list(dataset_id, method))
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

#' A fit's loadings (genes x factors). `n` keeps only the first n factors --
#' how the app views PCA's single max-rank fit at a smaller rank (its
#' components are nested); NULL keeps all.
load_loadings <- function(con, fit_id, n = NULL) {
  f <- get_fit(con, fit_id)
  if (nrow(f) == 0 || is.na(f$loadings_file)) return(NULL)
  path <- resolve_artifact(f$loadings_file)
  if (!file.exists(path)) return(NULL)
  first_n(readRDS(path), n)
}

first_n <- function(m, n) {
  if (is.null(n) || is.null(m) || is.na(n)) return(m)
  m[, seq_len(min(n, ncol(m))), drop = FALSE]
}

#' The dataset's single PCA fit (fits row), or a 0-row frame.
pca_fit_row <- function(con, dataset_id) {
  DBI::dbGetQuery(con,
    "SELECT * FROM fits WHERE dataset_id = ? AND method = 'pca' AND status = 'ok'
     ORDER BY rank DESC LIMIT 1",
    params = list(dataset_id))
}

#' PCA's diagnostics bundle -- list(sdev, center, mse_by_rank), see
#' R/lib/ingest/extract.R's pca branch.
load_pca_diag <- function(con, fit_id) {
  f <- get_fit(con, fit_id)
  if (nrow(f) == 0 || is.na(f$pca_diag_file)) return(NULL)
  path <- resolve_artifact(f$pca_diag_file)
  if (!file.exists(path)) return(NULL)
  readRDS(path)
}

#' Per-rank view of the single PCA fit: data.frame(rank, mse, pev) for
#' rank = 1..max -- reconstruction MSE of the first `rank` components
#' (mse_by_rank) and the cumulative proportion of variance they explain
#' (from the full sdev spectrum). NULL if the dataset has no PCA fit.
pca_rank_curve <- function(con, dataset_id) {
  f <- pca_fit_row(con, dataset_id)
  if (nrow(f) == 0) return(NULL)
  d <- load_pca_diag(con, f$fit_id)
  if (is.null(d) || is.null(d$mse_by_rank)) return(NULL)
  k <- seq_along(d$mse_by_rank)
  data.frame(rank = k, mse = d$mse_by_rank, pev = (cumsum(d$sdev^2) / sum(d$sdev^2))[k],
             fit_id = f$fit_id)
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

#' This fit's fastICA() diagnostics -- list(W, K, prewhiten_sdev), see
#' R/lib/ingest/db.R's ica_diag_file column comment. NULL for every ICA
#' fit currently in the DB (confirmed directly, 2026-09-23: 0/1,600 real
#' fits have this populated -- the capturing code is correct but every
#' currently-ingested ica_grid result predates it, and unlike sPCA/CoGAPS
#' there's no raw_result_file to backfill W/K from after the fact; see
#' R/methods/ica.R's header). Will return real data once ica_grid is
#' re-run with current code and re-ingested.
load_ica_diag <- function(con, fit_id) {
  f <- get_fit(con, fit_id)
  if (nrow(f) == 0 || is.na(f$ica_diag_file)) return(NULL)
  path <- resolve_artifact(f$ica_diag_file)
  if (!file.exists(path)) return(NULL)
  readRDS(path)
}

#' fastICA()'s own unmixing matrix W has no self-reported convergence
#' diagnostic -- a genuinely converged W should be close to orthonormal
#' (W %*% t(W) ~= I), the "cheap real check" R/methods/ica.R's own
#' header comment calls out. Returns the Frobenius-norm residual (0 =
#' perfectly orthonormal); NULL if this fit has no ica_diag_file yet.
ica_orthonormality_residual <- function(con, fit_id) {
  diag <- load_ica_diag(con, fit_id)
  if (is.null(diag) || is.null(diag$W)) return(NULL)
  W <- as.matrix(diag$W)
  norm(W %*% t(W) - diag(nrow(W)), type = "F")
}

#' How much of the pre-ICA PCA-prewhitening reduction's total variance
#' (pcfit$sdev, R/methods/ica.R:36-38) was actually retained by keeping
#' only n_pcs = n.comp + 5 components -- the literature (Lee & Batzoglou
#' 2003) found that reducing dimensionality via PCA before ICA measurably
#' hurts result quality vs. full-rank ICA; this lets a user see exactly
#' how much variance that tradeoff costs for a REAL fit, without this
#' function deciding anything about whether that cost is acceptable.
#' NULL if this fit has no ica_diag_file yet.
ica_prewhiten_variance_retained <- function(con, fit_id) {
  diag <- load_ica_diag(con, fit_id)
  if (is.null(diag) || is.null(diag$prewhiten_sdev)) return(NULL)
  sdev <- diag$prewhiten_sdev
  # K (fastICA's own pre-whitening matrix) is n_pcs x n.comp -- nrow(K) is
  # the PCA-prewhitening dimension (n.comp + 5 in R/methods/ica.R) fed
  # INTO fastICA, i.e. exactly how many of pcfit$sdev's components were
  # kept. Confirmed directly against a real fastICA() call: dim(K) =
  # c(n_pcs_in, n.comp), dim(W) = c(n.comp, n.comp) -- W is square and
  # does NOT carry n_pcs, unlike an earlier draft of this function assumed.
  n_pcs <- if (!is.null(diag$K)) nrow(diag$K) else length(sdev)
  sum(sdev[seq_len(n_pcs)]^2) / sum(sdev^2)
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

#' Per-component non-Gaussianity for one ICA fit (Lee & Batzoglou 2003 --
#' see R/lib/ingest/db.R's ica_component_kurtosis comment for the full
#' rationale). Plain query wrapper, not an artifact read -- unlike
#' load_ica_diag()'s W/K/prewhiten_sdev, this is small enough to live as
#' its own queryable table.
ica_component_kurtosis <- function(con, fit_id) {
  DBI::dbGetQuery(con,
    "SELECT component, kurtosis, excess_kurtosis FROM ica_component_kurtosis
     WHERE fit_id = ? ORDER BY component",
    params = list(fit_id))
}

#' ICASSO (Himberg, Hyvärinen & Esposito 2004) cluster-level summary for
#' one (dataset, rank) -- see R/lib/ingest/icasso.R. One row per cluster;
#' `iq` >= ~0.7 is the literature's usual "trustworthy" threshold.
icasso_clusters <- function(con, dataset_id, rank) {
  DBI::dbGetQuery(con,
    "SELECT cluster_id, iq, n_members, centrotype_fit_id, centrotype_factor_index
     FROM icasso_clusters WHERE dataset_id = ? AND method = 'ica' AND rank = ?
     ORDER BY cluster_id",
    params = list(dataset_id, rank))
}

#' ICASSO per-run membership for one (dataset, rank) -- every (fit,
#' factor)'s cluster assignment, needed to color a rendered dendrogram or
#' recompute a cluster's composition. See R/lib/ingest/icasso.R.
icasso_membership <- function(con, dataset_id, rank) {
  DBI::dbGetQuery(con,
    "SELECT cluster_id, fit_id, factor_index, intra_sim
     FROM icasso_membership WHERE dataset_id = ? AND method = 'ica' AND rank = ?",
    params = list(dataset_id, rank))
}

#' ICASSO cluster quality (Iq) across every rank this dataset's ICA grid
#' swept -- long frame (rank, cluster_id, iq, n_members), one row per
#' cluster, for the Level-1 "Iq vs n.comp" overview plot.
icasso_iq_by_rank <- function(con, dataset_id) {
  DBI::dbGetQuery(con,
    "SELECT rank, cluster_id, iq, n_members FROM icasso_clusters
     WHERE dataset_id = ? AND method = 'ica' ORDER BY rank, cluster_id",
    params = list(dataset_id))
}

#' The cross-run ICASSO clustering tree for one (dataset, rank),
#' reconstructed as a real `hclust`-classed object -- same shape/idiom as
#' load_wgcna_dendro(). `hc$labels` are "fit_id:factor_index" strings (see
#' R/lib/ingest/icasso.R) -- join back to icasso_membership() to color
#' leaves by cluster. NULL if this rank's ICASSO clustering hasn't been
#' computed yet (e.g. fewer than 2 ok fits still exist at this rank).
load_icasso_dendro <- function(con, dataset_id, rank) {
  f <- DBI::dbGetQuery(con,
    "SELECT dendro_file FROM icasso_dendrograms WHERE dataset_id = ? AND method = 'ica' AND rank = ?",
    params = list(dataset_id, rank))
  if (nrow(f) == 0 || is.na(f$dendro_file[1])) return(NULL)
  path <- resolve_artifact(f$dendro_file[1])
  if (!file.exists(path)) return(NULL)
  tree <- readRDS(path)
  hc <- list(merge = tree$merge, height = tree$height, order = tree$order,
             labels = tree$labels, method = "average", dist.method = "1 - |cosine|")
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

#' Module sizes at every power for a dataset, with modules ranked
#' largest-to-smallest WITHIN each power (module 0 / "unassigned" kept as
#' its own unranked category). Ranking rather than raw module number is
#' the robust key for a size-profile-across-power plot: WGCNA numbers
#' modules by descending size within a fit already, but the profile
#' shouldn't silently depend on that holding exactly.
wgcna_module_size_profile <- function(con, dataset_id) {
  fits <- wgcna_fits(con, dataset_id)
  if (nrow(fits) == 0) return(data.frame(power = integer(), module = integer(),
                                          n_genes = integer(), rank = integer()))
  rows <- lapply(seq_len(nrow(fits)), function(i) {
    sizes <- wgcna_module_sizes(con, fits$fit_id[i])
    sizes$power <- fits$power[i]
    sizes
  })
  d <- do.call(rbind, rows)
  d$rank <- NA_integer_
  for (p in unique(d$power)) {
    idx <- which(d$power == p & d$module != 0)
    d$rank[idx] <- rank(-d$n_genes[idx], ties.method = "first")
  }
  d[, c("power", "module", "n_genes", "rank")]
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

#' This fit pair's Hungarian-BEST-matched module correspondences only
#' (matched = 1), in both directions -- the building block for tracing a
#' module's identity forward one power step at a time.
wgcna_module_best_match_pairs <- function(con, fit_a, fit_b) {
  DBI::dbGetQuery(con,
    "SELECT module_a, module_b, jaccard FROM wgcna_module_pairs
     WHERE fit_a = ?1 AND fit_b = ?2 AND matched = 1
     UNION ALL
     SELECT module_b AS module_a, module_a AS module_b, jaccard
     FROM wgcna_module_pairs WHERE fit_a = ?2 AND fit_b = ?1 AND matched = 1",
    params = list(fit_a, fit_b))
}

#' Traces each module's identity across the dataset's full tested power
#' sequence by chaining the already-computed Hungarian best-matches
#' (wgcna_module_pairs) between ADJACENT powers only -- adjacent in the
#' tested grid (e.g. 10 -> 12 -> 14, respecting real gaps), not adjacent
#' by arithmetic power value. A lineage starts at the lowest power's
#' modules and is extended forward as long as a matched=1 row carries it
#' to the next power; when no such row exists (its genes were absorbed
#' elsewhere without it being anyone's best match), the lineage
#' terminates. A module at a later power that is nobody's match target
#' starts a new lineage there, so genuinely new/re-split modules show up
#' rather than being silently dropped. Returns long data.frame:
#' lineage_id, power, module, n_genes, jaccard_to_prev (NA at a lineage's
#' first power).
wgcna_module_lineages <- function(con, dataset_id) {
  fits <- wgcna_fits(con, dataset_id)
  empty <- data.frame(lineage_id = character(), power = integer(), module = integer(),
                       n_genes = integer(), jaccard_to_prev = double())
  if (nrow(fits) == 0) return(empty)

  sizes_by_fit <- lapply(fits$fit_id, function(fid) wgcna_module_sizes(con, fid))
  names(sizes_by_fit) <- as.character(fits$fit_id)
  size_of <- function(fit_id, module) {
    s <- sizes_by_fit[[as.character(fit_id)]]
    v <- s$n_genes[s$module == module]
    if (length(v) == 0) NA_integer_ else v[1]
  }

  first_sizes <- sizes_by_fit[[as.character(fits$fit_id[1])]]
  active <- first_sizes$module[first_sizes$module != 0]
  next_id <- 1L
  lineage_of <- setNames(paste0("L", seq_along(active)), active)
  next_id <- length(active) + 1L

  out <- list()
  for (m in active) {
    out[[length(out) + 1]] <- data.frame(
      lineage_id = lineage_of[[as.character(m)]], power = fits$power[1], module = m,
      n_genes = size_of(fits$fit_id[1], m), jaccard_to_prev = NA_real_)
  }

  for (i in seq_len(nrow(fits) - 1)) {
    matches <- wgcna_module_best_match_pairs(con, fits$fit_id[i], fits$fit_id[i + 1])
    new_lineage_of <- list()
    claimed_targets <- character()
    for (m_char in names(lineage_of)) {
      m <- as.integer(m_char)
      row <- matches[matches$module_a == m, , drop = FALSE]
      if (nrow(row) == 0) next
      target <- row$module_b[1]
      lid <- lineage_of[[m_char]]
      new_lineage_of[[as.character(target)]] <- lid
      claimed_targets <- c(claimed_targets, as.character(target))
      out[[length(out) + 1]] <- data.frame(
        lineage_id = lid, power = fits$power[i + 1], module = target,
        n_genes = size_of(fits$fit_id[i + 1], target), jaccard_to_prev = row$jaccard[1])
    }
    next_sizes <- sizes_by_fit[[as.character(fits$fit_id[i + 1])]]
    unclaimed <- setdiff(as.character(next_sizes$module[next_sizes$module != 0]), claimed_targets)
    for (m_char in unclaimed) {
      lid <- paste0("L", next_id); next_id <- next_id + 1L
      new_lineage_of[[m_char]] <- lid
      out[[length(out) + 1]] <- data.frame(
        lineage_id = lid, power = fits$power[i + 1], module = as.integer(m_char),
        n_genes = size_of(fits$fit_id[i + 1], as.integer(m_char)), jaccard_to_prev = NA_real_)
    }
    lineage_of <- new_lineage_of
  }
  if (length(out) == 0) return(empty)
  do.call(rbind, out)
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

## ---- projectR projections ------------------------------------------------------
## Read-only from this app's side -- all projection computation happens in
## the cluster ingest pipeline (R/ingest_jobs/projectr_job.R +
## R/ingest_projectr_results.R), which writes the `projections` table +
## each row's projection_file artifact directly.

#' Every OTHER dataset in this fit's own timecourse family that its basis
#' has been projectR-projected onto -- i.e. real data for a trajectory
#' view (see app.R's "Trajectory across time" tab). `projection_type =
#' 'within_dataset'` is exactly this project's definition of "same
#' family" (config/dataset_families.yml, e.g. ANEMONES's own 4
#' timepoint-split sub-datasets) -- see R/lib/ingest/projectr_pairs.R's
#' header. Empty for a fit whose dataset has no family siblings (most
#' datasets), or whose projectr_within_grid hasn't been run/ingested yet.
within_family_projections <- function(con, fit_id) {
  DBI::dbGetQuery(con,
    "SELECT target_dataset_id, projection_file FROM projections
     WHERE source_fit_id = ? AND projection_type = 'within_dataset'",
    params = list(fit_id))
}

#' Every dataset that has EVER been the source of a real projectR
#' projection -- the "Compare bases" screen's source-dataset choices.
projection_source_datasets <- function(con) {
  DBI::dbGetQuery(con, "SELECT DISTINCT source_dataset_id FROM projections ORDER BY source_dataset_id")$source_dataset_id
}

#' Every dataset a given source dataset's bases have actually been
#' projected onto (within OR cross family -- both are real, already-
#' computed projectR runs; see projectr_pairs.R for how the two grids
#' differ). Only pairs with this real, already-computed data can be
#' compared -- nothing here ever triggers a new full-dataset projection.
projection_targets_for_dataset <- function(con, source_dataset_id) {
  DBI::dbGetQuery(con,
    "SELECT DISTINCT target_dataset_id FROM projections WHERE source_dataset_id = ? ORDER BY target_dataset_id",
    params = list(source_dataset_id))$target_dataset_id
}

#' Every real projectR projection from one dataset onto another --
#' one row per (method, fit) already run, joined to `fits` for real
#' rank/seed/alpha labels via fit_descriptor(). This is the data behind
#' the "Compare bases" screen's R²-comparison and sample-space-alignment
#' tabs -- both entirely free reads, no new computation.
projections_for_pair <- function(con, source_dataset_id, target_dataset_id) {
  DBI::dbGetQuery(con,
    "SELECT p.source_fit_id, p.method, p.mean_r_squared, p.median_r_squared,
            p.n_genes_matched, p.n_samples, p.projection_file, f.rank, f.seed, f.alpha
     FROM projections p JOIN fits f ON f.fit_id = p.source_fit_id
     WHERE p.source_dataset_id = ? AND p.target_dataset_id = ?
     ORDER BY p.method, f.rank",
    params = list(source_dataset_id, target_dataset_id))
}

## ---- Overview: cross-dataset "which method/rank is best" -------------------
## Backs the new "Overview" screen (app.R, nav$mode == "overview"), answering
## four questions in order: (1) optimal hyperparameter per (dataset, method),
## (2) best method per dataset, (3) cross-method agreement within a dataset,
## (4) cross-dataset agreement for one method.

#' All currently-supported methods' names, used to drive Overview loops
#' generically -- kept here as a plain hardcoded vector, not read from
#' app.R's method-family constants (STABILITY_METHODS etc.), since
#' db_helpers.R is sourced before those are defined and this file's
#' existing convention is to hardcode method names directly (see e.g.
#' method_overview()'s wgcna branch) rather than depend on app.R's load
#' order. CP/Tucker removed 2026-09-29 (see FACTORIZATION_METHODS's header
#' in app.R) -- their historical fits remain in the DB but are excluded
#' here deliberately.
OVERVIEW_METHODS <- c("pca", "spca", "nmf", "cogaps", "ica", "wgcna")

#' crossrank_matrix()-shaped data.frame (rank_a, rank_b, cosine, ...) ->
#' data.frame(rank, cosine), symmetrized. `factor_pairs` stores each pair
#' once with fit_a always the lower-fit_id side (see
#' compute_factor_pairs_from_universe()'s `i < j` loop, R/lib/ingest/
#' pairs.R) -- for a seed-sweep/param-grid family where fit_id roughly
#' tracks rank, that means the HIGHEST rank tested only ever appears as
#' rank_b, never rank_a. Aggregating on rank_a alone (as crossrank_matrix()
#' itself returns it) silently gives that highest rank NA stability --
#' confirmed live as a real bug once optimal_fit_for_method() started
#' picking the single highest-quality value directly (almost always the
#' highest rank, since quality tends to rise monotonically with it).
#' Combining both directions before aggregating fixes this for every
#' crossrank_matrix() consumer at once.
symmetric_rank_stability <- function(cr) {
  if (nrow(cr) == 0) return(data.frame(rank = numeric(0), cosine = numeric(0)))
  long <- data.frame(rank = c(cr$rank_a, cr$rank_b), cosine = c(cr$cosine, cr$cosine))
  aggregate(cosine ~ rank, data = long, FUN = mean)
}

#' Per-hyperparameter (rank/power) in-sample R² + stability curve for one
#' (dataset, method) -- the data behind Overview Tab 1 and the input to
#' optimal_fit_for_method() below. One row per distinct hyperparameter
#' value, with `fit_id` a REPRESENTATIVE fit at that value (for methods
#' with a seed dimension, whichever seed has the highest in-sample R² --
#' same "collapse to one representative" idiom as best_fit_per_rank() in
#' app.R's Compare-bases screen). `stability` means three genuinely
#' different things depending on method family (seed-stability for
#' nmf/cogaps/ica; cross-rank persistence for pca/spca; cross-power ARI
#' for wgcna) -- always "how much does the answer change when this
#' hyperparameter's neighbors are perturbed", never blended across
#' families, and `rank_label` always says which hyperparameter it is so
#' this is never silently ambiguous. Returns
#' data.frame(rank_label, rank_key, fit_id, in_sample_r2, stability), or an
#' empty frame (not an error) if this method has no ok fits for this
#' dataset -- a real, common case (NMF/CoGAPS exist for only 10/32
#' datasets).
method_rank_curve <- function(con, dataset_id, method) {
  empty <- data.frame(rank_label = character(0), rank_key = numeric(0), fit_id = integer(0),
                       in_sample_r2 = numeric(0), stability = numeric(0))
  val_sd <- DBI::dbGetQuery(con, "SELECT val_sd FROM matrix_diagnostics WHERE dataset_id = ?",
                             params = list(dataset_id))$val_sd
  val_sd <- if (length(val_sd) == 1 && !is.na(val_sd)) val_sd else NA_real_

  if (method == "pca") {
    # One nested fit: rank n = its first n components. Cross-rank
    # persistence is 1 by construction (the rank-n components ARE the
    # first n of every larger rank), so stability is reported as 1.
    curve <- pca_rank_curve(con, dataset_id)
    if (is.null(curve)) return(empty)
    data.frame(rank_label = paste("rank", curve$rank), rank_key = curve$rank, fit_id = curve$fit_id,
               in_sample_r2 = 1 - curve$mse / val_sd^2, stability = 1)

  } else if (method == "spca") {
    fits <- DBI::dbGetQuery(con,
      "SELECT fit_id, rank, alpha, mse FROM fits
       WHERE dataset_id = ? AND method = ? AND status = 'ok' AND mse IS NOT NULL",
      params = list(dataset_id, method))
    if (nrow(fits) == 0) return(empty)
    # Multiple para (alpha) per rank -- keep the lowest-mse (highest
    # quality) alpha as that rank's representative.
    agg <- do.call(rbind, lapply(split(fits, fits$rank), function(g) g[which.min(g$mse), ]))
    stab <- symmetric_rank_stability(crossrank_matrix(con, dataset_id, method))
    agg$stability <- stab$cosine[match(agg$rank, stab$rank)]
    agg$in_sample_r2 <- 1 - agg$mse
    data.frame(rank_label = paste("rank", agg$rank), rank_key = agg$rank, fit_id = agg$fit_id,
               in_sample_r2 = agg$in_sample_r2, stability = agg$stability)

  } else if (method %in% c("nmf", "cogaps", "ica")) {
    q <- seed_sweep_mse_by_rank(con, dataset_id, method)   # mean_mse per rank, already averaged over seeds
    if (nrow(q) == 0) return(empty)
    fits <- DBI::dbGetQuery(con,
      "SELECT fit_id, rank, mse FROM fits
       WHERE dataset_id = ? AND method = ? AND family = 'seed_sweep' AND status = 'ok' AND mse IS NOT NULL
         AND (bootstrap IS NULL OR bootstrap = 0)",
      params = list(dataset_id, method))
    rep_fit <- do.call(rbind, lapply(split(fits, fits$rank), function(g) g[which.min(g$mse), ]))
    ss <- seed_stability_by_rank(con, dataset_id, method)
    stab <- if (nrow(ss) > 0) aggregate(cosine ~ rank, data = ss, FUN = mean) else
      data.frame(rank = numeric(0), cosine = numeric(0))
    q$in_sample_r2 <- 1 - q$mean_mse / val_sd^2
    q$fit_id <- rep_fit$fit_id[match(q$rank, rep_fit$rank)]
    q$stability <- stab$cosine[match(q$rank, stab$rank)]
    data.frame(rank_label = paste("rank", q$rank), rank_key = q$rank, fit_id = q$fit_id,
               in_sample_r2 = q$in_sample_r2, stability = q$stability)

  } else if (method == "wgcna") {
    f <- wgcna_fits(con, dataset_id)
    if (nrow(f) == 0) return(empty)
    kme_r2 <- DBI::dbGetQuery(con, sprintf(
      "SELECT fit_id, AVG(kme*kme) AS r2 FROM wgcna_kme WHERE module != 0 AND fit_id IN (%s) GROUP BY fit_id",
      paste(f$fit_id, collapse = ",")))
    ari <- wgcna_ari_matrix(con, dataset_id)
    stab <- if (nrow(ari) > 0) aggregate(ari ~ power_a, data = ari, FUN = mean) else
      data.frame(power_a = numeric(0), ari = numeric(0))
    f$in_sample_r2 <- kme_r2$r2[match(f$fit_id, kme_r2$fit_id)]
    f$stability <- stab$ari[match(f$power, stab$power_a)]
    data.frame(rank_label = paste("power", f$power), rank_key = f$power, fit_id = f$fit_id,
               in_sample_r2 = f$in_sample_r2, stability = f$stability)

  } else {
    empty
  }
}

#' The "optimal" hyperparameter setting for one (dataset, method), used as
#' that method's representative fit on Overview Tabs 2-3 -- simply the
#' hyperparameter value with the highest in-sample R². Quality tends to
#' rise with more components for every method here, so in practice this
#' picks (at or near) the largest hyperparameter value tested; Tab 1 shows
#' the full curve (quality AND stability) so that's always visible, not
#' hidden behind a single pick. Returns NULL if the method has no ok fits
#' for this dataset.
optimal_fit_for_method <- function(con, dataset_id, method) {
  d <- method_rank_curve(con, dataset_id, method)
  d <- d[!is.na(d$in_sample_r2), ]
  if (nrow(d) == 0) return(NULL)
  d[which.max(d$in_sample_r2), ]
}

#' One row per method PRESENT for this dataset, combining (a) the
#' optimal-fit in-sample quality + stability from optimal_fit_for_method(),
#' and (b) out-of-sample R² -- median `mean_r_squared` across every target
#' this dataset's fits have been projected onto, split by within/cross
#' family (`projections.projection_type`), since the investigation found
#' rankings differ between the two. WGCNA has no out-of-sample columns
#' (never projectable) -- NA, not zero. This is Overview Tab 2's data --
#' one method-comparison scatter per DATASET, not aggregated across the
#' whole collection (a method can win on one dataset and lose on another).
method_dataset_summary <- function(con, dataset_id) {
  rows <- lapply(OVERVIEW_METHODS, function(m) {
    opt <- optimal_fit_for_method(con, dataset_id, m)
    if (is.null(opt)) return(NULL)
    data.frame(method = m, fit_id = opt$fit_id, rank_label = opt$rank_label,
               in_sample_r2 = opt$in_sample_r2, stability = opt$stability)
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) return(data.frame())
  out <- do.call(rbind, rows)

  oos <- DBI::dbGetQuery(con,
    "SELECT method, projection_type, mean_r_squared FROM projections WHERE source_dataset_id = ?",
    params = list(dataset_id))
  if (nrow(oos) > 0) {
    agg <- aggregate(mean_r_squared ~ method + projection_type, data = oos, FUN = median)
    within <- agg[agg$projection_type == "within_dataset", ]
    cross  <- agg[agg$projection_type == "cross_dataset", ]
    out$oos_r2_within <- within$mean_r_squared[match(out$method, within$method)]
    out$oos_r2_cross  <- cross$mean_r_squared[match(out$method, cross$method)]
  } else {
    out$oos_r2_within <- NA_real_
    out$oos_r2_cross <- NA_real_
  }
  out
}

#' Method x method agreement matrix for ONE dataset -- Overview Tab 3.
#' Runs the same pairwise comparisons the "Compare methods" screen already
#' offers one-pair-at-a-time (comparator_gene_sets()/jaccard_matrix() for
#' gene-set overlap, cluster_and_ari() for sample-space agreement) across
#' every pair of methods PRESENT for this dataset, each represented by its
#' optimal_fit_for_method() pick. Returns
#' list(methods = character vector actually compared,
#'      gene_jaccard = methods x methods matrix of mean best-Hungarian-
#'        matched |Jaccard| (diagonal = 1),
#'      sample_ari = methods x methods matrix of Adjusted Rand Index
#'        between k=4 sample/eigengene clusters (diagonal = 1)).
#' A method is silently dropped from both matrices if it has no optimal
#' fit for this dataset (e.g. NMF/CoGAPS on 22/32 datasets) -- never
#' forced into an NA-filled row that would look like a real comparison.
method_pairwise_agreement <- function(con, dataset_id, methods = OVERVIEW_METHODS, top_n = 50) {
  opts <- lapply(methods, function(m) optimal_fit_for_method(con, dataset_id, m))
  names(opts) <- methods
  present <- methods[!vapply(opts, is.null, logical(1))]
  if (length(present) < 2) return(list(methods = present, gene_jaccard = NULL, sample_ari = NULL))

  gj <- matrix(NA_real_, length(present), length(present), dimnames = list(present, present))
  sa <- matrix(NA_real_, length(present), length(present), dimnames = list(present, present))
  diag(gj) <- 1; diag(sa) <- 1

  for (i in seq_along(present)) {
    for (j in seq_along(present)) {
      if (i >= j) next
      m_a <- present[i]; m_b <- present[j]
      fit_a <- opts[[m_a]]$fit_id; fit_b <- opts[[m_b]]$fit_id

      sets_a <- comparator_gene_sets(con, m_a, fit_a, top_n)
      sets_b <- comparator_gene_sets(con, m_b, fit_b, top_n)
      if (!is.null(sets_a) && !is.null(sets_b)) {
        jm <- jaccard_matrix(sets_a$sets, sets_b$sets)
        match_idx <- hungarian_match_abs(jm)
        val <- mean(abs(jm[match_idx]), na.rm = TRUE)
        gj[i, j] <- val; gj[j, i] <- val
      }

      ca <- tryCatch(cluster_and_ari(con, fit_a, fit_b), error = function(e) NULL)
      if (!is.null(ca)) { sa[i, j] <- ca$ari; sa[j, i] <- ca$ari }
    }
  }
  list(methods = present, gene_jaccard = gj, sample_ari = sa)
}

#' Dataset x dataset out-of-sample R², aggregated (median across every
#' fit/rank/seed of this method) per (source, target) pair -- Overview Tab
#' 4's data. Long format; coverage is exactly whatever `projections`
#' already has for this method (comprehensive for pca/spca/ica; sparse and
#' source-gated for nmf/cogaps -- see the plan's coverage-map finding).
#' Empty for wgcna (never projected).
dataset_similarity_matrix <- function(con, method) {
  d <- DBI::dbGetQuery(con,
    "SELECT source_dataset_id, target_dataset_id, mean_r_squared FROM projections WHERE method = ?",
    params = list(method))
  empty <- data.frame(source_dataset_id = character(0), target_dataset_id = character(0), median_r2 = numeric(0))
  if (nrow(d) == 0) return(empty)
  agg <- aggregate(mean_r_squared ~ source_dataset_id + target_dataset_id, data = d, FUN = median)
  names(agg)[3] <- "median_r2"
  agg
}

#' Dataset x dataset gene-space agreement for one method -- Overview Tab
#' 4's second view, reading `gene_space_agreement` (populated offline by
#' R/backfill_gene_space_agreement.R, NOT computed live -- see that
#' script's header for why). Deliberately sparser than
#' dataset_similarity_matrix()'s sample-space R²: only pairs the backfill
#' has actually been run for exist here at all (by design, scoped to
#' within-family pairs first -- see the plan). Empty data.frame (not an
#' error) if the backfill hasn't been run yet for this method.
gene_space_similarity_matrix <- function(con, method) {
  DBI::dbGetQuery(con,
    "SELECT source_dataset_id, target_dataset_id, mean_abs_diagonal, n_genes_matched
     FROM gene_space_agreement WHERE method = ?",
    params = list(method))
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

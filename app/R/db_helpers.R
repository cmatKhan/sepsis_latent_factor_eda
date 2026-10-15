# Queries against the stability DB (schema: R/db/schema.sql, docs/database.qmd).
# The app only reads; every function takes the open connection first.

METRICS <- c("cosine", "pearson", "spearman")

# Methods in display order. WGCNA has modules, not loadings.
METHOD_ORDER <- c("pca", "spca", "ica", "nmf", "cogaps", "wgcna")

#' Known methods in METHOD_ORDER, then any other method found in the database
#' (alphabetical), so a new method appears without an app edit.
order_methods <- function(m) {
  m <- unique(m[!is.na(m)])
  extra <- sort(setdiff(m, METHOD_ORDER))
  if (length(extra) > 0) note_missing("a display name and colour for method(s)", paste(extra, collapse = ", "))
  c(METHOD_ORDER[METHOD_ORDER %in% m], extra)
}

check_metric <- function(metric) {
  if (!metric %in% METRICS) stop("metric must be one of: ", paste(METRICS, collapse = ", "))
  metric
}

## ---- datasets and coverage ---------------------------------------------------

list_datasets <- function(con) {
  DBI::dbGetQuery(con, "SELECT dataset_id FROM datasets ORDER BY dataset_id")$dataset_id
}

#' Methods with at least one ok fit for a dataset, in METHOD_ORDER then any others.
methods_for_dataset <- function(con, dataset_id) {
  m <- DBI::dbGetQuery(con,
    "SELECT DISTINCT method FROM fits WHERE dataset_id = ? AND status = 'ok'",
    params = list(dataset_id))$method
  order_methods(m)
}

#' Fit counts per method for one dataset (v_ingest_summary).
fit_coverage <- function(con, dataset_id) {
  DBI::dbGetQuery(con,
    "SELECT method, family, n_fits, n_ok, n_failed, n_representative
     FROM v_ingest_summary WHERE dataset_id = ?",
    params = list(dataset_id))
}

dataset_val_sd <- function(con, dataset_id) {
  v <- DBI::dbGetQuery(con, "SELECT val_sd FROM matrix_diagnostics WHERE dataset_id = ?",
                       params = list(dataset_id))$val_sd
  if (length(v) == 1 && !is.na(v)) v else NA_real_
}

## ---- artifacts ---------------------------------------------------------------

#' Artifact paths are stored relative to the DB file's directory.
resolve_artifact <- function(path) {
  if (is.na(path) || !nzchar(path)) return(NA_character_)
  if (startsWith(path, "/")) return(path)
  db_dir <- getOption("stability.db_dir")
  stopifnot(!is.null(db_dir))
  file.path(db_dir, path)
}

load_fit_artifact <- function(con, fit_id, kind) {
  p <- DBI::dbGetQuery(con, "SELECT path FROM fit_artifacts WHERE fit_id = ? AND kind = ?",
                       params = list(fit_id, kind))$path
  if (length(p) != 1) return(NULL)
  path <- resolve_artifact(p)
  if (!file.exists(path)) return(NULL)
  readRDS(path)
}

#' A fit's loadings (genes x factors); NULL for WGCNA. `n` keeps the first n
#' factors, which is how PCA's single max-rank fit is viewed at a smaller rank.
load_loadings <- function(con, fit_id, n = NULL) {
  L <- load_fit_artifact(con, fit_id, "loadings")
  if (is.null(L) || is.null(n) || is.na(n)) return(L)
  L[, seq_len(min(n, ncol(L))), drop = FALSE]
}

#' A fit's sample scores (samples x factors; module eigengenes for WGCNA).
load_scores <- function(con, fit_id, n = NULL) {
  S <- load_fit_artifact(con, fit_id, "scores")
  if (is.null(S) || is.null(n) || is.na(n)) return(S)
  S[, seq_len(min(n, ncol(S))), drop = FALSE]
}

## ---- Q1/Q2: hyperparameter curves and the optimal fit --------------------------

CURVE_COLS <- c("method", "rank_label", "rank_key", "fit_id", "n_factors",
                "in_sample_r2", "stability")

empty_curve <- function() {
  data.frame(method = character(0), rank_label = character(0), rank_key = numeric(0),
             fit_id = integer(0), n_factors = integer(0), in_sample_r2 = numeric(0),
             stability = numeric(0), stringsAsFactors = FALSE)
}

#' One method's quality/stability curve over its hyperparameter, one row per
#' representative fit (`fits.representative = 1`). `rank_key` is the rank, or
#' the soft-threshold power for WGCNA; `rank_label` names which. Columns:
#'   in_sample_r2  1 - mse / val_sd^2 (PCA, NMF, CoGAPS, ICA); sPCA has no mse,
#'                 so its sparse PEV; WGCNA has no reconstruction, so the mean
#'                 member kME. Comparable within a method, not across methods.
#'   stability     mean over factors of `factors.stability_<metric>` (median
#'                 matched |similarity| across same-rank fit pairs); PCA's
#'                 nested components are stable by construction (1); WGCNA uses
#'                 the mean ARI against the other powers.
#' An empty frame (not an error) when the method has no ok fits here.
method_curve <- function(con, dataset_id, method, metric = "cosine", val_sd = NA_real_) {
  stab_col <- paste0("stability_", check_metric(metric))

  if (method == "pca") {
    f <- DBI::dbGetQuery(con,
      "SELECT fit_id, rank FROM fits WHERE dataset_id = ? AND method = 'pca' AND status = 'ok'
       ORDER BY rank DESC LIMIT 1", params = list(dataset_id))
    if (nrow(f) == 0) return(empty_curve())
    d <- load_fit_artifact(con, f$fit_id, "diag")
    if (is.null(d) || is.null(d$mse_by_rank)) return(empty_curve())
    k <- seq_along(d$mse_by_rank)
    return(data.frame(method = "pca", rank_label = paste("rank", k), rank_key = k, fit_id = f$fit_id,
                      n_factors = k, in_sample_r2 = 1 - d$mse_by_rank / val_sd^2, stability = 1,
                      stringsAsFactors = FALSE))
  }

  f <- DBI::dbGetQuery(con,
    "SELECT fit_id, rank, power, mse, n_factors FROM v_fits
     WHERE dataset_id = ? AND method = ? AND status = 'ok' AND representative = 1",
    params = list(dataset_id, method))
  if (nrow(f) == 0) return(empty_curve())
  ids <- paste(f$fit_id, collapse = ",")

  if (method == "wgcna") {
    kme <- DBI::dbGetQuery(con, sprintf(
      "SELECT fit_id, SUM(n * mean) / SUM(n) AS q FROM module_kme_summary
       WHERE grp = 'member' AND module != 0 AND fit_id IN (%s) GROUP BY fit_id", ids))
    ari <- DBI::dbGetQuery(con, sprintf(
      "SELECT fit_a, fit_b, ari FROM fit_pairs WHERE fit_a IN (%s) AND fit_b IN (%s)", ids, ids))
    ari_long <- rbind(data.frame(fit_id = ari$fit_a, ari = ari$ari), data.frame(fit_id = ari$fit_b, ari = ari$ari))
    stab <- if (nrow(ari_long) > 0) aggregate(ari ~ fit_id, ari_long, mean) else data.frame(fit_id = integer(0), ari = numeric(0))
    f <- f[order(f$power), ]
    return(data.frame(method = "wgcna", rank_label = paste("power", f$power), rank_key = f$power,
                      fit_id = f$fit_id, n_factors = f$n_factors,
                      in_sample_r2 = kme$q[match(f$fit_id, kme$fit_id)],
                      stability = stab$ari[match(f$fit_id, stab$fit_id)], stringsAsFactors = FALSE))
  }

  stab <- DBI::dbGetQuery(con, sprintf(
    "SELECT fit_id, AVG(%s) AS s FROM factors WHERE fit_id IN (%s) GROUP BY fit_id", stab_col, ids))
  if (method == "spca") {
    q <- DBI::dbGetQuery(con, sprintf(
      "SELECT fit_id, value FROM fit_metrics WHERE name = 'pev_sparse' AND fit_id IN (%s)", ids))
    r2 <- q$value[match(f$fit_id, q$fit_id)]
  } else {
    # mse averaged over the seed sweep at this rank would need every fit; the
    # representative is the lowest-mse seed, so this is the best-case quality.
    r2 <- 1 - f$mse / val_sd^2
  }
  f <- cbind(f, r2 = r2, s = stab$s[match(f$fit_id, stab$fit_id)])
  f <- f[order(f$rank), ]
  data.frame(method = method, rank_label = paste("rank", f$rank), rank_key = f$rank, fit_id = f$fit_id,
             n_factors = f$n_factors, in_sample_r2 = f$r2, stability = f$s, stringsAsFactors = FALSE)
}

#' Curves for every present method of one dataset, stacked.
dataset_curves <- function(con, dataset_id, methods = methods_for_dataset(con, dataset_id),
                           metric = "cosine") {
  val_sd <- dataset_val_sd(con, dataset_id)
  out <- lapply(methods, function(m) method_curve(con, dataset_id, m, metric, val_sd))
  do.call(rbind, c(list(empty_curve()), out))
}

#' The optimal row of one method's curve.
#'   "stable": among hyperparameters whose stability is at least `threshold`,
#'             the one with the highest in-sample quality; if none qualifies,
#'             the most stable one. Quality rises with rank, so this is the
#'             largest rank that is still reproducible across seeds.
#'   "max_r2": the highest in-sample quality, ignoring stability.
#' Returns a one-row data.frame, or NULL for an empty curve.
pick_optimal <- function(curve, rule = c("stable", "max_r2"), threshold = 0.9) {
  rule <- match.arg(rule)
  d <- curve[!is.na(curve$in_sample_r2), , drop = FALSE]
  if (nrow(d) == 0) return(NULL)
  if (rule == "max_r2") return(d[which.max(d$in_sample_r2), , drop = FALSE])
  ok <- !is.na(d$stability) & d$stability >= threshold
  if (any(ok)) return(d[ok, , drop = FALSE][which.max(d$in_sample_r2[ok]), , drop = FALSE])
  if (all(is.na(d$stability))) return(d[which.max(d$in_sample_r2), , drop = FALSE])
  d[which.max(d$stability), , drop = FALSE]
}

#' The optimal row for every method in `curves` (dataset_curves() output).
optimal_fits <- function(curves, rule = "stable", threshold = 0.9) {
  rows <- lapply(split(curves, curves$method), pick_optimal, rule = rule, threshold = threshold)
  out <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(out)) return(empty_curve())
  rownames(out) <- NULL
  out[order(match(out$method, METHOD_ORDER)), ]
}

#' Median out-of-sample R2 (projectR) of each method's optimal fit-family
#' sources, split by within-/cross-dataset targets. Medians are over every
#' projection whose source fit belongs to this dataset and method.
oos_r2_by_method <- function(con, dataset_id) {
  d <- DBI::dbGetQuery(con,
    "SELECT f.method, p.projection_type, p.mean_r_squared
     FROM projections p JOIN fits f ON f.fit_id = p.source_fit_id
     WHERE f.dataset_id = ?", params = list(dataset_id))
  if (nrow(d) == 0) return(data.frame(method = character(0), oos_r2_within = numeric(0), oos_r2_cross = numeric(0)))
  agg <- aggregate(mean_r_squared ~ method + projection_type, d, median)
  methods <- unique(agg$method)
  pick <- function(type) {
    a <- agg[agg$projection_type == type, ]
    a$mean_r_squared[match(methods, a$method)]
  }
  data.frame(method = methods, oos_r2_within = pick("within_dataset"), oos_r2_cross = pick("cross_dataset"))
}

## ---- Q4: cross-dataset similarity ----------------------------------------------

#' Long dataset x dataset similarity for one method.
#'   "r2":         median projectR out-of-sample R2 of source-dataset fits on
#'                 the target dataset's samples.
#'   "gene_space": mean |diagonal| of the loadings similarity between the two
#'                 datasets' representative fits, on shared genes.
#' Columns: source, target, value, n (projections or matched genes).
cross_dataset_similarity <- function(con, method, kind = c("r2", "gene_space")) {
  kind <- match.arg(kind)
  if (kind == "r2") {
    d <- DBI::dbGetQuery(con,
      "SELECT f.dataset_id AS source, p.target_dataset_id AS target, p.mean_r_squared AS value
       FROM projections p JOIN fits f ON f.fit_id = p.source_fit_id
       WHERE f.method = ? AND p.mean_r_squared IS NOT NULL", params = list(method))
    if (nrow(d) == 0) return(data.frame(source = character(0), target = character(0), value = numeric(0), n = integer(0)))
    agg <- aggregate(value ~ source + target, d, median)
    agg$n <- as.integer(aggregate(value ~ source + target, d, length)$value)
    return(agg)
  }
  d <- DBI::dbGetQuery(con,
    "SELECT a.dataset_id AS source, b.dataset_id AS target, g.mean_abs_diagonal AS value,
            g.n_genes_matched AS n
     FROM gene_space_agreement g
     JOIN fits a ON a.fit_id = g.source_fit_id JOIN fits b ON b.fit_id = g.target_fit_id
     WHERE a.method = ?", params = list(method))
  if (nrow(d) == 0) return(data.frame(source = character(0), target = character(0), value = numeric(0), n = integer(0)))
  agg <- aggregate(value ~ source + target, d, mean)
  agg$n <- aggregate(n ~ source + target, d, max)$n
  agg
}

#' Square, symmetrised dataset x dataset matrix from cross_dataset_similarity()
#' (the two directions are averaged; missing pairs stay NA), with rows ordered
#' by hierarchical clustering (NA treated as 0 for the ordering only).
similarity_matrix <- function(long) {
  if (nrow(long) == 0) return(NULL)
  ds <- sort(unique(c(long$source, long$target)))
  m <- matrix(NA_real_, length(ds), length(ds), dimnames = list(ds, ds))
  m[cbind(match(long$source, ds), match(long$target, ds))] <- long$value
  sym <- (m + t(m)) / 2
  sym[is.na(sym)] <- pmax(m, t(m), na.rm = TRUE)[is.na(sym)]   # one direction only: keep it
  diag(sym) <- NA
  if (length(ds) > 2) {
    z <- sym; z[is.na(z)] <- 0
    sym <- sym[hclust(dist(z))$order, hclust(dist(z))$order]
  }
  sym
}

## ---- Explore fits: one method on one dataset ------------------------------------

#' Every fit of a method on a dataset (ok or failed), with hyperparameters.
method_fits <- function(con, dataset_id, method) {
  DBI::dbGetQuery(con,
    "SELECT fit_id, rank, seed, bootstrap, representative, mse, n_factors, status, power, para
     FROM v_fits WHERE dataset_id = ? AND method = ?", params = list(dataset_id, method))
}

#' Matched |similarity| of every factor across same-rank fit pairs (the
#' values `factors.stability_*` summarises). `kind` says whether the pair is
#' two seeds or involves a bootstrap resample (ICA).
same_rank_matches <- function(con, dataset_id, method, metric = "cosine") {
  DBI::dbGetQuery(con,
    "SELECT fa.rank AS rank, p.fit_pair_id, ABS(s.value) AS value,
            CASE WHEN COALESCE(fa.bootstrap, 0) + COALESCE(fb.bootstrap, 0) > 0
                 THEN 'with bootstrap' ELSE 'seeds only' END AS kind
     FROM fit_pairs p
     JOIN fits fa ON fa.fit_id = p.fit_a JOIN fits fb ON fb.fit_id = p.fit_b
     JOIN factor_matches m ON m.fit_pair_id = p.fit_pair_id
     JOIN match_scores s ON s.match_id = m.match_id AND s.metric = ?
     WHERE fa.dataset_id = ? AND fa.method = ? AND p.same_rank = 1",
    params = list(check_metric(metric), dataset_id, method))
}

#' Mean 95th percentile of the unmatched-cell |similarity| per rank: how large
#' a match has to be to beat chance for that pair of fits.
same_rank_null <- function(con, dataset_id, method, metric = "cosine") {
  DBI::dbGetQuery(con,
    "SELECT fa.rank AS rank, AVG(n.p95) AS p95
     FROM fit_pairs p JOIN fits fa ON fa.fit_id = p.fit_a
     JOIN fit_pair_null n ON n.fit_pair_id = p.fit_pair_id AND n.metric = ?
     WHERE fa.dataset_id = ? AND fa.method = ? AND p.same_rank = 1 GROUP BY fa.rank",
    params = list(check_metric(metric), dataset_id, method))
}

#' Seed x seed matrix of mean matched |similarity| at one rank (non-bootstrap fits).
seed_pair_matrix <- function(con, dataset_id, method, rank, metric = "cosine") {
  d <- DBI::dbGetQuery(con,
    "SELECT fa.seed AS seed_a, fb.seed AS seed_b, AVG(ABS(s.value)) AS value
     FROM fit_pairs p
     JOIN fits fa ON fa.fit_id = p.fit_a JOIN fits fb ON fb.fit_id = p.fit_b
     JOIN factor_matches m ON m.fit_pair_id = p.fit_pair_id
     JOIN match_scores s ON s.match_id = m.match_id AND s.metric = ?
     WHERE fa.dataset_id = ? AND fa.method = ? AND fa.rank = ? AND p.same_rank = 1
       AND COALESCE(fa.bootstrap, 0) = 0 AND COALESCE(fb.bootstrap, 0) = 0
     GROUP BY p.fit_pair_id",
    params = list(check_metric(metric), dataset_id, method, rank))
  if (nrow(d) == 0) return(NULL)
  seeds <- as.character(sort(unique(c(d$seed_a, d$seed_b))))
  m <- matrix(NA_real_, length(seeds), length(seeds), dimnames = list(seeds, seeds))
  m[cbind(match(d$seed_a, seeds), match(d$seed_b, seeds))] <- d$value
  m[cbind(match(d$seed_b, seeds), match(d$seed_a, seeds))] <- d$value
  m
}

#' Per-factor stability of one fit.
fit_factor_stability <- function(con, fit_id, metric = "cosine") {
  col <- paste0("stability_", check_metric(metric))
  DBI::dbGetQuery(con, sprintf(
    "SELECT factor_index, label, %s AS stability, n_genes FROM factors WHERE fit_id = ? ORDER BY factor_index", col),
    params = list(fit_id))
}

#' PCA's per-rank reconstruction error and cumulative variance explained,
#' from the diag artifact of its single max-rank fit.
pca_scree <- function(con, dataset_id) {
  f <- DBI::dbGetQuery(con,
    "SELECT fit_id FROM fits WHERE dataset_id = ? AND method = 'pca' AND status = 'ok' ORDER BY rank DESC LIMIT 1",
    params = list(dataset_id))
  if (nrow(f) == 0) return(NULL)
  d <- load_fit_artifact(con, f$fit_id, "diag")
  if (is.null(d)) return(NULL)
  k <- seq_along(d$mse_by_rank)
  data.frame(rank = k, mse = d$mse_by_rank, pev = (cumsum(d$sdev^2) / sum(d$sdev^2))[k])
}

#' sPCA's Index-of-Sparseness inputs for every fit: rank, para, PEV of the
#' sparse and ordinary components, proportion sparse, IS, and whether chosen.
spca_is_table <- function(con, dataset_id) {
  d <- DBI::dbGetQuery(con,
    "SELECT f.fit_id, f.rank, f.para, m.name, m.value FROM v_fits f
     JOIN fit_metrics m ON m.fit_id = f.fit_id
     WHERE f.dataset_id = ? AND f.method = 'spca' AND f.status = 'ok'
       AND m.name IN ('pev_sparse', 'pev_pca', 'prop_sparse', 'index_of_sparseness', 'is_selected')",
    params = list(dataset_id))
  if (nrow(d) == 0) return(NULL)
  w <- reshape(d, idvar = c("fit_id", "rank", "para"), timevar = "name", direction = "wide")
  names(w) <- sub("^value\\.", "", names(w))
  w
}

#' ICASSO cluster quality (Iq) per rank, one row per cluster.
icasso_iq <- function(con, dataset_id) {
  DBI::dbGetQuery(con,
    "SELECT rank, cluster_id, iq, n_members FROM icasso_clusters WHERE dataset_id = ? ORDER BY rank, cluster_id",
    params = list(dataset_id))
}

#' WGCNA: ARI between every pair of powers (symmetric matrix).
wgcna_ari_matrix <- function(con, dataset_id) {
  d <- DBI::dbGetQuery(con,
    "SELECT ra.power AS power_a, rb.power AS power_b, p.ari
     FROM fit_pairs p JOIN v_fits ra ON ra.fit_id = p.fit_a JOIN v_fits rb ON rb.fit_id = p.fit_b
     WHERE ra.dataset_id = ? AND ra.method = 'wgcna'", params = list(dataset_id))
  if (nrow(d) == 0) return(NULL)
  pw <- as.character(sort(unique(c(d$power_a, d$power_b))))
  m <- matrix(NA_real_, length(pw), length(pw), dimnames = list(pw, pw))
  m[cbind(match(d$power_a, pw), match(d$power_b, pw))] <- d$ari
  m[cbind(match(d$power_b, pw), match(d$power_a, pw))] <- d$ari
  diag(m) <- 1
  m
}

#' WGCNA per-power summary: scale-free fit, modules, unassigned genes.
wgcna_power_summary <- function(con, dataset_id) {
  DBI::dbGetQuery(con,
    "SELECT f.power, f.n_factors AS n_modules, s.sft_r_sq, u.value AS n_unassigned
     FROM v_fits f
     LEFT JOIN wgcna_sft s ON s.dataset_id = f.dataset_id AND s.power = f.power
     LEFT JOIN fit_metrics u ON u.fit_id = f.fit_id AND u.name = 'n_unassigned'
     WHERE f.dataset_id = ? AND f.method = 'wgcna' AND f.status = 'ok' ORDER BY f.power",
    params = list(dataset_id))
}

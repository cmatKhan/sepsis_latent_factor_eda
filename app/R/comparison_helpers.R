# PCA-as-baseline comparison helpers: used from NMF/CoGAPS's and WGCNA's
# Level 2 "Compare to PCA" tab. Everything here is computed live from
# already-ingested loadings_file/scores_file artifacts -- no schema/ingest
# changes, matching the rest of this app's "read what's there" philosophy.
#
# Requires prepare_loadings()/pair_similarities() from
# R/lib/ingest/similarity.R to be sourced first (see app.R).

#' Best (min) masking-CV MSE at a given rank, across whatever grid that
#' method's masking-CV sweeps (CoGAPS: min over alpha; PCA/NMF: the single
#' value). NA if no masking-CV fit exists at that rank.
maskcv_best_mse <- function(con, dataset_id, method, rank) {
  d <- maskcv_curve(con, dataset_id, method)
  d <- d[!is.na(d$rank) & d$rank == rank, ]
  if (nrow(d) == 0 || all(is.na(d$mse))) return(NA_real_)
  min(d$mse, na.rm = TRUE)
}

#' Like hungarian_match() (R/lib/ingest/similarity.R) but costed on
#' |similarity| rather than raw similarity -- needed whenever a strongly
#' NEGATIVE similarity is just as informative a match as a strongly
#' positive one (e.g. a non-negative CoGAPS/NMF factor aligning with the
#' negative tail of a signed PCA component). Plain hungarian_match() would
#' systematically avoid such matches, since it maximizes raw (signed)
#' similarity.
hungarian_match_abs <- function(sim) {
  transposed <- FALSE
  if (nrow(sim) > ncol(sim)) { sim <- t(sim); transposed <- TRUE }
  a_sim <- abs(sim)
  cost <- max(a_sim) - a_sim
  assignment <- clue::solve_LSAP(cost)
  a <- seq_len(nrow(sim))
  b <- as.integer(assignment)
  if (transposed) cbind(a = b, b = a) else cbind(a = a, b = b)
}

#' Per-factor comparison of `fit_id` (method `method`, at `rank`) against
#' the PCA fit at the SAME rank, in gene-loadings space. Returns NULL if no
#' PCA fit exists at that rank. `component`/`pca_component` are 1-based
#' factor indices; `cosine`/`pearson`/`spearman` are SIGNED (a negative
#' value already tells you the match is to the opposite-signed side of the
#' PCA component -- no separate pos/neg bookkeeping needed).
compare_loadings_to_pca <- function(con, dataset_id, fit_id, rank) {
  pca_fits <- fits_at_rank(con, dataset_id, "pca", rank)
  if (nrow(pca_fits) == 0) return(NULL)
  pca_fit_id <- pca_fits$fit_id[1]

  L <- load_loadings(con, fit_id); if (is.null(L)) return(NULL)
  Lp <- load_loadings(con, pca_fit_id); if (is.null(Lp)) return(NULL)

  prep_a <- prepare_loadings(L)
  prep_b <- prepare_loadings(Lp)
  sims <- pair_similarities(prep_a, prep_b)
  if (is.null(sims)) return(NULL)

  match_idx <- hungarian_match_abs(sims$cosine)
  best <- data.frame(
    component = match_idx[, "a"],
    pca_component = match_idx[, "b"],
    cosine   = sims$cosine[match_idx],
    pearson  = sims$pearson[match_idx],
    spearman = sims$spearman[match_idx]
  )
  best <- best[order(best$component), ]
  attr(best, "pca_fit_id") <- pca_fit_id
  best
}

#' WGCNA-eigengenes-vs-PCA-scores analog of compare_loadings_to_pca(): plain
#' Pearson/Spearman correlation on shared samples (correlation is already
#' signed, so "negative side of a PC" shows up directly as a negative
#' value -- no absolute-value bookkeeping needed for the similarity itself,
#' only for finding each module's best match).
compare_scores_to_pca <- function(con, other_fit_id, pca_fit_id) {
  Lo <- load_scores(con, other_fit_id); if (is.null(Lo)) return(NULL)
  Lp <- load_scores(con, pca_fit_id);   if (is.null(Lp)) return(NULL)
  shared <- intersect(rownames(Lo), rownames(Lp))
  if (length(shared) < 3) return(NULL)
  Lo <- Lo[shared, , drop = FALSE]; Lp <- Lp[shared, , drop = FALSE]

  cos_p <- suppressWarnings(cor(Lo, Lp, method = "pearson"))
  cos_s <- suppressWarnings(cor(Lo, Lp, method = "spearman"))
  match_idx <- hungarian_match_abs(cos_p)
  best <- data.frame(
    component = match_idx[, "a"],
    pca_component = match_idx[, "b"],
    pearson  = cos_p[match_idx],
    spearman = cos_s[match_idx]
  )
  best[order(best$component), ]
}

#' Everything already cached (via Level 3) for one fit's factors -- per
#' factor/direction/query_type, count of significant terms + the smallest
#' p-value seen. Purely a live read of what's already been queried; never
#' triggers new g:Profiler calls.
cached_enrichment_summary <- function(con, fit_id) {
  DBI::dbGetQuery(con,
    "SELECT f.factor_index, q.query_type, q.direction,
            COUNT(c.term_id) AS n_significant_terms,
            MIN(c.p_value) AS min_p_value
     FROM factors f
     JOIN enrichment_queried q ON q.factor_id = f.factor_id
     LEFT JOIN enrichment_cache c
       ON c.factor_id = q.factor_id AND c.query_type = q.query_type AND c.direction = q.direction
     WHERE f.fit_id = ?
     GROUP BY f.factor_index, q.query_type, q.direction
     ORDER BY f.factor_index",
    params = list(fit_id))
}

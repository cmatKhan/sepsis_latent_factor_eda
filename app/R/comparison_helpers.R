# Cross-method comparison helpers, used by the standalone "Compare
# methods" screen (see app.R) -- everything here is computed live from
# already-ingested loadings_file/scores_file artifacts, no schema/ingest
# changes needed for any of it (WGCNA module-enrichment is the one
# exception, handled separately in R/lib/ingest/ -- see its own comment).
#
# Requires prepare_loadings()/pair_similarities() from
# R/lib/ingest/similarity.R to be sourced first (see app.R).

#' The rank with the lowest in-sample reconstruction MSE (best over the
#' alpha/para grid too, for sPCA) -- used as the DEFAULT (always
#' overridable) rank selection for PCA/sPCA (SCREE_RANK_METHODS, app.R).
#' NA if no ok fit exists for this method.
scree_best_rank <- function(con, dataset_id, method) {
  d <- scree_mse_by_rank(con, dataset_id, method)
  d <- d[!is.na(d$rank) & !is.na(d$mse), ]
  if (nrow(d) == 0) return(NA_integer_)
  best_by_rank <- aggregate(mse ~ rank, data = d, FUN = min)
  best_by_rank$rank[which.min(best_by_rank$mse)]
}

#' Per-K sparsity/quality tradeoff curve for sPCA, following the "Index
#' of Sparseness" (IS) approach (Trendafilov 2014, shown by Gu et al.
#' 2019 to outperform cross-validation/BIC for choosing sPCA's sparsity
#' penalty -- see Guerra-Urzola et al. 2021, "A Guide for Sparse PCA:
#' Model Comparison and Applications," Figs 6-7, PMC8636462). Replaces
#' the old rank x alpha MSE heatmap -- that view was both hard to act on
#' AND, until 2026-09-22, built on a buggy `mse` (see R/methods/spca.R's
#' comment) that reported ~0% variance explained for every fit.
#'
#' For a FIXED rank K, across every para/alpha this dataset's spca_grid
#' swept:
#'   PEV_sparse = sum(pev)  -- cumulative variance explained by the
#'     K-component sparse solution. spca_diag_file$pev is PER-COMPONENT
#'     (each entry is that one component's own incremental contribution),
#'     not cumulative -- summing it is what R/methods/spca.R's `mse` fix
#'     does too, for the same underlying reason.
#'   PS = 1 - sum(n_nonzero) / (n_genes * K)  -- proportion of sparsity
#'     (fraction of zero loadings across the whole K-component solution).
#'   PEV_pca = 1 - (mse_pca(K) * n_genes * n_samples) / var_all  -- the
#'     UNCONSTRAINED PCA reference at the same K, derived from this
#'     dataset's own PCA fit (`fits.mse`, a raw per-element MSE -- PCA's
#'     mse has no analogous bug, only sPCA's does) and any sPCA fit's
#'     `var_all` (mathematically the same total sum-of-squares
#'     `elasticnet::spca()` and `prcomp()` both compute from the
#'     identically-centered matrix -- verified numerically against real
#'     data, 2026-09-22: monotonically increasing, sane 0-1 values).
#'     Constant across the para sweep for this K -- no PCA re-fit needed.
#'   IS = PEV_sparse * PEV_pca * PS -- peaks at the paper's recommended
#'     "not too sparse, not too dense" sweet spot; a sparse solution can
#'     never explain more variance than unconstrained PCA at the same K,
#'     so PEV_sparse <= PEV_pca always holds by construction.
#'
#' Returns data.frame(fit_id, alpha, PS, PEV_sparse, IS) sorted by alpha
#' (IS is NA, not the whole row dropped, when no PCA fit exists at this
#' exact rank for this dataset -- callers should still plot PEV_sparse
#' vs PS in that case and just skip/note the IS series).
spca_sparsity_curve <- function(con, dataset_id, rank) {
  empty <- data.frame(fit_id = integer(0), alpha = numeric(0), PS = numeric(0),
                       PEV_sparse = numeric(0), IS = numeric(0))

  fits <- DBI::dbGetQuery(con,
    "SELECT fit_id, alpha FROM fits
     WHERE dataset_id = ? AND method = 'spca' AND rank = ? AND status = 'ok'
     ORDER BY alpha",
    params = list(dataset_id, rank))
  if (nrow(fits) == 0) return(empty)

  pca_fit <- DBI::dbGetQuery(con,
    "SELECT fit_id, mse FROM fits
     WHERE dataset_id = ? AND method = 'pca' AND rank = ? AND status = 'ok' LIMIT 1",
    params = list(dataset_id, rank))
  have_pca_ref <- nrow(pca_fit) == 1 && !is.na(pca_fit$mse)
  pca_n_samples <- if (have_pca_ref) {
    s <- load_scores(con, pca_fit$fit_id[1])
    if (!is.null(s)) nrow(s) else NA_integer_
  } else NA_integer_
  have_pca_ref <- have_pca_ref && !is.na(pca_n_samples)

  rows <- lapply(seq_len(nrow(fits)), function(i) {
    fit_id <- fits$fit_id[i]
    diag <- load_spca_diag(con, fit_id)
    if (is.null(diag) || is.null(diag$pev) || all(is.na(diag$pev)) ||
        is.null(diag$n_nonzero) || is.null(diag$var_all) || is.na(diag$var_all)) return(NULL)

    L <- load_loadings(con, fit_id)
    if (is.null(L)) return(NULL)
    n_genes <- nrow(L)

    pev_sparse <- sum(diag$pev)
    ps <- 1 - sum(diag$n_nonzero) / (n_genes * rank)

    is_val <- NA_real_
    if (have_pca_ref) {
      pev_pca <- 1 - (pca_fit$mse[1] * n_genes * pca_n_samples) / diag$var_all
      is_val <- pev_sparse * pev_pca * ps
    }
    data.frame(fit_id = fit_id, alpha = fits$alpha[i], PS = ps, PEV_sparse = pev_sparse, IS = is_val)
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) return(empty)
  do.call(rbind, rows)
}

#' Like hungarian_match() (R/lib/ingest/similarity.R) but costed on
#' |similarity| rather than raw similarity -- needed whenever a strongly
#' NEGATIVE similarity is just as informative a match as a strongly
#' positive one (e.g. a non-negative CoGAPS/NMF factor aligning with the
#' negative tail of a signed PCA component). Plain hungarian_match() would
#' systematically avoid such matches, since it maximizes raw (signed)
#' similarity. Also used on Jaccard matrices (always non-negative, so this
#' degenerates to ordinary Hungarian matching there).
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

#' Per-factor comparison between any two fits' gene loadings (continuous
#' cosine/Pearson/Spearman on the shared gene set). `component`/
#' `other_component` are 1-based factor indices; similarities are SIGNED
#' (a negative value already tells you the match is to the opposite-signed
#' side of the other fit's factor -- no separate pos/neg bookkeeping
#' needed for THIS metric; see top_genes()/jaccard_matrix() below for the
#' sign-aware gene-SET comparison you asked for specifically).
compare_loadings <- function(con, fit_a, fit_b) {
  L <- load_loadings(con, fit_a); if (is.null(L)) return(NULL)
  Lb <- load_loadings(con, fit_b); if (is.null(Lb)) return(NULL)

  prep_a <- prepare_loadings(L)
  prep_b <- prepare_loadings(Lb)
  sims <- pair_similarities(prep_a, prep_b)
  if (is.null(sims)) return(NULL)

  match_idx <- hungarian_match_abs(sims$cosine)
  best <- data.frame(
    component = match_idx[, "a"],
    other_component = match_idx[, "b"],
    cosine   = sims$cosine[match_idx],
    pearson  = sims$pearson[match_idx],
    spearman = sims$spearman[match_idx]
  )
  best[order(best$component), ]
}

#' Full factor x factor signed-cosine matrix between two fits (for the
#' heatmap -- compare_loadings() above only returns the best match per
#' factor).
loadings_similarity_matrix <- function(con, fit_a, fit_b) {
  L <- load_loadings(con, fit_a); if (is.null(L)) return(NULL)
  Lb <- load_loadings(con, fit_b); if (is.null(Lb)) return(NULL)
  sims <- pair_similarities(prepare_loadings(L), prepare_loadings(Lb))
  if (is.null(sims)) return(NULL)
  sims$cosine
}

#' Top (or bottom) `n` gene ids of one loading vector by `sign`:
#' "both" = top n by |loading|; "positive"/"negative" = top n by the
#' loading's own largest/smallest (most negative) values. This is what
#' lets you ask for e.g. "PCA's negative side vs CoGAPS's positive side"
#' specifically, unlike the signed-cosine comparison above.
top_genes <- function(loadings_vector, n, sign = c("both", "positive", "negative")) {
  sign <- match.arg(sign)
  v <- loadings_vector
  ord <- switch(sign,
    both     = order(abs(v), decreasing = TRUE),
    positive = order(v, decreasing = TRUE),
    negative = order(v, decreasing = FALSE)
  )
  names(v)[ord[seq_len(min(n, length(v)))]]
}

#' Named list of gene sets, one per factor (column) of a loadings matrix,
#' via top_genes(). Column names become the list's names.
gene_sets_from_loadings <- function(L, n, sign = "both") {
  sets <- lapply(seq_len(ncol(L)), function(j) top_genes(L[, j], n, sign))
  names(sets) <- colnames(L) %||% paste0("factor_", seq_len(ncol(L)))
  sets
}

#' WGCNA's analog of gene_sets_from_loadings(): one gene set per module,
#' its FULL membership (no n/sign concept -- module assignment is
#' categorical, not a continuous loading to threshold). List names are the
#' bare module number (as a string) -- deliberately NOT prefixed, so
#' `as.integer(names(...))` recovers the real module number cleanly (this
#' is what `factors.factor_index` is set to for WGCNA fits at ingest time).
gene_sets_from_wgcna <- function(con, fit_id) {
  d <- DBI::dbGetQuery(con, "SELECT gene, module FROM wgcna_modules WHERE fit_id = ? AND module != 0",
                        params = list(fit_id))
  if (nrow(d) == 0) return(list())
  split(d$gene, d$module)
}

#' One "side" of the Compare-methods screen: gene sets + the REAL
#' `factor_index` each one corresponds to (module number for WGCNA,
#' 1-based loading-column position for every other method -- exactly what
#' `factors.factor_index` holds, so a selection here can jump straight
#' into Level 3/the WGCNA module view) + a display label. Centralizing
#' this avoids fragile regex-parsing of column-name conventions that
#' differ by method (PCA's "PC1", NMF's "Pattern_1", sPCA/CP/Tucker's
#' "Component_1", ...).
comparator_gene_sets <- function(con, method, fit_id, n = 50, sign = "both") {
  if (method == "wgcna") {
    raw <- gene_sets_from_wgcna(con, fit_id)
    if (length(raw) == 0) return(NULL)
    list(sets = raw, factor_index = as.integer(names(raw)), label = paste0("module ", names(raw)))
  } else {
    L <- load_loadings(con, fit_id)
    if (is.null(L)) return(NULL)
    raw <- gene_sets_from_loadings(L, n, sign)
    list(sets = raw, factor_index = seq_along(raw), label = names(raw) %||% paste0("factor_", seq_along(raw)))
  }
}

jaccard <- function(a, b) {
  u <- union(a, b)
  if (length(u) == 0) return(NA_real_)
  length(intersect(a, b)) / length(u)
}

#' Jaccard overlap for every (gene set in `sets_a`) x (gene set in
#' `sets_b`) pair -- the core of the top-N gene-set comparison. Works
#' identically whether the sets came from top_genes() (loadings-based
#' methods) or gene_sets_from_wgcna() (module membership).
jaccard_matrix <- function(sets_a, sets_b) {
  m <- matrix(NA_real_, length(sets_a), length(sets_b),
              dimnames = list(names(sets_a), names(sets_b)))
  for (i in seq_along(sets_a)) {
    for (j in seq_along(sets_b)) m[i, j] <- jaccard(sets_a[[i]], sets_b[[j]])
  }
  m
}

#' Pearson/Spearman correlation between any two fits' sample (or
#' subject-mode, for CP/Tucker -- excluded from the app's comparator, but
#' this function itself doesn't care) score matrices, aligned on shared
#' rownames. Correlation is already signed, so "negative side of a
#' component" shows up directly as a negative value here -- no
#' absolute-value bookkeeping needed for the similarity itself, only for
#' finding each factor's best match.
compare_scores <- function(con, fit_a, fit_b) {
  La <- load_scores(con, fit_a); if (is.null(La)) return(NULL)
  Lb <- load_scores(con, fit_b); if (is.null(Lb)) return(NULL)
  shared <- intersect(rownames(La), rownames(Lb))
  if (length(shared) < 3) return(NULL)
  La <- La[shared, , drop = FALSE]; Lb <- Lb[shared, , drop = FALSE]

  cos_p <- suppressWarnings(cor(La, Lb, method = "pearson"))
  cos_s <- suppressWarnings(cor(La, Lb, method = "spearman"))
  match_idx <- hungarian_match_abs(cos_p)
  best <- data.frame(
    component = match_idx[, "a"],
    other_component = match_idx[, "b"],
    pearson  = cos_p[match_idx],
    spearman = cos_s[match_idx]
  )
  list(matrix = cos_p, best = best[order(best$component), ])
}

#' Cluster each of two fits' samples independently (hierarchical, `k`
#' clusters) and compare the two cluster assignments on their shared
#' samples via Adjusted Rand Index -- "do these two methods' factors
#' induce similar sample groupings?" NULL if fewer than `k` shared samples.
cluster_and_ari <- function(con, fit_a, fit_b, k = 4) {
  La <- load_scores(con, fit_a); if (is.null(La)) return(NULL)
  Lb <- load_scores(con, fit_b); if (is.null(Lb)) return(NULL)
  shared <- intersect(rownames(La), rownames(Lb))
  if (length(shared) < max(k, 3)) return(NULL)
  La <- La[shared, , drop = FALSE]; Lb <- Lb[shared, , drop = FALSE]

  cl_a <- cutree(hclust(dist(La)), k = k)
  cl_b <- cutree(hclust(dist(Lb)), k = k)
  list(
    ari = mclust::adjustedRandIndex(cl_a, cl_b),
    crosstab = as.data.frame.matrix(table(cluster_a = cl_a, cluster_b = cl_b)),
    n_shared = length(shared)
  )
}

#' Everything already cached (via Level 3, or the WGCNA module view) for
#' one fit's factors -- per factor/direction/query_type, count of
#' significant terms + the smallest p-value seen. Purely a live read of
#' what's already been queried; never triggers new enrichment computation.
#'
#' Resolves each factor through resolve_enrichment_factor_id() (see
#' db_helpers.R) rather than joining directly on this fit_id's own
#' factor_id -- a non-representative fit (e.g. any PCA rank below the
#' max, or an NMF/CoGAPS/ICA/sPCA seed that happens to match another
#' seed's loadings exactly) never has its OWN enrichment_queried rows
#' (representative_fit_ids() only runs fgsea_grid once per equivalence
#' class), but its factors' enrichment is still available via whichever
#' sibling fit IS representative -- confirmed directly (2026-09-22) this
#' was previously the reason the Compare tab's "Enrichment cross-
#' reference" showed "Nothing computed yet" for e.g. any non-max-rank
#' PCA fit even though real enrichment existed for that exact factor.
cached_enrichment_summary <- function(con, fit_id) {
  empty <- data.frame(factor_index = integer(0), query_type = character(0), direction = character(0),
                       n_significant_terms = integer(0), min_p_value = numeric(0))
  fac_idx <- DBI::dbGetQuery(con, "SELECT factor_index FROM factors WHERE fit_id = ? ORDER BY factor_index",
                              params = list(fit_id))$factor_index
  if (length(fac_idx) == 0) return(empty)

  rows <- lapply(fac_idx, function(fi) {
    resolved_id <- resolve_enrichment_factor_id(con, fit_id, fi)
    if (is.null(resolved_id)) return(NULL)
    d <- DBI::dbGetQuery(con,
      "SELECT q.query_type, q.direction,
              COUNT(c.term_id) AS n_significant_terms,
              MIN(c.p_value) AS min_p_value
       FROM enrichment_queried q
       LEFT JOIN enrichment_cache c
         ON c.factor_id = q.factor_id AND c.query_type = q.query_type AND c.direction = q.direction
       WHERE q.factor_id = ?
       GROUP BY q.query_type, q.direction",
      params = list(resolved_id))
    if (nrow(d) == 0) return(NULL)
    cbind(factor_index = fi, d)
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) return(empty)
  do.call(rbind, rows)
}

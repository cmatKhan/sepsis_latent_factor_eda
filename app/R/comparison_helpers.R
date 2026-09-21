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
#' what's already been queried; never triggers new g:Profiler calls.
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

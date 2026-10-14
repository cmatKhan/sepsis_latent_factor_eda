# Queries and computations behind the Factor detail, Enrichment, Metadata and
# WGCNA pages. Same conventions as db_helpers.R: read-only, connection first.

## ---- metadata and gene labels -----------------------------------------------------

# arrow must be loaded before reading the metadata artifacts: they were saved
# with arrow-backed string columns, which readRDS() cannot unserialize
# without it (the columns come back empty).
if (requireNamespace("arrow", quietly = TRUE)) suppressMessages(library(arrow))

dataset_artifact <- function(con, dataset_id, kind) {
  p <- DBI::dbGetQuery(con, "SELECT path FROM dataset_artifacts WHERE dataset_id = ? AND kind = ?",
                       params = list(dataset_id, kind))$path
  if (length(p) != 1) return(NULL)
  path <- resolve_artifact(p)
  if (!file.exists(path)) return(NULL)
  readRDS(path)
}

#' Sample metadata with the id column renamed `sample_id`; NULL if absent.
dataset_sample_metadata <- function(con, dataset_id) {
  m <- dataset_artifact(con, dataset_id, "sample_metadata")
  id <- DBI::dbGetQuery(con, "SELECT sample_id_col FROM datasets WHERE dataset_id = ?",
                        params = list(dataset_id))$sample_id_col
  if (is.null(m) || length(id) != 1 || !id %in% names(m)) return(NULL)
  m <- as.data.frame(m)
  names(m)[names(m) == id] <- "sample_id"
  m$sample_id <- as.character(m$sample_id)
  m
}

#' Named vector feature id -> display label (the dataset's symbol column, else
#' the id itself). Cached for the session.
.gene_label_cache <- new.env(parent = emptyenv())
gene_labels <- function(con, dataset_id) {
  if (!is.null(hit <- .gene_label_cache[[dataset_id]])) return(hit)
  out <- NULL
  d <- DBI::dbGetQuery(con, "SELECT feature_id_col, symbol_col FROM datasets WHERE dataset_id = ?",
                       params = list(dataset_id))
  fm <- tryCatch(dataset_artifact(con, dataset_id, "feature_metadata"), error = function(e) NULL)
  if (!is.null(fm) && nrow(d) == 1 && all(c(d$feature_id_col, d$symbol_col) %in% names(fm))) {
    ids <- as.character(fm[[d$feature_id_col]]); sym <- as.character(fm[[d$symbol_col]])
    ok <- !is.na(ids) & length(sym) == length(ids)
    if (all(ok)) out <- setNames(ifelse(is.na(sym) | !nzchar(sym), ids, sym), ids)
  }
  .gene_label_cache[[dataset_id]] <- if (is.null(out)) character(0) else out
  .gene_label_cache[[dataset_id]]
}

label_genes <- function(labels, ids) {
  out <- unname(labels[ids])
  ifelse(is.na(out), ids, out)
}

#' Association of every score column with every metadata field: Spearman for
#' numeric fields, Kruskal-Wallis for everything else (factors, logicals,
#' text). Fields with fewer than two usable values or levels are skipped.
#' `padj` is BH over the whole returned grid. Columns: component, field, test,
#' statistic (rho, or the H statistic), p_value, padj.
association_scan <- function(scores, meta, id_col = "sample_id") {
  empty <- data.frame(component = character(0), field = character(0), test = character(0),
                      statistic = numeric(0), p_value = numeric(0), padj = numeric(0))
  if (is.null(scores) || is.null(meta) || !id_col %in% names(meta)) return(empty)
  shared <- intersect(rownames(scores), meta[[id_col]])
  if (length(shared) < 3) return(empty)
  scores <- scores[shared, , drop = FALSE]
  meta <- meta[match(shared, meta[[id_col]]), , drop = FALSE]
  comps <- colnames(scores) %||% as.character(seq_len(ncol(scores)))
  rows <- list()
  for (field in setdiff(names(meta), id_col)) {
    v <- meta[[field]]
    if ((is.numeric(v) || is.integer(v)) && !is.factor(v)) {
      ok <- !is.na(v)
      if (sum(ok) < 3 || length(unique(v[ok])) < 2) next
      r <- suppressWarnings(stats::cor(scores[ok, , drop = FALSE], v[ok], method = "spearman",
                                       use = "pairwise.complete.obs"))[, 1]
      n <- colSums(!is.na(scores[ok, , drop = FALSE]))
      keep <- !is.na(r) & n >= 3
      if (!any(keep)) next
      t <- r[keep] * sqrt((n[keep] - 2) / pmax(1 - r[keep]^2, 1e-12))
      rows[[length(rows) + 1]] <- data.frame(component = comps[keep], field = field, test = "spearman",
                                             statistic = unname(r[keep]),
                                             p_value = 2 * stats::pt(-abs(t), n[keep] - 2))
    } else {
      f <- droplevels(as.factor(v))
      for (j in seq_along(comps)) {
        ok <- !is.na(f) & !is.na(scores[, j])
        if (sum(ok) < 3 || nlevels(droplevels(f[ok])) < 2) next
        kt <- tryCatch(stats::kruskal.test(scores[ok, j] ~ droplevels(f[ok])), error = function(e) NULL)
        if (is.null(kt)) next
        rows[[length(rows) + 1]] <- data.frame(component = comps[j], field = field, test = "kruskal",
                                               statistic = unname(kt$statistic), p_value = kt$p.value)
      }
    }
  }
  if (length(rows) == 0) return(empty)
  out <- do.call(rbind, rows)
  out$padj <- p.adjust(out$p_value, method = "BH")
  out
}

## ---- factors of one fit -----------------------------------------------------------

#' Factors of a fit; `n` keeps the first n (PCA's nested fit at a smaller rank).
fit_factors <- function(con, fit_id, n = NULL) {
  d <- DBI::dbGetQuery(con,
    "SELECT factor_id, factor_index, label, n_genes, stability_cosine, stability_pearson,
            stability_spearman, kurtosis, excess_kurtosis
     FROM factors WHERE fit_id = ? ORDER BY factor_index", params = list(fit_id))
  if (!is.null(n) && !is.na(n)) d <- d[d$factor_index <= n, , drop = FALSE]
  d
}

#' Top genes of one loadings column on each side. Returns data.frame(gene,
#' loading, side). Non-negative methods (all loadings >= 0) only have "positive".
top_loading_genes <- function(L, j, n = 20) {
  v <- L[, j]
  pos <- head(order(v, decreasing = TRUE), n)
  d <- data.frame(gene = rownames(L)[pos], loading = v[pos], side = "positive", stringsAsFactors = FALSE)
  if (any(v < 0)) {
    neg <- head(order(v), n)
    d <- rbind(d, data.frame(gene = rownames(L)[neg], loading = v[neg], side = "negative", stringsAsFactors = FALSE))
  }
  d
}

#' One factor's match in every other fit of its method (both orientations of
#' `v_factor_matches`): other fit's rank/seed, whether it is the same rank,
#' and the signed value of `metric` (WGCNA modules are matched by gene
#' Jaccard only, across powers).
factor_matches_everywhere <- function(con, factor_id, metric = "cosine") {
  DBI::dbGetQuery(con,
    "SELECT m.other_rank, f.power AS other_power, m.other_seed, m.other_fit_id, m.other_factor_index,
            m.same_rank, m.value, m.runner_up, m.margin
     FROM v_factor_matches m JOIN v_fits f ON f.fit_id = m.other_fit_id
     WHERE m.factor_id = ? AND m.metric = ?",
    params = list(factor_id, if (metric == "jaccard") metric else check_metric(metric)))
}

## ---- enrichment ---------------------------------------------------------------------

QUERY_LABELS <- c(fgsea = "fgsea (ranked loadings)", gsea = "GSEA", ora = "over-representation",
                  cogaps_fora = "CoGAPS marker over-representation")

#' query_type x direction combinations with results for any factor of a fit.
enrichment_queries_for_fit <- function(con, fit_id, n = NULL) {
  d <- DBI::dbGetQuery(con,
    "SELECT DISTINCT q.query_type, q.direction FROM enrichment_queries q
     JOIN factors x ON x.factor_id = q.factor_id WHERE x.fit_id = ? AND x.factor_index <= ?",
    params = list(fit_id, n %||% 1e9))
  d[order(d$query_type, d$direction), , drop = FALSE]
}

#' Stored enrichment terms for one factor (only BH-adjusted p < 0.05 are kept).
factor_enrichment <- function(con, factor_id, query_type, direction) {
  DBI::dbGetQuery(con,
    "SELECT r.source, r.term_id, r.term_name, r.p_value, r.nes, r.intersection_size, r.term_size, r.is_main_pathway
     FROM enrichment_results r JOIN enrichment_queries q ON q.query_id = r.query_id
     WHERE q.factor_id = ? AND q.query_type = ? AND q.direction = ? ORDER BY r.p_value",
    params = list(factor_id, query_type, direction))
}

#' Enrichment for every factor of a fit in one query (long format with factor_index/label).
fit_enrichment <- function(con, fit_id, query_type, direction, n = NULL) {
  DBI::dbGetQuery(con,
    "SELECT x.factor_index, x.label, r.source, r.term_name, r.p_value, r.nes, r.is_main_pathway
     FROM factors x JOIN enrichment_queries q ON q.factor_id = x.factor_id
     JOIN enrichment_results r ON r.query_id = q.query_id
     WHERE x.fit_id = ? AND x.factor_index <= ? AND q.query_type = ? AND q.direction = ?",
    params = list(fit_id, n %||% 1e9, query_type, direction))
}

## ---- WGCNA ---------------------------------------------------------------------------

wgcna_fits <- function(con, dataset_id) {
  DBI::dbGetQuery(con,
    "SELECT fit_id, power, n_factors FROM v_fits
     WHERE dataset_id = ? AND method = 'wgcna' AND status = 'ok' ORDER BY power", params = list(dataset_id))
}

#' Gene -> module assignment with kME for one fit.
wgcna_gene_modules <- function(con, fit_id) {
  DBI::dbGetQuery(con,
    "SELECT gene, module, kme_own, next_module, kme_next, margin, p_analytic, p_perm
     FROM gene_modules WHERE fit_id = ?", params = list(fit_id))
}

#' kME summary of members, non-members and the permutation null per module.
wgcna_kme_summary <- function(con, fit_id) {
  DBI::dbGetQuery(con,
    "SELECT module, grp, n, mean, median, p95 FROM module_kme_summary WHERE fit_id = ? AND module != 0",
    params = list(fit_id))
}

wgcna_kme_histograms <- function(con, fit_id, module) {
  DBI::dbGetQuery(con,
    "SELECT grp, bin_lo, bin_hi, count FROM module_kme_histograms WHERE fit_id = ? AND module = ?",
    params = list(fit_id, module))
}

#' hclust object for the gene dendrogram of a WGCNA fit (from the artifact).
wgcna_dendrogram <- function(con, fit_id) {
  d <- load_fit_artifact(con, fit_id, "dendrogram")
  if (is.null(d)) return(NULL)
  structure(list(merge = d$merge, height = d$height, order = d$order, labels = d$labels,
                 method = "average", dist.method = "unknown"), class = "hclust")
}

#' Fields with gene-significance tests for a dataset.
gene_significance_fields <- function(con, dataset_id) {
  DBI::dbGetQuery(con, "SELECT DISTINCT field, test FROM gene_significance WHERE dataset_id = ? ORDER BY field",
                  params = list(dataset_id))
}

gene_significance_for <- function(con, dataset_id, field) {
  DBI::dbGetQuery(con, "SELECT gene, statistic, p_value FROM gene_significance WHERE dataset_id = ? AND field = ?",
                  params = list(dataset_id, field))
}

#' Jaccard overlap of the modules of two WGCNA fits (rows: fit a, cols: fit b).
wgcna_module_overlap <- function(con, fit_a, fit_b) {
  ma <- wgcna_gene_modules(con, fit_a); mb <- wgcna_gene_modules(con, fit_b)
  sa <- split(ma$gene[ma$module != 0], ma$module[ma$module != 0])
  sb <- split(mb$gene[mb$module != 0], mb$module[mb$module != 0])
  if (length(sa) == 0 || length(sb) == 0) return(NULL)
  names(sa) <- paste0("ME", names(sa)); names(sb) <- paste0("ME", names(sb))
  jaccard_matrix(sa, sb)
}

## ---- biplots -------------------------------------------------------------------------------

#' Gene biplot data: all sample scores plus the `top_n` genes with the largest
#' combined loading on the two chosen components, scaled so the longest arrow
#' reaches 80% of the score cloud's radius (loadings and scores live on
#' different scales, so only direction and relative length carry meaning).
biplot_data <- function(L, S, cx, cy, top_n = 15) {
  scores <- data.frame(sample_id = rownames(S), x = S[, cx], y = S[, cy])
  lx <- L[, cx]; ly <- L[, cy]; mag <- sqrt(lx^2 + ly^2)
  keep <- order(-mag)[seq_len(min(top_n, sum(mag > 0)))]
  sr <- suppressWarnings(max(sqrt(scores$x^2 + scores$y^2), na.rm = TRUE)); lr <- max(mag[keep])
  k <- if (is.finite(lr) && lr > 0) 0.8 * sr / lr else 1
  list(scores = scores, arrows = data.frame(gene = rownames(L)[keep], x = lx[keep] * k, y = ly[keep] * k))
}

#' Metadata biplot data: every usable sample-metadata field on the same two
#' components. Numeric fields become a vector (Spearman correlation with each
#' axis, scaled like the gene arrows); categorical fields become one centroid
#' per level (levels with fewer than `min_n` samples are dropped; fields with
#' fewer than 2 or more than `max_levels` usable levels are skipped).
biplot_metadata_data <- function(S, meta, cx, cy, id_col = "sample_id", min_n = 3, max_levels = 8) {
  shared <- intersect(rownames(S), meta[[id_col]])
  if (length(shared) < min_n) return(NULL)
  S <- S[shared, , drop = FALSE]; meta <- meta[match(shared, meta[[id_col]]), , drop = FALSE]
  x <- S[, cx]; y <- S[, cy]
  scores <- data.frame(sample_id = shared, x = x, y = y)
  vec <- list(); cen <- list()
  for (field in setdiff(names(meta), id_col)) {
    v <- meta[[field]]
    if (is.numeric(v) && !is.factor(v)) {
      ok <- !is.na(v)
      if (sum(ok) < min_n) next
      dx <- suppressWarnings(cor(v[ok], x[ok], method = "spearman")); dy <- suppressWarnings(cor(v[ok], y[ok], method = "spearman"))
      if (!is.na(dx) && !is.na(dy)) vec[[length(vec) + 1]] <- data.frame(field = field, dx = dx, dy = dy)
    } else {
      f <- as.factor(v); tab <- table(f); levs <- names(tab)[tab >= min_n]
      if (length(levs) < 2 || length(levs) > max_levels) next
      for (lv in levs) {
        ok <- !is.na(f) & f == lv
        cen[[length(cen) + 1]] <- data.frame(field = field, level = lv, x = mean(x[ok]), y = mean(y[ok]))
      }
    }
  }
  vectors <- if (length(vec)) do.call(rbind, vec) else data.frame(field = character(0), dx = numeric(0), dy = numeric(0))
  centroids <- if (length(cen)) do.call(rbind, cen) else data.frame(field = character(0), level = character(0), x = numeric(0), y = numeric(0))
  if (nrow(vectors) > 0) {
    sr <- max(sqrt(scores$x^2 + scores$y^2)); mm <- max(sqrt(vectors$dx^2 + vectors$dy^2))
    k <- if (is.finite(mm) && mm > 0) 0.8 * sr / mm else 1
    vectors$x <- vectors$dx * k; vectors$y <- vectors$dy * k
  }
  list(scores = scores, vectors = vectors, centroids = centroids)
}

## ---- trajectory across time and projection comparison --------------------------------------------

#' Within-family projections of one source fit (its basis projected onto the
#' sibling per-timepoint datasets): target dataset and artifact path.
within_projections <- function(con, fit_id) {
  DBI::dbGetQuery(con,
    "SELECT target_dataset_id, path FROM projections
     WHERE source_fit_id = ? AND projection_type = 'within_dataset' ORDER BY target_dataset_id",
    params = list(fit_id))
}

#' Timepoint order for dataset ids of one family: the order listed in
#' config/dataset_families.yml when found, else by the last number in the id
#' (PRE first; DISCHARGE and POST last), else alphabetical. Returns a permutation.
order_timepoints <- function(ids) {
  for (p in c("config/dataset_families.yml", "../config/dataset_families.yml")) {
    if (!file.exists(p)) next
    for (fam in yaml::read_yaml(p)) {
      if (all(ids %in% fam)) return(order(match(ids, fam)))
    }
  }
  key <- rep(NA_real_, length(ids))
  key[grepl("PRE$", ids, ignore.case = TRUE)] <- -Inf
  key[is.na(key) & grepl("DISCHARGE$|POST", ids, ignore.case = TRUE)] <- Inf
  need <- is.na(key)
  num <- function(s) { m <- regmatches(s, gregexpr("\\d+", s))[[1]]; if (length(m)) as.numeric(m[length(m)]) else NA_real_ }
  key[need] <- vapply(ids[need], num, numeric(1))
  order(key, ids, na.last = TRUE)
}

#' Per-subject weight of one projected component at every timepoint dataset
#' (needs each target's subject id column). `component` is the row name of the
#' projection matrix. NULL when nothing is usable.
trajectory_data <- function(con, fit_id, component) {
  tg <- within_projections(con, fit_id)
  rows <- lapply(seq_len(nrow(tg)), function(i) {
    ds <- tg$target_dataset_id[i]
    P <- tryCatch(readRDS(resolve_artifact(tg$path[i]))$projection, error = function(e) NULL)
    if (is.null(P) || !component %in% rownames(P)) return(NULL)
    sc <- DBI::dbGetQuery(con, "SELECT subject_id_col FROM datasets WHERE dataset_id = ?", params = list(ds))$subject_id_col
    m <- dataset_sample_metadata(con, ds)
    if (is.null(m) || length(sc) != 1 || is.na(sc) || !sc %in% names(m)) return(NULL)
    w <- P[component, ]; subj <- m[[sc]][match(names(w), m$sample_id)]; ok <- !is.na(subj)
    if (!any(ok)) return(NULL)
    data.frame(dataset_id = ds, subject_id = as.character(subj[ok]), weight = unname(w[ok]), stringsAsFactors = FALSE)
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) NULL else do.call(rbind, rows)
}

#' Metadata fields of the source dataset that are constant within every
#' subject (so they can colour a subject's trajectory), one row per subject.
subject_level_fields <- function(con, dataset_id, max_levels = 8) {
  sc <- DBI::dbGetQuery(con, "SELECT subject_id_col FROM datasets WHERE dataset_id = ?", params = list(dataset_id))$subject_id_col
  m <- dataset_sample_metadata(con, dataset_id)
  if (is.null(m) || length(sc) != 1 || is.na(sc) || !sc %in% names(m)) return(NULL)
  subj <- as.character(m[[sc]])
  fields <- setdiff(names(m), c("sample_id", sc))
  const <- vapply(fields, function(f) all(vapply(split(m[[f]], subj), function(x) length(unique(x[!is.na(x)])) <= 1, logical(1))), logical(1))
  keep <- Filter(function(f) { v <- m[[f]]; is.numeric(v) || length(unique(na.omit(v))) <= max_levels }, fields[const])
  if (length(keep) == 0) return(NULL)
  out <- data.frame(subject_id = unique(subj), stringsAsFactors = FALSE)
  for (f in keep) out[[f]] <- m[[f]][match(out$subject_id, subj)]
  out
}

projection_targets <- function(con, dataset_id) {
  DBI::dbGetQuery(con,
    "SELECT DISTINCT p.target_dataset_id FROM projections p JOIN fits f ON f.fit_id = p.source_fit_id
     WHERE f.dataset_id = ? ORDER BY p.target_dataset_id", params = list(dataset_id))$target_dataset_id
}

#' Every projection of a dataset's fits onto one target, with the source fit's
#' method/rank/seed and whether it is the representative fit.
projections_for_pair <- function(con, dataset_id, target) {
  DBI::dbGetQuery(con,
    "SELECT p.source_fit_id, f.method, f.rank, f.seed, f.representative, p.projection_type, p.mean_r_squared, p.path
     FROM projections p JOIN fits f ON f.fit_id = p.source_fit_id
     WHERE f.dataset_id = ? AND p.target_dataset_id = ?", params = list(dataset_id, target))
}

#' Least-squares projection of fit A's loadings onto fit B's basis on shared
#' genes (the coefficients projectR computes): rows are B's components, columns
#' A's. NULL when fewer than 3 genes are shared.
loadings_projection <- function(La, Lb) {
  shared <- intersect(rownames(La), rownames(Lb))
  if (length(shared) < 3) return(NULL)
  A <- La[shared, , drop = FALSE]; B <- Lb[shared, , drop = FALSE]
  P <- solve(crossprod(B) + diag(1e-8, ncol(B)), crossprod(B, A))
  list(projection = P, n_genes = length(shared))
}

## ---- pattern drivers ------------------------------------------------------------------------

#' Stored projectionDriveR results for the factors of one fit.
fit_drivers <- function(con, fit_id, n = NULL) {
  d <- DBI::dbGetQuery(con,
    "SELECT d.driver_id, x.factor_index, x.label, d.grouping_col, d.group1_level, d.group2_level, d.mode,
            d.n_genes_considered, d.n_significant_shared, d.path
     FROM pattern_drivers d JOIN factors x ON x.factor_id = d.factor_id
     WHERE x.fit_id = ? ORDER BY x.factor_index, d.grouping_col, d.group1_level, d.group2_level, d.mode",
    params = list(fit_id))
  if (!is.null(n) && !is.na(n)) d <- d[d$factor_index <= n, , drop = FALSE]
  d
}

load_driver <- function(path) {
  p <- resolve_artifact(path)
  if (is.na(p) || !file.exists(p)) NULL else readRDS(p)
}

# Repeated-measures summaries built on the variance module (R/targets/variance.R)
# and kept whole for the app: how a dataset's expression variance splits into
# between-person, time and the rest, for every gene and for the genes the
# factor analysis actually uses; how repeatable a single-timepoint snapshot
# is; whether within-person change shares a direction or follows composition
# or technical changes; and (factor_variance_partition()) the same split for
# every fitted factor's sample scores.

#' Category of each variance-partition term
#'
#' subject; subject-level (variables fixed within a person: sex, age, arm);
#' time; group x time; group (cross-sectional); technical (the design's
#' `vp_covariates`: control-probe factor, sequencing QC); covariate (other
#' covariates); residual.
#'
#' @param design_obj prepare_design() output.
#' @param terms Term names (variance-partition columns).
#' @return Named character vector, one category per term.
variance_term_category <- function(design_obj, terms) {
  spec <- design_obj$spec
  subj_level <- subject_level_variables(design_obj)
  cat <- stats::setNames(rep("covariate", length(terms)), terms)
  cat[terms == "Residuals"] <- "residual"
  if (!is.null(spec[["subject"]])) cat[terms == spec[["subject"]]] <- "subject"
  cat[terms %in% subj_level] <- "subject-level"
  if (!is.null(spec[["time"]])) cat[terms == spec[["time"]]] <- "time"
  cat[terms %in% interaction_term(spec)] <- "group x time"
  if (!is.null(spec[["group"]]) && is.null(spec[["subject"]])) cat[terms == spec[["group"]]] <- "group"
  cat[terms %in% unlist(spec$vp_covariates)] <- "technical"
  cat
}

VARIANCE_CATEGORIES <- c("subject", "subject-level", "group", "time", "group x time",
                         "technical", "covariate", "residual")

#' Time + group x time fraction above which a gene counts as time-driven.
TIME_GENE_THRESHOLD <- 0.10

#' Genes in the top-variance view (the factor-analysis matrices' size).
TOP_N_GENES <- 8000

#' Between-subject and time totals per gene
#'
#' The one definition of the two derived quantities: between = subject +
#' subject-level variables; time = time + group x time.
#'
#' @param design_obj prepare_design() output.
#' @param frac Gene x term fractions (or variances).
#' @return List: `between_subject_total` (when the design has a subject) and
#'   `time_total` (when it has a time), each a per-gene vector.
category_totals <- function(design_obj, frac) {
  cats <- variance_term_category(design_obj, names(frac))
  total <- function(which) rowSums(frac[, names(cats)[cats %in% which], drop = FALSE])
  out <- list()
  if (!is.null(design_obj$spec[["subject"]])) out$between_subject_total <- total(c("subject", "subject-level"))
  if (!is.null(design_obj$spec[["time"]])) out$time_total <- total(c("time", "group x time"))
  out
}

#' Squared correlation of every row of a matrix with one vector
#'
#' @param m Matrix (rows = genes).
#' @param x Vector, one value per column.
#' @return Per-row R^2 (`NA` when `x` is constant or missing).
row_r2 <- function(m, x) {
  if (is.null(x) || stats::sd(x) == 0) return(rep(NA_real_, nrow(m)))
  mc <- m - rowMeans(m)
  xc <- x - mean(x)
  (drop(mc %*% xc) / (sqrt(rowSums(mc^2)) * sqrt(sum(xc^2))))^2
}

#' Which model genes the factor analysis uses
#'
#' The factor-analysis matrix (`mat_<dataset>`) and the variance models can
#' key genes differently (arrays: probe ids vs symbols); matched directly
#' first, else through gene symbols.
#'
#' @param mat The dataset's factor-analysis matrix.
#' @param meta The resolved dataset config.
#' @param genes The variance model's gene table (`id`, `symbol`).
#' @return Logical, one per row of `genes`.
factor_analysis_genes <- function(mat, meta, genes) {
  direct <- genes$id %in% rownames(mat)
  if (mean(rownames(mat) %in% genes$id) >= 0.5) return(direct)
  map <- build_symbol_map(meta)
  mat_symbols <- if (is.null(map)) rownames(mat) else unname(map[rownames(mat)])
  genes$symbol %in% stats::na.omit(mat_symbols) | direct
}

#' Within-person change, for subjects sampled more than once
#'
#' Each later sample minus the subject's first: whether that change shares a
#' direction across people (first principal component's share of it) and
#' how much of each gene's change follows the change in leukocyte
#' composition or in the technical covariates (median R^2 over genes).
#'
#' @param E Gene x sample expression the models used.
#' @param data The design's model data.
#' @param spec The design block.
#' @param composition analyze_composition() output.
#' @param genes_used Logical over rows of `E`: the factor-analysis genes.
#' @return One-row data frame, or `NULL` when fewer than 3 subjects repeat.
within_person_change <- function(E, data, spec, composition, genes_used) {
  if (is.null(spec[["subject"]]) || is.null(spec[["time"]])) return(NULL)
  subj <- as.character(data[[spec[["subject"]]]]); tp <- as.integer(data[[spec[["time"]]]])
  first <- tapply(seq_along(subj), subj, function(i) i[which.min(tp[i])])
  later <- which(seq_along(subj) != first[subj])
  if (length(unique(subj[later])) < 3) return(NULL)
  base <- first[subj[later]]
  dE <- E[, later, drop = FALSE] - E[, base, drop = FALSE]
  r2 <- function(x) row_r2(dE, x)
  cs <- composition$samples
  nl <- cs$nl_score[match(rownames(data), cs$sample)]
  d_nl <- if (is.null(nl)) NULL else nl[later] - nl[base]
  tech <- unlist(spec$vp_covariates)
  r2_tech <- if (length(tech)) {
    do.call(pmax, c(lapply(tech, function(v) r2(data[[v]][later] - data[[v]][base])), na.rm = TRUE))
  } else rep(NA_real_, nrow(dE))
  r2_nl <- r2(d_nl)
  pc <- stats::prcomp(t(dE[genes_used, , drop = FALSE]), center = TRUE)
  data.frame(
    n_subjects_repeated = length(unique(subj[later])), n_changes = length(later),
    change_pc1_share = pc$sdev[1]^2 / sum(pc$sdev^2),
    change_pc1_cor_nl = if (is.null(d_nl)) NA_real_ else stats::cor(pc$x[, 1], d_nl),
    median_r2_nl_all = stats::median(r2_nl, na.rm = TRUE),
    median_r2_nl_factor_genes = stats::median(r2_nl[genes_used], na.rm = TRUE),
    median_r2_technical_all = stats::median(r2_tech, na.rm = TRUE),
    median_r2_technical_factor_genes = stats::median(r2_tech[genes_used], na.rm = TRUE))
}

#' Repeated-measures summary of one dataset (the `vp_explore_<dataset>` target)
#'
#' @param dataset_id Dataset id.
#' @param meta The resolved dataset config.
#' @param design_obj prepare_design() output.
#' @param norm normalize_for_models() output.
#' @param vp run_variance_partition() output.
#' @param composition analyze_composition() output.
#' @param mat The dataset's factor-analysis matrix (`mat_<dataset>`).
#' @param n_top Genes in the top-variance view.
#' @return List:
#'   - `genes`: per gene, the category fractions, between (subject +
#'     subject-level), time_all (time + group x time), repeatability =
#'     between / (between + residual), time_within = time_all / (time_all +
#'     residual), total variance, mean expression, `factor_gene` and
#'     `top_variance` flags;
#'   - `categories`: per category, mean and median fraction and
#'     variance-weighted shares over all genes, the factor-analysis genes and
#'     the top `n_top`;
#'   - `headline`: one row of the numbers above per dataset;
#'   - `by_expression`, `by_variance_rank`: median between / time / residual by
#'     expression quartile and by variance-rank bin;
#'   - `change`: within_person_change().
explore_variance <- function(dataset_id, meta, design_obj, norm, vp, composition, mat, n_top = TOP_N_GENES) {
  frac <- vp$varPart; vars <- vp$variances
  cats <- variance_term_category(design_obj, names(frac))
  by_cat <- function(m) {
    out <- sapply(VARIANCE_CATEGORIES, function(cn) {
      cols <- names(cats)[cats == cn]
      if (length(cols)) rowSums(m[, cols, drop = FALSE]) else rep(0, nrow(m))
    })
    matrix(out, nrow = nrow(m), dimnames = list(rownames(m), VARIANCE_CATEGORIES))
  }
  fcat <- by_cat(frac); vcat <- by_cat(vars)
  total <- rowSums(vars)
  E <- model_expression(norm)[rownames(frac), rownames(design_obj$data), drop = FALSE]
  gene_tab <- vp$genes[match(rownames(frac), vp$genes$id), , drop = FALSE]
  used <- factor_analysis_genes(mat, meta, gene_tab)
  top <- rank(-total, ties.method = "first") <= n_top
  between <- fcat[, "subject"] + fcat[, "subject-level"]
  time_all <- fcat[, "time"] + fcat[, "group x time"]
  genes <- data.frame(gene = rownames(frac), symbol = gene_tab$symbol, fcat, check.names = FALSE,
                      between = between, time_all = time_all,
                      repeatability = between / (between + fcat[, "residual"]),
                      time_within = time_all / (time_all + fcat[, "residual"]),
                      total_variance = total, mean_expression = rowMeans(E),
                      factor_gene = used, top_variance = top, row.names = NULL)

  vw <- function(m, keep) colSums(m[keep, , drop = FALSE]) / sum(total[keep])
  categories <- data.frame(
    dataset_id = dataset_id, category = VARIANCE_CATEGORIES,
    present = VARIANCE_CATEGORIES %in% cats,
    mean = colMeans(fcat), median = apply(fcat, 2, stats::median),
    vw_all = vw(vcat, rep(TRUE, nrow(vcat))), vw_factor_genes = if (any(used)) vw(vcat, used) else NA_real_,
    vw_top = vw(vcat, top), row.names = NULL)

  wsum <- function(x, keep) sum((x * total)[keep]) / sum(total[keep])
  has_subject <- !is.null(design_obj$spec[["subject"]])
  headline <- data.frame(
    dataset_id = dataset_id, platform = design_obj$platform, repeated = has_subject,
    n_genes = nrow(frac), n_factor_genes = sum(used), n_top = sum(top),
    between_median = stats::median(between), time_median = stats::median(time_all),
    residual_median = stats::median(fcat[, "residual"]), technical_median = stats::median(fcat[, "technical"]),
    repeatability_median = if (has_subject) stats::median(genes$repeatability, na.rm = TRUE) else NA_real_,
    time_within_median = if (has_subject) stats::median(genes$time_within, na.rm = TRUE) else NA_real_,
    share_time_gt05 = mean(time_all > 0.05), share_time_gt10 = mean(time_all > 0.10),
    share_time_gt20 = mean(time_all > 0.20),
    between_vw_all = wsum(between, TRUE), time_vw_all = wsum(time_all, TRUE),
    residual_vw_all = wsum(fcat[, "residual"], TRUE),
    between_vw_factor_genes = if (any(used)) wsum(between, used) else NA_real_,
    time_vw_factor_genes = if (any(used)) wsum(time_all, used) else NA_real_,
    residual_vw_factor_genes = if (any(used)) wsum(fcat[, "residual"], used) else NA_real_,
    between_vw_top = wsum(between, top), time_vw_top = wsum(time_all, top),
    residual_vw_top = wsum(fcat[, "residual"], top))

  bin_summary <- function(bins, name) {
    out <- do.call(rbind, lapply(split(seq_along(bins), bins), function(i) data.frame(
      dataset_id = dataset_id, bin = as.character(bins[i[1]]), n_genes = length(i),
      between = stats::median(between[i]), time = stats::median(time_all[i]),
      residual = stats::median(fcat[i, "residual"]))))
    names(out)[names(out) == "bin"] <- name
    rownames(out) <- NULL
    out
  }
  ave_expr <- genes$mean_expression
  quart <- cut(ave_expr, unique(stats::quantile(ave_expr, 0:4 / 4)), include.lowest = TRUE,
               labels = FALSE)
  rank_bin <- cut(rank(-total, ties.method = "first"), c(0, 1000, 2000, 4000, n_top, Inf),
                  labels = c("1-1000", "1001-2000", "2001-4000", paste0("4001-", n_top), paste0(">", n_top)))

  list(genes = genes, categories = categories, headline = headline,
       by_expression = bin_summary(paste0("Q", quart), "expression_quartile"),
       by_variance_rank = bin_summary(rank_bin, "variance_rank"),
       change = within_person_change(E, design_obj$data, design_obj$spec, composition, used))
}

#' Repeated-measures summaries across datasets (the `variance_explore` target)
#'
#' @param parts Named list of explore_variance() outputs.
#' @param summaries Named list of summarize_variance() outputs (for the sex
#'   and composition checks shown next to the headline).
#' @param threshold Time + group x time fraction above which a gene counts
#'   as time-driven.
#' @return List: `headline`, `categories`, `by_expression`,
#'   `by_variance_rank`, `change` (stacked across datasets), `time_genes`
#'   (dataset_id, symbol, time_all, between for every time-driven gene) and
#'   `time_gene_recurrence` (symbol, number of datasets, which).
combine_variance_explore <- function(parts, summaries, threshold = TIME_GENE_THRESHOLD) {
  stack <- function(field) {
    rows <- Filter(Negate(is.null), lapply(names(parts), function(id) {
      x <- parts[[id]][[field]]
      if (is.null(x)) NULL else cbind(dataset_id = id, x[setdiff(names(x), "dataset_id")])
    }))
    if (length(rows)) do.call(rbind, rows) else NULL
  }
  checks <- do.call(rbind, lapply(summaries, function(s) s$dataset[, intersect(
    c("dataset_id", "sex_agreement", "n_subjects_sex_discordant", "rho_PC1_nl", "rho_mean_technical"),
    names(s$dataset))]))
  headline <- merge(stack("headline"), checks, by = "dataset_id", all.x = TRUE)
  time_genes <- do.call(rbind, lapply(names(parts), function(id) {
    g <- parts[[id]]$genes
    g <- g[g$time_all > threshold & !is.na(g$symbol), c("symbol", "time_all", "between")]
    if (nrow(g)) cbind(dataset_id = id, g) else NULL
  }))
  recurrence <- if (is.null(time_genes)) NULL else {
    r <- tapply(time_genes$dataset_id, time_genes$symbol, function(x) sort(unique(x)), simplify = FALSE)
    d <- data.frame(symbol = names(r), n_datasets = vapply(r, length, integer(1)),
                    datasets = vapply(r, paste, character(1), collapse = ","), row.names = NULL)
    d[order(-d$n_datasets, d$symbol), ]
  }
  list(headline = headline, categories = stack("categories"), by_expression = stack("by_expression"),
       by_variance_rank = stack("by_variance_rank"), change = stack("change"),
       time_genes = time_genes, time_gene_recurrence = recurrence, threshold = threshold)
}

#' The variance formula for a dataset's fits' scores
#'
#' A parent dataset uses its own design. A per-timepoint subset uses its
#' parent's design restricted to the subset's samples, without the subject,
#' time and group x time terms (one sample per person, one timepoint). Terms
#' whose variable no longer varies are dropped.
#'
#' @param design_obj The (parent) dataset's prepare_design() output.
#' @param samples The fitted samples.
#' @param subset Whether the fits are of a per-timepoint subset.
#' @return List: `formula` (or `NULL` when no term varies), `data` (rows =
#'   `samples` present in the design).
score_design <- function(design_obj, samples, subset) {
  spec <- design_obj$spec
  data <- design_obj$data[intersect(samples, rownames(design_obj$data)), , drop = FALSE]
  # One timepoint, one sample per person: no subject or time terms. Every
  # time field goes, not just `time`.
  if (subset) spec[c("subject", "time", "time_order", "group_at_first_time")] <- NULL
  numeric_vars <- unlist(spec$numeric)
  varies <- function(v) {
    x <- data[[v]]
    if (is.numeric(x)) stats::sd(x) > 0 else length(unique(as.character(x))) > 1
  }
  candidates <- c(spec[["subject"]], spec[["time"]], spec[["group"]], unlist(spec$covariates), unlist(spec$vp_covariates))
  keep <- Filter(varies, intersect(candidates, names(data)))
  terms <- vp_terms(spec)
  terms$random <- terms$random[vapply(strsplit(terms$random, ":"), function(v) all(v %in% keep), logical(1))]
  terms$fixed <- intersect(terms$fixed, keep)
  for (v in setdiff(keep, numeric_vars)) data[[v]] <- droplevels(factor(data[[v]]))
  list(formula = terms_formula(terms), data = data, spec = spec)
}

#' Variance partition of every fit's factor scores (the `factor_vp_<m>_<dataset>` target)
#'
#' Each factor's sample scores go through the same mixed model as the genes
#' (score_design()), so every factor gets its share of between-person,
#' time, group, technical and residual variance: whether it describes a
#' stable person-level trait, a change over time, a technical axis or
#' nothing structured. Fits without sample scores (ICA bootstrap fits) are
#' skipped; WGCNA's unassigned module (ME0) too.
#'
#' @param ingest The method's ingest_method_fits() output.
#' @param design_obj The parent dataset's prepare_design() output.
#' @param subset Whether `ingest` is a per-timepoint subset's.
#' @param db_path Path to the DB (artifact paths).
#' @param cpus Worker processes.
#' @return Data frame (fit_key, factor_index, term, category, fraction), or
#'   `NULL`; the formula used is the `formula` attribute.
factor_variance_partition <- function(ingest, design_obj, subset, db_path, cpus = 1L) {
  ok <- ingest$fits$fit_key[ingest$fits$status == "ok"]
  paths <- vapply(ok, fit_artifact, character(1), ingest = ingest, kind = "scores", db_path = db_path)
  paths <- paths[!is.na(paths)]
  if (!length(paths)) return(NULL)
  rows <- list(); index <- list()
  for (key in names(paths)) {
    s <- as.matrix(readRDS(paths[[key]]))
    idx <- if (ingest$method == "wgcna") suppressWarnings(as.integer(sub("^ME", "", colnames(s)))) else seq_len(ncol(s))
    keep <- !is.na(idx) & idx != 0L & apply(s, 2, function(x) stats::sd(x, na.rm = TRUE) > 0)
    if (!any(keep)) next
    rows[[key]] <- t(s[, keep, drop = FALSE])
    index[[key]] <- data.frame(fit_key = key, factor_index = idx[keep])
  }
  if (!length(rows)) return(NULL)
  samples <- Reduce(intersect, lapply(rows, colnames))
  sd_ <- score_design(design_obj, samples, subset)
  if (is.null(sd_$formula)) return(NULL)
  y <- do.call(rbind, lapply(rows, function(m) m[, rownames(sd_$data), drop = FALSE]))
  idx <- do.call(rbind, index)
  rownames(y) <- paste0("r", seq_len(nrow(y)))
  vp <- variancePartition::fitExtractVarPartModel(y, sd_$formula, sd_$data, REML = FALSE,
                                                  showWarnings = FALSE, BPPARAM = variance_bpparam(cpus))
  vp <- as.data.frame(vp)
  idx <- idx[match(rownames(vp), rownames(y)), , drop = FALSE]
  cats <- variance_term_category(c(design_obj["data"], list(spec = sd_$spec)), names(vp))
  out <- do.call(rbind, lapply(names(vp), function(t) data.frame(
    fit_key = idx$fit_key, factor_index = idx$factor_index, term = t, category = unname(cats[t]),
    fraction = vp[[t]])))
  rownames(out) <- NULL
  attr(out, "formula") <- formula_text(sd_$formula)
  out
}

# Shared helpers for the per-dataset differential-expression scripts under
# R/de/. Each dataset script is a thin wrapper: load its own raw matrix +
# metadata, do whatever small amount of dataset-specific cleanup its real
# metadata needs, then call run_de_timecourse() (the 8 repeated-measures
# datasets) or run_de_case_control() (the other datasets, no time axis) and
# write_de_results(). See docs/methods.qmd ("Differential expression") for the full per-dataset rationale
# (platform, design, contrasts) and why this lives outside the main
# ingest/app pipeline.
#
# Deliberately NOT using R/lib/method_registry.R, R/create_slurm_bundle.R,
# or any ingest/app machinery -- this is a standalone LOCAL analysis (no
# dataset here needs the cluster) that reads the sample/feature metadata
# parquet trio directly, the same way R/preprocessing/*.R scripts do,
# rather than through the app's SQLite DB (populated only after a dataset
# is ingested -- a separate, unrelated pipeline).
#
# Statistical approach: limma (array, log-scale data) / limma+voom
# (RNA-seq, raw counts), blocked on subject via duplicateCorrelation() for
# the repeated-measures datasets -- the standard way to handle a single
# repeated-measures blocking factor without spending a degree of freedom
# per subject (which a subject fixed effect would need, and which breaks
# entirely for an unbalanced/ragged design) or needing a full mixed-model
# package for what is, here, one random effect. See the limma User's Guide
# ch. 15 (voom) / ch. 18 (time series) for the underlying methodology.
#
# Design: rather than an intercept + main-effects + interaction
# parameterization (which requires deciding what the reference level even
# MEANS once a group covariate is added), every dataset uses a "cell-means"
# design -- one coefficient per (group x time) combination actually present
# in the data (`~ 0 + cell`, cell = interaction(group, time)). This is a
# single, uniform recipe that (a) degenerates cleanly to `~ 0 + time` when
# a dataset has no group covariate, (b) handles ragged group x time cells
# (e.g. a healthy-control arm only sampled at two of a sepsis arm's five
# timepoints) with no special-casing, and (c) lets every contrast of
# interest -- the time effect WITHIN one group -- be written directly as a
# difference of cell means via limma::makeContrasts(), rather than fighting
# interaction-term coefficient interpretation.

suppressMessages({
  library(limma)
})

#' Null-or-NA default
#'
#' @param a A value, possibly `NULL` or a single `NA`.
#' @param b The fallback.
#' @return `a`, or `b` when `a` is `NULL` or a single `NA`.
`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a

#' Collapse a probe x sample matrix to one row per gene symbol (microarray only)
#'
#' Keeps, per symbol, the probe with the highest mean expression -- the
#' convention the preprocessing scripts use -- applied to the raw matrix (no
#' top-variance filter, which would bias DE).
#'
#' @param mat Probe x sample matrix.
#' @param feature_meta Feature-metadata data frame.
#' @param feature_id_col Its probe-id column.
#' @param symbol_col Its gene-symbol column.
#' @return Symbol x sample matrix.
collapse_to_symbol <- function(mat, feature_meta, feature_id_col, symbol_col) {
  fm <- feature_meta[!is.na(feature_meta[[symbol_col]]) & nzchar(feature_meta[[symbol_col]]), ]
  fm <- fm[!duplicated(fm[[feature_id_col]]), ]
  fm <- fm[fm[[feature_id_col]] %in% rownames(mat), ]
  mat <- mat[fm[[feature_id_col]], , drop = FALSE]

  mean_expr <- rowMeans(mat, na.rm = TRUE)
  o <- order(fm[[symbol_col]], -mean_expr)
  fm <- fm[o, ]; mat <- mat[o, , drop = FALSE]
  keep <- !duplicated(fm[[symbol_col]])
  mat <- mat[keep, , drop = FALSE]
  rownames(mat) <- fm[[symbol_col]][keep]
  mat
}

#' write_de_results()'s `symbol_map` argument (a feature id -> symbol
#' lookup, for RNA-seq datasets whose native feature id isn't already a
#' readable symbol) is built via R/lib/ingest/symbol_mapping.R's own
#' build_symbol_map(list(dataset = ds_meta), fm = feature_meta) -- see
#' e.g. R/de/GSE110487_de.R. No separate DE-specific version of this
#' lookup exists (a redundant one used to live here -- see git history --
#' and its name collided with the ingest-side function when both got
#' sourced together, one silently shadowing the other with no warning;
#' consolidated onto the single canonical implementation instead of
#' re-diverging).

#' Drop the least variable genes
#'
#' Filters on total variance across every sample, not on the effect being
#' tested, so it isn't circular; it removes probes indistinguishable from noise
#' everywhere, for speed and multiple-testing power (array analogue of
#' edgeR::filterByExpr()).
#'
#' @param mat Gene x sample matrix.
#' @param q Fraction of genes to drop, by variance.
#' @return `mat` without its bottom-`q` genes.
filter_low_variance <- function(mat, q = 0.25) {
  v <- matrixStats::rowVars(mat)
  mat[v > stats::quantile(v, q, na.rm = TRUE), , drop = FALSE]
}

#' Cell-means design
#'
#' One column per (group x time) combination present (just time when `group`
#' is `NULL`). Column names are make.names()-sanitized, since real levels such
#' as "Day: 1" or "Post(24h)" break makeContrasts().
#'
#' @param time Factor of timepoints, one per sample.
#' @param group Factor of groups, or `NULL`.
#' @return `list(design, cell, level_map)`: the design matrix, the cell factor
#'   (readable levels), and readable level -> design column name.
build_cell_means_design <- function(time, group = NULL) {
  cell <- if (is.null(group)) droplevels(factor(time)) else
    droplevels(interaction(group, time, sep = ".", lex.order = TRUE))
  raw_levels <- levels(cell)
  safe_levels <- make.names(raw_levels, unique = TRUE)
  design <- model.matrix(~ 0 + cell)
  colnames(design) <- safe_levels
  list(design = design, cell = cell, level_map = setNames(safe_levels, raw_levels))
}

#' Repeated-measures differential expression
#'
#' `~ 0 + cell` (cell = group x time, or time alone), blocked on subject with
#' duplicateCorrelation() (two-pass for voom). For every group with at least
#' `min_pairs` subjects sampled at 2+ of its timepoints: an F-test across the
#' group's timepoint contrasts (`<group>_time_omnibus`) and each timepoint
#' against the group's earliest. Groups below `min_pairs` stay in the fit but
#' get no contrast.
#'
#' @param mat Feature x sample matrix (collapsed/filtered for the platform).
#' @param sample_meta Sample metadata, matched to `mat` by `sample_id_col`.
#' @param sample_id_col,subject_col,time_col Column names in `sample_meta`.
#' @param group_col Optional group column (e.g. disease arm).
#' @param platform `"array"` (limma) or `"rnaseq"` (edgeR + voom).
#' @param min_pairs Minimum subjects with repeated samples per reported group.
#' @return List with the fit, contrasts, topTables and the inputs
#'   write_de_results() needs.
run_de_timecourse <- function(mat, sample_meta, sample_id_col, subject_col, time_col,
                               group_col = NULL, platform = c("array", "rnaseq"),
                               min_pairs = 3) {
  platform <- match.arg(platform)

  sample_meta <- sample_meta[match(colnames(mat), sample_meta[[sample_id_col]]), ]
  if (anyNA(sample_meta[[sample_id_col]])) {
    stop("mat has columns not present in sample_meta[[sample_id_col]] -- check sample id alignment")
  }

  time <- factor(sample_meta[[time_col]])
  subject <- factor(sample_meta[[subject_col]])
  group <- if (!is.null(group_col)) factor(sample_meta[[group_col]]) else NULL

  cm <- build_cell_means_design(time, group)
  design <- cm$design
  level_map <- cm$level_map

  if (platform == "rnaseq") {
    dge <- edgeR::DGEList(counts = mat)
    keep <- edgeR::filterByExpr(dge, design)
    message("  filterByExpr: keeping ", sum(keep), "/", length(keep), " genes")
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    dge <- edgeR::calcNormFactors(dge, method = "TMM")

    v0 <- voom(dge, design, plot = FALSE)
    dupcor <- duplicateCorrelation(v0, design, block = subject)
    message("  duplicateCorrelation consensus (pass 1): ", round(dupcor$consensus, 3))
    v <- voom(dge, design, block = subject, correlation = dupcor$consensus, plot = FALSE)
    dupcor <- duplicateCorrelation(v, design, block = subject)
    message("  duplicateCorrelation consensus (pass 2): ", round(dupcor$consensus, 3))
    fit <- lmFit(v, design, block = subject, correlation = dupcor$consensus)
    expr_for_diagnostics <- v$E
  } else {
    dupcor <- duplicateCorrelation(mat, design, block = subject)
    message("  duplicateCorrelation consensus: ", round(dupcor$consensus, 3))
    fit <- lmFit(mat, design, block = subject, correlation = dupcor$consensus)
    expr_for_diagnostics <- mat
  }

  group_levels <- if (is.null(group)) "all" else levels(group)
  contrast_exprs <- character(0)
  omnibus_specs <- list()   # name -> character vector of contrast names spanning that group's time effect

  for (g in group_levels) {
    in_g <- if (is.null(group)) rep(TRUE, length(time)) else group == g
    time_g <- droplevels(time[in_g])
    subject_g <- droplevels(subject[in_g])
    time_levels_g <- levels(time_g)
    if (length(time_levels_g) < 2) next

    n_paired <- sum(tapply(as.character(time_g), subject_g, function(x) length(unique(x))) >= 2)
    if (n_paired < min_pairs) {
      message("  [", g, "] only ", n_paired, " subjects with >=2 timepoints (< min_pairs = ",
              min_pairs, ") -- included in the fit, no contrast reported")
      next
    }

    baseline <- time_levels_g[1]
    raw_prefix <- if (is.null(group)) "" else paste0(g, ".")
    # cname (used as both the results-list key and the output CSV
    # filename) is built from make.names()-sanitized labels -- separate
    # from level_map's design-column sanitization, since a level's
    # STANDALONE make.names() (e.g. "Day..1") can differ from how it looks
    # once combined into "group.time" and re-sanitized as a whole.
    safe <- function(x) make.names(x)
    this_contrasts <- character(0)
    for (t in time_levels_g[-1]) {
      cname <- if (is.null(group)) paste0(safe(t), "vs", safe(baseline))
                else paste0(safe(g), "_", safe(t), "vs", safe(baseline))
      cexpr <- paste0(level_map[[paste0(raw_prefix, t)]], " - ", level_map[[paste0(raw_prefix, baseline)]])
      contrast_exprs[cname] <- cexpr
      this_contrasts <- c(this_contrasts, cname)
    }
    omnibus_name <- if (is.null(group)) "time_omnibus" else paste0(safe(g), "_time_omnibus")
    omnibus_specs[[omnibus_name]] <- this_contrasts
  }

  if (length(contrast_exprs) == 0) {
    stop("No group had >= min_pairs (", min_pairs, ") subjects with repeated timepoints -- ",
         "nothing to test. Lower min_pairs or check subject/time column alignment.")
  }

  cmat <- makeContrasts(contrasts = unname(contrast_exprs), levels = design)
  colnames(cmat) <- names(contrast_exprs)
  fit2 <- contrasts.fit(fit, cmat)
  fit2 <- eBayes(fit2, robust = TRUE)

  results <- list()
  for (nm in names(contrast_exprs)) {
    results[[nm]] <- topTable(fit2, coef = nm, number = Inf, sort.by = "P")
  }
  for (nm in names(omnibus_specs)) {
    coefs <- omnibus_specs[[nm]]
    results[[nm]] <- if (length(coefs) == 1) topTable(fit2, coef = coefs, number = Inf, sort.by = "P")
                      else topTable(fit2, coef = coefs, number = Inf, sort.by = "F")
  }

  list(fit = fit, fit2 = fit2, design = design, contrasts = cmat, results = results,
       diagnostics = list(expr = expr_for_diagnostics, group = group, time = time,
                           subject = subject, dupcor = dupcor, platform = platform))
}

#' Cross-sectional differential expression
#'
#' `~ 0 + group`, no blocking.
#'
#' @param mat Feature x sample matrix.
#' @param sample_meta Sample metadata.
#' @param sample_id_col,group_col Column names in `sample_meta`.
#' @param platform `"array"` or `"rnaseq"`.
#' @param contrasts_spec Named makeContrasts()-style expressions (e.g.
#'   `c(SepsisVsHealthy = "Sepsis - Healthy")`); default every pairwise comparison.
#' @return Same shape as run_de_timecourse().
run_de_case_control <- function(mat, sample_meta, sample_id_col, group_col,
                                 platform = c("array", "rnaseq"), contrasts_spec = NULL) {
  platform <- match.arg(platform)
  sample_meta <- sample_meta[match(colnames(mat), sample_meta[[sample_id_col]]), ]
  if (anyNA(sample_meta[[sample_id_col]])) {
    stop("mat has columns not present in sample_meta[[sample_id_col]] -- check sample id alignment")
  }

  group <- droplevels(factor(sample_meta[[group_col]]))
  raw_levels <- levels(group)
  safe_levels <- make.names(raw_levels, unique = TRUE)
  level_map <- setNames(safe_levels, raw_levels)
  design <- model.matrix(~ 0 + group)
  colnames(design) <- safe_levels

  if (platform == "rnaseq") {
    dge <- edgeR::DGEList(counts = mat)
    keep <- edgeR::filterByExpr(dge, design)
    message("  filterByExpr: keeping ", sum(keep), "/", length(keep), " genes")
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    v <- voom(dge, design, plot = FALSE)
    fit <- lmFit(v, design)
    expr_for_diagnostics <- v$E
  } else {
    fit <- lmFit(mat, design)
    expr_for_diagnostics <- mat
  }

  if (is.null(contrasts_spec)) {
    combos <- utils::combn(raw_levels, 2, simplify = FALSE)
    nm <- vapply(combos, function(p) paste0(make.names(p[2]), "vs", make.names(p[1])), character(1))
    contrasts_spec <- setNames(
      vapply(combos, function(p) paste(level_map[[p[2]]], "-", level_map[[p[1]]]), character(1)), nm)
  }
  cmat <- makeContrasts(contrasts = unname(contrasts_spec), levels = design)
  colnames(cmat) <- names(contrasts_spec)
  fit2 <- contrasts.fit(fit, cmat)
  fit2 <- eBayes(fit2, robust = TRUE)

  results <- setNames(
    lapply(names(contrasts_spec), function(nm) topTable(fit2, coef = nm, number = Inf, sort.by = "P")),
    names(contrasts_spec))

  list(fit = fit, fit2 = fit2, design = design, contrasts = cmat, results = results,
       diagnostics = list(expr = expr_for_diagnostics, group = group, platform = platform))
}

#' Write a DE result to disk
#'
#' Every contrast's topTable as CSV, the fit as RDS, `de_bundle.rds`, and
#' diagnostic plots (sample PCA, p-value histograms, voom mean-variance trend
#' for RNA-seq).
#'
#' @param de A run_de_timecourse() or run_de_case_control() result.
#' @param out_dir Output directory (`results/de/<DATASET>/`).
#' @param color_by Sample factor for the PCA plot.
#' @param symbol_map Optional feature id -> symbol map, adding a `symbol`
#'   column (RNA-seq datasets).
#' @return `out_dir`, invisibly.
write_de_results <- function(de, out_dir, color_by = NULL, symbol_map = NULL) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  for (nm in names(de$results)) {
    tt <- de$results[[nm]]
    tt <- data.frame(gene = rownames(tt), tt, row.names = NULL)
    if (!is.null(symbol_map)) {
      tt <- data.frame(gene = tt$gene, symbol = unname(symbol_map[tt$gene]),
                        tt[, setdiff(names(tt), "gene"), drop = FALSE])
    }
    write.csv(tt, file.path(out_dir, paste0("topTable_", nm, ".csv")), row.names = FALSE)
  }
  saveRDS(de$fit2, file.path(out_dir, "fit.rds"))

  # Bundle for later reuse (gene/sample clustering, factor correlation) --
  # `fit` (pre-contrast cell-means fit)'s $coefficients is a genes x
  # (group.time) matrix of per-cell mean expression, exactly what
  # clustering genes by temporal SHAPE needs with no recomputation;
  # `expr` is the per-sample matrix the fit was run on, needed for
  # sample-level clustering. `time`/`subject` are NULL for
  # run_de_case_control() results (no time axis) -- expected, not an
  # error, for the 6 non-time-course datasets.
  saveRDS(list(fit = de$fit, expr = de$diagnostics$expr, time = de$diagnostics$time,
               subject = de$diagnostics$subject, group = de$diagnostics$group),
          file.path(out_dir, "de_bundle.rds"))

  expr <- de$diagnostics$expr
  grDevices::png(file.path(out_dir, "diag_pca.png"), width = 900, height = 700)
  pca <- prcomp(t(expr), scale. = FALSE)
  pct <- round(100 * pca$sdev[1:2]^2 / sum(pca$sdev^2), 1)
  col_factor <- if (!is.null(color_by)) color_by else if (!is.null(de$diagnostics$group)) de$diagnostics$group else de$diagnostics$time
  col_factor <- if (is.null(col_factor)) factor("all") else factor(col_factor)
  plot(pca$x[, 1], pca$x[, 2], col = as.integer(col_factor), pch = 19,
       xlab = paste0("PC1 (", pct[1], "%)"), ylab = paste0("PC2 (", pct[2], "%)"),
       main = "Sample PCA")
  legend("topright", legend = levels(col_factor), col = seq_along(levels(col_factor)), pch = 19)
  grDevices::dev.off()

  grDevices::png(file.path(out_dir, "diag_pval_hist.png"),
                 width = 1000, height = 250 * ceiling(length(de$results) / 3))
  graphics::par(mfrow = c(ceiling(length(de$results) / 3), min(3, length(de$results))))
  for (nm in names(de$results)) {
    graphics::hist(de$results[[nm]]$P.Value, breaks = 40, main = nm, xlab = "p-value", col = "steelblue")
  }
  grDevices::dev.off()

  if (identical(de$diagnostics$platform, "rnaseq")) {
    grDevices::png(file.path(out_dir, "diag_voom_mean_var.png"), width = 800, height = 600)
    limma::plotSA(de$fit2, main = "Post-eBayes mean-variance trend")
    grDevices::dev.off()
  }

  invisible(out_dir)
}

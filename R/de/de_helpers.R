# Shared helpers for the per-dataset differential-expression scripts under
# R/de/. Each dataset script is a thin wrapper: load its own raw matrix +
# metadata, do whatever small amount of dataset-specific cleanup its real
# metadata needs, then call run_de_timecourse() (the 8 repeated-measures
# datasets) or run_de_case_control() (the other datasets, no time axis) and
# write_de_results(). See R/de/README.md for the full per-dataset rationale
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

`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a

#' Collapse a probe/feature x sample matrix to one row per gene SYMBOL --
#' microarray-only (this project's RNA-seq datasets already use gene-level
#' Ensembl feature ids, so collapsing is skipped for those). Keeps, per
#' symbol, the probe with the highest mean expression across all samples --
#' the same convention every existing R/preprocessing/<dataset>_preprocessing.R
#' script already uses for its own probe collapse, applied here to the RAW
#' matrix (no upstream top-variance filtering, which would bias DE toward
#' genes pre-selected for high variance -- see R/de/README.md).
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

#' Drop the bottom `q` quantile of genes by overall variance (across ALL
#' samples, not per-contrast) -- array-only, the standard low-variance
#' pre-filter (akin in spirit to edgeR::filterByExpr() for RNA-seq, e.g.
#' genefilter::varFilter()'s default). This is NOT the same failure mode
#' as the factor-analysis pipeline's top-8000-by-variance filter: filtering
#' on TOTAL variance across every sample regardless of group/time doesn't
#' select for the specific effect this script then tests, it just drops
#' probes/symbols that are indistinguishable from noise everywhere. Mainly
#' matters for speed (duplicateCorrelation()/lmFit() cost scales with gene
#' count) and power (fewer wasted tests -> better FDR) on the larger
#' microarray platforms here (30k+ genes after symbol collapse).
filter_low_variance <- function(mat, q = 0.25) {
  v <- matrixStats::rowVars(mat)
  mat[v > stats::quantile(v, q, na.rm = TRUE), , drop = FALSE]
}

#' Build the cell-means design: one column per (group x time) combination
#' actually present. `group` may be NULL (no condition to include), in
#' which case `cell` is just `time`. Returns list(design, cell, level_map)
#' -- `cell`'s LEVELS are the human-readable "<group>.<time>" (or bare
#' "<time>") labels, but real-world time/group values are frequently not
#' syntactically valid R names (e.g. ANEMONES's "Day: 1", CORTICUS's
#' "Post(24h)" -- confirmed directly, both break makeContrasts() if used
#' as design column names as-is), so `design`'s actual column names are
#' `make.names()`-sanitized. `level_map` (named character vector, names =
#' raw human-readable level, values = sanitized design column name) is how
#' callers build contrast EXPRESSIONS against the real column names while
#' still working with the readable labels everywhere else.
build_cell_means_design <- function(time, group = NULL) {
  cell <- if (is.null(group)) droplevels(factor(time)) else
    droplevels(interaction(group, time, sep = ".", lex.order = TRUE))
  raw_levels <- levels(cell)
  safe_levels <- make.names(raw_levels, unique = TRUE)
  design <- model.matrix(~ 0 + cell)
  colnames(design) <- safe_levels
  list(design = design, cell = cell, level_map = setNames(safe_levels, raw_levels))
}

#' Repeated-measures differential expression: expression ~ 0 + group:time
#' (or ~ 0 + time if `group_col` is NULL), blocked on subject via
#' duplicateCorrelation(). For every group level with at least `min_pairs`
#' subjects sampled at >=2 of that group's timepoints, reports an overall
#' F-test across all of that group's timepoint contrasts (the headline
#' "does expression change over time within this group" result -- an
#' ANOVA-style omnibus test built from CONTRASTS between cell means, not a
#' raw test of the cell means themselves) plus each individual pairwise
#' timepoint contrast against that group's earliest timepoint (for
#' interpretability). A group level with too few paired subjects still
#' contributes its data to the fit (nothing is dropped), it just gets no
#' reported contrast.
#'
#' @param mat feature x sample matrix, already collapsed/filtered as
#'   appropriate for this platform (see collapse_to_symbol()).
#' @param sample_meta one row per sample, matched to mat's columns
#'   internally by `sample_id_col` (any row order is fine).
#' @param platform "array" (plain limma) or "rnaseq" (edgeR + voom).
#' @param min_pairs see above; default 3 is deliberately low (this is a
#'   first pass, not a well-powered study for every dataset -- see
#'   R/de/README.md) but still rules out a "pair" of exactly one subject.
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

#' Simple (no repeated measures) differential expression for cross-sectional
#' datasets: expression ~ 0 + group, no blocking. `contrasts_spec` is a
#' named character vector of makeContrasts()-style expressions (e.g.
#' c(SepsisVsHealthy = "Sepsis - Healthy")); defaults to every pairwise
#' comparison among `group_col`'s levels if not supplied.
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

#' Write every contrast's topTable to CSV, the full eBayes fit to RDS, and
#' a handful of diagnostic PNGs (PCA colored by the main grouping factor,
#' a p-value histogram per contrast, and -- RNA-seq only -- the voom
#' mean-variance trend) -- diagnostics are not optional, see R/de/README.md.
#'
#' `symbol_map` (optional, named character vector: names = mat's rownames
#' i.e. the feature id the fit was run on, values = gene symbol) adds a
#' `symbol` column to every CSV -- purely a display convenience for the
#' RNA-seq datasets (whose feature id is Entrez/Ensembl, not already a
#' readable symbol; the microarray datasets already collapsed to symbol
#' via collapse_to_symbol(), so this is typically NULL for those).
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

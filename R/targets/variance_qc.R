# Quality control and introspection for the variance module
# (R/targets/variance.R), following the variancePartition vignettes: the
# dream-vs-duplicateCorrelation comparison (dream vignette, "Comparing
# p-values"), and per-dataset plots and tables -- variance-fraction violins and
# per-gene bars, expression stratified by the top gene of each term,
# canonical correlations between the design variables, the correlation
# structure of one gene's fit, contrast matrices, voom mean-variance trend,
# p-value histograms, and leukocyte composition versus technical scale
# (R/targets/composition.R). Written to results/variance/<dataset>/.

#' limma duplicateCorrelation() fit on the same data and contrasts as dream
#'
#' The dream vignette's comparison: one genome-wide subject correlation
#' (two voom passes for RNA-seq), lmFit(), contrasts.fit() with dream's
#' contrasts, eBayes(). The `dupcor_<dataset>` target; repeated-measures
#' datasets only.
#'
#' @param expr prepare_expression() output.
#' @param norm normalize_for_models() output (arrays: the matrix fitted).
#' @param design_obj prepare_design() output.
#' @return List: `consensus` (the genome-wide correlation), `p_value` (gene x
#'   contrast matrix), `fit` (eBayes-moderated MArrayLM), `seconds`.
run_duplicate_correlation <- function(expr, norm, design_obj) {
  t0 <- proc.time()[["elapsed"]]
  data <- design_obj$data
  block <- data[[design_obj$spec[["subject"]]]]
  rhs <- paste(c("0", "cell", unlist(design_obj$spec$de_covariates)), collapse = " + ")
  X <- stats::model.matrix(stats::as.formula(paste("~", rhs)), data)
  if (expr$platform == "rnaseq") {
    v <- limma::voom(expr$values, X)
    dc <- limma::duplicateCorrelation(v, X, block = block)
    v <- limma::voom(expr$values, X, block = block, correlation = dc$consensus)
    dc <- limma::duplicateCorrelation(v, X, block = block)
    y <- v
  } else {
    y <- norm$values
    dc <- limma::duplicateCorrelation(y, X, block = block)
  }
  fit <- limma::lmFit(y, X, block = block, correlation = dc$consensus)
  L <- design_obj$contrasts$L
  fit <- limma::eBayes(limma::contrasts.fit(fit, L[colnames(X), , drop = FALSE]))
  list(consensus = dc$consensus, p_value = fit$p.value, fit = fit,
       seconds = proc.time()[["elapsed"]] - t0)
}

#' Expression matrix the models used
#'
#' @param norm normalize_for_models() output.
#' @return Gene x sample matrix (the EList's `E` for RNA-seq).
model_expression <- function(norm) {
  if (inherits(norm$values, "EList")) norm$values$E else norm$values
}

#' Write one ggplot or base-graphics figure
#'
#' @param path Output `.png` path.
#' @param plot A ggplot, or a function drawing a base-graphics plot.
#' @param width,height Size in inches.
#' @return `path`.
save_qc_plot <- function(path, plot, width = 7, height = 5) {
  if (inherits(plot, "ggplot")) {
    ggplot2::ggsave(path, plot, width = width, height = height, dpi = 120)
  } else {
    # ragg: the compute nodes' R has neither X11 nor cairo for grDevices::png().
    ragg::agg_png(path, width = width, height = height, units = "in", res = 120)
    on.exit(grDevices::dev.off())
    plot()
  }
  path
}

#' QC plots and tables for one dataset (the `vp_qc_<dataset>` file target)
#'
#' Every figure is wrapped so one failing plot doesn't lose the rest; failures
#' are listed in `qc_errors.csv`.
#'
#' @param dataset_id Dataset id.
#' @param design_obj prepare_design() output.
#' @param norm normalize_for_models() output.
#' @param vp run_variance_partition() output.
#' @param de dream_tables() output.
#' @param summary summarize_variance() output.
#' @param dupcor run_duplicate_correlation() output, or `NULL` (cross-sectional).
#' @param composition analyze_composition() output, or `NULL`.
#' @param out_root Output root (`results/variance`).
#' @param n_top Genes per term in the per-gene bar plots.
#' @return Paths of the files written.
write_variance_qc <- function(dataset_id, design_obj, norm, vp, de, summary, dupcor = NULL,
                              composition = NULL, out_root = "results/variance", n_top = 20) {
  out <- file.path(out_root, dataset_id)
  unlink(out, recursive = TRUE)
  dir.create(out, recursive = TRUE)
  files <- character(0)
  errors <- character(0)
  attempt <- function(name, expr_) {
    res <- tryCatch(expr_, error = function(e) { errors[name] <<- conditionMessage(e); NULL })
    files <<- c(files, unlist(res))
  }
  path <- function(...) file.path(out, paste0(...))
  spec <- design_obj$spec
  data <- design_obj$data
  vpm <- vp$varPart
  E <- model_expression(norm)
  symbol <- stats::setNames(vp$genes$symbol, vp$genes$id)
  label <- function(g) ifelse(is.na(symbol[g]) | symbol[g] == g, g, paste0(symbol[g], " (", g, ")"))

  # Variance fractions: genome-wide violins and the top genes of each term.
  attempt("violin", save_qc_plot(path("vp_violin.png"),
    variancePartition::plotVarPart(variancePartition::sortCols(vpm), label.angle = 30,
                                   main = paste(dataset_id, formula_text(vp$formula))), width = 8))
  terms <- setdiff(names(vpm), "Residuals")
  for (term in terms) {
    top <- utils::head(order(vpm[[term]], decreasing = TRUE), n_top)
    bars <- vpm[top, , drop = FALSE]
    rownames(bars) <- label(rownames(bars))
    attempt(paste0("bars_", term), save_qc_plot(path("vp_top_", make.names(term), ".png"),
      variancePartition::plotPercentBars(variancePartition::sortCols(bars)) +
        ggplot2::ggtitle(paste("Top", n_top, "genes by", term)), width = 8, height = 6))
  }

  # Expression of the top gene of each term, stratified by that term.
  for (term in terms) {
    g <- rownames(vpm)[which.max(vpm[[term]])]
    strat <- if (grepl(":", term)) "cell" else term
    if (!strat %in% names(data) || is.numeric(data[[strat]])) next
    GE <- data.frame(Expression = E[g, rownames(data)], x = data[[strat]])
    is_subject <- identical(strat, spec[["subject"]])
    args <- list(Expression ~ x, GE, xlab = strat, main = label(g), x.labels = !is_subject,
                 text = paste0(term, ": ", format(100 * vpm[g, term], digits = 3), "%"))
    if (is_subject) args["colorBy"] <- list(NULL)
    attempt(paste0("stratify_", term), save_qc_plot(path("stratify_", make.names(term), ".png"),
      do.call(variancePartition::plotStratify, args), width = if (is_subject) 10 else 7))
  }

  # Correlation between design variables (canonical correlation).
  if (!is.null(design_obj$cancor$matrix)) {
    attempt("cancor", save_qc_plot(path("cancor_pairs.png"),
      function() variancePartition::plotCorrMatrix(design_obj$cancor$matrix), width = 6, height = 6))
  }

  # Correlation structure between samples implied by one gene's fit (the
  # top subject gene, refitted: the per-gene fits aren't stored).
  if (!is.null(spec[["subject"]])) {
    g <- rownames(vpm)[which.max(vpm[[spec[["subject"]]]])]
    attempt("corr_structure", {
      one <- if (inherits(norm$values, "EList")) norm$values[g, ] else norm$values[g, , drop = FALSE]
      fit <- variancePartition::fitVarPartModel(one, vp$formula, data, REML = FALSE,
                                                showWarnings = FALSE)[[1]]
      # As in the vignette, by named variables: the default (every random
      # term) fails on the group:time interaction.
      c(save_qc_plot(path("corr_structure_subject_", make.names(label(g)), ".png"),
                     function() variancePartition::plotCorrStructure(fit, spec[["subject"]]), width = 7, height = 7),
        if (!is.null(spec[["time"]])) save_qc_plot(path("corr_structure_subject_time_", make.names(label(g)), ".png"),
          function() variancePartition::plotCorrStructure(fit, c(spec[["subject"]], spec[["time"]])), width = 7, height = 7))
    })
  }

  # dream: contrast matrix, p-value histograms, dream vs duplicateCorrelation.
  L <- design_obj$contrasts$L
  if (!is.null(L)) {
    attempt("contrasts", save_qc_plot(path("contrast_matrix.png"),
      variancePartition::plotContrasts(L), width = 9, height = 1 + 0.4 * ncol(L)))
    attempt("pvalue_hist", save_qc_plot(path("dream_pvalue_hist.png"), function() {
      cn <- colnames(L)
      graphics::par(mfrow = c(ceiling(length(cn) / 3), min(3, length(cn))))
      for (c in cn) graphics::hist(de$tables[[c]]$P.Value, breaks = 40, main = c, xlab = "p-value")
    }, width = 9, height = 3 * ceiling(ncol(L) / 3)))
  }
  if (!is.null(dupcor)) {
    for (c in colnames(dupcor$p_value)) {
      dream_p <- stats::setNames(de$tables[[c]]$P.Value, rownames(de$tables[[c]]))
      genes <- intersect(names(dream_p), rownames(dupcor$p_value))
      attempt(paste0("compareP_", c), save_qc_plot(path("compareP_", make.names(c), ".png"),
        variancePartition::plotCompareP(dupcor$p_value[genes, c], dream_p[genes],
                                        vpm[genes, spec[["subject"]]], dupcor$consensus) +
          ggplot2::ggtitle(c), width = 6, height = 6))
    }
  }
  if (inherits(norm$values, "EList") && !is.null(norm$values$voom.xy)) {
    attempt("voom_trend", save_qc_plot(path("voom_mean_variance.png"), function() {
      xy <- norm$values$voom.xy
      graphics::plot(xy$x, xy$y, pch = 16, cex = 0.25, xlab = xy$xlab, ylab = xy$ylab,
                     main = "voomWithDreamWeights mean-variance trend")
      graphics::lines(norm$values$voom.line, col = "red")
    }))
  }

  # Composition versus technical scale (R/targets/composition.R).
  if (!is.null(composition)) {
    cs <- composition$samples
    colour <- if (!is.null(spec[["time"]])) spec[["time"]] else spec[["group"]]
    technical <- intersect(c("control_mean", "log2_lib_size"), names(cs))
    technical <- technical[vapply(cs[technical], function(x) any(!is.na(x)), logical(1))]
    scatter <- function(x, y) {
      ggplot2::ggplot(cs, ggplot2::aes(.data[[x]], .data[[y]], colour = .data[[colour]])) +
        ggplot2::geom_point(alpha = 0.7) +
        ggplot2::labs(title = sprintf("%s vs %s (Spearman rho %.2f)", y, x, composition_rho(composition, x, y)))
    }
    for (xy in list(c("nl_score", "PC1"), c("nl_score", "expr_sd"), c("nl_score", "expr_mean"),
                    if (length(technical)) c(technical[1], "expr_sd"),
                    if (length(technical)) c(technical[1], "expr_mean"))) {
      if (is.null(xy) || !all(xy %in% names(cs))) next
      attempt(paste0("composition_", xy[2], "_", xy[1]),
              save_qc_plot(path("composition_", xy[2], "_vs_", xy[1], ".png"), scatter(xy[1], xy[2])))
    }
    if ("nl_score" %in% names(cs)) {
      attempt("composition_nl_by_cell", save_qc_plot(path("composition_nl_score_by_cell.png"),
        ggplot2::ggplot(cs, ggplot2::aes(cell, nl_score)) + ggplot2::geom_boxplot(outlier.shape = NA) +
          ggplot2::geom_jitter(width = 0.2, alpha = 0.4, size = 0.8) +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)) +
          ggplot2::labs(title = "Neutrophil - lymphocyte score by cell"), width = 8))
    }
    if (any(!is.na(cs$chrY_score))) {
      sc <- composition$sex_check
      # Without XIST (e.g. EARLI's expression has none), chrY score by recorded sex.
      x <- if (any(!is.na(cs$xist))) "xist" else "recorded_sex"
      attempt("sex_check", save_qc_plot(path("sex_check.png"),
        ggplot2::ggplot(cs, ggplot2::aes(.data[[x]], chrY_score, colour = recorded_sex, shape = inferred_sex)) +
          ggplot2::geom_jitter(alpha = 0.8, width = if (x == "xist") 0 else 0.2, height = 0) +
          ggplot2::labs(title = sprintf("Sex check: agreement %s, separation %.1f SD, discordant subjects %s",
                                        format(round(sc$agreement, 2)), sc$separation,
                                        format(sc$n_subjects_discordant)))))
    }
    tv <- composition$trait_varpart
    if (!"error" %in% names(tv)) {
      tvm <- tv[, setdiff(names(tv), "trait"), drop = FALSE]
      rownames(tvm) <- tv$trait
      attempt("composition_trait_varpart", save_qc_plot(path("composition_trait_varpart.png"),
        variancePartition::plotPercentBars(variancePartition::sortCols(tvm)) +
          ggplot2::ggtitle("Sample-level traits: variance fractions"), width = 8, height = 5))
    }
  }

  # Tables.
  tables <- list(
    terms = summary$terms,
    gene_flags = merge(vp$genes, vp$flags, by.x = "id", by.y = "gene"),
    top_genes = do.call(rbind, lapply(terms, function(term) {
      top <- utils::head(order(vpm[[term]], decreasing = TRUE), n_top)
      data.frame(term = term, gene = rownames(vpm)[top], symbol = unname(symbol[rownames(vpm)[top]]),
                 fraction = vpm[top, term])
    })),
    dropped = design_obj$dropped
  )
  if (!is.null(composition)) {
    tables$composition_samples <- composition$samples
    tables$composition_correlations <- composition$correlations
    tables$composition_trait_varpart <- composition$trait_varpart
    tables$composition_markers <- composition$markers
    tables$composition_pc_variance <- composition$pc_variance
    tables$sex_check <- composition$sex_check
  }
  if (length(vp$errors)) tables$vp_errors <- data.frame(gene = names(vp$errors), error = unname(vp$errors))
  if (length(de$errors)) tables$dream_errors <- data.frame(gene = names(de$errors), error = unname(de$errors))
  # Every contrast and joint test (F.std is comparable across genes, per the
  # dream vignette), most significant first.
  for (c in names(de$tables)) {
    tt <- de$tables[[c]]
    tt <- tt[order(tt$P.Value), , drop = FALSE]
    tables[[paste0("dream_", make.names(c))]] <- data.frame(gene = rownames(tt), symbol = unname(symbol[rownames(tt)]), tt)
  }
  for (nm in names(tables)) {
    f <- path(nm, ".csv")
    utils::write.csv(tables[[nm]], f, row.names = FALSE)
    files <- c(files, f)
  }
  if (length(errors)) {
    f <- path("qc_errors.csv")
    utils::write.csv(data.frame(item = names(errors), error = unname(errors)), f, row.names = FALSE)
    files <- c(files, f)
  }
  files
}

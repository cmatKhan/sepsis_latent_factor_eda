# Leukocyte composition versus technical scale, per dataset (the
# `composition_<dataset>` target). Sepsis blood shifts strongly between
# neutrophils and lymphocytes, which moves many genes at once and can change a
# sample's whole expression distribution -- the same signature a technical
# scale artifact leaves. This separates the two: per-sample neutrophil and
# lymphocyte marker scores (config/composition_markers.yml) against the
# sample's expression spread, location and leading principal components, and
# against a technical reference (the array's control probes,
# control_probe_matrix(); RNA-seq library size and TMM factor). The
# sample-level traits also go through the dataset's variance-partition
# formula, splitting composition and spread into subject, time and group
# components. A sex check (chrY genes vs recorded sex) flags metadata that
# doesn't belong to its samples. Kept whole for ingestion and the app; plots
# are drawn by write_variance_qc().

#' Read the composition marker lists
#'
#' @param path Path to `config/composition_markers.yml`.
#' @return Named list of character vectors (one per cell type).
read_composition_markers <- function(path) yaml::read_yaml(path)

#' Marker score per sample
#'
#' @param E Gene x sample log-expression matrix with symbol row names.
#' @param genes Marker symbols.
#' @return List: `score` (mean over present markers of each marker's
#'   across-sample z-score), `found` (markers present).
marker_score <- function(E, genes) {
  found <- intersect(genes, rownames(E))
  if (length(found) == 0) return(list(score = rep(NA_real_, ncol(E)), found = found))
  z <- t(scale(t(E[found, , drop = FALSE])))
  list(score = colMeans(z, na.rm = TRUE), found = found)
}

#' Sex inferred from expression, against the recorded sex
#'
#' chrY genes against XIST, split into two groups (k-means on the chrY
#' score; the higher group is male). Catches sample metadata that doesn't
#' belong to the arrays/libraries it is attached to: recorded sex that
#' disagrees with the inferred one, or one subject's samples disagreeing.
#'
#' @param E_sym Gene x sample log-expression matrix, symbol row names.
#' @param data The design's model data frame.
#' @param spec The design block.
#' @return List: `samples` (chrY score, XIST, inferred and recorded sex),
#'   `summary` (one row: separation of the two groups in within-group SDs,
#'   agreement with recorded sex, subjects whose samples disagree).
sex_check <- function(E_sym, data, spec) {
  y_genes <- intersect(c("RPS4Y1", "DDX3Y", "KDM5D", "EIF1AY", "UTY", "ZFY", "USP9Y"), rownames(E_sym))
  n <- ncol(E_sym)
  empty <- list(samples = data.frame(chrY_score = rep(NA_real_, n), xist = NA_real_, inferred_sex = NA_character_,
                                     recorded_sex = NA_character_),
                summary = data.frame(n_chrY_genes = length(y_genes), separation = NA_real_, n_recorded = NA_integer_,
                                     agreement = NA_real_, n_subjects_discordant = NA_integer_))
  if (length(y_genes) == 0) return(empty)
  chrY <- colMeans(E_sym[y_genes, , drop = FALSE])
  xist <- if ("XIST" %in% rownames(E_sym)) E_sym["XIST", ] else rep(NA_real_, n)
  set.seed(1)
  km <- stats::kmeans(chrY, centers = range(chrY))
  male <- which.max(km$centers[, 1])
  inferred <- ifelse(km$cluster == male, "M", "F")
  within_sd <- sqrt(sum(km$withinss) / max(n - 2, 1))
  sex_col <- intersect(design_variables(spec), c("sex", "gender"))[1]
  recorded <- if (is.na(sex_col)) rep(NA_character_, n) else {
    r <- toupper(substr(as.character(data[[sex_col]]), 1, 1))
    ifelse(r %in% c("F", "M"), r, NA_character_)
  }
  n_discordant <- if (is.null(spec[["subject"]])) NA_integer_ else
    sum(tapply(inferred, data[[spec[["subject"]]]], function(x) length(unique(x)) > 1))
  list(
    samples = data.frame(chrY_score = unname(chrY), xist = unname(xist), inferred_sex = inferred, recorded_sex = recorded),
    summary = data.frame(n_chrY_genes = length(y_genes), separation = abs(diff(km$centers[, 1])) / within_sd,
                         n_recorded = sum(!is.na(recorded)),
                         agreement = if (all(is.na(recorded))) NA_real_ else mean(inferred == recorded, na.rm = TRUE),
                         n_subjects_discordant = n_discordant)
  )
}

#' Composition versus technical scale for one dataset
#'
#' @param meta The resolved dataset config.
#' @param design_obj prepare_design() output.
#' @param expr prepare_expression() output.
#' @param norm normalize_for_models() output.
#' @param markers read_composition_markers() output.
#' @param controls control_probe_matrix() output, or `NULL` (no control probes).
#' @param n_pcs Principal components kept.
#' @return List: `samples` (per-sample design variables, distribution stats,
#'   PCs, marker scores, technical reference), `markers` (cell type, gene,
#'   present), `pc_variance`, `correlations` (Spearman, every trait pair of
#'   interest), `trait_varpart` (variance fractions of the sample-level traits
#'   under the dataset's variance formula), `control` (probe count and label),
#'   `sex_check` (sex_check()'s summary; per-sample values are in `samples`).
analyze_composition <- function(meta, design_obj, expr, norm, markers, controls = NULL, n_pcs = 5) {
  data <- design_obj$data
  samples <- rownames(data)
  E <- model_expression(norm)[, samples, drop = FALSE]
  sym <- stats::setNames(expr$genes$symbol, expr$genes$id)[rownames(E)]
  # RNA-seq rows are feature ids: score on symbols, first row per symbol.
  E_sym <- E[!is.na(sym) & !duplicated(sym), , drop = FALSE]
  rownames(E_sym) <- sym[!is.na(sym) & !duplicated(sym)]

  scores <- lapply(markers, function(g) marker_score(E_sym, g))
  pcs <- stats::prcomp(t(E[matrixStats::rowVars(E) > 0, , drop = FALSE]), center = TRUE, rank. = n_pcs)
  pc_var <- pcs$sdev^2 / sum(pcs$sdev^2)

  out <- data.frame(sample = samples, data[, setdiff(names(data), "cell"), drop = FALSE], cell = data$cell,
                    row.names = NULL, check.names = FALSE)
  out$expr_mean <- colMeans(E)
  out$expr_sd <- matrixStats::colSds(E)
  q <- matrixStats::colQuantiles(E, probs = c(0.1, 0.5, 0.9))
  out$expr_q10 <- q[, 1]; out$expr_median <- q[, 2]; out$expr_q90 <- q[, 3]
  for (k in seq_len(ncol(pcs$x))) out[[paste0("PC", k)]] <- pcs$x[samples, k]
  for (ct in names(scores)) out[[paste0(ct, "_score")]] <- scores[[ct]]$score
  if (all(c("neutrophil", "lymphocyte") %in% names(scores))) {
    out$nl_score <- out$neutrophil_score - out$lymphocyte_score
  }
  sex <- sex_check(E_sym, data, design_obj$spec)
  out <- cbind(out, sex$samples)
  control <- list(n_probes = 0L, label = NA_character_)
  if (expr$platform == "array") {
    has <- !is.null(controls$values)
    out$control_mean <- if (has) unname(colMeans(controls$values[, samples, drop = FALSE], na.rm = TRUE)) else NA_real_
    control <- list(n_probes = controls$n_probes %||% 0L, label = controls$label %||% NA_character_,
                    covariate = isTRUE(meta$design$control_probes$covariate))
  } else {
    out$log2_lib_size <- log2(expr$values$samples[samples, "lib.size"])
    out$norm_factor <- expr$values$samples[samples, "norm.factors"]
  }

  traits <- intersect(c("expr_mean", "expr_sd", "expr_q90", "PC1", "PC2", "PC3", "neutrophil_score",
                        "lymphocyte_score", "nl_score", "control_mean", "log2_lib_size", "norm_factor"),
                      names(out))
  traits <- traits[vapply(out[traits], function(x) sum(!is.na(x)) > 2 && stats::sd(x, na.rm = TRUE) > 0, logical(1))]
  pairs <- utils::combn(traits, 2, simplify = FALSE)
  correlations <- do.call(rbind, lapply(pairs, function(p) {
    ct <- suppressWarnings(stats::cor.test(out[[p[1]]], out[[p[2]]], method = "spearman", exact = FALSE))
    data.frame(x = p[1], y = p[2], rho = unname(ct$estimate), p_value = ct$p.value)
  }))

  trait_mat <- t(as.matrix(out[, traits]))
  colnames(trait_mat) <- samples
  trait_varpart <- tryCatch({
    vp <- variancePartition::fitExtractVarPartModel(trait_mat, design_obj$vp_formula, data, REML = FALSE,
                                                    showWarnings = FALSE)
    data.frame(trait = rownames(vp), as.data.frame(vp), row.names = NULL, check.names = FALSE)
  }, error = function(e) data.frame(trait = traits, error = conditionMessage(e)))

  list(
    samples = out,
    markers = do.call(rbind, lapply(names(markers), function(ct) {
      data.frame(cell_type = ct, gene = markers[[ct]], present = markers[[ct]] %in% scores[[ct]]$found)
    })),
    pc_variance = data.frame(pc = paste0("PC", seq_along(pc_var))[seq_len(n_pcs)], fraction = pc_var[seq_len(n_pcs)]),
    correlations = correlations,
    trait_varpart = trait_varpart,
    control = control,
    sex_check = sex$summary
  )
}

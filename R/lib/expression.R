# Expression input for variancePartition and dream (R/targets/variance.R).
# Built from the raw values (pivot_expression_long()), not the factor-analysis
# matrix (`mat_<dataset>`, filtered to the top-variance genes and
# transformed): each platform gets the processing the variancePartition
# vignette recommends. Arrays: the log-scale values, one probe per gene
# symbol, no variance filter (a technical factor from the array's control
# probes can enter the models as a covariate, compute_control_factor()).
# RNA-seq: counts, filterByExpr() and TMM; voom weights come later
# (normalize_for_models()).

#' Collapse a probe x sample matrix to one row per gene symbol
#'
#' Keeps, per symbol, the probe with the highest mean expression (the
#' convention the preprocessing scripts use).
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

#' A dataset's control probes
#'
#' From the design's `control_probes`: feature ids matching `id_pattern`, or
#' features whose feature-metadata `column` is in `values`.
#'
#' @param cp The design's `control_probes` block.
#' @param ids Feature ids present in the expression data.
#' @param fm Feature-metadata data frame.
#' @param fid Its feature-id column.
#' @return Control feature ids present in `ids`.
control_probe_ids <- function(cp, ids, fm, fid) {
  if (!is.null(cp$id_pattern)) return(grep(cp$id_pattern, ids, value = TRUE))
  intersect(fm[[fid]][fm[[cp$column]] %in% unlist(cp$values)], ids)
}

#' The array's control probes, raw (the `controls_<dataset>` target)
#'
#' The probes named by the design's `control_probes`, for every sample in the
#' expression file. Read once and shared by compute_control_factor() and
#' analyze_composition().
#'
#' @param meta The resolved dataset config.
#' @return List: `values` (control probe x sample matrix, `NULL` when none
#'   match), `n_probes`, `label` (how the probes were chosen).
control_probe_matrix <- function(meta) {
  ds <- meta$dataset
  cp <- meta$design$control_probes
  raw <- pivot_expression_long(ds)
  fm <- read_parquet_plain(ds$feature_metadata_path)
  ids <- control_probe_ids(cp, rownames(raw), fm, ds$feature_id_col %||% "feature_id")
  label <- if (!is.null(cp$id_pattern)) cp$id_pattern else
    paste0(cp$column, " in [", paste(unlist(cp$values), collapse = ", "), "]")
  list(values = if (length(ids)) raw[ids, , drop = FALSE], n_probes = length(ids), label = label)
}

#' Technical factor from the array's control probes (the `control_factor_<dataset>` target)
#'
#' RUVg-style (Risso et al. 2014): the first principal component of the
#' control probes, oriented to follow their mean, used as a per-sample
#' covariate (design `control_probes: {covariate: true}`) so each gene gets its
#' own coefficient on it. Shifting every probe by a control offset instead
#' over-corrected ANEMONES (its controls move ~2x as much as genes).
#'
#' @param controls control_probe_matrix() output.
#' @return List: `factor` (named per-sample scores), `fraction` (control-probe
#'   variance it explains), `n_probes`.
compute_control_factor <- function(controls) {
  if (controls$n_probes < 10) {
    stop("control_probes.covariate needs at least 10 control probes, found ", controls$n_probes)
  }
  ctl <- controls$values
  ctl <- ctl[stats::complete.cases(ctl) & matrixStats::rowVars(ctl) > 0, , drop = FALSE]
  pc <- stats::prcomp(t(ctl), center = TRUE)
  score <- pc$x[, 1]
  if (stats::cor(score, colMeans(ctl)) < 0) score <- -score
  list(factor = score, fraction = pc$sdev[1]^2 / sum(pc$sdev^2), n_probes = nrow(ctl))
}

#' Expression prepared for the variance models (the `expr_<dataset>` target)
#'
#' Samples are those of the design (in its order). Arrays: probes collapsed
#' to symbols, genes with any missing value dropped. RNA-seq:
#' integer counts in an edgeR DGEList, filterByExpr() against the cell-means
#' design, TMM factors.
#'
#' @param meta The resolved dataset config.
#' @param design_obj prepare_design() output.
#' @return List: `platform`, `values` (gene x sample matrix for arrays, DGEList
#'   for RNA-seq), `genes` (data frame `id`, `symbol`), `n_features_raw`,
#'   `n_genes`, `notes` (what was dropped).
prepare_expression <- function(meta, design_obj) {
  ds <- meta$dataset
  raw <- pivot_expression_long(ds)
  samples <- rownames(design_obj$data)
  absent <- setdiff(samples, colnames(raw))
  if (length(absent) > 0) {
    stop(ds$id, ": ", length(absent), " design samples have no expression column (e.g. ",
         paste(utils::head(absent, 3), collapse = ", "), ")")
  }
  notes <- c(samples_without_metadata = sum(!colnames(raw) %in% samples))
  raw <- raw[, samples, drop = FALSE]
  fm <- read_parquet_plain(ds$feature_metadata_path)
  fid <- ds$feature_id_col %||% "feature_id"

  if (design_obj$platform == "array") {
    mat <- collapse_to_symbol(raw, fm, fid, ds$symbol_col)
    incomplete <- !stats::complete.cases(mat)
    notes["genes_with_missing_values"] <- sum(incomplete)
    mat <- mat[!incomplete, , drop = FALSE]
    genes <- data.frame(id = rownames(mat), symbol = rownames(mat))
    return(list(platform = "array", values = mat, genes = genes, n_features_raw = nrow(raw),
                n_genes = nrow(mat), notes = notes))
  }

  if (anyNA(raw) || any(raw != round(raw))) {
    stop(ds$id, ": design.platform is rnaseq but the expression values aren't whole-number counts")
  }
  storage.mode(raw) <- "integer"
  dge <- edgeR::DGEList(counts = raw)
  keep <- edgeR::filterByExpr(dge, stats::model.matrix(~ 0 + cell, design_obj$data))
  notes["genes_failing_filterByExpr"] <- sum(!keep)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::normLibSizes(dge, method = "TMM")
  symbols <- build_symbol_map(list(dataset = ds), fm = fm)
  genes <- data.frame(id = rownames(dge), symbol = unname(symbols[rownames(dge)]))
  list(platform = "rnaseq", values = dge, genes = genes, n_features_raw = nrow(raw),
       n_genes = nrow(dge), notes = notes)
}

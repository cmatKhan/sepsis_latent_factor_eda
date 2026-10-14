# From the variance module's target values to DB rows and dataset artifacts:
# write_dataset_db()'s `variance` argument (variance_db_rows()), the
# `variance_files_<dataset>` and `factor_gene_files_<dataset>` file targets,
# and the cross-dataset time-gene recurrence.

#' JSON text of a character vector (or `NA`)
#'
#' @param x Character vector or `NULL`.
#' @return A JSON array string, or `NA` when `x` is empty.
json_text <- function(x) {
  x <- unlist(x)
  if (!length(x)) NA_character_ else as.character(jsonlite::toJSON(as.character(x)))
}

#' A column of a data frame, or `NA` when absent
#'
#' @param df A data frame.
#' @param col Column name.
#' @param as Coercion function.
#' @return The column, or `NA`s of the right length.
column_or_na <- function(df, col, as = as.numeric) {
  if (!is.null(col) && col %in% names(df)) as(df[[col]]) else as(rep(NA, nrow(df)))
}

#' The variance inputs of one dataset's DB write
#'
#' @param dataset_id Dataset id.
#' @param parent_id The dataset whose design applies (itself, or the parent
#'   of a per-timepoint subset).
#' @param design prepare_design() output (the parent's for a subset).
#' @param composition analyze_composition() output (the parent's for a subset).
#' @param samples The dataset's samples (`colnames(mat_<dataset>)`).
#' @param summary summarize_variance() output, or `NULL` (subsets).
#' @param explore explore_variance() output, or `NULL` (subsets).
#' @return The list write_dataset_db() takes as `variance`.
variance_inputs <- function(dataset_id, parent_id, design, composition, samples, summary = NULL,
                            explore = NULL) {
  list(dataset_id = dataset_id, parent_id = parent_id, design = design, composition = composition,
       samples = samples, summary = summary, explore = explore)
}

#' Rows for the repeated-measures tables
#'
#' @param v variance_inputs() output.
#' @param now Timestamp.
#' @return Named list of data frames, one per table (`dataset_design`,
#'   `samples`, and for datasets with their own design `variance_summary`,
#'   `variance_terms`, `variance_categories`, `variance_bins`,
#'   `trait_variance`, `composition_correlations`, `de_contrasts`,
#'   `variance_time_genes`).
variance_db_rows <- function(v, now) {
  id <- v$dataset_id
  d <- v$design; spec <- d$spec; comp <- v$composition
  is_subset <- !identical(v$parent_id, id)
  in_design <- intersect(v$samples, rownames(d$data))
  data <- d$data[in_design, , drop = FALSE]
  out <- list()
  out$dataset_design <- data.frame(
    dataset_id = id, parent_dataset_id = v$parent_id, is_subset = as.integer(is_subset),
    platform = d$platform, subject_col = spec[["subject"]] %||% NA_character_,
    time_col = spec[["time"]] %||% NA_character_, time_order = json_text(spec$time_order),
    group_col = spec[["group"]] %||% NA_character_, covariates = json_text(spec$covariates),
    technical_covariates = json_text(spec$vp_covariates),
    vp_formula = formula_text(d$vp_formula), dream_formula = formula_text(d$dream_formula),
    n_samples = length(in_design),
    n_subjects = if (is.null(spec[["subject"]])) NA_integer_ else length(unique(as.character(data[[spec[["subject"]]]]))),
    n_dropped = if (is_subset) NA_integer_ else as.integer(sum(d$dropped$n)))

  cs <- comp$samples[comp$samples$sample %in% in_design, , drop = FALSE]
  chr <- function(col) column_or_na(cs, col, as.character)
  num <- function(col) column_or_na(cs, col)
  out$samples <- data.frame(
    dataset_id = id, sample_id = cs$sample, subject = chr(spec[["subject"]]), timepoint = chr(spec[["time"]]),
    grp = chr(spec[["group"]]), cell = chr("cell"), expr_mean = num("expr_mean"), expr_sd = num("expr_sd"),
    expr_q10 = num("expr_q10"), expr_median = num("expr_median"), expr_q90 = num("expr_q90"),
    pc1 = num("PC1"), pc2 = num("PC2"), pc3 = num("PC3"), neutrophil_score = num("neutrophil_score"),
    lymphocyte_score = num("lymphocyte_score"), nl_score = num("nl_score"),
    control_mean = num("control_mean"), control_factor = num("control_factor"),
    log2_lib_size = num("log2_lib_size"), norm_factor = num("norm_factor"),
    chrY_score = num("chrY_score"), xist = num("xist"), inferred_sex = chr("inferred_sex"),
    recorded_sex = chr("recorded_sex"))
  if (is_subset || is.null(v$summary)) return(out)

  # Columns copied as named, from the explore headline, the within-person
  # change and the summary row; a few summary columns are renamed.
  h <- v$explore$headline; ch <- v$explore$change; s <- v$summary$dataset
  summary_cols <- c(n_genes = "n_genes_fit", n_genes_failed = "n_genes_failed",
                    share_singular = "share_singular", max_colinearity = "max_colinearity",
                    dupcor_consensus = "dupcor_consensus", rho_pc1_nl = "rho_PC1_nl", rho_sd_nl = "rho_sd_nl",
                    rho_mean_nl = "rho_mean_nl", rho_sd_technical = "rho_sd_technical",
                    rho_mean_technical = "rho_mean_technical", n_control_probes = "n_control_probes",
                    sex_separation = "sex_separation", sex_agreement = "sex_agreement",
                    n_subjects_sex_discordant = "n_subjects_sex_discordant",
                    nl_subject_fraction = "nl_subject_fraction", nl_time_fraction = "nl_time_fraction",
                    nl_group_fraction = "nl_group_fraction", vp_seconds = "vp_seconds",
                    dream_seconds = "dream_seconds")
  headline_cols <- c("n_factor_genes", "between_median", "time_median", "residual_median", "technical_median",
                     "repeatability_median", "time_within_median", "share_time_gt05", "share_time_gt10",
                     "share_time_gt20", "between_vw_all", "time_vw_all", "residual_vw_all",
                     "between_vw_factor_genes", "time_vw_factor_genes", "residual_vw_factor_genes",
                     "between_vw_top", "time_vw_top", "residual_vw_top")
  change_cols <- c("n_subjects_repeated", "n_changes", "change_pc1_share", "change_pc1_cor_nl",
                   "median_r2_nl_all", "median_r2_nl_factor_genes", "median_r2_technical_all",
                   "median_r2_technical_factor_genes")
  first_or_na <- function(df, cols) lapply(stats::setNames(cols, names(cols) %||% cols),
                                           function(c) if (!is.null(df) && c %in% names(df)) df[[c]][1] else NA)
  out$variance_summary <- data.frame(
    dataset_id = id, first_or_na(s, summary_cols), first_or_na(h, headline_cols), first_or_na(ch, change_cols),
    control_probes = comp$control$label %||% NA_character_,
    control_factor_fraction = d$control_factor_fraction %||% NA_real_, fitted_at = now)

  terms <- v$summary$terms
  cats <- variance_term_category(d, terms$term)
  out$variance_terms <- data.frame(
    dataset_id = id, term = terms$term,
    category = ifelse(terms$kind == "derived", "derived", unname(cats[terms$term])),
    n_genes = terms$n_genes, mean = terms$mean, q10 = terms$q10, median = terms$median, q90 = terms$q90,
    share_gt05 = terms$share_gt05, share_gt10 = terms$share_gt10)
  out$variance_categories <- transform(v$explore$categories, present = as.integer(present))
  bins <- function(df, binning, col) data.frame(dataset_id = id, binning = binning, bin = df[[col]],
                                                n_genes = df$n_genes, between_subject = df$between, time = df$time,
                                                residual = df$residual)
  out$variance_bins <- rbind(bins(v$explore$by_expression, "expression_quartile", "expression_quartile"),
                             bins(v$explore$by_variance_rank, "variance_rank", "variance_rank"))
  tv <- comp$trait_varpart
  if (!"error" %in% names(tv)) {
    tcat <- variance_term_category(d, setdiff(names(tv), "trait"))
    out$trait_variance <- do.call(rbind, lapply(setdiff(names(tv), "trait"), function(t) data.frame(
      dataset_id = id, trait = tv$trait, term = t, category = unname(tcat[t]), fraction = tv[[t]])))
  }
  out$composition_correlations <- cbind(data.frame(dataset_id = id), comp$correlations)
  ct <- v$summary$contrasts
  if (nrow(ct)) out$de_contrasts <- data.frame(dataset_id = id, contrast = ct$contrast, kind = ct$kind,
                                               readable = ct$readable, n_genes = ct$n_genes,
                                               n_fdr05 = as.integer(ct$n_fdr05))
  g <- v$explore$genes
  g <- g[g$time_all > TIME_GENE_THRESHOLD & !is.na(g$symbol), , drop = FALSE]
  g <- g[order(-g$time_all), , drop = FALSE]
  g <- g[!duplicated(g$symbol), , drop = FALSE]
  if (nrow(g)) out$variance_time_genes <- data.frame(dataset_id = id, symbol = g$symbol,
                                                     time_all = g$time_all, between_subject = g$between)
  out
}

#' Rows for factor_variance and factor_variance_models
#'
#' @param methods write_dataset_db()'s per-method parts (`factor_variance`).
#' @param dataset_id Dataset id.
#' @return `list(models, rows)`: per-method formula and factor count, and the
#'   per-factor term fractions keyed by fit_key and factor_index.
factor_variance_rows <- function(methods, dataset_id) {
  parts <- Filter(function(p) !is.null(p$factor_variance) && nrow(p$factor_variance), methods)
  if (!length(parts)) return(NULL)
  models <- do.call(rbind, lapply(names(parts), function(m) {
    fv <- parts[[m]]$factor_variance
    data.frame(dataset_id = dataset_id, method = m, formula = attr(fv, "formula") %||% NA_character_,
               n_factors = nrow(unique(fv[, c("fit_key", "factor_index")])))
  }))
  rows <- do.call(rbind, lapply(parts, function(p) as.data.frame(p$factor_variance)))
  list(models = models, rows = rows)
}

#' Write a dataset's variance artifacts (the `variance_files_<dataset>` file target)
#'
#' `variance_partition.rds`: per-gene fractions, variances, fit flags, the
#' gene table and explore_variance()'s per-gene table; `de_results.rds`: every
#' dream contrast and joint test as a full topTable (dream_tables());
#' `composition.rds`: analyze_composition() whole. File names are the
#' `dataset_artifacts` kinds.
#'
#' @param dataset_id Dataset id.
#' @param db_path Path to the DB.
#' @param vp run_variance_partition() output.
#' @param de dream_tables() output.
#' @param composition analyze_composition() output.
#' @param explore explore_variance() output.
#' @return The written paths.
write_variance_artifacts <- function(dataset_id, db_path, vp, de, composition, explore) {
  dir <- artifacts_dir(db_path, dataset_id)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  files <- c(
    variance_partition = file.path(dir, "variance_partition.rds"),
    de_results = file.path(dir, "de_results.rds"),
    composition = file.path(dir, "composition.rds"))
  saveRDS(list(fractions = vp$varPart, variances = vp$variances, flags = vp$flags, genes = vp$genes,
               formula = formula_text(vp$formula), gene_summary = explore$genes),
          files[["variance_partition"]])
  saveRDS(de[c("tables", "contrasts")], files[["de_results"]])
  saveRDS(composition, files[["composition"]])
  unname(files)
}

#' Write a dataset's per-factor gene lists (the `factor_gene_files_<dataset>` file target)
#'
#' `enrichment_genes.rds`: each enrichment result's leading-edge or overlap
#' genes (fit_key, factor_index, query_type, direction, term_id, genes);
#' `pattern_markers.rds`: each representative fit's marker genes and scores
#' (fit_key, factor_index, gene, score). Kept out of the DB tables for size.
#'
#' @param dataset_id Dataset id.
#' @param db_path Path to the DB.
#' @param enrichment enrichment_rows() output.
#' @param redundancy Named list of compute_method_redundancy() outputs.
#' @return The written paths.
write_factor_gene_artifacts <- function(dataset_id, db_path, enrichment, redundancy) {
  dir <- artifacts_dir(db_path, dataset_id)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  er <- enrichment$results
  genes <- if (is.null(er) || !nrow(er)) NULL else
    er[, intersect(c("fit_key", "factor_index", "query_type", "direction", "term_id", "genes"), names(er))]
  markers <- do.call(rbind, lapply(Filter(Negate(is.null), redundancy), `[[`, "pattern_markers"))
  files <- c(file.path(dir, "enrichment_genes.rds"), file.path(dir, "pattern_markers.rds"))
  saveRDS(genes, files[1])
  saveRDS(markers, files[2])
  files
}

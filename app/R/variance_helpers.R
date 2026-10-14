# Queries behind the first view (the repeated-measures tables of schema 2,
# docs/database.qmd) -- every one tolerates a missing table or artifact.

VARIANCE_COLORS <- c(subject = "#2C5F8A", `subject-level` = "#8DB9DC", time = "#E07B00",
                     `group x time` = "#F2B36B", technical = "#7B5EA7", group = "#4E9F6B",
                     covariate = "#A9D6B5", residual = "#BDBDBD")
VARIANCE_LABELS <- c(subject = "patient", `subject-level` = "patient-level (group, sex, age)",
                     time = "time", `group x time` = "group x time", technical = "technical",
                     group = "group", covariate = "covariate", residual = "residual (within-person + noise)")
# The summary column each choice reads from variance_categories.
VARIANCE_GENESETS <- c("Factor-analysis genes, variance-weighted" = "vw_factor_genes",
                       "All modelled genes, variance-weighted" = "vw_all",
                       "Top 8,000 by variance, variance-weighted" = "vw_top",
                       "Typical gene (median per gene)" = "median")

variance_ready <- function(con) {
  need <- c("dataset_design", "variance_summary", "variance_categories", "variance_terms")
  miss <- setdiff(need, DBI::dbListTables(con))
  if (length(miss) > 0) {
    note_missing("repeated-measures tables", paste("absent:", paste(miss, collapse = ", "),
                 "-- the database predates schema 2, or the variance module has not run"))
    return(FALSE)
  }
  TRUE
}

design_table <- function(con) safe_query(con, "SELECT * FROM dataset_design", what = "dataset_design")

#' The design that applies to a dataset: a subset uses its parent's. NULL (logged) when the
#' dataset has none (e.g. no sample metadata).
design_for <- function(con, dataset_id) {
  d <- design_table(con)
  if (is.null(d)) return(NULL)
  row <- d[d$dataset_id == dataset_id, , drop = FALSE]
  if (nrow(row) == 0) { note_missing("a repeated-measures design", paste("dataset", dataset_id)); return(NULL) }
  list(dataset = dataset_id, parent = row$parent_dataset_id, is_subset = row$is_subset == 1, row = row)
}

time_course_parents <- function(design) {
  design$dataset_id[design$is_subset == 0 & !is.na(design$subject_col) & !is.na(design$time_col)]
}
cross_sectional_parents <- function(design) {
  design$dataset_id[design$is_subset == 0 & (is.na(design$subject_col) | is.na(design$time_col))]
}

variance_summary_table <- function(con) safe_query(con, "SELECT * FROM variance_summary", what = "variance_summary")

#' Share of variance per category and dataset for the chosen gene set / weighting
#' (columns of `variance_categories`); only categories in the dataset's model.
variance_long <- function(con, datasets, column) {
  column <- match.arg(column, unname(VARIANCE_GENESETS))
  d <- safe_query(con, sprintf("SELECT dataset_id, category, present, %s AS share FROM variance_categories", column),
                  what = "variance_categories")
  if (is.null(d)) return(NULL)
  d <- d[d$present == 1 & d$dataset_id %in% datasets, , drop = FALSE]
  d$category <- factor(d$category, levels = names(VARIANCE_COLORS))
  d
}

#' One row per dataset: between patients, time, technical, residual (summed categories).
headline_table <- function(long, summary) {
  part <- function(cats) tapply(ifelse(long$category %in% cats, long$share, 0), long$dataset_id, sum)
  out <- data.frame(dataset = names(part("subject")), between = as.numeric(part(c("subject", "subject-level"))),
                    time = as.numeric(part(c("time", "group x time"))), technical = as.numeric(part("technical")),
                    residual = as.numeric(part("residual")), stringsAsFactors = FALSE)
  s <- summary[match(out$dataset, summary$dataset_id), ]
  out$repeatability <- s$repeatability_median
  out$time_within <- s$time_within_median
  out$genes_time_gt10 <- s$share_time_gt10
  out$subjects_repeated <- s$n_subjects_repeated
  out$n_genes <- s$n_genes
  out$flags <- dataset_flags(s)
  out
}

dataset_flags <- function(s) {
  vapply(seq_len(nrow(s)), function(i) {
    f <- character(0)
    if (!is.na(s$sex_agreement[i]) && s$sex_agreement[i] < 0.9) f <- c(f, sprintf("sex agreement %.2f", s$sex_agreement[i]))
    if (!is.na(s$share_singular[i]) && s$share_singular[i] > 0.9) f <- c(f, sprintf("%.0f%% singular fits", 100 * s$share_singular[i]))
    if (!is.na(s$n_genes_failed[i]) && s$n_genes_failed[i] > 0) f <- c(f, sprintf("%d genes failed", s$n_genes_failed[i]))
    paste(f, collapse = "; ")
  }, character(1))
}

variance_terms_for <- function(con, dataset_id)
  safe_query(con, "SELECT * FROM variance_terms WHERE dataset_id = ?", list(dataset_id), what = "variance_terms")
variance_bins_for <- function(con, dataset_id)
  safe_query(con, "SELECT * FROM variance_bins WHERE dataset_id = ?", list(dataset_id), what = "variance_bins")
time_genes_for <- function(con, dataset_id, n = 25)
  safe_query(con, "SELECT symbol, time_all, between_subject FROM variance_time_genes WHERE dataset_id = ?
                   ORDER BY time_all DESC LIMIT ?", list(dataset_id, n), what = "variance_time_genes")
time_gene_recurrence <- function(con, min_datasets = 2)
  safe_query(con, "SELECT symbol, n_datasets, datasets FROM time_gene_recurrence WHERE n_datasets >= ?
                   ORDER BY n_datasets DESC, symbol", list(min_datasets), what = "time_gene_recurrence")
samples_for <- function(con, dataset_id)
  safe_query(con, "SELECT * FROM samples WHERE dataset_id = ?", list(dataset_id), what = "samples")
trait_variance_for <- function(con, dataset_id)
  safe_query(con, "SELECT * FROM trait_variance WHERE dataset_id = ?", list(dataset_id), what = "trait_variance")
composition_correlations_for <- function(con, dataset_id)
  safe_query(con, "SELECT * FROM composition_correlations WHERE dataset_id = ?", list(dataset_id), what = "composition_correlations")
de_contrasts_for <- function(con, dataset_id)
  safe_query(con, "SELECT contrast, kind, readable, n_genes, n_fdr05 FROM de_contrasts WHERE dataset_id = ?",
             list(dataset_id), what = "de_contrasts")

#' Per-gene variance shares of one dataset from the `variance_partition` artifact
#' (symbol, between, time_all, residual, technical, repeatability, total_variance,
#' mean_expression, factor_gene, top_variance), or NULL when the file is absent.
gene_variance_for <- function(con, dataset_id) {
  p <- safe_query(con, "SELECT path FROM dataset_artifacts WHERE dataset_id = ? AND kind = 'variance_partition'",
                  list(dataset_id), what = "variance_partition artifact")
  if (is.null(p) || nrow(p) == 0) { note_missing("variance_partition artifact", paste("dataset", dataset_id)); return(NULL) }
  path <- resolve_artifact(p$path[1])
  if (!file.exists(path)) { note_missing("variance_partition artifact file", path); return(NULL) }
  tryCatch(readRDS(path)$gene_summary, error = function(e) { note_missing("variance_partition artifact", conditionMessage(e)); NULL })
}

#' PC1's Spearman correlation with each leukocyte marker score, one row per dataset
#' (stored once per pair, in either orientation).
pc1_composition_table <- function(con) {
  d <- safe_query(con, "SELECT dataset_id, x, y, rho, p_value FROM composition_correlations
                        WHERE (x = 'PC1' AND y IN ('neutrophil_score', 'lymphocyte_score', 'nl_score'))
                           OR (y = 'PC1' AND x IN ('neutrophil_score', 'lymphocyte_score', 'nl_score'))",
                  what = "composition_correlations")
  if (is.null(d) || nrow(d) == 0) return(NULL)
  d$score <- ifelse(d$x == "PC1", d$y, d$x)
  w <- reshape(d[, c("dataset_id", "score", "rho", "p_value")], idvar = "dataset_id", timevar = "score", direction = "wide")
  col <- function(nm, v) if (paste0(v, ".", nm) %in% names(w)) w[[paste0(v, ".", nm)]] else rep(NA_real_, nrow(w))
  data.frame(dataset = w$dataset_id, neutrophil = col("neutrophil_score", "rho"), lymphocyte = col("lymphocyte_score", "rho"),
             nl = col("nl_score", "rho"), p_nl = col("nl_score", "p_value"), stringsAsFactors = FALSE)
}

#' The marker genes behind the composition scores (config/composition_markers.yml), as
#' list(neutrophil = ..., lymphocyte = ...), or NULL (logged) when the file isn't found.
composition_markers <- function() {
  for (p in c("config/composition_markers.yml", "../config/composition_markers.yml")) {
    if (file.exists(p)) return(yaml::read_yaml(p))
  }
  note_missing("composition marker list", "config/composition_markers.yml not found")
  NULL
}

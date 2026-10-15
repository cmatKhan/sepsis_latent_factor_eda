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

## ---- explaining a dataset's models ---------------------------------------------------------

DOCS_URL <- "https://cmatkhan.github.io/sepsis_latent_factor_eda/"
PUBLIC_COLLECTION <- "cmatkhan/Sepsis transcriptomics"

# One line per kind of term the variance model can contain (variance_terms.category);
# an unknown category is shown as its own name.
TERM_MEANINGS <- c(
  subject = "stable differences between people",
  `subject-level` = "differences between people explained by a person-level variable (group, sex, age)",
  time = "the average change over time that every person shares",
  `group x time` = "time trajectories that differ between groups",
  group = "differences between groups",
  technical = "a technical axis (for example array controls or library size)",
  covariate = "a further covariate in the model",
  residual = "everything else: change within a person that is not shared by time, plus measurement noise",
  derived = "a sum of the terms above")

term_meaning <- function(category) {
  out <- unname(TERM_MEANINGS[category])
  ifelse(is.na(out), category, out)
}

#' The terms of a model formula as written in `dataset_design`, split on top-level `+`
#' and marked random (`(1 | x)`) or fixed, with the underlying variable name.
formula_terms <- function(formula_text) {
  if (is.null(formula_text) || is.na(formula_text)) return(NULL)
  txt <- gsub("\\s+", " ", trimws(sub("^\\s*~", "", formula_text)))
  parts <- character(0); depth <- 0; cur <- ""
  for (ch in strsplit(txt, "")[[1]]) {
    if (ch == "(") depth <- depth + 1
    if (ch == ")") depth <- depth - 1
    if (ch == "+" && depth == 0) { parts <- c(parts, trimws(cur)); cur <- "" } else cur <- paste0(cur, ch)
  }
  parts <- c(parts, trimws(cur))
  parts <- parts[nzchar(parts) & parts != "0" & parts != "1"]
  random <- grepl("^\\(.*\\|.*\\)$", parts)
  data.frame(term = parts, random = random,
             variable = ifelse(random, trimws(sub("^\\(\\s*1\\s*\\|\\s*(.*)\\)$", "\\1", parts)), parts),
             stringsAsFactors = FALSE)
}

#' Samples and subjects per group x time cell, in the design's time order (NULL when
#' the dataset has no time or no samples). Empty cells (a group not sampled then) are kept.
cell_table <- function(con, dataset_id, time_order = NULL) {
  s <- samples_for(con, dataset_id)
  if (is.null(s) || nrow(s) == 0 || all(is.na(s$timepoint))) return(NULL)
  s$grp <- ifelse(is.na(s$grp), "all samples", s$grp)
  tp <- unique(s$timepoint)
  lv <- if (length(time_order)) c(intersect(time_order, tp), setdiff(tp, time_order)) else sort(tp)
  n <- table(factor(s$grp), factor(s$timepoint, levels = lv))
  nsub <- tapply(s$subject, list(factor(s$grp), factor(s$timepoint, levels = lv)), function(x) length(unique(x)))
  nsub[is.na(nsub)] <- 0
  out <- matrix(sprintf("%d / %d", as.integer(n), as.integer(nsub)), nrow(n), ncol(n), dimnames = dimnames(n))
  out[n == 0] <- "-"
  as.data.frame.matrix(out)
}

#' Map dataset id -> hf_collection from config/*_config.yml (read only; configs aren't named
#' by id, so each file's declared dataset id is used). Empty when config/ isn't found.
dataset_collections <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    dir <- Find(dir.exists, c("config", "../config"))
    if (is.null(dir)) { note_missing("dataset configs", "config/ not found; no documentation links"); cache <<- character(0); return(cache) }
    files <- list.files(dir, pattern = "_config\\.yml$", full.names = TRUE)
    out <- character(0)
    for (f in files) {
      d <- tryCatch(yaml::read_yaml(f)$dataset, error = function(e) NULL)
      if (!is.null(d$id) && !is.null(d$hf_collection)) out[[d$id]] <- d$hf_collection
    }
    cache <<- out
    cache
  }
})

#' The docs page of a dataset's family, or NULL when its configured collection isn't the public one.
dataset_docs_url <- function(dataset_id, parent_id = dataset_id) {
  col <- dataset_collections()
  if (!isTRUE(unname(col[dataset_id]) == PUBLIC_COLLECTION) && !isTRUE(unname(col[parent_id]) == PUBLIC_COLLECTION)) return(NULL)
  paste0(DOCS_URL, "datasets.html#", tolower(gsub("[^A-Za-z0-9]+", "-", parent_id)))
}

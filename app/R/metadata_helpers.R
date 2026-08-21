# Generic sample/feature metadata access + sample-score helpers for the app.
#
# Sample/feature metadata are stored in the DB as a POINTER (path + id
# column) via R/ingest_results.R, never copied in -- read_metadata_table()
# is the only place that resolves that pointer, so the app always reflects
# the metadata file's *current* contents, and edits to the source file are
# picked up on next read with no re-ingest needed.
#
# Field names beyond the canonical sample_id/feature_id id column are never
# hardcoded here or anywhere in the app: whatever other columns a metadata
# table happens to contain are what's offered to the user generically.

#' Look up the registered pointer for a dataset's sample/feature metadata,
#' or NULL if none was registered at ingest time (e.g. the dataset config
#' didn't set that path/id_col).
metadata_source <- function(con, dataset_id, kind = c("sample", "feature")) {
  kind <- match.arg(kind)
  d <- DBI::dbGetQuery(con,
    "SELECT path, id_col FROM dataset_metadata_sources WHERE dataset_id = ? AND kind = ?",
    params = list(dataset_id, kind))
  if (nrow(d) == 0) return(NULL)
  list(path = d$path[1], id_col = d$id_col[1])
}

#' Read a metadata table from its registered path (local file or http(s)
#' URL; parquet/csv/tsv by extension), renaming `id_col` to
#' `canonical_id` ("sample_id" or "feature_id"). Returns NULL (with a
#' notification, when running inside a reactive context) if unreachable or
#' malformed, rather than erroring the whole view.
read_metadata_table <- function(path, id_col, canonical_id) {
  tryCatch({
    ext <- tolower(tools::file_ext(path))
    is_url <- grepl("^https?://", path)
    local_path <- path
    if (is_url && ext %in% c("csv", "tsv", "txt")) {
      local_path <- tempfile(fileext = paste0(".", ext))
      utils::download.file(path, local_path, quiet = TRUE, mode = "wb")
    }
    df <- if (ext == "parquet") {
      arrow::read_parquet(local_path)
    } else if (ext == "tsv") {
      utils::read.delim(local_path, stringsAsFactors = FALSE, check.names = FALSE)
    } else {
      utils::read.csv(local_path, stringsAsFactors = FALSE, check.names = FALSE)
    }
    df <- as.data.frame(df)
    if (!(id_col %in% names(df))) {
      stop("id column '", id_col, "' not found in ", path,
           " (columns: ", paste(names(df), collapse = ", "), ")")
    }
    names(df)[names(df) == id_col] <- canonical_id
    df[[canonical_id]] <- as.character(df[[canonical_id]])
    df
  }, error = function(e) {
    if (!is.null(shiny::getDefaultReactiveDomain())) {
      showNotification(paste("Could not read metadata:", conditionMessage(e)), type = "error")
    } else {
      warning(conditionMessage(e))
    }
    NULL
  })
}

#' Convenience: resolve + read a dataset's sample or feature metadata table
#' in one call; NULL if no pointer registered or the read failed.
dataset_metadata <- function(con, dataset_id, kind = c("sample", "feature")) {
  kind <- match.arg(kind)
  src <- metadata_source(con, dataset_id, kind)
  if (is.null(src)) return(NULL)
  canonical <- if (kind == "sample") "sample_id" else "feature_id"
  read_metadata_table(src$path, src$id_col, canonical)
}

#' scores_file paths use the same relative-to-DB-dir convention as
#' loadings_file (see resolve_artifact() in app/R/db_helpers.R).
load_scores <- function(con, fit_id) {
  f <- get_fit(con, fit_id)
  if (nrow(f) == 0 || is.na(f$scores_file)) return(NULL)
  path <- resolve_artifact(f$scores_file)
  if (!file.exists(path)) return(NULL)
  readRDS(path)
}

#' Generic per-column association scan: for every column of `meta_df`
#' except `id_col`, tests each column of `scores_mat` (samples x
#' factor/module) against it -- Spearman correlation for numeric fields,
#' Kruskal-Wallis for everything else. No field names are hardcoded; skips
#' fields with too little signal (all-NA, a single level, or < 3 samples
#' after alignment) rather than erroring. `padj` is BH-adjusted across every
#' row returned by this call (i.e. across the whole factor x field grid
#' shown together in one heatmap).
generic_association_scan <- function(scores_mat, meta_df, id_col = "sample_id") {
  empty <- data.frame(component = character(0), field = character(0),
                       test = character(0), statistic = numeric(0),
                       p_value = numeric(0))
  if (is.null(scores_mat) || is.null(meta_df) || !(id_col %in% names(meta_df))) return(empty)

  shared <- intersect(rownames(scores_mat), meta_df[[id_col]])
  if (length(shared) < 3) return(empty)
  scores_mat <- scores_mat[shared, , drop = FALSE]
  meta_df <- meta_df[match(shared, meta_df[[id_col]]), , drop = FALSE]

  fields <- setdiff(names(meta_df), id_col)
  comps <- colnames(scores_mat)
  if (is.null(comps)) comps <- as.character(seq_len(ncol(scores_mat)))

  rows <- list()
  for (field in fields) {
    v <- meta_df[[field]]
    is_numeric <- is.numeric(v)
    for (j in seq_along(comps)) {
      y <- scores_mat[, j]
      ok <- !is.na(v) & !is.na(y)
      if (sum(ok) < 3) next
      res <- tryCatch({
        if (is_numeric) {
          if (length(unique(v[ok])) < 2) return(NULL)
          ct <- suppressWarnings(cor.test(y[ok], v[ok], method = "spearman", exact = FALSE))
          list(test = "spearman", statistic = unname(ct$estimate), p_value = ct$p.value)
        } else {
          f <- as.factor(v[ok])
          if (nlevels(droplevels(f)) < 2) return(NULL)
          kt <- kruskal.test(y[ok] ~ droplevels(f))
          list(test = "kruskal", statistic = unname(kt$statistic), p_value = kt$p.value)
        }
      }, error = function(e) NULL)
      if (is.null(res)) next
      rows[[length(rows) + 1]] <- data.frame(
        component = comps[j], field = field, test = res$test,
        statistic = res$statistic, p_value = res$p_value)
    }
  }
  if (length(rows) == 0) return(empty)
  out <- do.call(rbind, rows)
  out$padj <- p.adjust(out$p_value, method = "BH")
  out
}

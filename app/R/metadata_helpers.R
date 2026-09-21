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
#' didn't set that path/id_col). `ensembl_col`/`symbol_col` (kind =
#' "feature" only; NA otherwise) are NA for datasets ingested before those
#' columns existed -- re-run --stage core (or R/ingest_results.R) for that
#' dataset to populate them.
metadata_source <- function(con, dataset_id, kind = c("sample", "feature")) {
  kind <- match.arg(kind)
  d <- DBI::dbGetQuery(con,
    "SELECT path, id_col, ensembl_col, symbol_col FROM dataset_metadata_sources WHERE dataset_id = ? AND kind = ?",
    params = list(dataset_id, kind))
  if (nrow(d) == 0) return(NULL)
  list(path = d$path[1], id_col = d$id_col[1], ensembl_col = d$ensembl_col[1], symbol_col = d$symbol_col[1])
}

#' Build (and cache for the life of this Shiny session) a dataset's
#' (feature_id -> Ensembl gene id) map -- used by build_display_map()
#' below for gene-label display (DISPLAY ONLY; all enrichment computation
#' now happens in the cluster ingest pipeline, R/ingest_jobs/fgsea_job.R,
#' which builds its own ensembl maps independently -- see that file's
#' header). THE canonical cross-dataset identifier used everywhere else in
#' this pipeline (see R/lib/ingest/symbol_mapping.R's header).
#'
#' Deliberately reuses build_ensembl_map() (sourced from R/lib/ingest/
#' symbol_mapping.R -- see app.R's top) by constructing a minimal
#' config-shaped list from dataset_metadata_sources instead of reading
#' config/*.yml directly -- this app is "fully decoupled from ingest" (see
#' this file's header) and must get everything through the DB.
#'
#' Returns NULL (logged once via a message, not a per-call notification)
#' if this dataset predates the ensembl_col column, has no feature
#' metadata registered at all, or the map fails to build for any reason
#' -- callers should treat that as "no remap available, fall back to the
#' feature_id itself" exactly like remap_to_ensembl() itself does for an
#' empty/NULL map.
.ensembl_map_cache <- new.env(parent = emptyenv())
ensembl_map_for_dataset <- function(con, dataset_id) {
  if (exists(dataset_id, envir = .ensembl_map_cache, inherits = FALSE)) {
    return(get(dataset_id, envir = .ensembl_map_cache, inherits = FALSE))
  }
  meta <- metadata_source(con, dataset_id, "feature")
  map <- NULL
  if (!is.null(meta) && !is.na(meta$ensembl_col) && nzchar(meta$ensembl_col)) {
    pseudo_yaml <- list(dataset = list(
      id = dataset_id,
      feature_metadata_path = meta$path,
      feature_id_col = meta$id_col,
      ensembl_col = meta$ensembl_col
    ))
    map <- tryCatch(build_ensembl_map(pseudo_yaml), error = function(e) NULL)
  }
  if (is.null(map)) {
    message("ensembl_map_for_dataset('", dataset_id, "'): no usable Ensembl map -- ",
            "gene labels for this dataset will fall back to the feature_id itself ",
            "(re-run --stage core / R/ingest_results.R for this dataset if it predates ",
            "the ensembl_col column)")
  }
  assign(dataset_id, map, envir = .ensembl_map_cache)
  map
}

#' DISPLAY ONLY -- never used for cross-dataset matching/enrichment (see
#' this file's/symbol_mapping.R's headers; ensembl_map_for_dataset()/
#' remap_to_ensembl() are the computational path). Which feature_metadata
#' column to show gene names by, absent an explicit user pick: the
#' dataset's registered `dataset.symbol_col` (set explicitly in every
#' dataset's config -- see config/dataset_metadata.example.yml -- and
#' registered at ingest time via register_metadata_source()) if it names
#' a real column, else NA. Deliberately NOT auto-detected/guessed from a
#' fixed list of likely column names (an earlier version of this function
#' did that) -- explicit config is more reliable than a guess, and this
#' project's configs all set it explicitly now anyway.
resolve_symbol_col <- function(con, dataset_id, available_cols) {
  meta <- metadata_source(con, dataset_id, "feature")
  if (!is.null(meta) && !is.na(meta$symbol_col) && meta$symbol_col %in% available_cols) {
    return(meta$symbol_col)
  }
  NA_character_
}

#' Build a (feature_id -> display string) map for showing recognizable
#' gene names to users -- DISPLAY ONLY (see this function's/
#' resolve_symbol_col()'s headers; never used for cross-dataset matching --
#' enrichment computation happens entirely in the cluster ingest pipeline
#' now, not in this app).
#'
#' Fallback chain per feature, first non-blank wins: (1) `label_col` (an
#' arbitrary feature_metadata column name -- defaults to
#' resolve_symbol_col()'s pick, i.e. gene symbol, when NULL); (2) this
#' dataset's canonical Ensembl gene id (via ensembl_map_for_dataset(),
#' already version-stripped/first-of-multi-mapping -- see
#' R/lib/ingest/symbol_mapping.R); (3) the feature_id itself (always
#' present by definition, so this is the guaranteed final fallback -- a
#' displayed gene "name" should never be blank).
#'
#' Reads feature_metadata via dataset_metadata() (this file's existing
#' live-read helper -- reflects the source file's *current* contents, per
#' this file's header), not a fresh parquet read, for consistency with
#' every other metadata display in the app.
#'
#' Returns NULL if this dataset has no feature metadata registered at all.
build_display_map <- function(con, dataset_id, label_col = NULL) {
  m <- dataset_metadata(con, dataset_id, "feature")
  if (is.null(m)) return(NULL)

  resolved_label_col <- label_col %||% resolve_symbol_col(con, dataset_id, names(m))
  chosen <- if (!is.null(resolved_label_col) && !is.na(resolved_label_col) &&
                resolved_label_col %in% names(m)) {
    as.character(m[[resolved_label_col]])
  } else {
    rep(NA_character_, nrow(m))
  }

  ensembl_map <- ensembl_map_for_dataset(con, dataset_id)
  ensembl_vals <- if (!is.null(ensembl_map)) unname(ensembl_map[m$feature_id]) else rep(NA_character_, nrow(m))

  display <- ifelse(!is.na(chosen) & nzchar(chosen), chosen,
             ifelse(!is.na(ensembl_vals) & nzchar(ensembl_vals), ensembl_vals,
                    m$feature_id))
  setNames(display, m$feature_id)
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

    if (is_numeric) {
      # Vectorized across every component at once via cor()'s own
      # `use = "pairwise.complete.obs"` (per-column NA handling, exactly
      # matching the per-component `ok <- !is.na(v) & !is.na(y)` mask the
      # old per-gene loop computed one column at a time) + WGCNA's own
      # corPvalueStudent() (confirmed to reproduce cor.test(method=
      # "spearman", exact=FALSE)'s p-value exactly, real-data check
      # 2026-09-21). ~1000x faster than the equivalent per-column
      # cor.test() loop at real scale (8000 genes: ~2.0s loop vs ~0.06s
      # vectorized) -- this branch is what made WGCNA gene-significance's
      # numeric-field half a genuine minutes-per-dataset cost; the
      # categorical/Kruskal-Wallis branch below has no comparable
      # vectorized equivalent available in this project's dependencies
      # and is left as a per-component loop.
      ok_v <- !is.na(v)
      if (sum(ok_v) >= 3 && length(unique(v[ok_v])) >= 2) {
        sub_mat <- scores_mat[ok_v, , drop = FALSE]
        v_sub <- v[ok_v]
        cor_vals <- suppressWarnings(stats::cor(sub_mat, v_sub, method = "spearman",
                                                 use = "pairwise.complete.obs"))[, 1]
        n_vals <- colSums(!is.na(sub_mat))
        keep <- !is.na(cor_vals) & n_vals >= 3
        if (any(keep)) {
          p_vals <- WGCNA::corPvalueStudent(cor_vals[keep], n_vals[keep])
          rows[[length(rows) + 1]] <- data.frame(
            component = comps[keep], field = field, test = "spearman",
            statistic = unname(cor_vals[keep]), p_value = p_vals)
        }
      }
      next
    }

    for (j in seq_along(comps)) {
      y <- scores_mat[, j]
      ok <- !is.na(v) & !is.na(y)
      if (sum(ok) < 3) next
      res <- tryCatch({
        f <- as.factor(v[ok])
        if (nlevels(droplevels(f)) < 2) {
          NULL   # skip just this (field, component) cell -- NOT return(),
                  # which would abort the whole scan (see this function's
                  # header; a bug fixed 2026-09-21 after it surfaced via
                  # WGCNA gene-significance's much larger field x gene grid)
        } else {
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

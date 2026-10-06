# The stability DB (R/db/schema.sql, documented in docs/database.qmd) as a
# derived output of the targets store. Two writers, both deployment =
# "main" so writes never overlap:
#   write_dataset_db()       db_<dataset>: everything about one dataset,
#                            in one transaction
#   write_cross_dataset_db() db_cross_dataset: projections and gene-space
#                            agreement (they link two datasets' fits), then
#                            removal of artifact files no row references
# Other targets produce rows keyed by fit_key (+ factor_index); the writers
# resolve those keys to integer ids. Fits are upserted on fit_key, so a fit
# keeps its fit_id across rewrites; everything that hangs off a fit is
# deleted and re-inserted.

#' Current time as an ISO-8601 UTC string
#'
#' @return e.g. `"2026-10-05T16:13:15Z"`.
utc_now <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")

#' Run code in a transaction
#'
#' @param con DB connection.
#' @param code Code to run; rolled back if it errors.
#' @return `NULL`, invisibly.
with_transaction <- function(con, code) {
  DBI::dbExecute(con, "BEGIN")
  ok <- FALSE
  on.exit(if (!ok) DBI::dbExecute(con, "ROLLBACK"), add = TRUE)
  force(code)
  DBI::dbExecute(con, "COMMIT")
  ok <- TRUE
  invisible(NULL)
}

#' Append rows to a table
#'
#' @param con DB connection.
#' @param table Table name.
#' @param df Rows; `NULL` or empty is skipped.
#' @return Number of rows appended, invisibly.
append_rows <- function(con, table, df) {
  if (is.null(df) || nrow(df) == 0) return(invisible(0L))
  DBI::dbAppendTable(con, table, df)
  invisible(nrow(df))
}

#' Map keys to ids, failing on unmapped keys
#'
#' @param keys Keys (fit_key, or `"<fit_key> <factor_index>"`).
#' @param ids Named integer vector (names = keys).
#' @param what Label for the error message.
#' @return Integer ids; errors naming the keys that have none.
map_ids <- function(keys, ids, what) {
  out <- unname(ids[keys])
  bad <- unique(keys[is.na(out)])
  if (length(bad)) {
    stop(length(bad), " ", what, " key(s) have no id in the DB: ",
         paste(utils::head(bad, 5), collapse = ", "), call. = FALSE)
  }
  out
}

#' The next unused id of a table
#'
#' @param con DB connection.
#' @param table Table name.
#' @param col Id column.
#' @return `MAX(col) + 1`.
next_id <- function(con, table, col) {
  as.integer(DBI::dbGetQuery(con, sprintf("SELECT COALESCE(MAX(%s), 0) AS m FROM %s", col, table))$m) + 1L
}

#' Row-bind one component across methods
#'
#' @param methods Named list (by method) of the writer's per-method inputs.
#' @param getter `function(parts, method)` returning a data frame or `NULL`.
#' @return The bound data frame, or `NULL`.
rbind_parts <- function(methods, getter) {
  rows <- Filter(Negate(is.null), lapply(names(methods), function(m) getter(methods[[m]], m)))
  if (length(rows)) do.call(rbind, rows) else NULL
}

#' Write one dataset's rows (the `db_<dataset>` target)
#'
#' In one transaction: upserts the dataset and methods; deletes everything
#' hanging off the dataset's fits and its dataset-level rows; upserts fits on
#' fit_key (a fit keeps its fit_id) and deletes fits no longer in any grid;
#' then inserts every child table, resolving keys to ids.
#'
#' @param db_path Path to the DB.
#' @param schema The schema file target (a dependency only).
#' @param meta The resolved dataset config.
#' @param dataset_files Paths of the dataset's cached files (matrix.rds,
#'   sample_metadata.rds, feature_metadata.rds); the kind is the file name.
#' @param methods Named list (by method) of `list(ingest, pairs,
#'   representatives, redundancy, icasso, selection, wgcna_fits)`.
#' @param method_table method_table() output.
#' @param enrichment enrichment_rows() output.
#' @param matrix_diagnostics matrix_diagnostics_row() output.
#' @param wgcna wgcna_dataset_diagnostics() output, or `NULL`.
#' @param drivers pattern_driver_rows() output.
#' @return `list(dataset_id, fits, written_at)`: fit counts by status.
write_dataset_db <- function(db_path, schema, meta, dataset_files, methods, method_table,
                             enrichment = NULL, matrix_diagnostics = NULL, wgcna = NULL,
                             drivers = NULL) {
  ds <- meta$dataset
  id <- ds$id
  now <- utc_now()
  con <- open_stability_db(db_path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  q <- function(sql, ...) DBI::dbExecute(con, sql, params = list(...))

  with_transaction(con, {
    # ---- dataset + methods ----
    q("INSERT INTO datasets (dataset_id, description, n_samples, n_genes, sample_id_col, feature_id_col,
         ensembl_col, symbol_col, subject_id_col, sample_metadata_source, feature_metadata_source, written_at)
       VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
       ON CONFLICT(dataset_id) DO UPDATE SET description = excluded.description,
         n_samples = excluded.n_samples, n_genes = excluded.n_genes,
         sample_id_col = excluded.sample_id_col, feature_id_col = excluded.feature_id_col,
         ensembl_col = excluded.ensembl_col, symbol_col = excluded.symbol_col,
         subject_id_col = excluded.subject_id_col,
         sample_metadata_source = excluded.sample_metadata_source,
         feature_metadata_source = excluded.feature_metadata_source, written_at = excluded.written_at",
      id, trimws(ds$description %||% NA_character_), matrix_diagnostics$n_samples, matrix_diagnostics$n_genes,
      ds$sample_id_col %||% NA_character_, ds$feature_id_col %||% NA_character_,
      ds$ensembl_col %||% NA_character_, ds$symbol_col %||% NA_character_,
      ds$subject_id_col %||% NA_character_, ds$sample_metadata_path %||% NA_character_,
      ds$feature_metadata_path %||% NA_character_, now)
    for (i in seq_len(nrow(method_table))) {
      r <- method_table[i, ]
      q("INSERT INTO methods (method, family, sign_ambiguous, has_loadings) VALUES (?, ?, ?, ?)
         ON CONFLICT(method) DO UPDATE SET family = excluded.family,
           sign_ambiguous = excluded.sign_ambiguous, has_loadings = excluded.has_loadings",
        r$method, r$family, as.integer(r$sign_ambiguous), as.integer(r$has_loadings))
    }

    # ---- clear what hangs off this dataset's fits, and its dataset-level rows ----
    in_ds <- "(SELECT fit_id FROM fits WHERE dataset_id = ?)"
    q(sprintf("DELETE FROM factors WHERE fit_id IN %s", in_ds), id)    # cascades factor-level rows
    for (tab in c("fit_params", "fit_metrics", "fit_artifacts", "gene_modules", "module_kme_summary",
                  "module_kme_histograms", "fit_redundancy")) {
      q(sprintf("DELETE FROM %s WHERE fit_id IN %s", tab, in_ds), id)
    }
    q(sprintf("DELETE FROM fit_pairs WHERE fit_a IN %s OR fit_b IN %s", in_ds, in_ds), id, id)
    for (tab in c("dataset_artifacts", "matrix_diagnostics", "similarity_histograms", "wgcna_sft",
                  "gene_significance", "icasso_clusters")) {
      q(sprintf("DELETE FROM %s WHERE dataset_id = ?", tab), id)
    }

    # ---- fits: upsert by fit_key, drop fits no longer in any grid ----
    fits <- rbind_parts(methods, function(p, m) {
      f <- p$ingest$fits
      cbind(f, dataset_id = id, method = m, representative = as.integer(f$fit_key %in% p$representatives))
    })
    fits$written_at <- now
    DBI::dbWriteTable(con, "new_fits", fits, temporary = TRUE, overwrite = TRUE)
    q("DELETE FROM fits WHERE dataset_id = ? AND fit_key NOT IN (SELECT fit_key FROM temp.new_fits)", id)
    cols <- c("fit_key", "dataset_id", "method", "rank", "seed", "bootstrap", "representative", "mse",
              "n_factors", "status", "error", "written_at")
    DBI::dbExecute(con, sprintf(
      "INSERT INTO fits (%s) SELECT %s FROM temp.new_fits WHERE true
       ON CONFLICT(fit_key) DO UPDATE SET %s",
      paste(cols, collapse = ", "), paste(cols, collapse = ", "),
      paste(sprintf("%s = excluded.%s", cols[-1], cols[-1]), collapse = ", ")))
    DBI::dbExecute(con, "DROP TABLE temp.new_fits")
    fit_ids <- with(DBI::dbGetQuery(con, "SELECT fit_id, fit_key FROM fits WHERE dataset_id = ?",
                                    params = list(id)), stats::setNames(fit_id, fit_key))
    fid <- function(keys) map_ids(keys, fit_ids, "fit")
    by_fit <- function(df) {
      if (is.null(df) || nrow(df) == 0) return(NULL)
      df$fit_id <- fid(df$fit_key)
      df[, c("fit_id", setdiff(names(df), c("fit_id", "fit_key"))), drop = FALSE]
    }

    # ---- per-fit long tables ----
    append_rows(con, "fit_params", by_fit(rbind_parts(methods, function(p, m) p$ingest$params)))
    append_rows(con, "fit_metrics", by_fit(rbind_parts(methods, function(p, m) {
      rbind(p$ingest$metrics, selection_metric_rows(p$ingest, p$selection))
    })))
    append_rows(con, "fit_artifacts", by_fit(rbind_parts(methods, function(p, m) {
      rbind(p$ingest$artifacts, p$redundancy$artifacts,
            do.call(rbind, lapply(p$wgcna_fits, `[[`, "artifacts")))
    })))

    # ---- factors (with per-factor stability) ----
    factors <- rbind_parts(methods, function(p, m) {
      f <- p$ingest$factors
      if (is.null(f)) return(NULL)
      st <- p$pairs$stability
      hit <- if (is.null(st)) rep(NA_integer_, nrow(f)) else
        match(paste(f$fit_key, f$factor_index), paste(st$fit_key, st$factor_index))
      f$stability_cosine <- if (is.null(st)) NA_real_ else st$cosine[hit]
      f$stability_pearson <- if (is.null(st)) NA_real_ else st$pearson[hit]
      f$stability_spearman <- if (is.null(st)) NA_real_ else st$spearman[hit]
      f
    })
    append_rows(con, "factors", by_fit(factors))
    factor_ids <- with(DBI::dbGetQuery(con,
      "SELECT x.factor_id, f.fit_key, x.factor_index FROM factors x JOIN fits f ON f.fit_id = x.fit_id
       WHERE f.dataset_id = ?", params = list(id)), stats::setNames(factor_id, paste(fit_key, factor_index)))
    xid <- function(keys, idx) map_ids(paste(keys, idx), factor_ids, "factor")

    # ---- pair stability ----
    next_pair <- next_id(con, "fit_pairs", "fit_pair_id")
    next_match <- next_id(con, "factor_matches", "match_id")
    for (m in names(methods)) {
      pr <- methods[[m]]$pairs
      if (is.null(pr) || is.null(pr$fit_pairs)) next
      fp <- pr$fit_pairs
      ia <- fid(fp$fit_a); ib <- fid(fp$fit_b)
      swap <- ia > ib                                  # schema stores fit_a < fit_b
      fp$fit_pair_id <- seq.int(next_pair, length.out = nrow(fp)); next_pair <- next_pair + nrow(fp)
      pair_id <- stats::setNames(fp$fit_pair_id, paste(fp$fit_a, fp$fit_b))
      append_rows(con, "fit_pairs", data.frame(
        fit_pair_id = fp$fit_pair_id, fit_a = ifelse(swap, ib, ia), fit_b = ifelse(swap, ia, ib),
        same_rank = fp$same_rank, n_factors_a = ifelse(swap, fp$n_factors_b, fp$n_factors_a),
        n_factors_b = ifelse(swap, fp$n_factors_a, fp$n_factors_b), ari = fp$ari))
      if (!is.null(pr$matches) && nrow(pr$matches)) {
        mt <- pr$matches
        xa <- xid(mt$fit_a, mt$factor_a); xb <- xid(mt$fit_b, mt$factor_b)
        sw <- fid(mt$fit_a) > fid(mt$fit_b)
        match_id <- stats::setNames(seq.int(next_match, length.out = nrow(mt)), mt$match)
        next_match <- next_match + nrow(mt)
        append_rows(con, "factor_matches", data.frame(
          match_id = unname(match_id), fit_pair_id = unname(pair_id[paste(mt$fit_a, mt$fit_b)]),
          factor_a = ifelse(sw, xb, xa), factor_b = ifelse(sw, xa, xb)))
        sc <- pr$scores
        append_rows(con, "match_scores", data.frame(
          match_id = unname(match_id[as.character(sc$match)]), metric = sc$metric, value = sc$value,
          runner_up = sc$runner_up, margin = sc$margin))
      }
      nl <- pr$null
      if (!is.null(nl) && nrow(nl)) {
        append_rows(con, "fit_pair_null", cbind(
          data.frame(fit_pair_id = unname(pair_id[paste(nl$fit_a, nl$fit_b)])),
          nl[, c("metric", "n", "mean", "sd", "median", "p95", "max")]))
      }
      if (!is.null(pr$histograms)) {
        append_rows(con, "similarity_histograms", cbind(data.frame(dataset_id = id, method = m), pr$histograms))
      }
    }

    # ---- WGCNA per-fit kME ----
    for (m in names(methods)) {
      wf <- Filter(Negate(is.null), methods[[m]]$wgcna_fits)
      if (length(wf) == 0) next
      append_rows(con, "gene_modules", by_fit(do.call(rbind, lapply(wf, `[[`, "gene_modules"))))
      append_rows(con, "module_kme_summary", by_fit(do.call(rbind, lapply(wf, `[[`, "summary"))))
      append_rows(con, "module_kme_histograms", by_fit(do.call(rbind, lapply(wf, `[[`, "histograms"))))
    }

    # ---- representative-fit analyses ----
    mk <- rbind_parts(methods, function(p, m) p$redundancy$pattern_markers)
    if (!is.null(mk)) append_rows(con, "pattern_markers",
                                  data.frame(factor_id = xid(mk$fit_key, mk$factor_index),
                                             gene = mk$gene, score = mk$score))
    append_rows(con, "fit_redundancy", by_fit(rbind_parts(methods, function(p, m) p$redundancy$fit_redundancy)))
    if (!is.null(drivers) && nrow(drivers)) {
      append_rows(con, "pattern_drivers", cbind(
        data.frame(factor_id = xid(drivers$fit_key, drivers$factor_index)),
        drivers[, setdiff(names(drivers), c("fit_key", "factor_index"))]))
    }

    # ---- ICASSO ----
    for (m in names(methods)) {
      for (ic in methods[[m]]$icasso) {
        cl <- ic$clusters
        append_rows(con, "icasso_clusters", data.frame(
          dataset_id = id, rank = ic$rank, cluster_id = cl$cluster_id, iq = cl$iq, n_members = cl$n_members,
          centrotype_factor_id = xid(cl$centrotype_fit_key, cl$centrotype_factor_index),
          dendro_file = cl$dendro_file))
        mb <- ic$membership
        append_rows(con, "icasso_membership", data.frame(
          factor_id = xid(mb$fit_key, mb$factor_index), dataset_id = id, rank = ic$rank,
          cluster_id = mb$cluster_id, intra_sim = mb$intra_sim))
      }
    }

    # ---- enrichment ----
    eq <- enrichment$queries
    if (!is.null(eq) && nrow(eq)) {
      known <- paste(eq$fit_key, eq$factor_index) %in% names(factor_ids)
      if (any(!known)) warning(sum(!known), " enrichment queries for factors not in the DB were dropped")
      eq <- eq[known, , drop = FALSE]
      qkey <- paste(eq$fit_key, eq$factor_index, eq$query_type, eq$direction)
      eq$query_id <- seq.int(next_id(con, "enrichment_queries", "query_id"), length.out = nrow(eq))
      append_rows(con, "enrichment_queries", data.frame(
        query_id = eq$query_id, factor_id = xid(eq$fit_key, eq$factor_index), query_type = eq$query_type,
        direction = eq$direction, query_size = eq$query_size, queried_at = now))
      er <- enrichment$results
      if (!is.null(er) && nrow(er)) {
        rq <- stats::setNames(eq$query_id, qkey)[paste(er$fit_key, er$factor_index, er$query_type, er$direction)]
        er <- er[!is.na(rq), , drop = FALSE]
        append_rows(con, "enrichment_results", cbind(
          data.frame(query_id = unname(rq[!is.na(rq)])),
          er[, c("term_id", "source", "term_name", "p_value", "intersection_size", "term_size",
                 "nes", "es", "log2err", "is_main_pathway", "genes")]))
      }
    }

    # ---- dataset-level ----
    if (length(dataset_files)) {
      append_rows(con, "dataset_artifacts", data.frame(
        dataset_id = id, kind = sub("\\.rds$", "", basename(dataset_files)),
        path = artifact_rel(id, basename(dataset_files))))
    }
    if (!is.null(matrix_diagnostics)) {
      append_rows(con, "matrix_diagnostics", cbind(
        data.frame(dataset_id = id),
        matrix_diagnostics[, setdiff(names(matrix_diagnostics), c("n_genes", "n_samples"))],
        data.frame(computed_at = now)))
    }
    if (!is.null(wgcna)) {
      if (!is.null(wgcna$sft)) append_rows(con, "wgcna_sft", cbind(data.frame(dataset_id = id), wgcna$sft))
      if (!is.null(wgcna$significance)) {
        append_rows(con, "gene_significance", cbind(data.frame(dataset_id = id), wgcna$significance))
      }
    }
  })
  counts <- DBI::dbGetQuery(con, "SELECT status, COUNT(*) AS n FROM fits WHERE dataset_id = ? GROUP BY status",
                            params = list(id))
  list(dataset_id = id, fits = stats::setNames(counts$n, counts$status), written_at = now)
}

#' Write projections and gene-space agreement (the `db_cross_dataset` target)
#'
#' Runs after every `db_<dataset>`, so fits on both sides resolve; replaces
#' both tables, then removes artifact files no row references.
#'
#' @param db_path Path to the DB.
#' @param written The `db_<dataset>` values (a dependency only).
#' @param projections Flat list of projection_rows() data frames.
#' @param gene_space List of gene_space_rows() data frames.
#' @return `list(n_projections, n_gene_space, orphans_removed, written_at)`.
write_cross_dataset_db <- function(db_path, written, projections = list(), gene_space = list()) {
  con <- open_stability_db(db_path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  fit_ids <- with(DBI::dbGetQuery(con, "SELECT fit_id, fit_key FROM fits"), stats::setNames(fit_id, fit_key))
  pr <- do.call(rbind, Filter(Negate(is.null), projections))
  gs <- do.call(rbind, Filter(Negate(is.null), gene_space))
  with_transaction(con, {
    DBI::dbExecute(con, "DELETE FROM projections")
    DBI::dbExecute(con, "DELETE FROM gene_space_agreement")
    if (!is.null(pr) && nrow(pr)) {
      append_rows(con, "projections", cbind(
        data.frame(source_fit_id = map_ids(pr$source_fit_key, fit_ids, "fit")),
        pr[, setdiff(names(pr), "source_fit_key")]))
    }
    if (!is.null(gs) && nrow(gs)) {
      append_rows(con, "gene_space_agreement", data.frame(
        source_fit_id = map_ids(gs$source_fit_key, fit_ids, "fit"),
        target_fit_id = map_ids(gs$target_fit_key, fit_ids, "fit"),
        n_genes_matched = gs$n_genes_matched, mean_abs_diagonal = gs$mean_abs_diagonal))
    }
  })
  removed <- vapply(DBI::dbGetQuery(con, "SELECT dataset_id FROM datasets")$dataset_id,
                    function(d) remove_orphan_artifacts(con, db_path, d), integer(1))
  list(n_projections = if (is.null(pr)) 0L else nrow(pr), n_gene_space = if (is.null(gs)) 0L else nrow(gs),
       orphans_removed = removed, written_at = utc_now())
}

#' Delete a dataset's unreferenced artifact files
#'
#' Files under the dataset's artifact directory not referenced by
#' `fit_artifacts`, `dataset_artifacts`, `pattern_drivers`, `projections` or
#' `icasso_clusters` (e.g. fits dropped from a grid, keys from an older run).
#'
#' @param con DB connection.
#' @param db_path Path to the DB.
#' @param dataset_id Dataset id.
#' @return Number of files deleted.
remove_orphan_artifacts <- function(con, db_path, dataset_id) {
  ds <- list(dataset_id)
  q <- function(sql) unlist(DBI::dbGetQuery(con, sql, params = ds), use.names = FALSE)
  in_ds <- "(SELECT fit_id FROM fits WHERE dataset_id = ?)"
  referenced <- c(
    q(sprintf("SELECT path FROM fit_artifacts WHERE fit_id IN %s", in_ds)),
    q("SELECT path FROM dataset_artifacts WHERE dataset_id = ?"),
    q(sprintf("SELECT d.path FROM pattern_drivers d JOIN factors x ON x.factor_id = d.factor_id
               WHERE x.fit_id IN %s", in_ds)),
    q(sprintf("SELECT path FROM projections WHERE source_fit_id IN %s", in_ds)),
    q("SELECT DISTINCT dendro_file FROM icasso_clusters WHERE dataset_id = ?"))
  prefix <- paste0("stability_artifacts/", dataset_id, "/")
  referenced <- sub(prefix, "", referenced[!is.na(referenced)], fixed = TRUE)
  dir <- artifacts_dir(db_path, dataset_id)
  orphans <- setdiff(list.files(dir, recursive = TRUE), referenced)
  unlink(file.path(dir, orphans))
  length(orphans)
}

# Pairwise stability computation, run after a family's fits are inserted.
#
# Incremental by design: given the fit_ids ADDED this run, pairs are
# computed as new x (new + existing) fits of the same method -- so adding
# a family to a DB that already holds others of the same method still
# yields the complete all-rank-pairs x all-seed-pairs set, without
# recomputing pairs between fits that were both already present.
#
# Every "pure" function below (suffix `_from_universe`) takes a plain
# `fits`/`universe` data.frame instead of a live DB connection, so it can
# run during a PARALLEL, no-DB-writer compute phase (see
# R/lib/ingest/ingest_dataset.R::compute_ingest_bundle()) -- a fit not
# yet inserted into the DB has no real `fit_id` yet, so these functions
# use a NEGATIVE id (`-local_id`) to mean "a fit from THIS bundle, not yet
# merged" versus a positive id, which always means a real, already-merged
# `fits.fit_id`. The thin wrapper below each pure function preserves the
# original DB-querying, con-taking API for the single-dataset direct path
# (R/ingest_results.R) and for the merge phase's use once real fit_ids
# exist for everything.

#' All-factor-combination similarity rows for every (new x all) fit pair
#' of one method, from an in-memory universe (no DB access). `universe`
#' must have columns `id` (positive = already-merged real fit_id,
#' negative = this bundle's local id), `rank`, `loadings_file_abs`
#' (already resolved to an absolute path). `new_ids` is the subset of
#' `universe$id` to treat as new (at least one side of a pair must be
#' new). Loading matrices are prepared (see similarity.R) once per fit
#' and cached in memory for the duration of the call.
compute_factor_pairs_from_universe <- function(universe, new_ids) {
  if (nrow(universe) < 2 || length(new_ids) == 0) return(NULL)

  is_new <- universe$id %in% new_ids
  idx <- which(upper.tri(diag(nrow(universe))), arr.ind = TRUE)
  keep <- is_new[idx[, 1]] | is_new[idx[, 2]]
  idx <- idx[keep, , drop = FALSE]
  if (nrow(idx) == 0) return(NULL)

  prep_cache <- new.env(parent = emptyenv())
  get_prep <- function(i) {
    key <- as.character(universe$id[i])
    if (!exists(key, envir = prep_cache)) {
      mat <- readRDS(universe$loadings_file_abs[i])
      assign(key, prepare_loadings(mat), envir = prep_cache)
    }
    get(key, envir = prep_cache)
  }

  rows <- list()
  for (p in seq_len(nrow(idx))) {
    i <- idx[p, 1]; j <- idx[p, 2]
    sims <- pair_similarities(get_prep(i), get_prep(j))
    if (is.null(sims)) next
    match_idx <- hungarian_match(sims$cosine)
    matched_flag <- matrix(0L, nrow(sims$cosine), ncol(sims$cosine))
    matched_flag[match_idx] <- 1L

    grid_idx <- expand.grid(factor_a = seq_len(nrow(sims$cosine)),
                            factor_b = seq_len(ncol(sims$cosine)))
    rows[[length(rows) + 1]] <- data.frame(
      fit_a    = universe$id[i],
      fit_b    = universe$id[j],
      factor_a = grid_idx$factor_a,
      factor_b = grid_idx$factor_b,
      cosine   = as.vector(sims$cosine),
      pearson  = as.vector(sims$pearson),
      spearman = as.vector(sims$spearman),
      matched  = as.integer(matched_flag[cbind(grid_idx$factor_a, grid_idx$factor_b)]),
      same_rank = as.integer(!is.na(universe$rank[i]) && !is.na(universe$rank[j]) &&
                               universe$rank[i] == universe$rank[j])
    )
  }
  if (length(rows) == 0) return(NULL)
  do.call(rbind, rows)
}

#' Thin DB-querying wrapper around compute_factor_pairs_from_universe() --
#' unchanged public behavior for the single-dataset direct path
#' (R/ingest_results.R) and for merge-time use once every fit (old and
#' newly-merged) already has a real fit_id in the DB.
compute_factor_pairs <- function(con, db_path, dataset_id, method, new_fit_ids) {
  fits <- DBI::dbGetQuery(con,
    "SELECT fit_id AS id, rank, loadings_file FROM fits
     WHERE dataset_id = ? AND method = ? AND status = 'ok'
       AND loadings_file IS NOT NULL",
    params = list(dataset_id, method))
  if (nrow(fits) < 2 || length(new_fit_ids) == 0) return(invisible(0L))
  fits$loadings_file_abs <- vapply(fits$loadings_file, resolve_artifact,
                                    character(1), db_path = db_path)

  rows <- compute_factor_pairs_from_universe(
    fits[, c("id", "rank", "loadings_file_abs")], new_fit_ids)
  if (is.null(rows)) return(invisible(0L))

  for (start in seq(1, nrow(rows), by = 200)) {
    chunk <- rows[start:min(start + 199, nrow(rows)), , drop = FALSE]
    DBI::dbWriteTable(con, "factor_pairs", chunk, append = TRUE)
  }
  invisible(nrow(rows))
}

#' Per-factor stability summaries: median matched similarity across all
#' SAME-RANK pairs the factor participates in, from an in-memory
#' `factor_pairs`-shaped data.frame (no DB access) -- `pairs` may mix
#' negative (this bundle's local id) and positive (already-merged) `fit_a`/
#' `fit_b` values, exactly as produced by compute_factor_pairs_from_universe().
#' Returns a data.frame(fit_id, factor_index, cosine, pearson, spearman)
#' -- one row per (fit, factor) needing a stability update; the caller
#' applies it (remapping negative ids to real ones first, if needed).
update_factor_stability_from_pairs <- function(pairs) {
  pairs <- pairs[pairs$matched == 1 & pairs$same_rank == 1, , drop = FALSE]
  if (nrow(pairs) == 0) return(NULL)

  long <- rbind(
    data.frame(fit_id = pairs$fit_a, factor_index = pairs$factor_a,
               cosine = pairs$cosine, pearson = pairs$pearson, spearman = pairs$spearman),
    data.frame(fit_id = pairs$fit_b, factor_index = pairs$factor_b,
               cosine = pairs$cosine, pearson = pairs$pearson, spearman = pairs$spearman)
  )
  aggregate(cbind(cosine, pearson, spearman) ~ fit_id + factor_index,
            data = long, FUN = median)
}

#' Thin DB-querying wrapper: reads this dataset+method's factor_pairs
#' (every fit already has a real fit_id at this point) and applies the
#' resulting stability UPDATEs. No BEGIN/COMMIT here -- this function's
#' callers (R/lib/ingest/ingest_dataset.R) always run it inside their own
#' wrapping transaction; a nested BEGIN here would error ("cannot start a
#' transaction within a transaction").
update_factor_stability <- function(con, dataset_id, method) {
  pairs <- DBI::dbGetQuery(con,
    "SELECT fp.fit_a, fp.fit_b, fp.factor_a, fp.factor_b,
            fp.cosine, fp.pearson, fp.spearman, fp.matched, fp.same_rank
     FROM factor_pairs fp
     JOIN fits fa ON fa.fit_id = fp.fit_a
     WHERE fa.dataset_id = ? AND fa.method = ?",
    params = list(dataset_id, method))
  agg <- update_factor_stability_from_pairs(pairs)
  if (is.null(agg)) return(invisible(NULL))

  for (r in seq_len(nrow(agg))) {
    DBI::dbExecute(con,
      "UPDATE factors SET stability_cosine = ?, stability_pearson = ?, stability_spearman = ?
       WHERE fit_id = ? AND factor_index = ?",
      params = list(agg$cosine[r], agg$pearson[r], agg$spearman[r],
                    agg$fit_id[r], agg$factor_index[r]))
  }
  invisible(NULL)
}

#' WGCNA: ARI between module assignments for every (new x all) fit pair,
#' plus per-module Jaccard overlaps with Hungarian matching -- pure,
#' in-memory version (no DB access). `mod_list` is a named list (keys =
#' as.character(id), matching `ids`) of data.frames with columns
#' `gene`/`module`, for every fit (old + new) of this dataset. `ids`
#' holds every fit's id (positive = already-merged, negative = this
#' bundle's local id); `new_ids` the subset to treat as new. Returns
#' list(fit_pairs = data.frame(fit_a, fit_b, ari, n_modules_a,
#' n_modules_b), module_pairs = data.frame(fit_a, module_a, fit_b,
#' module_b, jaccard, matched)) or NULL if nothing to compute.
compute_wgcna_pairs_from_universe <- function(mod_list, ids, new_ids) {
  if (length(ids) < 2 || length(new_ids) == 0) return(NULL)
  is_new <- ids %in% new_ids

  fit_pair_rows <- list()
  module_pair_rows <- list()
  for (i in seq_along(ids)) {
    for (j in seq_along(ids)) {
      if (i >= j) next
      if (!is_new[i] && !is_new[j]) next
      a <- mod_list[[as.character(ids[i])]]
      b <- mod_list[[as.character(ids[j])]]
      shared <- intersect(a$gene, b$gene)
      la <- a$module[match(shared, a$gene)]
      lb <- b$module[match(shared, b$gene)]
      ari <- mclust::adjustedRandIndex(la, lb)

      fit_pair_rows[[length(fit_pair_rows) + 1]] <- data.frame(
        fit_a = ids[i], fit_b = ids[j], ari = ari,
        n_modules_a = length(setdiff(unique(la), 0L)),
        n_modules_b = length(setdiff(unique(lb), 0L)))

      mods_a <- setdiff(sort(unique(la)), 0L)
      mods_b <- setdiff(sort(unique(lb)), 0L)
      if (length(mods_a) == 0 || length(mods_b) == 0) next
      jac <- matrix(0, length(mods_a), length(mods_b))
      genes_a <- lapply(mods_a, function(m) shared[la == m])
      genes_b <- lapply(mods_b, function(m) shared[lb == m])
      for (x in seq_along(mods_a)) {
        for (y in seq_along(mods_b)) {
          inter <- length(intersect(genes_a[[x]], genes_b[[y]]))
          uni   <- length(genes_a[[x]]) + length(genes_b[[y]]) - inter
          jac[x, y] <- if (uni > 0) inter / uni else 0
        }
      }
      match_idx <- hungarian_match(jac)
      matched_flag <- matrix(0L, nrow(jac), ncol(jac))
      matched_flag[match_idx] <- 1L
      grid_idx <- expand.grid(x = seq_along(mods_a), y = seq_along(mods_b))
      module_pair_rows[[length(module_pair_rows) + 1]] <- data.frame(
        fit_a = ids[i], module_a = mods_a[grid_idx$x],
        fit_b = ids[j], module_b = mods_b[grid_idx$y],
        jaccard = jac[cbind(grid_idx$x, grid_idx$y)],
        matched = as.integer(matched_flag[cbind(grid_idx$x, grid_idx$y)]))
    }
  }
  if (length(fit_pair_rows) == 0) return(NULL)
  list(fit_pairs = do.call(rbind, fit_pair_rows),
       module_pairs = if (length(module_pair_rows) == 0) NULL else do.call(rbind, module_pair_rows))
}

#' Thin DB-querying wrapper around compute_wgcna_pairs_from_universe() --
#' unchanged public behavior for the single-dataset direct path and for
#' merge-time use once every fit already has a real fit_id.
compute_wgcna_pairs <- function(con, dataset_id, new_fit_ids) {
  fits <- DBI::dbGetQuery(con,
    "SELECT fit_id FROM fits
     WHERE dataset_id = ? AND method = 'wgcna' AND status = 'ok'",
    params = list(dataset_id))
  if (nrow(fits) < 2 || length(new_fit_ids) == 0) return(invisible(NULL))

  mods <- DBI::dbGetQuery(con, sprintf(
    "SELECT fit_id, gene, module FROM wgcna_modules WHERE fit_id IN (%s)",
    paste(fits$fit_id, collapse = ",")))
  mod_list <- split(mods[, c("gene", "module")], mods$fit_id)

  out <- compute_wgcna_pairs_from_universe(mod_list, fits$fit_id, new_fit_ids)
  if (is.null(out)) return(invisible(NULL))
  DBI::dbWriteTable(con, "wgcna_fit_pairs", out$fit_pairs, append = TRUE)
  if (!is.null(out$module_pairs)) {
    DBI::dbWriteTable(con, "wgcna_module_pairs", out$module_pairs, append = TRUE)
  }
  invisible(NULL)
}

# Pairwise stability computation, run after a family's fits are inserted.
#
# Incremental by design: given the fit_ids ADDED this run, pairs are
# computed as new x (new + existing) fits of the same method -- so adding
# a family to a DB that already holds others of the same method still
# yields the complete all-rank-pairs x all-seed-pairs set, without
# recomputing pairs between fits that were both already present.

#' All-factor-combination similarity rows for every (new x all) fit pair
#' of one method. Loading matrices are prepared (see similarity.R) once
#' per fit and cached in memory for the duration of the call.
compute_factor_pairs <- function(con, db_path, dataset_id, method, new_fit_ids) {
  fits <- DBI::dbGetQuery(con,
    "SELECT fit_id, rank, loadings_file FROM fits
     WHERE dataset_id = ? AND method = ? AND status = 'ok'
       AND loadings_file IS NOT NULL",
    params = list(dataset_id, method))
  if (nrow(fits) < 2 || length(new_fit_ids) == 0) return(invisible(0L))
  fits$loadings_file <- vapply(fits$loadings_file, resolve_artifact,
                               character(1), db_path = db_path)

  is_new <- fits$fit_id %in% new_fit_ids
  # candidate pairs: canonical order i < j, at least one side new (pairs
  # where both fits were already in the DB were computed on a prior ingest)
  idx <- which(upper.tri(diag(nrow(fits))), arr.ind = TRUE)
  keep <- is_new[idx[, 1]] | is_new[idx[, 2]]
  idx <- idx[keep, , drop = FALSE]
  if (nrow(idx) == 0) return(invisible(0L))

  prep_cache <- new.env(parent = emptyenv())
  get_prep <- function(i) {
    key <- as.character(fits$fit_id[i])
    if (!exists(key, envir = prep_cache)) {
      mat <- readRDS(fits$loadings_file[i])
      assign(key, prepare_loadings(mat), envir = prep_cache)
    }
    get(key, envir = prep_cache)
  }

  n_rows_written <- 0L
  batch <- list()
  flush <- function() {
    if (length(batch) > 0) {
      DBI::dbWriteTable(con, "factor_pairs", do.call(rbind, batch), append = TRUE)
    }
    length(batch) > 0
  }

  for (p in seq_len(nrow(idx))) {
    i <- idx[p, 1]; j <- idx[p, 2]
    sims <- pair_similarities(get_prep(i), get_prep(j))
    if (is.null(sims)) next
    match_idx <- hungarian_match(sims$cosine)
    matched_flag <- matrix(0L, nrow(sims$cosine), ncol(sims$cosine))
    matched_flag[match_idx] <- 1L

    grid_idx <- expand.grid(factor_a = seq_len(nrow(sims$cosine)),
                            factor_b = seq_len(ncol(sims$cosine)))
    batch[[length(batch) + 1]] <- data.frame(
      fit_a    = fits$fit_id[i],
      fit_b    = fits$fit_id[j],
      factor_a = grid_idx$factor_a,
      factor_b = grid_idx$factor_b,
      cosine   = as.vector(sims$cosine),
      pearson  = as.vector(sims$pearson),
      spearman = as.vector(sims$spearman),
      matched  = as.integer(matched_flag[cbind(grid_idx$factor_a, grid_idx$factor_b)]),
      same_rank = as.integer(!is.na(fits$rank[i]) && !is.na(fits$rank[j]) &&
                               fits$rank[i] == fits$rank[j])
    )
    n_rows_written <- n_rows_written + nrow(grid_idx)

    if (length(batch) >= 200) {
      flush(); batch <- list()
    }
  }
  flush()
  invisible(n_rows_written)
}

#' Per-factor stability summaries: median matched similarity across all
#' SAME-RANK pairs the factor participates in. Recomputed for every fit of
#' the method (new pairs change existing factors' summaries too).
update_factor_stability <- function(con, dataset_id, method) {
  pairs <- DBI::dbGetQuery(con,
    "SELECT fp.fit_a, fp.fit_b, fp.factor_a, fp.factor_b,
            fp.cosine, fp.pearson, fp.spearman
     FROM factor_pairs fp
     JOIN fits fa ON fa.fit_id = fp.fit_a
     WHERE fp.matched = 1 AND fp.same_rank = 1
       AND fa.dataset_id = ? AND fa.method = ?",
    params = list(dataset_id, method))
  if (nrow(pairs) == 0) return(invisible(NULL))

  long <- rbind(
    data.frame(fit_id = pairs$fit_a, factor_index = pairs$factor_a,
               cosine = pairs$cosine, pearson = pairs$pearson, spearman = pairs$spearman),
    data.frame(fit_id = pairs$fit_b, factor_index = pairs$factor_b,
               cosine = pairs$cosine, pearson = pairs$pearson, spearman = pairs$spearman)
  )
  agg <- aggregate(cbind(cosine, pearson, spearman) ~ fit_id + factor_index,
                   data = long, FUN = median)

  DBI::dbExecute(con, "BEGIN")
  for (r in seq_len(nrow(agg))) {
    DBI::dbExecute(con,
      "UPDATE factors SET stability_cosine = ?, stability_pearson = ?, stability_spearman = ?
       WHERE fit_id = ? AND factor_index = ?",
      params = list(agg$cosine[r], agg$pearson[r], agg$spearman[r],
                    agg$fit_id[r], agg$factor_index[r]))
  }
  DBI::dbExecute(con, "COMMIT")
  invisible(NULL)
}

#' WGCNA: ARI between module assignments for every (new x all) fit pair,
#' plus per-module Jaccard overlaps with Hungarian matching.
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

  ids <- fits$fit_id
  is_new <- ids %in% new_fit_ids
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

      DBI::dbExecute(con,
        "INSERT INTO wgcna_fit_pairs (fit_a, fit_b, ari, n_modules_a, n_modules_b)
         VALUES (?, ?, ?, ?, ?)",
        params = list(ids[i], ids[j], ari,
                      length(setdiff(unique(la), 0L)),
                      length(setdiff(unique(lb), 0L))))

      # module x module Jaccard (excluding module 0 = unassigned)
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
      DBI::dbWriteTable(con, "wgcna_module_pairs", data.frame(
        fit_a = ids[i], module_a = mods_a[grid_idx$x],
        fit_b = ids[j], module_b = mods_b[grid_idx$y],
        jaccard = jac[cbind(grid_idx$x, grid_idx$y)],
        matched = as.integer(matched_flag[cbind(grid_idx$x, grid_idx$y)])
      ), append = TRUE)
    }
  }
  invisible(NULL)
}

#' wTO: per-edge value correlation + significant-edge Jaccard for every
#' (new x all) run pair. Edge tables are parquet artifacts (path stored in
#' fits.loadings_file), read via arrow.
compute_wto_pairs <- function(con, db_path, dataset_id, new_fit_ids, padj_cutoff = 0.05) {
  fits <- DBI::dbGetQuery(con,
    "SELECT fit_id, loadings_file FROM fits
     WHERE dataset_id = ? AND method = 'wto' AND status = 'ok'
       AND loadings_file IS NOT NULL",
    params = list(dataset_id))
  if (nrow(fits) < 2 || length(new_fit_ids) == 0) return(invisible(NULL))
  fits$loadings_file <- vapply(fits$loadings_file, resolve_artifact,
                               character(1), db_path = db_path)

  edge_cache <- new.env(parent = emptyenv())
  get_edges <- function(i) {
    key <- as.character(fits$fit_id[i])
    if (!exists(key, envir = edge_cache)) {
      e <- arrow::read_parquet(fits$loadings_file[i])
      e$key <- paste(pmin(e$node1, e$node2), pmax(e$node1, e$node2), sep = "|")
      assign(key, e, envir = edge_cache)
    }
    get(key, envir = edge_cache)
  }

  is_new <- fits$fit_id %in% new_fit_ids
  for (i in seq_len(nrow(fits))) {
    for (j in seq_len(nrow(fits))) {
      if (i >= j) next
      if (!is_new[i] && !is_new[j]) next
      a <- get_edges(i); b <- get_edges(j)
      m <- match(a$key, b$key)
      ok <- !is.na(m)
      if (sum(ok) < 2) next
      wa <- a$wto[ok]; wb <- b$wto[m[ok]]
      sig_a <- a$key[!is.na(a$padj) & a$padj < padj_cutoff]
      sig_b <- b$key[!is.na(b$padj) & b$padj < padj_cutoff]
      uni <- length(union(sig_a, sig_b))
      DBI::dbExecute(con,
        "INSERT INTO wto_fit_pairs (fit_a, fit_b, pearson, spearman, jaccard_sig, padj_cutoff, n_edges_common)
         VALUES (?, ?, ?, ?, ?, ?, ?)",
        params = list(fits$fit_id[i], fits$fit_id[j],
                      cor(wa, wb, method = "pearson"),
                      cor(wa, wb, method = "spearman"),
                      if (uni > 0) length(intersect(sig_a, sig_b)) / uni else NA_real_,
                      padj_cutoff, sum(ok)))
    }
  }
  invisible(NULL)
}

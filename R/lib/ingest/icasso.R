# ICASSO (Himberg, Hyvärinen & Esposito 2004, "Validating the independent
# components of neuroimaging time series via clustering and
# visualization," NeuroImage 22(3):1214-1222): cluster every recovered ICA
# component across ALL runs at one rank -- both seed-sweep (reinit-only)
# and bootstrap (bootstrap+reinit) fits together, the pooled randomization
# set the method calls for -- by |cosine| similarity, then score each
# cluster's robustness via Iq = mean intra-cluster similarity - mean
# similarity to everything outside the cluster.
#
# Run at ingest time (see R/lib/ingest/ingest_dataset.R), one full
# delete-and-recompute per (dataset_id, rank) whenever that rank's ICA fit
# set changes -- a cluster analysis isn't meaningfully "incremental" the
# way pairwise factor_pairs comparisons are (adding one more run changes
# where every existing component's cluster boundary falls, not just one
# new row).

#' Pure computation: given every ok ICA fit's ALREADY-LOADED loadings
#' matrix (genes x n.comp) at one rank, cluster all (fit, factor) pairs
#' and compute each cluster's Iq + centrotype.
#'
#' @param loadings_list named list (names = as.character(fit_id)) of
#'   genes x n.comp loading matrices, one per fit at this rank -- assumed
#'   to share an identical gene set/order (true by construction: every
#'   ICA fit for one dataset is computed from the same cached matrix).
#' @param rank expected component count -- clusters are cut to exactly
#'   this many (Himberg et al.'s convention: cut to the number of sources
#'   the model was told to look for).
#'
#' Returns list(clusters = data.frame(cluster_id, iq, n_members,
#' centrotype_fit_id, centrotype_factor_index), membership =
#' data.frame(cluster_id, fit_id, factor_index, intra_sim), hclust = the
#' real hclust object, comp_fit/comp_factor = parallel vectors identifying
#' the hclust object's leaves) or NULL if fewer than 2 fits are given.
compute_icasso_clusters_from_loadings <- function(loadings_list, rank) {
  if (length(loadings_list) < 2) return(NULL)

  fit_ids <- names(loadings_list)
  ref_genes <- rownames(loadings_list[[1]])

  cols <- list()
  comp_fit <- integer(0)
  comp_factor <- integer(0)
  for (fid in fit_ids) {
    L <- as.matrix(loadings_list[[fid]])[ref_genes, , drop = FALSE]
    nrm <- sqrt(colSums(L^2)); nrm[nrm == 0] <- 1
    cols[[fid]] <- sweep(L, 2, nrm, "/")
    comp_fit <- c(comp_fit, rep(as.integer(fid), ncol(cols[[fid]])))
    comp_factor <- c(comp_factor, seq_len(ncol(cols[[fid]])))
  }
  M <- do.call(cbind, cols)   # genes x N (N = total components across every fit)
  a_sim <- abs(crossprod(M))  # N x N |cosine| -- sign-robust, same reasoning as hungarian_match_abs()
  diag(a_sim) <- 1
  n <- ncol(a_sim)

  hc <- hclust(as.dist(1 - a_sim), method = "average")
  cl <- cutree(hc, k = rank)

  clusters <- list(); membership <- list()
  for (k in seq_len(rank)) {
    members <- which(cl == k)
    if (length(members) == 0) next
    outside <- setdiff(seq_len(n), members)

    if (length(members) > 1) {
      intra <- a_sim[members, members, drop = FALSE]
      diag(intra) <- NA
      intra_per_member <- rowMeans(intra, na.rm = TRUE)
      mean_intra <- mean(intra[upper.tri(intra)])
    } else {
      # A singleton "cluster" has no internal agreement to speak of --
      # Iq degenerates to -mean_inter, appropriately penalizing an
      # isolated, unreplicated component.
      intra_per_member <- 0
      mean_intra <- 0
    }
    mean_inter <- if (length(outside) > 0) mean(a_sim[members, outside, drop = FALSE]) else 0
    iq <- mean_intra - mean_inter
    centrotype_idx <- members[which.max(intra_per_member)]

    clusters[[length(clusters) + 1]] <- data.frame(
      cluster_id = k, iq = iq, n_members = length(members),
      centrotype_fit_id = comp_fit[centrotype_idx], centrotype_factor_index = comp_factor[centrotype_idx])
    membership[[length(membership) + 1]] <- data.frame(
      cluster_id = k, fit_id = comp_fit[members], factor_index = comp_factor[members],
      intra_sim = intra_per_member)
  }

  list(clusters = do.call(rbind, clusters), membership = do.call(rbind, membership),
       hclust = hc, comp_fit = comp_fit, comp_factor = comp_factor)
}

#' DB-querying wrapper: full delete-and-recompute for one (dataset_id,
#' rank)'s ICA fits. `method` is hardcoded to "ica" for now -- the
#' icasso_* tables' `method` column is forward-looking schema, not yet
#' meaningful for another method. Loads every ok fit's loadings_file
#' (seed-sweep AND bootstrap both -- see R/methods/ica.R's header),
#' computes clusters, and replaces whatever icasso_clusters/
#' icasso_membership/icasso_dendrograms rows already existed for this
#' (dataset_id, rank). No-op (returns NULL) if fewer than 2 ok fits exist
#' at this rank yet.
compute_icasso_clusters <- function(con, db_path, dataset_id, rank) {
  fits <- DBI::dbGetQuery(con,
    "SELECT fit_id, loadings_file FROM fits
     WHERE dataset_id = ? AND method = 'ica' AND rank = ? AND status = 'ok'
       AND loadings_file IS NOT NULL",
    params = list(dataset_id, rank))
  if (nrow(fits) < 2) return(invisible(NULL))

  loadings_list <- setNames(
    lapply(fits$loadings_file, function(p) as.matrix(readRDS(resolve_artifact(p, db_path)))),
    as.character(fits$fit_id))

  out <- compute_icasso_clusters_from_loadings(loadings_list, rank)
  if (is.null(out)) return(invisible(NULL))

  now <- as.character(Sys.time())

  DBI::dbExecute(con, "DELETE FROM icasso_clusters WHERE dataset_id = ? AND method = 'ica' AND rank = ?",
                 params = list(dataset_id, rank))
  DBI::dbExecute(con, "DELETE FROM icasso_membership WHERE dataset_id = ? AND method = 'ica' AND rank = ?",
                 params = list(dataset_id, rank))
  DBI::dbExecute(con, "DELETE FROM icasso_dendrograms WHERE dataset_id = ? AND method = 'ica' AND rank = ?",
                 params = list(dataset_id, rank))

  clusters <- out$clusters
  clusters$dataset_id <- dataset_id; clusters$method <- "ica"; clusters$rank <- rank; clusters$computed_at <- now
  DBI::dbWriteTable(con, "icasso_clusters", clusters, append = TRUE)

  membership <- out$membership
  membership$dataset_id <- dataset_id; membership$method <- "ica"; membership$rank <- rank
  DBI::dbWriteTable(con, "icasso_membership", membership, append = TRUE)

  art_dir <- artifacts_dir(db_path, dataset_id)
  dir.create(art_dir, recursive = TRUE, showWarnings = FALSE)
  fname <- sprintf("icasso_dendro_ica_rank%d.rds", rank)
  # Same small hclust-reconstructable shape as fits.wgcna_dendro_file (see
  # R/lib/ingest/ingest_dataset.R::compute_wgcna_dendro()) -- `labels`
  # encodes "fit_id:factor_index" per leaf so a rendered dendrogram can be
  # joined back to icasso_membership for cluster-colored display.
  tree <- list(merge = out$hclust$merge, height = out$hclust$height, order = out$hclust$order,
               labels = paste0(out$comp_fit, ":", out$comp_factor))
  saveRDS(tree, file.path(art_dir, fname))
  DBI::dbExecute(con,
    "INSERT INTO icasso_dendrograms (dataset_id, method, rank, dendro_file, computed_at) VALUES (?, ?, ?, ?, ?)",
    params = list(dataset_id, "ica", rank, file.path("stability_artifacts", dataset_id, fname), now))

  invisible(nrow(clusters))
}

# ICASSO (Himberg, Hyvärinen & Esposito 2004, "Validating the independent
# components of neuroimaging time series via clustering and visualization",
# NeuroImage 22(3):1214-1222): cluster every ICA component from every run at
# one rank -- seed sweep and bootstrap resamples together -- by |cosine|, and
# score each cluster with Iq = mean similarity within the cluster - mean
# similarity to everything outside it. Run per rank by
# compute_method_icasso() (R/targets/ingest.R).

#' Cluster the ICA components of every fit at one rank
#'
#' Average-linkage clustering on 1 - |cosine|, cut to `rank` clusters (one per
#' source the model was asked for), with each cluster's Iq and centrotype (the
#' member most similar to the rest of its cluster).
#'
#' @param loadings_list List of genes x n.comp loading matrices, one per fit
#'   at this rank, sharing one gene set; names identify the fits.
#' @param rank Number of clusters to cut to.
#' @return `list(clusters, membership, hclust, comp_fit, comp_factor)`:
#'   clusters (cluster_id, iq, n_members, centrotype_fit_id,
#'   centrotype_factor_index); membership (cluster_id, fit_id, factor_index,
#'   intra_sim); the `hclust` object; and the fit name and factor index of
#'   each leaf. `NULL` for fewer than 2 fits.
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


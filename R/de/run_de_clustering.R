# Cluster the DE results for each of the 8 time-course datasets, and
# check specifically whether there's clustering by TIME. Two analyses per
# (dataset, group) -- "group" being the disease/condition label baked
# into that dataset's cell-means design where one exists (see
# docs/methods.qmd ("Differential expression")), or the whole dataset when there isn't one:
#
#   1. Gene clustering by temporal SHAPE: significant genes (from that
#      group's `*_time_omnibus` contrast) pulled directly from the
#      cell-means fit's $coefficients (genes x group.time -- already
#      fitted, no recomputation), z-scored per gene, hclust()'d on
#      correlation distance -- groups genes into "early responder,"
#      "late responder," "monotonic," etc. modules.
#   2. Does DE-gene expression actually recover time structure? Samples
#      (within that group) clustered on the SAME significant genes' raw
#      per-sample expression, cut to the number of real timepoint levels,
#      cross-tabbed against the TRUE time labels via
#      mclust::adjustedRandIndex() (same function
#      app/R/comparison_helpers.R::cluster_and_ari() already uses
#      elsewhere in this repo) -- directly answers "is there clustering
#      by time."
#
# Needs R/de/de_bundle.rds (fit, expr, time, subject, group) -- see
# R/de/de_helpers.R::write_de_results(). Local, one-time run.

library(here)

TIME_DATASETS <- c("ANEMONES", "CORTICUS", "GSE110487", "GSE13904", "GSE273700",
                    "GSE54514", "GSE95233", "ROSE")
K_GENE_CLUSTERS <- 6
SIG_THRESHOLD <- 0.05

summary_rows <- list()

for (dataset_id in TIME_DATASETS) {
  out_dir <- here("results/de", dataset_id)
  bundle_path <- file.path(out_dir, "de_bundle.rds")
  if (!file.exists(bundle_path)) {
    message("[", dataset_id, "] no de_bundle.rds -- skipping (re-run R/de/", dataset_id, "_de.R first)")
    next
  }
  bundle <- readRDS(bundle_path)
  coefs <- bundle$fit$coefficients

  omnibus_files <- list.files(out_dir, pattern = "^topTable_.*omnibus.*\\.csv$", full.names = TRUE)
  if (length(omnibus_files) == 0) {
    message("[", dataset_id, "] no omnibus contrast files -- skipping")
    next
  }

  for (f in omnibus_files) {
    contrast_name <- sub("^topTable_", "", sub("\\.csv$", "", basename(f)))
    group_name <- if (identical(contrast_name, "time_omnibus")) NULL else sub("_time_omnibus$", "", contrast_name)
    tag <- paste0(dataset_id, "/", contrast_name)

    tt <- read.csv(f, stringsAsFactors = FALSE)
    sig_genes <- tt$gene[!is.na(tt$adj.P.Val) & tt$adj.P.Val < SIG_THRESHOLD]
    if (length(sig_genes) < 3) {
      message("[", tag, "] only ", length(sig_genes), " significant genes (< 3) -- skipping clustering")
      next
    }

    # Resolve this contrast's group back to the raw factor level (design
    # column names/sample masking both use the RAW group string, not the
    # make.names()-sanitized one the contrast filename carries) and this
    # group's sample mask.
    if (is.null(group_name)) {
      cell_cols <- colnames(coefs)
      sample_mask <- rep(TRUE, length(bundle$time))
    } else {
      group_levels_raw <- levels(factor(bundle$group))
      matched <- group_levels_raw[make.names(group_levels_raw) == group_name]
      if (length(matched) != 1) {
        message("[", tag, "] could not resolve group name back to a raw factor level -- skipping")
        next
      }
      prefix <- paste0(make.names(matched), ".")
      cell_cols <- colnames(coefs)[startsWith(colnames(coefs), prefix)]
      sample_mask <- as.character(bundle$group) == matched
    }
    sig_genes <- intersect(sig_genes, rownames(coefs))

    ## ---- 1. gene clustering by temporal shape ----
    cell_means <- coefs[sig_genes, cell_cols, drop = FALSE]
    z <- t(scale(t(cell_means)))
    z <- z[stats::complete.cases(z), , drop = FALSE]   # drop any zero-variance gene (flat across all cells)
    k <- min(K_GENE_CLUSTERS, nrow(z) - 1)
    if (k < 2) {
      message("[", tag, "] too few genes with non-zero variance for clustering -- skipping")
    } else {
      d <- stats::as.dist(1 - stats::cor(t(z)))
      hc <- stats::hclust(d, method = "average")
      cl <- stats::cutree(hc, k = k)
      write.csv(data.frame(gene = rownames(z), cluster = cl),
                file.path(out_dir, paste0("gene_clusters_", contrast_name, ".csv")), row.names = FALSE)

      grDevices::png(file.path(out_dir, paste0("gene_cluster_heatmap_", contrast_name, ".png")),
                     width = 900, height = 900)
      stats::heatmap(z[hc$order, , drop = FALSE], Rowv = NA, Colv = NA, scale = "none",
                     col = grDevices::colorRampPalette(c("steelblue", "white", "firebrick"))(50),
                     main = paste0(tag, " -- ", nrow(z), " sig. genes, ", k, " temporal-shape clusters"))
      grDevices::dev.off()
      message("[", tag, "] gene clustering: ", nrow(z), " genes -> ", k, " clusters")
    }

    ## ---- 2. does DE-gene expression recover time structure? ----
    expr_sub <- bundle$expr[sig_genes, sample_mask, drop = FALSE]
    time_sub <- droplevels(factor(bundle$time[sample_mask]))
    n_time_levels <- length(levels(time_sub))
    if (n_time_levels < 2 || ncol(expr_sub) < n_time_levels + 1) {
      message("[", tag, "] too few samples/timepoints for sample-clustering check -- skipping")
      next
    }
    d_samp <- stats::dist(t(expr_sub))
    hc_samp <- stats::hclust(d_samp, method = "average")
    cl_samp <- stats::cutree(hc_samp, k = n_time_levels)
    ari <- mclust::adjustedRandIndex(cl_samp, as.character(time_sub))
    message("[", tag, "] sample clustering vs. true time: ARI = ", round(ari, 3),
            " (", n_time_levels, " time levels, ", ncol(expr_sub), " samples)")

    pca <- stats::prcomp(t(expr_sub), scale. = FALSE)
    pct <- round(100 * pca$sdev[1:2]^2 / sum(pca$sdev^2), 1)
    grDevices::png(file.path(out_dir, paste0("sample_time_clustering_", contrast_name, ".png")),
                   width = 1200, height = 600)
    graphics::par(mfrow = c(1, 2))
    graphics::plot(pca$x[, 1], pca$x[, 2], col = as.integer(time_sub), pch = 19,
                   xlab = paste0("PC1 (", pct[1], "%)"), ylab = paste0("PC2 (", pct[2], "%)"),
                   main = paste0(tag, " -- colored by TRUE time"))
    graphics::legend("topright", legend = levels(time_sub), col = seq_along(levels(time_sub)), pch = 19)
    graphics::plot(pca$x[, 1], pca$x[, 2], col = cl_samp, pch = 19,
                   xlab = paste0("PC1 (", pct[1], "%)"), ylab = paste0("PC2 (", pct[2], "%)"),
                   main = paste0("colored by cluster (ARI = ", round(ari, 3), ")"))
    grDevices::dev.off()

    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      dataset_id = dataset_id, contrast = contrast_name, n_sig_genes = length(sig_genes),
      n_time_levels = n_time_levels, n_samples = ncol(expr_sub), ari = ari)
  }
}

if (length(summary_rows) > 0) {
  summary_df <- do.call(rbind, summary_rows)
  write.csv(summary_df, here("results/de/sample_time_clustering_summary.csv"), row.names = FALSE)
  message("\nSummary written to results/de/sample_time_clustering_summary.csv")
  print(summary_df)
}
message("done.")

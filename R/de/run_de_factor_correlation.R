# Correlate each dataset's DE results with its factor-analysis latent
# factors (PCA/NMF/CoGAPS/sPCA/ICA), two ways -- see R/de/de_factor_helpers.R
# and docs/methods.qmd ("Differential expression") for the full rationale:
#   (a) overlap_test(): Fisher's exact test between a factor's gene set
#       and a DE contrast's significant genes -- run for EVERY contrast
#       (including the time-course datasets' `*_time_omnibus` ones, which
#       have no sign but don't need one for a binary overlap test).
#   (b) factor_gsea_test(): rank all genes by a DE contrast's signed `t`
#       stat, test whether each factor's gene set is enriched at either
#       end -- PAIRWISE contrasts only (needs a sign).
#
# One fixed rank per method (10, or K=10 for spca -- see
# R/de/de_factor_helpers.R::pick_fit_at_rank()) -- deliberately sidesteps
# model-selection for this first pass. Local, one-time run -- reads
# results/stability.sqlite READ-ONLY (never writes to it).

library(here)
library(yaml)
library(arrow)
library(DBI)
library(RSQLite)
source(here("R/lib/data_paths.R"))
source(here("R/db/connect.R"))                 # open_stability_db(), resolve_artifact()
source(here("R/lib/ingest/symbol_mapping.R"))  # build_symbol_map(), remap_to_symbol()
source(here("R/lib/ingest/redundancy.R"))      # representative_fit_ids()
source(here("R/de/de_helpers.R"))
source(here("R/de/de_factor_helpers.R"))
source(here("R/de/de_dataset_manifest.R"))

FACTOR_METHODS <- c("pca", "nmf", "cogaps", "spca", "ica")
SIG_THRESHOLD <- 0.05
TARGET_RANK <- 10

db_path <- here(yaml::read_yaml(here("config/pipeline.yml"))$db_path %||% "results/targets/stability.sqlite")
con <- open_stability_db(db_path, schema_dir = here("R/db"))
on.exit(DBI::dbDisconnect(con), add = TRUE)

for (dataset_id in names(DE_DATASETS)) {
  spec <- DE_DATASETS[[dataset_id]]
  out_dir <- here("results/de", dataset_id)
  contrast_files <- list.files(out_dir, pattern = "^topTable_.*\\.csv$", full.names = TRUE)
  if (length(contrast_files) == 0) {
    message("[", dataset_id, "] no DE results -- skipping")
    next
  }

  avail_methods <- DBI::dbGetQuery(con,
    "SELECT DISTINCT method FROM fits WHERE dataset_id = ? AND status = 'ok'",
    params = list(dataset_id))$method
  avail_methods <- intersect(FACTOR_METHODS, avail_methods)
  if (length(avail_methods) == 0) {
    message("[", dataset_id, "] no ingested factor-analysis fits -- skipping")
    next
  }

  ds_meta <- read_dataset_yaml(here(spec$config))$dataset
  feature_meta <- arrow::read_parquet(ds_meta$feature_metadata_path)

  overlap_rows <- list()
  gsea_rows <- list()

  for (method in avail_methods) {
    f <- pick_fit_at_rank(con, dataset_id, method, TARGET_RANK)
    if (is.null(f)) {
      message("[", dataset_id, "/", method, "] no representative fit found -- skipping")
      next
    }
    loadings <- as.matrix(readRDS(resolve_artifact(f$loadings_file, db_path)))
    loadings <- reconcile_loadings_to_de_ids(loadings, ds_meta, feature_meta, platform = spec$platform)
    message("[", dataset_id, "/", method, "] rank ", f$rank, " (fit_id ", f$fit_id, "), ",
            ncol(loadings), " factor(s), ", nrow(loadings), " genes post-reconciliation")

    factor_sets <- setNames(
      lapply(seq_len(ncol(loadings)), function(i) factor_gene_set(loadings[, i], method, n = 100)),
      colnames(loadings))

    for (cf in contrast_files) {
      contrast_name <- sub("^topTable_", "", sub("\\.csv$", "", basename(cf)))
      tt <- read.csv(cf, stringsAsFactors = FALSE)
      de_sig <- tt$gene[!is.na(tt$adj.P.Val) & tt$adj.P.Val < SIG_THRESHOLD]
      universe <- intersect(rownames(loadings), tt$gene)
      if (length(universe) < 20) next   # too little shared identifier overlap to trust a test

      # (a) overlap test -- every contrast, every factor
      for (fi in colnames(loadings)) {
        ov <- overlap_test(factor_sets[[fi]], de_sig, universe)
        overlap_rows[[length(overlap_rows) + 1]] <- cbind(
          dataset_id = dataset_id, method = method, rank = f$rank, factor = fi,
          contrast = contrast_name, ov)
      }

      # (b) GSEA-style test -- pairwise contrasts only (needs a signed stat)
      if ("t" %in% names(tt)) {
        rank_vector <- sort(setNames(tt$t, tt$gene)[universe], decreasing = TRUE)
        gs <- factor_gsea_test(rank_vector, factor_sets)
        if (!is.null(gs) && nrow(gs) > 0) {
          gs$leadingEdge <- vapply(gs$leadingEdge, paste, character(1), collapse = ",")
          gsea_rows[[length(gsea_rows) + 1]] <- cbind(
            dataset_id = dataset_id, method = method, rank = f$rank,
            contrast = contrast_name, as.data.frame(gs))
        }
      }
    }
  }

  if (length(overlap_rows) > 0) {
    df <- do.call(rbind, overlap_rows)
    df$adj_p_value <- p.adjust(df$p_value, method = "BH")
    df <- df[order(df$p_value), ]
    write.csv(df, file.path(out_dir, "factor_overlap.csv"), row.names = FALSE)
    message("[", dataset_id, "] factor_overlap.csv: ", nrow(df), " (factor x contrast) tests")
  }
  if (length(gsea_rows) > 0) {
    df <- do.call(rbind, gsea_rows)
    df <- df[order(df$padj), ]
    write.csv(df, file.path(out_dir, "factor_gsea.csv"), row.names = FALSE)
    message("[", dataset_id, "] factor_gsea.csv: ", nrow(df), " factor-hit(s) across contrasts")
  }
}

message("done.")

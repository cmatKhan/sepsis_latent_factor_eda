# fgsea enrichment for every DE contrast, mirroring R/ingest_jobs/
# fgsea_job.R's conventions (same 6 msigdbr collections, same fgsea()
# parameters, same canonical Ensembl identifier space) -- see
# R/de/README.md and R/de/de_fgsea_helpers.R for the full rationale.
#
# Only PAIRWISE contrasts get fgsea (files with "omnibus" in the name are
# skipped) -- GSEA needs a signed ranking statistic (up in A, down in B),
# and the omnibus F-test has no sign. Each pairwise contrast's `t` column
# is the rank vector.
#
# Writes results/de/<DATASET>/fgsea_<contrast>.csv (all 6 sources stacked,
# `source` column) for every dataset that has DE results (13 -- excludes
# SHIP-TREND, which has none). Local, one-time run -- no cluster needed,
# same as the rest of R/de/.

library(here)
library(yaml)
library(arrow)
source(here("R/lib/ingest/symbol_mapping.R"))
source(here("R/de/de_helpers.R"))
source(here("R/de/de_fgsea_helpers.R"))
source(here("R/de/de_dataset_manifest.R"))
DATASETS <- DE_DATASETS

message("Fetching msigdbr pathways (once for this whole run)...")
pathways_by_source <- fetch_de_pathways()
message("  ", paste(sprintf("%s: %d gene sets", names(pathways_by_source),
                             lengths(pathways_by_source)), collapse = "; "))

for (dataset_id in names(DATASETS)) {
  spec <- DATASETS[[dataset_id]]
  out_dir <- here("results/de", dataset_id)
  contrast_files <- list.files(out_dir, pattern = "^topTable_.*\\.csv$", full.names = TRUE)
  contrast_files <- contrast_files[!grepl("omnibus", contrast_files)]
  if (length(contrast_files) == 0) {
    message("[", dataset_id, "] no pairwise contrast files found -- skipping")
    next
  }

  ds_meta <- yaml::read_yaml(here(spec$config))$dataset
  feature_meta <- arrow::read_parquet(ds_meta$feature_metadata_path)
  ens_map <- build_de_ensembl_map(ds_meta, feature_meta, platform = spec$platform)
  if (is.null(ens_map)) {
    message("[", dataset_id, "] no usable Ensembl map -- skipping (check ensembl_col/symbol_col in config)")
    next
  }

  for (f in contrast_files) {
    contrast_name <- sub("^topTable_", "", sub("\\.csv$", "", basename(f)))
    tt <- read.csv(f, stringsAsFactors = FALSE)
    if (!("t" %in% names(tt))) {
      message("[", dataset_id, "/", contrast_name, "] no `t` column (not a pairwise contrast?) -- skipping")
      next
    }
    stats_raw <- setNames(tt$t, tt$gene)
    stats_raw <- stats_raw[!is.na(stats_raw)]
    stats_ens <- remap_named_vector_to_ensembl(stats_raw, ens_map)
    pct_mapped <- round(100 * length(stats_ens) / length(stats_raw), 1)
    message("[", dataset_id, "/", contrast_name, "] ", length(stats_raw), " genes -> ",
            length(stats_ens), " Ensembl-mapped (", pct_mapped, "%)")

    all_res <- list()
    for (src in names(pathways_by_source)) {
      res <- run_de_fgsea(stats_ens, pathways_by_source[[src]])
      if (is.null(res) || nrow(res) == 0) next
      res <- res[!is.na(res$padj) & res$padj < 0.05, ]
      if (nrow(res) == 0) next
      res$leadingEdge <- vapply(res$leadingEdge, paste, character(1), collapse = ",")
      res$source <- src
      all_res[[src]] <- as.data.frame(res)
    }

    out_file <- file.path(out_dir, paste0("fgsea_", contrast_name, ".csv"))
    if (length(all_res) == 0) {
      message("  no significant pathways across any source")
      write.csv(data.frame(pathway = character(0), source = character(0)), out_file, row.names = FALSE)
    } else {
      combined <- do.call(rbind, all_res)
      combined <- combined[order(combined$padj), ]
      write.csv(combined, out_file, row.names = FALSE)
      message("  ", nrow(combined), " significant pathway-hit(s) written to ", basename(out_file))
    }
  }
}

message("done.")

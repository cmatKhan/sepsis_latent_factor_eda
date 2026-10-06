# Differential expression for hfgp-500fg (GSE134080) -- CROSS-SECTIONAL
# (confirmed: 100 samples, 100 distinct geo_accession, no repeat/subject
# column). See docs/methods.qmd ("Differential expression").
#
# ** This cohort has NO sepsis or disease-relevant covariate ** -- same
# situation as dilgom: a healthy population functional-genomics reference
# cohort (real metadata: geo_accession, age, sex, accession -- nothing
# else). `sex` is used below purely for pipeline completeness, NOT as a
# sepsis-relevant result -- see docs/methods.qmd ("Differential expression") and dilgom_de.R's matching
# caveat.
#
# Platform: RNA-seq -- R/preprocessing/hfgp-500fg_preprocessing.R confirms
# `as.integer(value)` -> DESeq2 VST (same pattern as EARLI/GSE110487/
# GSE273700/ROSE), raw counts. limma+voom.

library(here)
library(yaml)
library(arrow)
source(here("R/lib/data_paths.R"))
source(here("R/lib/matrices.R"))
source(here("R/lib/ingest/symbol_mapping.R"))
source(here("R/de/de_helpers.R"))

ds_meta <- read_dataset_yaml(here("config/hfgp-500fg_config.yml"))$dataset

message("[hfgp-500fg] loading raw count matrix...")
mat <- pivot_expression_long(ds_meta)
storage.mode(mat) <- "integer"

sample_meta <- arrow::read_parquet(ds_meta$sample_metadata_path)
feature_meta <- arrow::read_parquet(ds_meta$feature_metadata_path)
symbol_map <- build_symbol_map(list(dataset = ds_meta), fm = feature_meta)
message("[hfgp-500fg] ", nrow(mat), " genes x ", ncol(mat), " samples (raw, pre-filterByExpr)")

de <- run_de_case_control(
  mat, sample_meta,
  sample_id_col = ds_meta$sample_id_col,
  group_col     = "sex",
  platform      = "rnaseq"
)

out_dir <- here("results/de/hfgp-500fg")
write_de_results(de, out_dir, symbol_map = symbol_map)
message("[hfgp-500fg] done -- results in ", out_dir, " (NOTE: no sepsis-relevant covariate exists for this cohort -- see file header)")

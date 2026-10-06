# Differential expression for EARLI -- CROSS-SECTIONAL (confirmed
# directly: 189 samples, 189 distinct sample_id, no repeat/subject column
# at all), so this is a simple case/control comparison, not a time-course
# model. See docs/methods.qmd ("Differential expression").
#
# Real metadata: sample_id, geo_accession, age, sex, lca_label,
# imputed_age, age_scaled. `lca_label` (Hypo/Hyper) is a latent-class
# immune endotype label -- the closest thing this cohort has to a
# sepsis-relevant grouping (there is no explicit sepsis/healthy label),
# so it's the contrast used here.
#
# Platform: RNA-seq -- R/preprocessing/earli_preprocessing.R confirms
# `as.integer(value)` -> DESeqDataSetFromMatrix, raw counts. limma+voom.

library(here)
library(yaml)
library(arrow)
source(here("R/lib/data_paths.R"))
source(here("R/lib/matrices.R"))
source(here("R/lib/ingest/symbol_mapping.R"))
source(here("R/de/de_helpers.R"))

ds_meta <- read_dataset_yaml(here("config/EARLI_config.yml"))$dataset

message("[EARLI] loading raw count matrix...")
mat <- pivot_expression_long(ds_meta)
storage.mode(mat) <- "integer"

sample_meta <- arrow::read_parquet(ds_meta$sample_metadata_path)
feature_meta <- arrow::read_parquet(ds_meta$feature_metadata_path)
symbol_map <- build_symbol_map(list(dataset = ds_meta), fm = feature_meta)
message("[EARLI] ", nrow(mat), " genes x ", ncol(mat), " samples (raw, pre-filterByExpr)")

de <- run_de_case_control(
  mat, sample_meta,
  sample_id_col = ds_meta$sample_id_col,
  group_col     = "lca_label",
  platform      = "rnaseq"
)

out_dir <- here("results/de/EARLI")
write_de_results(de, out_dir, symbol_map = symbol_map)
message("[EARLI] done -- results in ", out_dir)

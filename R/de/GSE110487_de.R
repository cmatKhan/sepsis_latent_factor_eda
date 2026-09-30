# Differential expression over time for GSE110487, RNA-seq, 2-level
# timecourse (T1/T2). See R/de/README.md.
#
# Real metadata (62 samples, 31 patients, balanced): sample_id,
# patient_id, response_to_treatment, timepoint. No healthy/other-condition
# comparator arm exists in this cohort (all patients are the same sepsis
# population), so no `group_col` is used -- plain `~ 0 + timepoint`,
# blocked on patient. `response_to_treatment` is left out of the model
# deliberately, same reasoning as CORTICUS's `treatment` -- see
# R/de/README.md.
#
# Platform: RNA-seq -- R/preprocessing/gse110487_preprocessing.R confirms
# `as.integer(value)` -> DESeqDataSetFromMatrix(countData=...), and the
# raw expression.parquet's `value` column is itself integer with no
# negatives -- raw counts, read directly here (bypassing the
# preprocessing_script's own VST step, which is for the factor-analysis
# pipeline, not for voom). limma+voom, not plain limma.

library(here)
library(yaml)
library(arrow)
source(here("R/lib/matrices.R"))
source(here("R/lib/ingest/symbol_mapping.R"))
source(here("R/de/de_helpers.R"))

ds_meta <- yaml::read_yaml(here("config/GSE110487_config.yml"))$dataset

message("[GSE110487] loading raw count matrix...")
mat <- pivot_expression_long(ds_meta)
storage.mode(mat) <- "integer"

sample_meta <- arrow::read_parquet(ds_meta$sample_metadata_path)
feature_meta <- arrow::read_parquet(ds_meta$feature_metadata_path)
symbol_map <- build_symbol_map(list(dataset = ds_meta), fm = feature_meta)
message("[GSE110487] ", nrow(mat), " genes x ", ncol(mat), " samples (raw, pre-filterByExpr)")

de <- run_de_timecourse(
  mat, sample_meta,
  sample_id_col = ds_meta$sample_id_col,
  subject_col   = "patient_id",
  time_col      = "timepoint",
  group_col     = NULL,
  platform      = "rnaseq"
)

out_dir <- here("results/de/GSE110487")
write_de_results(de, out_dir, symbol_map = symbol_map)
message("[GSE110487] done -- results in ", out_dir)

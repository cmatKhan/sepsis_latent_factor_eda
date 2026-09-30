# Differential expression over time for GSE273700, RNA-seq, 2-level
# timecourse (Day1/Day8). See R/de/README.md.
#
# Real metadata (104 samples, 52 patients, balanced): sample_id,
# author_id, patient_id, thrombocytopenia, timepoint. No healthy/other-
# condition comparator arm -- no `group_col` used (plain `~ 0 + timepoint`,
# blocked on patient). `thrombocytopenia` left out of the model
# deliberately, same reasoning as CORTICUS/GSE110487 -- see R/de/README.md.
#
# Platform: RNA-seq -- R/preprocessing/gse273700_preprocessing.R confirms
# `as.integer(value)` -> DESeqDataSetFromMatrix, raw expression.parquet's
# `value` is integer, no negatives. limma+voom.

library(here)
library(yaml)
library(arrow)
source(here("R/lib/matrices.R"))
source(here("R/lib/ingest/symbol_mapping.R"))
source(here("R/de/de_helpers.R"))

ds_meta <- yaml::read_yaml(here("config/GSE273700_config.yml"))$dataset

message("[GSE273700] loading raw count matrix...")
mat <- pivot_expression_long(ds_meta)
storage.mode(mat) <- "integer"

sample_meta <- arrow::read_parquet(ds_meta$sample_metadata_path)
feature_meta <- arrow::read_parquet(ds_meta$feature_metadata_path)
symbol_map <- build_symbol_map(list(dataset = ds_meta), fm = feature_meta)
message("[GSE273700] ", nrow(mat), " genes x ", ncol(mat), " samples (raw, pre-filterByExpr)")

de <- run_de_timecourse(
  mat, sample_meta,
  sample_id_col = ds_meta$sample_id_col,
  subject_col   = "patient_id",
  time_col      = "timepoint",
  group_col     = NULL,
  platform      = "rnaseq"
)

out_dir <- here("results/de/GSE273700")
write_de_results(de, out_dir, symbol_map = symbol_map)
message("[GSE273700] done -- results in ", out_dir)

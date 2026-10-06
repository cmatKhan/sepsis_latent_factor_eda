# Differential expression over time for CORTICUS (GSE106878), array
# platform, 2-level timecourse (Pre / Post(24h)). See docs/methods.qmd ("Differential expression").
#
# Real metadata (94 samples, 47 patients, all balanced Pre+Post24h pairs):
# sample_id, patient_id, sex, age, responder, pre_treatment,
# survive_28_days, timepoint, treatment. Every patient here is a septic
# shock patient enrolled in this steroid RCT (hydrocortisone vs. placebo)
# -- there is no healthy/other-condition comparator arm to include, unlike
# ANEMONES/GSE13904/GSE54514, so no `group_col` is used (plain `~ 0 +
# timepoint`, blocked on patient). `treatment` is left out of the model
# deliberately -- see docs/methods.qmd ("Differential expression") for why this first pass
# keeps every dataset's model to the minimum the data structurally
# requires rather than chasing every available covariate.
#
# Platform: array (log-scale microarray intensity, no DESeq2/count logic
# in R/preprocessing/corticus_preprocessing.R) -- plain limma, no voom.

library(here)
library(yaml)
library(arrow)
source(here("R/lib/data_paths.R"))
source(here("R/lib/matrices.R"))
source(here("R/de/de_helpers.R"))

ds_meta <- read_dataset_yaml(here("config/CORTICUS_config.yml"))$dataset

message("[CORTICUS] loading raw expression matrix...")
mat <- pivot_expression_long(ds_meta)

sample_meta <- arrow::read_parquet(ds_meta$sample_metadata_path)
feature_meta <- arrow::read_parquet(ds_meta$feature_metadata_path)

message("[CORTICUS] collapsing probes to gene symbol...")
mat <- collapse_to_symbol(mat, feature_meta, ds_meta$feature_id_col, ds_meta$symbol_col)
mat <- filter_low_variance(mat)
message("[CORTICUS] ", nrow(mat), " genes x ", ncol(mat), " samples after collapse + low-variance filter")

de <- run_de_timecourse(
  mat, sample_meta,
  sample_id_col = ds_meta$sample_id_col,
  subject_col   = "patient_id",
  time_col      = "timepoint",
  group_col     = NULL,
  platform      = "array"
)

out_dir <- here("results/de/CORTICUS")
write_de_results(de, out_dir)
message("[CORTICUS] done -- results in ", out_dir)

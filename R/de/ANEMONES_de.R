# Differential expression over time for ANEMONES (GSE236713), array
# platform, 4-level timecourse (Day 1/2/5/Discharge). See docs/methods.qmd ("Differential expression")
# for the shared design rationale.
#
# Real metadata (verified directly against sample_metadata.parquet, 447
# samples): sample_id, patient_id, timepoint, outcome, disease,
# disease_type, sex. `disease` (Sepsis 324 / SIRS 93 / Control 30) is the
# condition included in the model per the user's direction -- Control
# patients are sampled ONLY at Day 1 (0/30 have a second timepoint), so
# they contribute to the fit's residual/noise estimate but get no
# reported time-omnibus contrast (run_de_timecourse() drops any group
# below `min_pairs` automatically -- no special-casing needed here).
# Sepsis (106/126 patients paired) and SIRS (33/38 paired) both get a
# real within-group time-omnibus F-test.
#
# Platform: array (log-scale, already-normalized microarray intensities;
# confirmed no DESeq2/count logic in R/preprocessing/anemones_preprocessing.R
# and the raw expression values include negatives) -- plain limma, no voom.

library(here)
library(yaml)
library(arrow)
source(here("R/lib/data_paths.R"))
source(here("R/lib/matrices.R"))       # pivot_expression_long()
source(here("R/de/de_helpers.R"))

ds_meta <- read_dataset_yaml(here("config/ANEMONES_config.yml"))$dataset

message("[ANEMONES] loading raw expression matrix...")
mat <- pivot_expression_long(ds_meta)

sample_meta <- arrow::read_parquet(ds_meta$sample_metadata_path)
feature_meta <- arrow::read_parquet(ds_meta$feature_metadata_path)

message("[ANEMONES] collapsing probes to gene symbol...")
mat <- collapse_to_symbol(mat, feature_meta, ds_meta$feature_id_col, ds_meta$symbol_col)
mat <- filter_low_variance(mat)
message("[ANEMONES] ", nrow(mat), " genes x ", ncol(mat), " samples after collapse + low-variance filter")

de <- run_de_timecourse(
  mat, sample_meta,
  sample_id_col = ds_meta$sample_id_col,
  subject_col   = "patient_id",
  time_col      = "timepoint",
  group_col     = "disease",
  platform      = "array"
)

out_dir <- here("results/de/ANEMONES")
write_de_results(de, out_dir)
message("[ANEMONES] done -- results in ", out_dir)

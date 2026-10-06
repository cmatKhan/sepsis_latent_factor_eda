# Differential expression over time for GSE95233, array platform, 3-level
# timecourse (D01/D02/D03). See docs/methods.qmd ("Differential expression").
#
# Real metadata (124 samples, 71 patients): sample_id, patient_id, age,
# sex, survival, timepoint. 22 samples are the study's healthy-control arm
# -- NOT true NA, but the literal STRING "NA" in both `timepoint` and
# `survival` (confirmed directly: every timepoint=="NA" row also has
# survival=="NA"). These 22 samples are excluded from the time model here
# -- not a subsetting choice like ANEMONES/GSE13904/GSE54514's healthy/
# other-condition arms (which DO have real timepoints and stay in via a
# `group_col`), but a structural necessity: "NA" isn't a real timepoint
# value, so there is no time axis for these samples to contribute to at
# all. Once they're dropped, every remaining sample is from the sepsis
# arm and there is no other disease/condition column left to include as
# `group_col` for this dataset specifically -- `survival` is an outcome
# of the sepsis course, not a baseline condition label like ANEMONES's
# `disease`. Paired-subject count among the remaining sepsis patients:
# 51/71 have >=2 timepoints.
#
# Platform: array (no DESeq2/count logic in
# R/preprocessing/gse95233_preprocessing.R) -- plain limma, no voom.

library(here)
library(yaml)
library(arrow)
source(here("R/lib/data_paths.R"))
source(here("R/lib/matrices.R"))
source(here("R/de/de_helpers.R"))

ds_meta <- read_dataset_yaml(here("config/GSE95233_config.yml"))$dataset

message("[GSE95233] loading raw expression matrix...")
mat <- pivot_expression_long(ds_meta)

sample_meta <- arrow::read_parquet(ds_meta$sample_metadata_path)
feature_meta <- arrow::read_parquet(ds_meta$feature_metadata_path)

n_before <- nrow(sample_meta)
sample_meta <- sample_meta[sample_meta$timepoint != "NA", ]
message("[GSE95233] dropped ", n_before - nrow(sample_meta),
        " healthy-control samples with no real timepoint")
mat <- mat[, sample_meta[[ds_meta$sample_id_col]], drop = FALSE]

message("[GSE95233] collapsing probes to gene symbol...")
mat <- collapse_to_symbol(mat, feature_meta, ds_meta$feature_id_col, ds_meta$symbol_col)
mat <- filter_low_variance(mat)
message("[GSE95233] ", nrow(mat), " genes x ", ncol(mat), " samples after collapse + low-variance filter")

de <- run_de_timecourse(
  mat, sample_meta,
  sample_id_col = ds_meta$sample_id_col,
  subject_col   = "patient_id",
  time_col      = "timepoint",
  group_col     = NULL,
  platform      = "array"
)

out_dir <- here("results/de/GSE95233")
write_de_results(de, out_dir)
message("[GSE95233] done -- results in ", out_dir)

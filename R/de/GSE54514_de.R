# Differential expression over time for GSE54514, array platform, up to
# 5-level timecourse (day 1-5). See R/de/README.md.
#
# Real metadata (163 samples, 54 subjects, subject col `group_id` not
# `patient_id`): sample_id, accession, disease_status, group_day, group_id,
# sex, age, apache_ii, neutrophil_proportion, site_of_infection, day.
# `disease_status` (healthy / sepsis survivor / sepsis nonsurvivor) is the
# condition included in the model per the user's direction. The healthy
# arm is sampled ONLY at day1/day5 (not the sepsis arms' full 1-5 range) --
# the cell-means design handles this without confounding the sepsis arms'
# 5-day trajectory with the healthy arm's 2-day one. Paired-subject counts
# per group: healthy 18/18, sepsis survivor 24/26, sepsis nonsurvivor
# 7/10 -- all three clear run_de_timecourse()'s default min_pairs = 3, so
# all three get a real within-group time-omnibus contrast.
#
# Platform: array (no DESeq2/count logic in
# R/preprocessing/gse54514_preprocessing.R) -- plain limma, no voom.

library(here)
library(yaml)
library(arrow)
source(here("R/lib/matrices.R"))
source(here("R/de/de_helpers.R"))

ds_meta <- yaml::read_yaml(here("config/GSE54514_config.yml"))$dataset

message("[GSE54514] loading raw expression matrix...")
mat <- pivot_expression_long(ds_meta)

sample_meta <- arrow::read_parquet(ds_meta$sample_metadata_path)
feature_meta <- arrow::read_parquet(ds_meta$feature_metadata_path)

message("[GSE54514] collapsing probes to gene symbol...")
mat <- collapse_to_symbol(mat, feature_meta, ds_meta$feature_id_col, ds_meta$symbol_col)
mat <- filter_low_variance(mat)
message("[GSE54514] ", nrow(mat), " genes x ", ncol(mat), " samples after collapse + low-variance filter")

de <- run_de_timecourse(
  mat, sample_meta,
  sample_id_col = ds_meta$sample_id_col,
  subject_col   = "group_id",
  time_col      = "day",
  group_col     = "disease_status",
  platform      = "array"
)

out_dir <- here("results/de/GSE54514")
write_de_results(de, out_dir)
message("[GSE54514] done -- results in ", out_dir)

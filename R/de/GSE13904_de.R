# Differential expression over time for GSE13904, array platform, 2-level
# timecourse (day1/day3). See R/de/README.md.
#
# Real metadata (227 samples, 186 patients): sample_id, patient_id,
# timepoint, clinical_status (Control 18 / SIRS 27 / SIRS resolved 24 /
# Sepsis 52 / Septic Shock 106). `clinical_status` is the condition
# included in the model per the user's direction. Most patients are
# single-timepoint (145/186) -- per-group paired-subject counts: Control 0,
# SIRS resolved 0 (sampled ONLY at day3, by definition of "resolved"),
# SIRS 2, Sepsis 3, Septic Shock 19. run_de_timecourse()'s default
# min_pairs = 3 means only Septic Shock (and borderline Sepsis) get a real
# within-group time-omnibus contrast; the other groups still contribute
# their samples to the fit, they just get no reported time contrast --
# no special-casing needed, this falls out of the shared engine directly.
#
# Platform: array (no DESeq2/count logic in
# R/preprocessing/gse13904_preprocessing.R) -- plain limma, no voom.

library(here)
library(yaml)
library(arrow)
source(here("R/lib/matrices.R"))
source(here("R/de/de_helpers.R"))

ds_meta <- yaml::read_yaml(here("config/GSE13904_config.yml"))$dataset

message("[GSE13904] loading raw expression matrix...")
mat <- pivot_expression_long(ds_meta)

sample_meta <- arrow::read_parquet(ds_meta$sample_metadata_path)
feature_meta <- arrow::read_parquet(ds_meta$feature_metadata_path)

message("[GSE13904] collapsing probes to gene symbol...")
mat <- collapse_to_symbol(mat, feature_meta, ds_meta$feature_id_col, ds_meta$symbol_col)
mat <- filter_low_variance(mat)
message("[GSE13904] ", nrow(mat), " genes x ", ncol(mat), " samples after collapse + low-variance filter")

de <- run_de_timecourse(
  mat, sample_meta,
  sample_id_col = ds_meta$sample_id_col,
  subject_col   = "patient_id",
  time_col      = "timepoint",
  group_col     = "clinical_status",
  platform      = "array"
)

out_dir <- here("results/de/GSE13904")
write_de_results(de, out_dir)
message("[GSE13904] done -- results in ", out_dir)

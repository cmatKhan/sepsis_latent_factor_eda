# Differential expression over time for ROSE, RNA-seq, 2-level timecourse
# (Day 0/Day 2). See docs/methods.qmd ("Differential expression").
#
# Real metadata (221 samples, 128 subjects, subject col `subject_number`):
# a large block of columns beyond the usual clinical ones, including
# pseudoalignment QC and CIBERSORT-style cell-composition fractions.
# `class` (Hyperinflammatory/Hypoinflammatory -- this cohort's own
# endotype label, the closest analogue here to ANEMONES's `disease`) is
# the condition included in the model per the user's direction. Both
# classes are well-paired (38/62 and 55/66 subjects with both timepoints).
# Cell-composition fractions (neutrophils, monocytes, etc.) are NOT added
# as covariates -- real compositional shifts over a sepsis time-course are
# arguably part of the biology being asked about here, not solely a
# confound to adjust away, and adding them would be exactly the kind of
# "as many covariates as possible" completeness the user said not to
# chase for this first pass. See docs/methods.qmd ("Differential expression").
#
# Platform: RNA-seq -- R/preprocessing/rose_preprocessing.R confirms
# `as.integer(value)` -> DESeq2 (explicit estimateSizeFactors/
# estimateDispersions, not just vst()), raw data/rose_counts.parquet's
# `value` is integer, no negatives. limma+voom.

library(here)
library(yaml)
library(arrow)
source(here("R/lib/data_paths.R"))
source(here("R/lib/matrices.R"))
source(here("R/lib/ingest/symbol_mapping.R"))
source(here("R/de/de_helpers.R"))

ds_meta <- read_dataset_yaml(here("config/ROSE_config.yml"))$dataset

message("[ROSE] loading raw count matrix...")
mat <- pivot_expression_long(ds_meta)
storage.mode(mat) <- "integer"

sample_meta <- arrow::read_parquet(ds_meta$sample_metadata_path)
feature_meta <- arrow::read_parquet(ds_meta$feature_metadata_path)
symbol_map <- build_symbol_map(list(dataset = ds_meta), fm = feature_meta)
message("[ROSE] ", nrow(mat), " genes x ", ncol(mat), " samples (raw, pre-filterByExpr)")

de <- run_de_timecourse(
  mat, sample_meta,
  sample_id_col = ds_meta$sample_id_col,
  subject_col   = "subject_number",
  time_col      = "timepoint",
  group_col     = "class",
  platform      = "rnaseq"
)

out_dir <- here("results/de/ROSE")
write_de_results(de, out_dir, symbol_map = symbol_map)
message("[ROSE] done -- results in ", out_dir)

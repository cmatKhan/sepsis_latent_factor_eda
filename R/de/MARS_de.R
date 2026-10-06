# Differential expression for MARS (GSE65682) -- CROSS-SECTIONAL (confirmed
# directly: 802 samples, 802 distinct sample_id; `title_identifier1/2` are
# NOT subject ids -- checked directly, rows sharing a value have different
# ages/sexes, so it's a site/batch code, not a repeat-subject marker), so
# this is a simple case/control comparison. See docs/methods.qmd ("Differential expression"). Note
# `preprocessing_script: R/preprocessing/mars_preprocessing.R` in the
# config does not exist on disk -- irrelevant here since this script reads
# `expression_path` directly via pivot_expression_long(), the same as
# every other R/de/*.R script.
#
# Real metadata (rich): sample_id, title_identifier1/2, healthy_control,
# sex, age, pneumonia_diagnosis, thrombocytopenia, endotype_cohort,
# endotype_class (Mars1-4, this cohort's own namesake sepsis endotypes),
# mortality_28d, time_to_event_28d, icu_acquired_infection(_paired),
# diabetes_mellitus, abdominal_sepsis_or_control. `healthy_control`
# (FALSE 760 / TRUE 42) is used as the primary contrast -- the most direct
# "relevant to sepsis" signal available. `endotype_class`/`mortality_28d`
# are natural alternative contrasts for a deeper follow-up pass, not run
# here (see docs/methods.qmd ("Differential expression") scope).
#
# Platform: array (log2-scale intensity, range 0.68-13.5, no DESeq2/count
# logic anywhere) -- plain limma, no voom.

library(here)
library(yaml)
library(arrow)
source(here("R/lib/data_paths.R"))
source(here("R/lib/matrices.R"))
source(here("R/de/de_helpers.R"))

ds_meta <- read_dataset_yaml(here("config/MARS_config.yml"))$dataset

message("[MARS] loading raw expression matrix...")
mat <- pivot_expression_long(ds_meta)

sample_meta <- arrow::read_parquet(ds_meta$sample_metadata_path)
feature_meta <- arrow::read_parquet(ds_meta$feature_metadata_path)

message("[MARS] collapsing probes to gene symbol...")
mat <- collapse_to_symbol(mat, feature_meta, ds_meta$feature_id_col, ds_meta$symbol_col)
mat <- filter_low_variance(mat)
message("[MARS] ", nrow(mat), " genes x ", ncol(mat), " samples after collapse + low-variance filter")

# healthy_control is logical (TRUE/FALSE) in the parquet -- make.names()
# on "TRUE"/"FALSE" is safe as-is, but spelling it out as a factor with
# readable labels makes the output contrast name self-explanatory
# ("SepsisvsHealthy" rather than "TRUEvsFALSE" or the reverse).
sample_meta$sepsis_group <- ifelse(sample_meta$healthy_control, "Healthy", "Sepsis")

de <- run_de_case_control(
  mat, sample_meta,
  sample_id_col = ds_meta$sample_id_col,
  group_col     = "sepsis_group",
  platform      = "array"
)

out_dir <- here("results/de/MARS")
write_de_results(de, out_dir)
message("[MARS] done -- results in ", out_dir)

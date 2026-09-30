# Differential expression for dilgom (E-TABM-1036) -- CROSS-SECTIONAL
# (confirmed: 518 samples, 518 distinct sample_id, no repeat/subject
# column). See R/de/README.md.
#
# ** This cohort has NO sepsis or disease-relevant covariate at all ** --
# it's a healthy population reference cohort (real metadata: sample_id,
# accession, age_band, sex -- nothing else). The `age_band` (oldest vs.
# youngest quintile) contrast below is included ONLY for pipeline
# completeness (so every configured dataset has a script, per the user's
# direction), NOT because it's a sepsis-relevant finding -- treat this
# script's output as a working example / QC check that the pipeline runs
# on this dataset, not as a result to interpret biologically.
#
# Platform: array (log2 intensity, range 6.7-16.0, no count/DESeq2 logic)
# -- plain limma, no voom.

library(here)
library(yaml)
library(arrow)
source(here("R/lib/matrices.R"))
source(here("R/de/de_helpers.R"))

ds_meta <- yaml::read_yaml(here("config/dilgom_config.yml"))$dataset

message("[dilgom] loading raw expression matrix...")
mat <- pivot_expression_long(ds_meta)

sample_meta <- arrow::read_parquet(ds_meta$sample_metadata_path)
feature_meta <- arrow::read_parquet(ds_meta$feature_metadata_path)

message("[dilgom] collapsing probes to gene symbol...")
mat <- collapse_to_symbol(mat, feature_meta, ds_meta$feature_id_col, ds_meta$symbol_col)
mat <- filter_low_variance(mat)
message("[dilgom] ", nrow(mat), " genes x ", ncol(mat), " samples after collapse + low-variance filter")

# Oldest vs. youngest quintile only -- keeps this an illustrative 2-group
# contrast rather than a 5-level age model (this dataset has no
# disease-relevant question to actually answer, see header).
keep <- sample_meta$age_band %in% c("25-34", "65-74")
sample_meta <- sample_meta[keep, ]
mat <- mat[, sample_meta[[ds_meta$sample_id_col]], drop = FALSE]

de <- run_de_case_control(
  mat, sample_meta,
  sample_id_col = ds_meta$sample_id_col,
  group_col     = "age_band",
  platform      = "array"
)

out_dir <- here("results/de/dilgom")
write_de_results(de, out_dir)
message("[dilgom] done -- results in ", out_dir, " (NOTE: no sepsis-relevant covariate exists for this cohort -- see file header)")

# Differential expression for "gains" -- CROSS-SECTIONAL per its current
# usable metadata (no timepoint/subject column). See docs/methods.qmd ("Differential expression") for
# why this is treated as cross-sectional-with-caveats rather than a true
# repeated-measures dataset: raw sample ids hint at undocumented repeated
# sampling (e.g. "CAP0173.B.5") that isn't safely reconstructable from the
# metadata as it stands today.
#
# `preprocessing_script: R/preprocessing/gains_preprocessing.R` in the
# config is an unfinished stub (never builds a matrix) -- this script
# builds the matrix itself. `expression_path` is NOT a single parquet file
# for this dataset -- it's a directory of 4 accession-partitioned files
# (expression/accession=<id>/part-0.parquet), so pivot_expression_long()
# (which reads one file via arrow::read_parquet()) doesn't apply;
# arrow::open_dataset() + dplyr::collect() reads all 4 at once instead.
# `sample_id` is only unique WITHIN one accession (confirmed: 56 sample_id
# values recur across different accessions in the raw metadata, almost
# certainly different physical samples reusing the same code in a
# different sub-cohort) -- pivoting on sample_id alone would silently
# merge unrelated columns, so every join/pivot here uses the composite key
# accession + "__" + sample_id instead (confirmed unique: 0 duplicated
# (accession, sample_id) pairs in sample_metadata.parquet).
#
# Real metadata (729 rows, 8 cols): sample_id, accession, age, sex,
# survived_28d, srs_group, disease_state, author_use. `srs_group`
# (SRS1/SRS2 -- the Sepsis Response Signature, a well-established
# sepsis severity/immune classification for exactly this kind of cohort)
# is the contrast used here.
#
# Platform: array (log-scale, range -13 to 16.2, no count/DESeq2 logic
# despite the stub script's `library(DESeq2)`) -- plain limma, no voom.

library(here)
library(yaml)
library(arrow)
library(dplyr)
library(tidyr)
source(here("R/lib/data_paths.R"))
source(here("R/de/de_helpers.R"))

ds_meta <- read_dataset_yaml(here("config/gains_config.yml"))$dataset

message("[gains] loading raw expression matrix (4 partitioned files)...")
expr_long <- arrow::open_dataset(ds_meta$expression_path) |> dplyr::collect()
expr_long$key <- paste(expr_long$accession, expr_long$sample_id, sep = "__")
wide <- tidyr::pivot_wider(expr_long[, c("IlluminaID", "key", "value")],
                            id_cols = "IlluminaID", names_from = "key", values_from = "value")
mat <- as.matrix(wide[, setdiff(names(wide), "IlluminaID")])
rownames(mat) <- wide$IlluminaID

# The 4 accessions don't all share the identical probe set (different chip
# batches) -- pivot_wider() leaves NA for a probe missing from a given
# accession. Only 3732/28538 probes are affected; dropping them (rather
# than imputing) is the simplest safe fix, confirmed cheap here.
complete <- rowSums(is.na(mat)) == 0
message("[gains] dropping ", sum(!complete), "/", nrow(mat), " probes not shared across all 4 accessions")
mat <- mat[complete, , drop = FALSE]

sample_meta <- arrow::read_parquet(ds_meta$sample_metadata_path)
sample_meta$key <- paste(sample_meta$accession, sample_meta$sample_id, sep = "__")
sample_meta <- sample_meta[sample_meta$key %in% colnames(mat) & !is.na(sample_meta$srs_group), ]
mat <- mat[, sample_meta$key, drop = FALSE]

feature_meta <- arrow::read_parquet(ds_meta$feature_metadata_path)
message("[gains] collapsing probes to gene symbol...")
mat <- collapse_to_symbol(mat, feature_meta, ds_meta$feature_id_col, ds_meta$symbol_col)
mat <- filter_low_variance(mat)
message("[gains] ", nrow(mat), " genes x ", ncol(mat), " samples after collapse + low-variance filter")

de <- run_de_case_control(
  mat, sample_meta,
  sample_id_col = "key",
  group_col     = "srs_group",
  platform      = "array"
)

out_dir <- here("results/de/gains")
write_de_results(de, out_dir)
message("[gains] done -- results in ", out_dir)

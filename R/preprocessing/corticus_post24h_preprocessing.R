library(tidyverse)
library(arrow)
library(here)

# Per-timepoint slice of CORTICUS (config/CORTICUS_post24h_config.yml) --
# filters to timepoint == "Post(24h)" from the SAME raw parquet trio
# the parent dataset config points at, then collapses to one probe per
# gene and keeps the top-variance genes, mirroring
# R/preprocessing/anemones_preprocessing.R's pattern but adapted to this
# Illumina annotation's own columns (ProbeQuality / GenomicLocation /
# SYMBOL rather than control_type / chromosomal_location / gene_symbol --
# CORTICUS has no base (non-timepoint) preprocessing script to mirror
# directly, see config/CORTICUS_config.yml's TODO).
framework_run <- exists("tmp_dir", inherits = FALSE)

if (!framework_run) {
    expression_path <- "~/projects/hf_sepsis_collection/CORTICUS/expression.parquet"
    sample_metadata_path <- "~/projects/hf_sepsis_collection/CORTICUS/sample_metadata.parquet"
    feature_metadata_path <- "~/projects/hf_sepsis_collection/CORTICUS/feature_metadata.parquet"
}

TIMEPOINT_COL <- "timepoint"
TIMEPOINT_VALUE <- "Post(24h)"
FEATURE_ID_COL <- "feature_id"

raw <- list(
    ex = dplyr::collect(arrow::open_dataset(expression_path)),
    meta = arrow::read_parquet(sample_metadata_path),
    feature = arrow::read_parquet(feature_metadata_path)
)

meta_tp <- raw$meta |> filter(.data[[TIMEPOINT_COL]] == TIMEPOINT_VALUE)
stopifnot(nrow(meta_tp) > 0)

ex_tp <- raw$ex |>
    dplyr::rename(feature_id = all_of(FEATURE_ID_COL)) |>
    filter(sample_id %in% meta_tp$sample_id)

feature <- raw$feature |> dplyr::rename(feature_id = all_of(FEATURE_ID_COL))

# Drop low-quality/unmatched probes ("Bad"/"No match" ProbeQuality) and
# sex-chromosome/mitochondrial probes (parsed from "chrN:start:end:strand"
# GenomicLocation strings) -- the closest equivalent here to
# anemones_preprocessing.R's control_type/chromosomal_location filter;
# this annotation has no explicit spike-in control_type column.
bad_or_sex_mito_probes <- feature |>
    filter(
        ProbeQuality %in% c("Bad", "No match") |
        str_detect(GenomicLocation, "^chr[XYM]:")
    ) |>
    pull(feature_id)

ann <- feature |>
    filter(!feature_id %in% bad_or_sex_mito_probes, SYMBOL != "", !is.na(SYMBOL)) |>
    select(feature_id, SYMBOL)

# One probe per gene, keeping the highest-mean probe (mean taken WITHIN
# this timepoint's samples only, not the full dataset).
best_probe_per_gene <- ex_tp |>
    filter(!feature_id %in% bad_or_sex_mito_probes) |>
    group_by(feature_id) |>
    summarize(mean_value = mean(value), .groups = "drop") |>
    inner_join(ann, by = "feature_id") |>
    group_by(SYMBOL) |>
    slice_max(mean_value, n = 1, with_ties = FALSE) |>
    ungroup() |>
    pull(feature_id)

mat <- ex_tp |>
    filter(feature_id %in% best_probe_per_gene) |>
    pivot_wider(id_cols = feature_id, names_from = sample_id, values_from = value) |>
    dplyr::select(feature_id, all_of(meta_tp$sample_id)) |>
    column_to_rownames("feature_id") |>
    as.matrix()

gene_var <- matrixStats::rowVars(mat)
top_n <- 8000
mat_filt <- mat[order(gene_var, decreasing = TRUE)[seq_len(min(top_n, nrow(mat)))], ]

if (framework_run) {
    saveRDS(mat_filt, file.path(tmp_dir, "matrix.rds"))
}


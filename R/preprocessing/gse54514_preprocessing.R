library(tidyverse)
library(arrow)
library(here)

# All samples/features from GSE54514 (config/GSE54514_config.yml) --
# collapses to one probe per gene and keeps the top-variance genes,
# mirroring R/preprocessing/anemones_preprocessing.R's pattern but adapted
# to this Illumina annotation's own columns (ProbeQuality / GenomicLocation
# / SYMBOL rather than control_type / chromosomal_location / gene_symbol).
# Feature id here is IlluminaID itself -- there is no separate literal
# "feature_id" column in this dataset's metadata (see
# config/GSE54514_config.yml's feature_id_col). See
# R/preprocessing/gse54514_day1_preprocessing.R (etc., one per day 1-5) for
# the same pattern applied to a single timepoint's samples.
framework_run <- exists("tmp_dir", inherits = FALSE)

if (!framework_run) {
    expression_path <- "~/projects/hf_sepsis_collection/GSE54514/expression.parquet"
    sample_metadata_path <- "~/projects/hf_sepsis_collection/GSE54514/sample_metadata.parquet"
    feature_metadata_path <- "~/projects/hf_sepsis_collection/GSE54514/feature_metadata.parquet"
}

raw <- list(
    ex = dplyr::collect(arrow::open_dataset(expression_path)) |>
        dplyr::rename(feature_id = IlluminaID),
    meta = arrow::read_parquet(sample_metadata_path),
    feature = arrow::read_parquet(feature_metadata_path) |>
        dplyr::rename(feature_id = IlluminaID)
)

# Drop low-quality/unmatched probes ("Bad"/"No match" ProbeQuality) and
# sex-chromosome/mitochondrial probes (parsed from "chrN:start:end:strand"
# GenomicLocation strings) -- the closest equivalent here to
# anemones_preprocessing.R's control_type/chromosomal_location filter;
# this annotation has no explicit spike-in control_type column.
bad_or_sex_mito_probes <- raw$feature |>
    filter(
        ProbeQuality %in% c("Bad", "No match") |
        str_detect(GenomicLocation, "^chr[XYM]:")
    ) |>
    pull(feature_id)

ann <- raw$feature |>
    filter(!feature_id %in% bad_or_sex_mito_probes, SYMBOL != "", !is.na(SYMBOL)) |>
    select(feature_id, SYMBOL)

# One probe per gene, keeping the highest-mean probe.
best_probe_per_gene <- raw$ex |>
    filter(!feature_id %in% bad_or_sex_mito_probes) |>
    group_by(feature_id) |>
    summarize(mean_value = mean(value), .groups = "drop") |>
    inner_join(ann, by = "feature_id") |>
    group_by(SYMBOL) |>
    slice_max(mean_value, n = 1, with_ties = FALSE) |>
    ungroup() |>
    pull(feature_id)

mat <- raw$ex |>
    filter(feature_id %in% best_probe_per_gene) |>
    pivot_wider(id_cols = feature_id, names_from = sample_id, values_from = value) |>
    dplyr::select(feature_id, all_of(raw$meta$sample_id)) |>
    column_to_rownames("feature_id") |>
    as.matrix()

gene_var <- matrixStats::rowVars(mat)
top_n <- 8000
mat_filt <- mat[order(gene_var, decreasing = TRUE)[seq_len(min(top_n, nrow(mat)))], ]

if (framework_run) {
    saveRDS(mat_filt, file.path(tmp_dir, "matrix.rds"))
}

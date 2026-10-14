library(tidyverse)
library(arrow)
library(here)

# All samples/features from MARS (GSE65682, config/MARS_config.yml): RMA
# log2 Affymetrix values as deposited. Collapses to one probe set per gene
# symbol (highest mean) and keeps the top-variance genes, the same pattern
# as gse13904_preprocessing.R. The deposited data has no AFFX- control probe
# sets; the filter is kept so the pattern matches the other Affymetrix
# datasets.
framework_run <- exists("tmp_dir", inherits = FALSE)

if (!framework_run) {
    expression_path <- "~/projects/hf_sepsis_collection/MARS/expression.parquet"
    sample_metadata_path <- "~/projects/hf_sepsis_collection/MARS/sample_metadata.parquet"
    feature_metadata_path <- "~/projects/hf_sepsis_collection/MARS/feature_metadata.parquet"
}

raw <- list(
    ex = dplyr::collect(arrow::open_dataset(expression_path)),
    meta = arrow::read_parquet(sample_metadata_path),
    feature = arrow::read_parquet(feature_metadata_path)
)

control_probes <- raw$feature |>
    filter(str_detect(feature_id, "^AFFX")) |>
    pull(feature_id)

ann <- raw$feature |>
    filter(!feature_id %in% control_probes, SYMBOL != "", !is.na(SYMBOL)) |>
    select(feature_id, symbol = SYMBOL)

# One probe set per gene, keeping the highest-mean probe set.
best_probe_per_gene <- raw$ex |>
    filter(!feature_id %in% control_probes) |>
    group_by(feature_id) |>
    summarize(mean_value = mean(value), .groups = "drop") |>
    inner_join(ann, by = "feature_id") |>
    group_by(symbol) |>
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

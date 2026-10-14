library(tidyverse)
library(arrow)
library(here)

# Per-timepoint slice of GSE95233 (config/GSE95233_d02_config.yml) --
# filters to timepoint == "D02" from the SAME raw parquet trio
# the parent dataset config points at, then collapses to one probe set
# per gene and keeps the top-variance genes, mirroring
# R/preprocessing/anemones_preprocessing.R's pattern but adapted to this
# Affymetrix annotation (AFFX-prefixed control probe sets / lowercase
# `symbol` column rather than control_type / gene_symbol -- GSE95233 has no
# base (non-timepoint) preprocessing script to mirror directly, see
# config/GSE95233_config.yml's TODO).
framework_run <- exists("tmp_dir", inherits = FALSE)

if (!framework_run) {
    expression_path <- "~/projects/hf_sepsis_collection/GSE95233/expression.parquet"
    sample_metadata_path <- "~/projects/hf_sepsis_collection/GSE95233/sample_metadata.parquet"
    feature_metadata_path <- "~/projects/hf_sepsis_collection/GSE95233/feature_metadata.parquet"
}

TIMEPOINT_COL <- "timepoint"
TIMEPOINT_VALUE <- "D02"

raw <- list(
    ex = dplyr::collect(arrow::open_dataset(expression_path)),
    meta = arrow::read_parquet(sample_metadata_path),
    feature = arrow::read_parquet(feature_metadata_path)
)

meta_tp <- raw$meta |> filter(.data[[TIMEPOINT_COL]] == TIMEPOINT_VALUE)
stopifnot(nrow(meta_tp) > 0)

ex_tp <- raw$ex |> filter(sample_id %in% meta_tp$sample_id)

# AFFX-* probe sets are Affymetrix spike-in/housekeeping controls, not
# real transcripts -- the closest equivalent here to
# anemones_preprocessing.R's control_type filter.
control_probes <- raw$feature |>
    filter(str_detect(feature_id, "^AFFX")) |>
    pull(feature_id)

ann <- raw$feature |>
    filter(!feature_id %in% control_probes, symbol != "", !is.na(symbol)) |>
    select(feature_id, symbol)

# One probe set per gene, keeping the highest-mean probe set (mean taken
# WITHIN this timepoint's samples only, not the full dataset).
best_probe_per_gene <- ex_tp |>
    filter(!feature_id %in% control_probes) |>
    group_by(feature_id) |>
    summarize(mean_value = mean(value), .groups = "drop") |>
    inner_join(ann, by = "feature_id") |>
    group_by(symbol) |>
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


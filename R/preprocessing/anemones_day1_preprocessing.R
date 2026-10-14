library(tidyverse)
library(arrow)
library(here)

# Per-timepoint slice of ANEMONES (config/ANEMONES_day1_config.yml) --
# filters to timepoint == "Day: 1" from the SAME raw parquet trio
# the parent dataset config points at, then runs the same probe-collapse +
# top-variance filter as R/preprocessing/anemones_preprocessing.R (the
# full-timecourse version of this pattern), on just this subset.
framework_run <- exists("tmp_dir", inherits = FALSE)

if (!framework_run) {
    expression_path <- "~/projects/hf_sepsis_collection/ANEMONES/expression.parquet"
    sample_metadata_path <- "~/projects/hf_sepsis_collection/ANEMONES/sample_metadata.parquet"
    feature_metadata_path <- "~/projects/hf_sepsis_collection/ANEMONES/feature_metadata.parquet"
}

TIMEPOINT_COL <- "timepoint"
TIMEPOINT_VALUE <- "Day: 1"

raw <- list(
    ex = dplyr::collect(arrow::open_dataset(expression_path)),
    meta = arrow::read_parquet(sample_metadata_path),
    feature = arrow::read_parquet(feature_metadata_path)
)

meta_tp <- raw$meta |> filter(.data[[TIMEPOINT_COL]] == TIMEPOINT_VALUE)
stopifnot(nrow(meta_tp) > 0)

ex_tp <- raw$ex |> filter(sample_id %in% meta_tp$sample_id)

control_probes_germ_and_mito <- raw$feature |>
    filter(control_type %in% c("pos", "neg") | str_detect(chromosomal_location, "chrX|chrY|chrM")) |>
    pull(feature_id)

ann <- raw$feature |>
    filter(!feature_id %in% control_probes_germ_and_mito, gene_symbol != "", !is.na(gene_symbol)) |>
    select(feature_id, gene_symbol)

# One probe per gene, keeping the highest-mean probe (mean taken WITHIN
# this timepoint's samples only, not the full dataset).
best_probe_per_gene <- ex_tp |>
    filter(!feature_id %in% control_probes_germ_and_mito) |>
    group_by(feature_id) |>
    summarize(mean_value = mean(value), .groups = "drop") |>
    inner_join(ann, by = "feature_id") |>
    group_by(gene_symbol) |>
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
# min(top_n, nrow(mat)): see R/preprocessing/anemones_preprocessing.R's
# comment -- a small enough gene pool after collapsing/annotation makes
# order(...)[seq_len(top_n)] index past the end and insert all-NA rows.
mat_filt <- mat[order(gene_var, decreasing = TRUE)[seq_len(min(top_n, nrow(mat)))], ]

if (framework_run) {
    saveRDS(mat_filt, file.path(tmp_dir, "matrix.rds"))
}


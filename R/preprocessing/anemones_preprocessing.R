library(tidyverse)
library(arrow)
library(tidyverse)
library(here)

# When the `mat_<dataset>` target runs this script
# (run_preprocessing_script(), R/lib/matrices.R), expression_path/
# sample_metadata_path/feature_metadata_path/tmp_dir are already defined in
# its environment; tmp_dir's presence is what distinguishes that from running
# it line by line in the console, which falls back to local paths and skips
# the final saveRDS().
framework_run <- exists("tmp_dir", inherits = FALSE)

if (!framework_run) {
    expression_path <- "~/projects/hf_sepsis_collection/ANEMONES/expression.parquet"
    sample_metadata_path <- "~/projects/hf_sepsis_collection/ANEMONES/sample_metadata.parquet"
    feature_metadata_path <- "~/projects/hf_sepsis_collection/ANEMONES/feature_metadata.parquet"
}

raw <- list(
    ex = dplyr::collect(arrow::open_dataset(expression_path)),
    meta = arrow::read_parquet(sample_metadata_path),
    feature = arrow::read_parquet(feature_metadata_path)
)

control_probes_germ_and_mito <- raw$feature |>
    filter(control_type %in% c("pos", "neg") | str_detect(chromosomal_location, "chrX|chrY|chrM")) |>
    pull(feature_id)

ann <- raw$feature |>
    filter(!feature_id %in% control_probes_germ_and_mito, gene_symbol != "", !is.na(gene_symbol)) |>
    select(feature_id, gene_symbol)

# One probe per gene, keeping the highest-mean probe.
# Also selects only autosomal loci
best_probe_per_gene <- raw$ex |>
    filter(!feature_id %in% control_probes_germ_and_mito) |>
    group_by(feature_id) |>
    summarize(mean_value = mean(value), .groups = "drop") |>
    inner_join(ann, by = "feature_id") |>
    group_by(gene_symbol) |>
    slice_max(mean_value, n = 1, with_ties = FALSE) |>
    ungroup() |>
    pull(feature_id)

mat <- raw$ex |>
    filter(feature_id %in% best_probe_per_gene) |>
    pivot_wider(id_cols = feature_id, names_from = sample_id, values_from = value) |>
    column_to_rownames("feature_id") |>
    as.matrix()

gene_var <- matrixStats::rowVars(mat)
top_n <- 8000
# min(top_n, nrow(mat)): if fewer genes survive collapsing/annotation
# than top_n, order(...)[seq_len(top_n)] indexes past the end and
# returns NA for the missing positions -- mat[NA_rows, ] then silently
# inserts all-NA rows, which blows up downstream (svd()/prcomp() error
# with "infinite or missing values") rather than erroring here where the
# actual cause is obvious. Confirmed this bites dilgom specifically
# (~6,071 genes survive there, below top_n) even though it's not
# dilgom-specific in cause -- any small enough gene pool triggers it.
mat_filt <- mat[order(gene_var, decreasing = TRUE)[seq_len(min(top_n, nrow(mat)))], ]

if (framework_run) {
    saveRDS(mat_filt, file.path(tmp_dir, "matrix.rds"))
}

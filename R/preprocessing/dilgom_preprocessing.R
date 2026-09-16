library(tidyverse)
library(arrow)
library(tidyverse)
library(here)

# When this script is sourced by the rslurm setup scripts (see
# R/lib/matrices.R::run_preprocessing_script()), expression_path/
# sample_metadata_path/feature_metadata_path/tmp_dir are already defined in
# this script's environment -- tmp_dir is unique to that call path, so its
# presence is what distinguishes "run by the framework" from "run directly
# (e.g. line-by-line in the console) for interactive exploration." In the
# interactive case, fall back to local paths and skip the final saveRDS()
# below -- exploring intermediate objects is the point, not writing output.
framework_run <- exists("tmp_dir", inherits = FALSE)

if (!framework_run) {
    expression_path <- "~/projects/hf_sepsis_collection/dilgom/expression.parquet"
    sample_metadata_path <- "~/projects/hf_sepsis_collection/dilgom/sample_metadata.parquet"
    feature_metadata_path <- "~/projects/hf_sepsis_collection/dilgom/feature_metadata.parquet"
}

raw <- list(
    ex = dplyr::collect(arrow::open_dataset(expression_path)),
    meta = arrow::read_parquet(sample_metadata_path),
    feature = arrow::read_parquet(feature_metadata_path)
)

# removing probequality == bad also removes intergenic probes
remove_probes <- raw$feature |>
    filter(ProbeQuality == "Bad" |
        str_detect(GenomicLocation, "chrX|chrY|chrM")) |>
    pull(IlluminaID)

ann <- raw$feature |>
    filter(!IlluminaID %in% remove_probes, SYMBOL != "", !is.na(SYMBOL)) |>
    select(IlluminaID, gene_symbol = SYMBOL)

# One probe per gene, keeping the highest-mean probe.
# Also selects only autosomal loci
best_probe_per_gene <- raw$ex |>
    filter(!IlluminaID %in% remove_probes) |>
    group_by(IlluminaID) |>
    summarize(mean_value = mean(value), .groups = "drop") |>
    inner_join(ann, by = "IlluminaID") |>
    group_by(gene_symbol) |>
    slice_max(mean_value, n = 1, with_ties = FALSE) |>
    ungroup() |>
    pull(IlluminaID)

mat <- raw$ex |>
    filter(IlluminaID %in% best_probe_per_gene) |>
    pivot_wider(id_cols = IlluminaID, names_from = sample_id, values_from = value) |>
    column_to_rownames("IlluminaID") |>
    as.matrix()

gene_var <- matrixStats::rowVars(mat)
top_n <- 8000
# min(top_n, nrow(mat)): dilgom's Illumina chip only has 8,765 probes
# total, and only ~6,071 genes survive quality/annotation
# filtering+collapsing above -- below top_n=8000. Without the min(),
# order(...)[seq_len(top_n)] indexes past the end of a length-6071
# vector and returns NA for the missing positions; mat[NA_rows, ] then
# silently inserts all-NA rows into mat_filt, which is exactly what made
# spca()'s svd() fail downstream with "infinite or missing values in
# 'x'" (the actual cause was here, not in spca()).
mat_filt <- mat[order(gene_var, decreasing = TRUE)[seq_len(min(top_n, nrow(mat)))], ]

if (framework_run) {
    saveRDS(mat_filt, file.path(tmp_dir, "matrix.rds"))
}

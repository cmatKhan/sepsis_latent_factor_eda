library(tidyverse)
library(arrow)
library(here)

# All samples/features from SHIP-TREND (GSE36382, config/SHIP-TREND_config.yml):
# a healthy population cohort, 991 samples, Illumina HumanHT-12 v3 values
# deposited quantile-normalized and log2. There is no sample metadata
# (phenotype data is access-restricted), so samples are taken from the
# expression file itself. Same pattern as dilgom_preprocessing.R (also
# Illumina v3): drop probes the reannotation marks "Bad" (which also removes
# intergenic probes) and those on chrX/chrY/chrM, collapse to one probe per
# gene symbol (highest mean), keep the top-variance genes. The expression
# file's extra `accession` column is ignored.
framework_run <- exists("tmp_dir", inherits = FALSE)

if (!framework_run) {
    expression_path <- "~/projects/hf_sepsis_collection/SHIP-TREND/expression.parquet"
    feature_metadata_path <- "~/projects/hf_sepsis_collection/SHIP-TREND/feature_metadata.parquet"
}

raw <- list(
    ex = dplyr::collect(dplyr::select(arrow::open_dataset(expression_path), IlluminaID, sample_id, value)),
    feature = arrow::read_parquet(feature_metadata_path)
)

remove_probes <- raw$feature |>
    filter(ProbeQuality == "Bad" |
        str_detect(GenomicLocation, "chrX|chrY|chrM")) |>
    pull(IlluminaID)

ann <- raw$feature |>
    filter(!IlluminaID %in% remove_probes, SYMBOL != "", !is.na(SYMBOL)) |>
    select(IlluminaID, gene_symbol = SYMBOL)

# One probe per gene, keeping the highest-mean probe.
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
mat_filt <- mat[order(gene_var, decreasing = TRUE)[seq_len(min(top_n, nrow(mat)))], ]

if (framework_run) {
    saveRDS(mat_filt, file.path(tmp_dir, "matrix.rds"))
}

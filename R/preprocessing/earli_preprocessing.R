library(tidyverse)
library(here)
library(DESeq2)

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
    expression_path <- "~/projects/hf_sepsis_collection/EARLI/expression.parquet"
    sample_metadata_path <- "~/projects/hf_sepsis_collection/EARLI/sample_metadata.parquet"
    feature_metadata_path <- "~/projects/hf_sepsis_collection/EARLI/feature_metadata.parquet"
}

raw <- list(
    ex = arrow::read_parquet(expression_path),
    meta = arrow::read_parquet(sample_metadata_path),
    feature = arrow::read_parquet(feature_metadata_path)
)



mat <- raw$ex |>
    mutate(value = as.integer(value)) |>
    pivot_wider(id_cols = feature_id, names_from = sample_id, values_from = value) |>
    dplyr::select(feature_id, all_of(raw$meta$sample_id)) |>
    column_to_rownames(var = "feature_id") |>
    as.matrix()

mat <- mat[raw$feature$feature_id, ]

stopifnot(identical(colnames(mat), raw$meta$sample_id))
stopifnot(identical(rownames(mat), raw$feature$feature_id))

dds <- DESeqDataSetFromMatrix(
    countData = mat,
    colData = raw$meta,
    rowRanges = GenomicRanges::GRanges(raw$feature),
    design = ~1
)


# Keep genes where at least 3 samples have a count of 10 or more
# AND remove hemoglobin genes and ribosomal RNA loci
smallest_group_size <- 189
keep <- rowSums(counts(dds) >= 10) >= smallest_group_size &
    rowRanges(dds)$autosomal_protein_coding &
    !rowRanges(dds)$hemoglobin_related

# Subset the dds object
dds_filt <- dds[keep, ]

dds_filt <- estimateSizeFactors(dds_filt)
# note: the local fit looked quite wonky, with the end behavior being pulled
# up by some outlier genes that appear to have high leverage. mean looks more
# appropriate, especially given the low expression filtering already
# performed above
dds_filt <- estimateDispersions(dds_filt, fitType = "mean")
plotDispEsts(dds_filt)

vsd <- vst(dds_filt, blind = TRUE)
vst_mat <- assay(vsd)

vst_mat <- assay(vsd)

gene_var <- rowVars(vst_mat)
gene_mean <- rowMeans(vst_mat)

hist(gene_var, breaks = 100)
sapply(c(0.1, 0.3, 0.5, 1, 1.5), function(v) sum(gene_var > v))

min_var <- 0.3
keep_var <- gene_var > min_var
vst_filt <- vst_mat[keep_var, ]

if (framework_run) {
    saveRDS(vst_filt, file.path(tmp_dir, "matrix.rds"))
}

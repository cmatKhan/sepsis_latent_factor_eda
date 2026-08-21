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
    expression_path <- "~/projects/hf_sepsis_collection/GSE273700/expression.parquet"
    sample_metadata_path <- "~/projects/hf_sepsis_collection/GSE273700/sample_metadata.parquet"
    feature_metadata_path <- "~/projects/hf_sepsis_collection/GSE273700/feature_metadata.parquet"
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
    rowData = raw$feature,
    design = ~1
)

# remove hemoglobin genes
hb_ids <- c(
    HBA2  = "ENSG00000188536",
    HBG2  = "ENSG00000196565",
    HBA1  = "ENSG00000206172",
    HBG1  = "ENSG00000213934",
    HBD   = "ENSG00000223609",
    HBBP1 = "ENSG00000229988",
    HBB   = "ENSG00000244734"
)


# Keep genes where at least 3 samples have a count of 10 or more
# AND remove hemoglobin genes and ribosomal RNA loci
smallest_group_size <- 62
keep <- rowSums(counts(dds) >= 10) >= smallest_group_size &
    !rownames(dds) %in% hb_ids

# Subset the dds object
dds_filt <- dds[keep, ]

# plotDispEsts(dds_filt)

vsd <- vst(dds_filt, blind = TRUE)
vst_mat <- assay(vsd)

vst_mat <- assay(vsd)

gene_var <- rowVars(vst_mat)
gene_mean <- rowMeans(vst_mat)

# hist(gene_var, breaks = 100)
# sapply(c(0.1, 0.3, 0.5, 1, 1.5), function(v) sum(gene_var > v))

min_var <- 0.3
keep_var <- gene_var > min_var
vst_filt <- vst_mat[keep_var, ]

if (framework_run) {
    saveRDS(vst_filt, file.path(tmp_dir, "matrix.rds"))
}

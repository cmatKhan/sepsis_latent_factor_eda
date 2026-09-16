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
    expression_path <- "~/projects/hf_sepsis_collection/GSE110487/expression.parquet"
    sample_metadata_path <- "~/projects/hf_sepsis_collection/GSE110487/sample_metadata.parquet"
    feature_metadata_path <- "~/projects/hf_sepsis_collection/GSE110487/feature_metadata.parquet"
}

raw <- list(
    ex = arrow::read_parquet(expression_path),
    meta = arrow::read_parquet(sample_metadata_path),
    feature = arrow::read_parquet(feature_metadata_path)
)

mat <- raw$ex |>
    mutate(value = as.integer(value)) |>
    pivot_wider(id_cols = GeneID, names_from = sample_id, values_from = value) |>
    dplyr::select(GeneID, all_of(raw$meta$sample_id)) |>
    column_to_rownames(var = "GeneID") |>
    as.matrix()

mat <- mat[as.character(raw$feature$GeneID), ]

stopifnot(identical(colnames(mat), raw$meta$sample_id))
stopifnot(identical(rownames(mat), as.character(raw$feature$GeneID)))

dds <- DESeqDataSetFromMatrix(
    countData = mat,
    colData = raw$meta,
    rowData = raw$feature,
    design = ~1
)

# remove hemoglobin genes
hb_entrez_ids <- c(
    HBA1 = "3039",
    HBA2 = "3040",
    HBB = "3043",
    HBG1 = "3047",
    HBBP1 = "3044",
    HBD = "3045",
    HBG2 = "3048"
)

# remove ribosomal RNA loci -- residual rRNA depletion carryover shows up as
# a handful of very-high-mean, very-high-dispersion genes that break
# DESeq2's parametric dispersion fit (see conversation history / dispersion
# diagnostics: these were among the top-15 highest-dispersion genes despite
# high mean expression, the opposite of the technical-noise trend the
# parametric fit assumes)
rrna_entrez_ids <- c(
    RNA18SN1 = "106631781",
    RNA18SN2 = "109864280",
    RNA18SN4 = "109864273",
    `RNA5-8SN3` = "109910381",
    `RNA5-8SN4` = "109864274",
    RNA28SN5 = "100008589"
)


# Keep genes where at least 3 samples have a count of 10 or more
# AND remove hemoglobin genes and ribosomal RNA loci
smallest_group_size <- 62
keep <- rowSums(counts(dds) >= 10) >= smallest_group_size &
    !rownames(dds) %in% hb_entrez_ids &
    !rownames(dds) %in% rrna_entrez_ids

# Subset the dds object
dds_filt <- dds[keep, ]

vsd <- vst(dds_filt, blind = TRUE, fitType = "local")
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

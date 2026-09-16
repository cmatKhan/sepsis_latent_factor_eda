library(tidyverse)
library(here)
library(DESeq2)

# All samples/features from hfgp-500fg (GSE134080) -- config/hfgp-500fg_config.yml.
# RNA-seq raw counts (verified directly: expression.parquet's `value` is
# integer, range 0-225809), same DESeq2 VST + gene-filter pattern as
# R/preprocessing/gse110487_preprocessing.R / gse273700_preprocessing.R.
#
# Two things genuinely differ from those two scripts (not just relabeled):
#   - this dataset's own sample-identifying column is `geo_accession`, NOT
#     `sample_id` (verified directly against both parquet files -- no
#     column literally named `sample_id` exists here at all), so this
#     script uses the framework-provided `sample_id_col` variable
#     (contract: R/README.md's "Preprocessing script contract") rather
#     than hardcoding a column name that doesn't exist for this dataset.
#   - only hemoglobin genes are removed, not ribosomal RNA loci --
#     GSE110487's rRNA-depletion-carryover dispersion issue was diagnosed
#     specifically against ITS OWN dispersion plot; nothing here confirms
#     the same issue exists in this dataset, so that filter isn't copied
#     over unverified. Revisit (add the rRNA filter, or force
#     fitType = "local") if plotDispEsts(dds_filt) shows the same pattern.
framework_run <- exists("tmp_dir", inherits = FALSE)

if (!framework_run) {
    expression_path <- "~/projects/hf_sepsis_collection/hfgp-500fg/expression.parquet"
    sample_metadata_path <- "~/projects/hf_sepsis_collection/hfgp-500fg/sample_metadata.parquet"
    feature_metadata_path <- "~/projects/hf_sepsis_collection/hfgp-500fg/feature_metadata.parquet"
    sample_id_col <- "geo_accession"
    feature_id_col <- "GeneID"
}

raw <- list(
    ex = arrow::read_parquet(expression_path),
    meta = arrow::read_parquet(sample_metadata_path),
    feature = arrow::read_parquet(feature_metadata_path)
)

mat <- raw$ex |>
    mutate(value = as.integer(value)) |>
    pivot_wider(id_cols = all_of(feature_id_col), names_from = all_of(sample_id_col), values_from = value) |>
    dplyr::select(all_of(feature_id_col), all_of(raw$meta[[sample_id_col]])) |>
    column_to_rownames(var = feature_id_col) |>
    as.matrix()

mat <- mat[as.character(raw$feature[[feature_id_col]]), ]

stopifnot(identical(colnames(mat), raw$meta[[sample_id_col]]))
stopifnot(identical(rownames(mat), as.character(raw$feature[[feature_id_col]])))

dds <- DESeqDataSetFromMatrix(
    countData = mat,
    colData = raw$meta,
    rowData = raw$feature,
    design = ~1
)

# remove hemoglobin genes (Entrez GeneIDs -- same global ids used in
# R/preprocessing/gse110487_preprocessing.R, this dataset's feature
# annotation uses the same NCBI Entrez GeneID scheme)
hb_entrez_ids <- c(
    HBA1 = "3039",
    HBA2 = "3040",
    HBB = "3043",
    HBG1 = "3047",
    HBBP1 = "3044",
    HBD = "3045",
    HBG2 = "3048"
)

# Keep genes where every sample has a count of 10 or more (design = ~1 is
# a single group, so "smallest group" is just every sample) AND remove
# hemoglobin genes.
smallest_group_size <- ncol(dds)
keep <- rowSums(counts(dds) >= 10) >= smallest_group_size &
    !rownames(dds) %in% hb_entrez_ids

dds_filt <- dds[keep, ]

vsd <- vst(dds_filt, blind = TRUE)
vst_mat <- assay(vsd)

gene_var <- rowVars(vst_mat)

# hist(gene_var, breaks = 100)
# sapply(c(0.1, 0.3, 0.5, 1, 1.5), function(v) sum(gene_var > v))

min_var <- 0.3
keep_var <- gene_var > min_var
vst_filt <- vst_mat[keep_var, ]

if (framework_run) {
    saveRDS(vst_filt, file.path(tmp_dir, "matrix.rds"))
}

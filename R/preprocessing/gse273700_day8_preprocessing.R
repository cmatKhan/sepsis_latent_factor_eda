library(tidyverse)
library(here)
library(DESeq2)

# Per-timepoint slice of GSE273700 (config/GSE273700_day8_config.yml) --
# filters to timepoint == "Day8" from the SAME raw parquet trio
# the parent dataset config points at, then runs the same DESeq2 VST +
# gene filter as R/preprocessing/gse273700_preprocessing.R (the full-timecourse
# version of this pattern), on just this subset.
framework_run <- exists("tmp_dir", inherits = FALSE)

if (!framework_run) {
    expression_path <- "~/projects/hf_sepsis_collection/GSE273700/expression.parquet"
    sample_metadata_path <- "~/projects/hf_sepsis_collection/GSE273700/sample_metadata.parquet"
    feature_metadata_path <- "~/projects/hf_sepsis_collection/GSE273700/feature_metadata.parquet"
}

TIMEPOINT_COL <- "timepoint"
TIMEPOINT_VALUE <- "Day8"

raw <- list(
    ex = arrow::read_parquet(expression_path),
    meta = arrow::read_parquet(sample_metadata_path),
    feature = arrow::read_parquet(feature_metadata_path)
)

meta_tp <- raw$meta |> filter(.data[[TIMEPOINT_COL]] == TIMEPOINT_VALUE)
stopifnot(nrow(meta_tp) > 0)

ex_tp <- raw$ex |> filter(sample_id %in% meta_tp$sample_id)

mat <- ex_tp |>
    mutate(value = as.integer(value)) |>
    pivot_wider(id_cols = feature_id, names_from = sample_id, values_from = value) |>
    dplyr::select(feature_id, all_of(meta_tp$sample_id)) |>
    column_to_rownames(var = "feature_id") |>
    as.matrix()

mat <- mat[as.character(raw$feature$feature_id), ]

stopifnot(identical(colnames(mat), meta_tp$sample_id))
stopifnot(identical(rownames(mat), as.character(raw$feature$feature_id)))

dds <- DESeqDataSetFromMatrix(
    countData = mat,
    colData = meta_tp,
    rowData = raw$feature,
    design = ~1
)

# remove hemoglobin genes
hb_entrez_ids <- c(
    HBA2  = "ENSG00000188536",
    HBG2  = "ENSG00000196565",
    HBA1  = "ENSG00000206172",
    HBG1  = "ENSG00000213934",
    HBD   = "ENSG00000223609",
    HBBP1 = "ENSG00000229988",
    HBB   = "ENSG00000244734"
)

# Keep genes where EVERY sample in this timepoint has a count of 10 or
# more (smallest_group_size = the number of samples in this subset,
# since design = ~1 is a single group) AND remove hemoglobin genes.
smallest_group_size <- ncol(dds)
keep <- rowSums(counts(dds) >= 10) >= smallest_group_size &
    !rownames(dds) %in% hb_entrez_ids

dds_filt <- dds[keep, ]

vsd <- vst(dds_filt, blind = TRUE)
vst_mat <- assay(vsd)

gene_var <- rowVars(vst_mat)

min_var <- 0.3
keep_var <- gene_var > min_var
vst_filt <- vst_mat[keep_var, ]

if (framework_run) {
    saveRDS(vst_filt, file.path(tmp_dir, "matrix.rds"))
}


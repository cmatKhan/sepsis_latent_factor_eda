library(tidyverse)
library(here)
library(DESeq2)

# Per-timepoint slice of GSE110487 (config/GSE110487_t2_config.yml) --
# filters to timepoint == "T2" from the SAME raw parquet trio
# the parent dataset config points at, then runs the same DESeq2 VST +
# gene filter as R/preprocessing/gse110487_preprocessing.R (the full-timecourse
# version of this pattern), on just this subset.
framework_run <- exists("tmp_dir", inherits = FALSE)

if (!framework_run) {
    expression_path <- "~/projects/hf_sepsis_collection/GSE110487/expression.parquet"
    sample_metadata_path <- "~/projects/hf_sepsis_collection/GSE110487/sample_metadata.parquet"
    feature_metadata_path <- "~/projects/hf_sepsis_collection/GSE110487/feature_metadata.parquet"
}

TIMEPOINT_COL <- "timepoint"
TIMEPOINT_VALUE <- "T2"

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
    pivot_wider(id_cols = GeneID, names_from = sample_id, values_from = value) |>
    dplyr::select(GeneID, all_of(meta_tp$sample_id)) |>
    column_to_rownames(var = "GeneID") |>
    as.matrix()

mat <- mat[as.character(raw$feature$GeneID), ]

stopifnot(identical(colnames(mat), meta_tp$sample_id))
stopifnot(identical(rownames(mat), as.character(raw$feature$GeneID)))

dds <- DESeqDataSetFromMatrix(
    countData = mat,
    colData = meta_tp,
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

# remove ribosomal RNA loci -- residual rRNA depletion carryover shows up
# as a handful of very-high-mean, very-high-dispersion genes that break
# DESeq2's parametric dispersion fit (see R/preprocessing/gse110487_preprocessing.R)
rrna_entrez_ids <- c(
    RNA18SN1 = "106631781",
    RNA18SN2 = "109864280",
    RNA18SN4 = "109864273",
    `RNA5-8SN3` = "109910381",
    `RNA5-8SN4` = "109864274",
    RNA28SN5 = "100008589"
)

# Keep genes where EVERY sample in this timepoint has a count of 10 or
# more (smallest_group_size = the number of samples in this subset,
# since design = ~1 is a single group) AND remove hemoglobin/rRNA genes.
smallest_group_size <- ncol(dds)
keep <- rowSums(counts(dds) >= 10) >= smallest_group_size &
    !rownames(dds) %in% hb_entrez_ids &
    !rownames(dds) %in% rrna_entrez_ids

dds_filt <- dds[keep, ]

vsd <- vst(dds_filt, blind = TRUE, fitType = "local")
vst_mat <- assay(vsd)

gene_var <- rowVars(vst_mat)

min_var <- 0.3
keep_var <- gene_var > min_var
vst_filt <- vst_mat[keep_var, ]

if (framework_run) {
    saveRDS(vst_filt, file.path(tmp_dir, "matrix.rds"))
}


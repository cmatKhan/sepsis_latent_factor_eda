library(tidyverse)
library(arrow)
library(tidyverse)
library(here)
library(DESeq2)

framework_run <- exists("tmp_dir", inherits = FALSE)

if (!framework_run) {
    expression_path <- "~/projects/hf_sepsis_collection/gains/expression"
    sample_metadata_path <- "~/projects/hf_sepsis_collection/gains/sample_metadata.parquet"
    feature_metadata_path <- "~/projects/hf_sepsis_collection/gains/feature_metadata.parquet"
}

raw <- list(
    ex = dplyr::collect(arrow::open_dataset(expression_path)),
    meta = arrow::read_parquet(sample_metadata_path),
    feature = arrow::read_parquet(feature_metadata_path)
)

# there are some difficulties in gains. there may be repeated measure samples,
# eg CAP0003.B.1, CAP0003.B.3, CAP0003.B.5 that are not documented in the
# metadata from GEO. Additionally, there are 53 subjects with the same id
# measured in more than 1 accession without exactly the same probe
# values, eg CAP0063

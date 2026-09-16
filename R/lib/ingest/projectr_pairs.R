# Builds the "reduced all-pairs" projectr job grid: every ordered
# (source dataset, target dataset) pair with source != target, crossed
# with every method's representative fits (representative_fit_ids(), see
# R/lib/ingest/redundancy.R -- one fit PER RANK for nmf/cogaps/ica, not
# collapsed across rank, only across seed; PCA's single max-rank fit still
# covers every smaller rank via column-slicing).
#
# include_intercept is fixed by whether a pair is within the SAME dataset
# family (config/dataset_families.yml) or across families -- this is
# baked in explicitly as a real column (and as which of
# projectr_within_grid/projectr_cross_grid the row goes into), never
# inferred implicitly downstream.

same_family <- function(families, ds_a, ds_b) {
  fam_a <- names(Filter(function(v) ds_a %in% v, families))
  fam_b <- names(Filter(function(v) ds_b %in% v, families))
  length(fam_a) == 1 && length(fam_b) == 1 && fam_a == fam_b
}

#' @param methods loadings-bearing methods to project -- cp/tucker's gene
#'   loadings are included by default despite their subject/time modes
#'   being tensor-specific (the gene x pattern loadings themselves are a
#'   perfectly valid plain-matrix projectR input).
build_projectr_pairs <- function(con, dataset_ids, families,
                                  methods = c("pca", "nmf", "cogaps", "spca", "ica", "cp", "tucker")) {
  rows <- list()
  for (ds_a in dataset_ids) {
    for (ds_b in setdiff(dataset_ids, ds_a)) {
      within <- same_family(families, ds_a, ds_b)
      for (m in methods) {
        fids <- representative_fit_ids(con, ds_a, m)
        if (length(fids) == 0) next
        rows[[length(rows) + 1]] <- data.frame(
          source_dataset_id = ds_a, target_dataset_id = ds_b, method = m,
          source_fit_id = fids,
          projection_type = if (within) "within_dataset" else "cross_dataset",
          include_intercept = !within,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(rows) == 0) return(data.frame(
    source_dataset_id = character(0), target_dataset_id = character(0), method = character(0),
    source_fit_id = integer(0), projection_type = character(0), include_intercept = logical(0)
  ))
  do.call(rbind, rows)
}

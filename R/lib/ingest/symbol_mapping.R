# Gene identifiers across datasets:
#   - build_ensembl_map() / remap_to_ensembl(): version-stripped Ensembl gene
#     ids, the identifier for all cross-dataset computation (enrichment,
#     projectR). Symbols are ambiguous and renamed across annotation
#     releases; Ensembl ids aren't.
#   - build_symbol_map(): gene symbols, for display only.
# The Ensembl map is built once per dataset (the ensembl_map_<dataset>
# target) from its feature metadata and its matrix, and is one-to-one: each
# gene is one feature of that dataset, so converting is a rename.

#' Strip an Ensembl id's version suffix
#'
#' @param x Character vector of ids (`ENSG00000001234.2` -> `ENSG00000001234`;
#'   ids without a suffix are unchanged).
#' @return `x` without version suffixes.
strip_ensembl_version <- function(x) sub("\\.[0-9]+$", "", x)

#' A dataset's feature id -> Ensembl gene id map, one feature per gene
#'
#' Reads `dataset.ensembl_col` (default `"ensembl"`) from the feature
#' metadata. A `;`-separated multi-gene entry (one probe, several genes) maps to
#' its first gene -- a simple, transparent choice; expanding one probe across
#' every listed gene is the alternative. The map is restricted to the rows of
#' `mat`, and where several of them map to one gene only the one with the
#' highest mean expression is kept (ties: the first in `mat`), so a gene is
#' the same feature in the matrix, in every fit's loadings and in every gene
#' list.
#'
#' @param dataset_yaml A resolved dataset config.
#' @param mat The dataset's feature x sample matrix.
#' @param fm Optional already-loaded feature-metadata data frame.
#' @return Named character vector (names = feature ids, in `mat`'s order), or
#'   `NULL` when the dataset has no such column or feature metadata -- its row
#'   names are then used as they are.
build_ensembl_map <- function(dataset_yaml, mat, fm = NULL) {
  ds <- dataset_yaml$dataset
  if (is.null(fm)) {
    if (is.null(ds$feature_metadata_path) || !file.exists(ds$feature_metadata_path)) return(NULL)
    fm <- arrow::read_parquet(ds$feature_metadata_path)
  }
  ens_col <- ds$ensembl_col %||% "ensembl"
  if (!(ens_col %in% names(fm))) return(NULL)
  id_col <- ds$feature_id_col %||% "feature_id"
  if (!(id_col %in% names(fm))) return(NULL)

  raw <- as.character(fm[[ens_col]])
  first_id <- vapply(strsplit(raw, ";", fixed = TRUE), function(x) {
    if (length(x) == 0) NA_character_ else x[1]
  }, character(1))

  map <- setNames(strip_ensembl_version(first_id), as.character(fm[[id_col]]))
  map <- map[!is.na(map) & nzchar(map)]
  feats <- intersect(rownames(mat), names(map))
  ens <- map[feats]
  ord <- order(ens, -rowMeans(mat[feats, , drop = FALSE]), seq_along(feats))
  keep <- feats[ord][!duplicated(ens[ord])]
  map[keep[order(match(keep, rownames(mat)))]]
}

#' A dataset's feature id -> gene symbol map (display only)
#'
#' Uses `dataset.symbol_col`, set explicitly in each config (no guessing);
#' warns and returns `NULL` when it's unset or names a missing column.
#'
#' @param dataset_yaml A resolved dataset config.
#' @param fm Optional already-loaded feature-metadata data frame.
#' @return Named character vector (names = feature ids), or `NULL`.
build_symbol_map <- function(dataset_yaml, fm = NULL) {
  ds <- dataset_yaml$dataset
  if (is.null(fm)) {
    if (is.null(ds$feature_metadata_path) || !file.exists(ds$feature_metadata_path)) return(NULL)
    fm <- arrow::read_parquet(ds$feature_metadata_path)
  }
  sym_col <- if (!is.null(ds$symbol_col) && ds$symbol_col %in% names(fm)) {
    ds$symbol_col
  } else {
    if (!is.null(ds$symbol_col)) {
      warning("dataset.symbol_col '", ds$symbol_col, "' (dataset '", ds$id %||% "?",
               "') not found in feature_metadata_path")
    } else {
      warning("dataset.symbol_col not set (dataset '", ds$id %||% "?", "')")
    }
    NA_character_
  }
  if (is.na(sym_col) || is.null(sym_col)) return(NULL)
  id_col <- ds$feature_id_col %||% "feature_id"
  if (!(id_col %in% names(fm))) return(NULL)
  map <- setNames(as.character(fm[[sym_col]]), as.character(fm[[id_col]]))
  map[!is.na(map) & nzchar(map)]
}

#' Remap row names to Ensembl gene ids
#'
#' Rows without an Ensembl id are dropped; the rest are renamed and sorted by
#' id. The map must be one-to-one (build_ensembl_map()), so no rows are
#' merged.
#'
#' @param mat Matrix with platform-native row names (features x samples, or
#'   features x factors).
#' @param ensembl_map build_ensembl_map() output.
#' @return The matrix with Ensembl gene id row names; `mat` unchanged when the
#'   map is empty or matches nothing.
remap_to_ensembl <- function(mat, ensembl_map) {
  if (length(ensembl_map) == 0) return(mat)
  ids <- ensembl_map[rownames(mat)]
  keep <- !is.na(ids)
  if (!any(keep)) return(mat)
  if (anyDuplicated(ids[keep])) {
    stop("remap_to_ensembl(): the map sends several rows to one gene; build it with ",
         "build_ensembl_map() from this dataset's matrix", call. = FALSE)
  }
  out <- mat[keep, , drop = FALSE]
  rownames(out) <- ids[keep]
  out[order(rownames(out)), , drop = FALSE]
}

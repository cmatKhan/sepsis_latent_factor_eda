# Gene identifiers across datasets (docs: Data, "Gene identifiers across
# datasets"):
#   - build_ensembl_map() / remap_to_ensembl(): version-stripped Ensembl gene
#     ids, the identifier for all cross-dataset computation (enrichment,
#     projectR). Symbols are ambiguous and renamed across annotation
#     releases; Ensembl ids aren't.
#   - build_symbol_map() / remap_to_symbol(): gene symbols, for display only.
# Maps are built once per dataset (the ensembl_map_<dataset> target) from
# its feature metadata.

#' Strip an Ensembl id's version suffix
#'
#' @param x Character vector of ids (`ENSG00000001234.2` -> `ENSG00000001234`;
#'   ids without a suffix are unchanged).
#' @return `x` without version suffixes.
strip_ensembl_version <- function(x) sub("\\.[0-9]+$", "", x)

#' A dataset's feature id -> Ensembl gene id map
#'
#' Reads `dataset.ensembl_col` (default `"ensembl"`) from the feature
#' metadata. A `;`-separated multi-gene entry (one probe, several genes) maps to
#' its first gene -- a simple, transparent choice; expanding one probe across
#' every listed gene is the alternative.
#'
#' @param dataset_yaml A resolved dataset config.
#' @param fm Optional already-loaded feature-metadata data frame.
#' @return Named character vector (names = feature ids), or `NULL` when the
#'   dataset has no such column or feature metadata -- its row names are then
#'   used as they are.
build_ensembl_map <- function(dataset_yaml, fm = NULL) {
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
  map[!is.na(map) & nzchar(map)]
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

#' Remap a matrix's row names through an id map
#'
#' Several rows mapping to one id are collapsed, keeping the row with the
#' largest |value| per column.
#'
#' @param mat Feature x sample (or feature x factor) matrix.
#' @param id_map Named character vector (names = current row names).
#' @return The remapped matrix; `mat` unchanged when `id_map` is empty or
#'   matches nothing.
.remap_ids <- function(mat, id_map) {
  if (is.null(id_map) || length(id_map) == 0) return(mat)
  ids <- id_map[rownames(mat)]
  keep <- !is.na(ids) & nzchar(ids)
  if (sum(keep) == 0) return(mat)
  m <- mat[keep, , drop = FALSE]
  ids <- ids[keep]
  collapsed <- apply(m, 2, function(col) tapply(col, ids, function(x) x[which.max(abs(x))]))
  as.matrix(collapsed)
}

#' Remap row names to Ensembl gene ids
#'
#' @param mat Matrix with platform-native row names.
#' @param ensembl_map build_ensembl_map() output.
#' @return The matrix with Ensembl gene id row names (see .remap_ids()).
remap_to_ensembl <- function(mat, ensembl_map) .remap_ids(mat, ensembl_map)

#' Remap row names to gene symbols (display only)
#'
#' @param mat Matrix with platform-native row names.
#' @param symbol_map build_symbol_map() output.
#' @return The matrix with gene-symbol row names (see .remap_ids()).
remap_to_symbol <- function(mat, symbol_map) .remap_ids(mat, symbol_map)

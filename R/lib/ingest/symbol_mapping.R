# Cross-dataset identifier mapping. Two DISTINCT purposes -- don't mix them:
#
#   - build_ensembl_map() / remap_to_ensembl(): the CANONICAL identifier
#     space for ALL cross-dataset COMPUTATIONAL work (fgsea_grid/
#     gprofiler_grid/projectr_within_grid/projectr_cross_grid). Ensembl
#     gene id, version suffix stripped (ENSG00000001234.2 ->
#     ENSG00000001234) -- chosen over gene symbol because symbols are
#     ambiguous/aliased/renamed across annotation releases and are
#     exact-string-matched here, both of which silently produce spurious
#     mismatches; Ensembl ids don't have that problem.
#   - build_symbol_map() / remap_to_symbol(): DISPLAY ONLY (e.g. showing
#     recognizable gene names alongside results in the app) -- NEVER used
#     for the actual cross-dataset matching/comparison itself.
#
# Built ONCE per dataset on the login node (from its feature_metadata.parquet)
# and baked into each job as plain named lists (`ensembl_maps`/`symbol_maps`,
# dataset_id -> named character vector) -- compute nodes never touch the DB
# or re-read feature_metadata themselves.

#' Strip an Ensembl id's version suffix (ENSG00000001234.2 ->
#' ENSG00000001234). A safe no-op on ids that already lack one.
strip_ensembl_version <- function(x) sub("\\.[0-9]+$", "", x)

#' One dataset's (feature_id -> Ensembl gene id) map, built from its
#' config -- THE canonical cross-dataset identifier for all computational
#' work (see this file's header). Uses `dataset.ensembl_col` (default
#' "ensembl", same convention as R/lib/matrices.R's preprocessing_script
#' contract). Returns NULL if there's no such column (or no
#' feature_metadata_path at all) -- callers should then treat the
#' matrix's own rownames as already Ensembl (or otherwise already a
#' shared id space).
#'
#' Some datasets' ensembl column holds ";"-delimited MULTI-gene mappings
#' (one probe/feature -> several Ensembl genes -- confirmed directly
#' against several of this project's array-platform datasets, e.g.
#' "ENSG00000204580;ENSG00000229767;..."). Rather than silently treat the
#' whole compound string as one (never-matching) id, this takes the FIRST
#' listed gene as a simple, transparent default -- a real analytical
#' choice, not obviously "correct" for every use case (an alternative
#' would be expanding one probe's value across every listed gene instead
#' of picking one) -- revisit if that matters for your analysis.
#'
#' @param fm optional pre-loaded feature_metadata data.frame (e.g. from
#'   the cached artifact -- see R/lib/ingest/ingest_dataset.R::
#'   cache_dataset_metadata()) -- skips reading feature_metadata_path
#'   entirely when supplied. Needed wherever this is called from a
#'   machine that doesn't have raw access to that path (e.g.
#'   R/create_ingest_slurm_bundle.R, typically invoked on the cluster
#'   login node, not wherever the raw HuggingFace data lives) --
#'   confirmed directly (2026-09-18) that omitting this on such a machine
#'   silently returns NULL for every dataset (file.exists() on the raw
#'   path is always FALSE there), not an error.
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

#' One dataset's (feature_id -> symbol) map, built from its config --
#' DISPLAY ONLY (see this file's header; never used for cross-dataset
#' computational matching -- see build_ensembl_map() for that).
#'
#' Column choice: `dataset.symbol_col`, set EXPLICITLY in every dataset's
#' config (see config/dataset_metadata.example.yml) -- naming exactly
#' which feature_metadata column holds the gene symbol. Deliberately NOT
#' auto-detected/guessed from a fixed list of likely column names (an
#' earlier version of this function did that) -- explicit config is more
#' reliable than a guess, and this project's configs all set it
#' explicitly now anyway. Warns (and returns NULL) if `symbol_col` is
#' unset or names a column that doesn't actually exist in
#' feature_metadata_path.
#'
#' Returns NULL (no map at all) if there's genuinely no usable symbol
#' column -- callers should then leave the matrix's own rownames as-is
#' for display.
#'
#' @param fm optional pre-loaded feature_metadata data.frame -- see
#'   build_ensembl_map()'s matching doc for why this exists.
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

#' Shared collapse-by-id logic for remap_to_ensembl()/remap_to_symbol():
#' collapses many-probes-to-one-id by keeping the max-|value| row per id
#' per column. Returns `mat` unchanged if `id_map` is NULL/empty or
#' nothing matches -- assumes rownames are already a usable shared id in
#' that case.
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

#' Remap a matrix's rownames (platform-native ids) to Ensembl gene id
#' (version-stripped) via a prebuilt map (see build_ensembl_map()) --
#' THE canonical remap for all cross-dataset computational work.
remap_to_ensembl <- function(mat, ensembl_map) .remap_ids(mat, ensembl_map)

#' Remap a matrix's rownames to gene symbol -- DISPLAY ONLY, never used
#' for cross-dataset computational matching (see remap_to_ensembl()).
remap_to_symbol <- function(mat, symbol_map) .remap_ids(mat, symbol_map)

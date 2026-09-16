# Gene-symbol mapping used by fgsea_job.R/projectr_job.R -- both need a
# COMMON gene identifier across datasets (loadings/matrices are otherwise
# keyed by each dataset's native platform id: Illumina/Affymetrix probe,
# Entrez GeneID, Ensembl feature_id) and msigdbr's gene sets are
# symbol-keyed. Built ONCE per dataset on the login node (from its
# feature_metadata.parquet) and baked into each job as a plain named list
# (`symbol_maps`, dataset_id -> named character vector) -- compute nodes
# never touch the DB or re-read feature_metadata themselves.

#' One dataset's (feature_id -> symbol) map, built from its config.
#' Tries a prioritized list of likely symbol-column names; if none exists
#' (e.g. GSE110487's GeneID is already Entrez), returns NULL and callers
#' should treat the matrix's own rownames as already-usable ids.
build_symbol_map <- function(dataset_yaml) {
  ds <- dataset_yaml$dataset
  if (is.null(ds$feature_metadata_path) || !file.exists(ds$feature_metadata_path)) return(NULL)
  fm <- arrow::read_parquet(ds$feature_metadata_path)
  sym_col <- intersect(c("symbol", "gene_symbol", "SYMBOL", "Symbol"), names(fm))[1]
  if (is.na(sym_col) || is.null(sym_col)) return(NULL)
  id_col <- ds$feature_id_col %||% "feature_id"
  if (!(id_col %in% names(fm))) return(NULL)
  map <- setNames(as.character(fm[[sym_col]]), as.character(fm[[id_col]]))
  map[!is.na(map) & nzchar(map)]
}

#' Remap a matrix's rownames (platform-native ids) to gene symbol via a
#' prebuilt map (see build_symbol_map()) -- collapses many-probes-to-one-
#' symbol by keeping the max-|value| row per symbol per column. Returns
#' `mat` unchanged (with a message) if `symbol_map` is NULL/empty or
#' nothing matches -- assumes rownames are already a usable id (Entrez/
#' Ensembl) in that case.
remap_to_symbol <- function(mat, symbol_map) {
  if (is.null(symbol_map) || length(symbol_map) == 0) return(mat)
  sym <- symbol_map[rownames(mat)]
  keep <- !is.na(sym) & nzchar(sym)
  if (sum(keep) == 0) return(mat)
  m <- mat[keep, , drop = FALSE]
  sym <- sym[keep]
  collapsed <- apply(m, 2, function(col) tapply(col, sym, function(x) x[which.max(abs(x))]))
  as.matrix(collapsed)
}

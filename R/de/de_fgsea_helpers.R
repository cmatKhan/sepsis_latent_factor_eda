# fgsea helpers for the DE contrasts under R/de/ -- mirrors
# R/ingest_jobs/fgsea_job.R's conventions (same msigdbr collections, same
# fgsea() parameters, same canonical Ensembl identifier space) rather than
# inventing new ones, but runs locally/standalone (no cluster, no DB
# writes) -- consistent with the rest of R/de/. Needs R/lib/ingest/
# symbol_mapping.R sourced first (build_ensembl_map()/remap_to_ensembl()/
# .remap_ids()).

suppressMessages({
  library(fgsea)
})

#' The MSigDB collections for DE enrichment
#'
#' The same six collections as the pipeline (fetch_msigdb_pathways()), Ensembl
#' keyed, fetched once per driver run.
#'
#' @return Named list by collection of gene sets.
fetch_de_pathways <- function() {
  if (!requireNamespace("msigdbr", quietly = TRUE)) stop("Package 'msigdbr' is required")
  fetch_msig <- function(collection, subcollection = NULL) {
    msig <- msigdbr::msigdbr(species = "Homo sapiens", collection = collection, subcollection = subcollection)
    msig <- msig[!is.na(msig$ensembl_gene) & nzchar(msig$ensembl_gene), ]
    split(msig$ensembl_gene, msig$gs_name)
  }
  hallmark <- fetch_msig("H")
  list(
    HALLMARK = hallmark,
    `GO:BP`  = fetch_msig("C5", "GO:BP"),
    `GO:MF`  = fetch_msig("C5", "GO:MF"),
    KEGG     = fetch_msig("C2", "CP:KEGG_LEGACY"),
    REAC     = fetch_msig("C2", "CP:REACTOME"),
    WP       = fetch_msig("C2", "CP:WIKIPATHWAYS")
  )
}

#' Gene id -> Ensembl map for a DE dataset's results
#'
#' RNA-seq results keep native feature ids, so build_ensembl_map() applies;
#' array results were collapsed to symbols, so this maps each symbol to the
#' first non-missing Ensembl id among its probes.
#'
#' @param ds_meta The dataset's `dataset:` block.
#' @param feature_meta Feature-metadata data frame.
#' @param platform `"array"` or `"rnaseq"`.
#' @return Named character vector, or `NULL`.
build_de_ensembl_map <- function(ds_meta, feature_meta, platform = c("array", "rnaseq")) {
  platform <- match.arg(platform)
  if (platform == "rnaseq") {
    return(build_ensembl_map(list(dataset = ds_meta), fm = feature_meta))
  }
  sym_col <- ds_meta$symbol_col
  ens_col <- ds_meta$ensembl_col %||% "ensembl"
  if (is.null(sym_col) || !(sym_col %in% names(feature_meta)) || !(ens_col %in% names(feature_meta))) {
    return(NULL)
  }
  fm <- feature_meta[!is.na(feature_meta[[sym_col]]) & nzchar(feature_meta[[sym_col]]) &
                        !is.na(feature_meta[[ens_col]]) & nzchar(feature_meta[[ens_col]]), ]
  first_id <- vapply(strsplit(as.character(fm[[ens_col]]), ";", fixed = TRUE),
                      function(x) if (length(x) == 0) NA_character_ else x[1], character(1))
  fm$.ens <- strip_ensembl_version(first_id)
  fm <- fm[!duplicated(fm[[sym_col]]), ]
  setNames(fm$.ens, as.character(fm[[sym_col]]))
}

#' Remap a named numeric vector to Ensembl ids
#'
#' Keeps the max-|value| entry per Ensembl id (the vector analogue of
#' .remap_ids()).
#'
#' @param v Named numeric vector.
#' @param ens_map Id -> Ensembl map, or `NULL` to leave names as they are.
#' @return The remapped vector, sorted decreasing.
remap_named_vector_to_ensembl <- function(v, ens_map) {
  if (is.null(ens_map) || length(ens_map) == 0) return(sort(v, decreasing = TRUE))
  ids <- ens_map[names(v)]
  keep <- !is.na(ids) & nzchar(ids)
  v <- v[keep]; ids <- ids[keep]
  v <- tapply(v, ids, function(x) x[which.max(abs(x))])
  sort(v, decreasing = TRUE)
}

#' fgsea of one contrast against one collection
#'
#' Same parameters as the pipeline (minSize 10, maxSize 500).
#'
#' @param rank_vector Ensembl-keyed signed statistics, sorted decreasing.
#' @param pathways One collection's gene sets.
#' @return fgsea() result table.
run_de_fgsea <- function(rank_vector, pathways) {
  tryCatch(
    fgsea::fgsea(pathways = pathways, stats = rank_vector, minSize = 10, maxSize = 500),
    error = function(e) NULL
  )
}

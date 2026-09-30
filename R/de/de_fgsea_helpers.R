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

#' Same 6 msigdbr collections as R/create_ingest_slurm_bundle.R's
#' `pathways_by_source` (species = "Homo sapiens", ensembl_gene-keyed) --
#' fetched ONCE per driver run, not per dataset/contrast (msigdbr's own
#' query is the slow part here).
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

#' Build a (gene id -> Ensembl gene id) map for one DE dataset's topTable
#' gene identifiers -- dispatches on platform, since R/de/*.R's gene
#' column convention differs by platform (see R/de/README.md):
#'   - RNA-seq datasets never collapse (R/de/de_helpers.R's
#'     collapse_to_symbol() is array-only) -- their `gene` column already
#'     equals that dataset's `feature_id_col`, so the existing
#'     R/lib/ingest/symbol_mapping.R::build_ensembl_map() (feature_id ->
#'     Ensembl, via `ensembl_col`) applies directly.
#'   - Array datasets collapsed to gene SYMBOL before fitting -- there is
#'     no existing symbol -> Ensembl map (build_ensembl_map() is
#'     feature_id/probe -keyed), so this builds one directly from
#'     feature_meta: for each symbol, the first non-missing `ensembl_col`
#'     value among the probes sharing that symbol. Deliberately simple --
#'     doesn't try to reconstruct which exact probe
#'     collapse_to_symbol() itself picked (picking a DIFFERENT probe's
#'     Ensembl annotation for the same symbol essentially never matters
#'     in practice; genuinely conflicting Ensembl annotations across
#'     probes for one symbol would indicate a feature_metadata quality
#'     issue this first pass isn't trying to solve).
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

#' Remap a NAMED numeric vector's names (platform-native gene ids) to
#' Ensembl via a prebuilt map -- the named-vector analogue of
#' R/lib/ingest/symbol_mapping.R's `.remap_ids()` (which operates on a
#' matrix's rownames instead). Same collapse rule for a many-to-one
#' remap: keep the max-|value| entry per Ensembl id. Returns a vector
#' sorted decreasing (fgsea()'s expected input shape), or the vector
#' unchanged (just re-sorted) if `ens_map` is NULL/empty.
remap_named_vector_to_ensembl <- function(v, ens_map) {
  if (is.null(ens_map) || length(ens_map) == 0) return(sort(v, decreasing = TRUE))
  ids <- ens_map[names(v)]
  keep <- !is.na(ids) & nzchar(ids)
  v <- v[keep]; ids <- ids[keep]
  v <- tapply(v, ids, function(x) x[which.max(abs(x))])
  sort(v, decreasing = TRUE)
}

#' fgsea::fgsea() with the same parameters as the factor pipeline
#' (minSize=10, maxSize=500) -- one call, one pathway source. `rank_vector`
#' must already be Ensembl-keyed and sorted decreasing (ties handled by
#' fgsea() itself, same as R/ingest_jobs/fgsea_job.R -- no pre-jittering).
run_de_fgsea <- function(rank_vector, pathways) {
  tryCatch(
    fgsea::fgsea(pathways = pathways, stats = rank_vector, minSize = 10, maxSize = 500),
    error = function(e) NULL
  )
}

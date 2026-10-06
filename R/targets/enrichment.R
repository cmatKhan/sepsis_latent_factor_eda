# Enrichment as targets: fgsea (Hallmark GSEA, CoGAPS marker ORA, local
# GSEA/ORA over MSigDB collections) per representative fit, and module ORA
# per WGCNA fit. The job functions are R/ingest_jobs/fgsea_job.R and
# wgcna_ora_job.R. Each branch keeps its job's result trimmed to the
# significant terms (trim_enrichment_result()); enrichment_rows() turns
# those into `enrichment_queries` / `enrichment_results` rows, so changing
# the DB layout never reruns the jobs.

#' MSigDB gene sets for enrichment, in Ensembl gene ids
#'
#' Built once on the main process (the `msigdb_pathways` target; msigdbr may
#' download its data).
#'
#' @return Named list by collection (`HALLMARK`, `GO:BP`, `GO:MF`, `KEGG`,
#'   `REAC`, `WP`), each a named list of gene sets.
fetch_msigdb_pathways <- function() {
  fetch <- function(collection, subcollection = NULL) {
    msig <- msigdbr::msigdbr(species = "Homo sapiens", collection = collection, subcollection = subcollection)
    msig <- msig[!is.na(msig$ensembl_gene) & nzchar(msig$ensembl_gene), ]
    split(msig$ensembl_gene, msig$gs_name)
  }
  list(
    HALLMARK = fetch("H"),
    `GO:BP` = fetch("C5", "GO:BP"),
    `GO:MF` = fetch("C5", "GO:MF"),
    KEGG = fetch("C2", "CP:KEGG_LEGACY"),
    REAC = fetch("C2", "CP:REACTOME"),
    WP = fetch("C2", "CP:WIKIPATHWAYS")
  )
}

#' Keep a job table branchable
#'
#' targets can't branch over an empty table, so an empty job table becomes
#' one placeholder row (`col` = `NA`); its branch returns `NULL` (has_job()),
#' which the row builders and writers skip.
#'
#' @param jobs Job table, or `NULL`.
#' @param col Key column to create in the placeholder.
#' @return `jobs`, or a one-row placeholder if it had no rows.
placeholder_jobs <- function(jobs, col) {
  if (!is.null(jobs) && nrow(jobs) > 0) return(jobs)
  stats::setNames(data.frame(NA_character_), col)
}

#' Whether a branch's job is real
#'
#' @param x The job's key value.
#' @return `TRUE` unless `x` is the placeholder (`NA`).
has_job <- function(x) length(x) == 1 && !is.na(x)

#' fgsea jobs for a method's representative fits
#'
#' @param ingest An ingest_method_fits() result.
#' @param representatives Representative fit_keys.
#' @param redundancy compute_method_redundancy() output (CoGAPS marker genes).
#' @param db_path Path to the DB.
#' @return Job table (fit_key, loadings_file, cogaps_marker_genes), one row
#'   per fit, or a placeholder.
fgsea_jobs <- function(ingest, representatives, redundancy, db_path) {
  if (length(representatives) == 0) return(placeholder_jobs(NULL, "fit_key"))
  files <- vapply(representatives, fit_artifact, character(1), ingest = ingest, kind = "loadings",
                  db_path = db_path)
  keys <- representatives[!is.na(files)]
  if (length(keys) == 0) return(placeholder_jobs(NULL, "fit_key"))
  jobs <- data.frame(fit_key = keys, loadings_file = unname(files[!is.na(files)]))
  jobs$cogaps_marker_genes <- lapply(keys, function(key) {
    mk <- redundancy$pattern_markers
    if (ingest$method != "cogaps" || is.null(mk)) return(NULL)
    mk <- mk[mk$fit_key == key, , drop = FALSE]
    if (nrow(mk) == 0) NULL else split(mk$gene, mk$factor_index)
  })
  jobs
}

#' Module-ORA jobs for every WGCNA fit with modules
#'
#' @param ingest The WGCNA ingest_method_fits() result.
#' @return Job table (fit_key, module_genes, universe_genes), or a placeholder.
wgcna_ora_jobs <- function(ingest) {
  mods <- ingest$modules
  if (is.null(mods)) return(placeholder_jobs(NULL, "fit_key"))
  keys <- unique(mods$fit_key)
  jobs <- data.frame(fit_key = keys)
  jobs$module_genes <- lapply(keys, function(k) {
    m <- mods[mods$fit_key == k & mods$module != 0, , drop = FALSE]
    split(m$gene, m$module)
  })
  jobs$universe_genes <- lapply(keys, function(k) unique(mods$gene[mods$fit_key == k]))
  placeholder_jobs(jobs[lengths(jobs$module_genes) > 0, , drop = FALSE], "fit_key")
}

ENRICHMENT_PADJ <- 0.05

#' Cut a job result to its significant terms
#'
#' What each enrichment branch stores (padj < `ENRICHMENT_PADJ`), so rows can
#' be rebuilt without rerunning the job.
#'
#' @param x A run_fgsea_job() or run_wgcna_ora_job() result.
#' @return `x` with every fgsea()/fora() table filtered.
trim_enrichment_result <- function(x) {
  cut <- function(r) {
    if (is.null(r$result) || nrow(r$result) == 0) return(r)
    r$result <- r$result[!is.na(r$result$padj) & r$result$padj < ENRICHMENT_PADJ, ]
    r
  }
  for (part in c("gsea", "local_gsea", "local_ora", "ora")) {
    if (!is.null(x[[part]])) x[[part]] <- lapply(x[[part]], cut)
  }
  if (!is.null(x$fora)) {
    x$fora <- lapply(x$fora, function(r) {
      if (is.null(r) || nrow(r) == 0) r else r[!is.na(r$padj) & r$padj < ENRICHMENT_PADJ, ]
    })
  }
  x
}

# ---- trimmed results -> DB rows ------------------------------------------
# `queries`: one row per (fit_key, factor_index, query_type, direction) that
# was run, even with no significant term (the app distinguishes "no hits"
# from "not run"). `results`: significant terms, linked to their query by
# the same four columns. GSEA queries are direction "both" -- a hit's sign
# is its NES. p_value holds fgsea's BH-adjusted p (padj).

#' Gene lists as JSON arrays
#'
#' @param x List of character vectors.
#' @return Character vector of JSON arrays (`enrichment_results.genes`).
json_genes <- function(x) vapply(x, function(g) as.character(jsonlite::toJSON(as.character(g))), character(1))

#' One `enrichment_queries` row, keyed
#'
#' @param key fit_key.
#' @param fi factor_index.
#' @param query_type Query type.
#' @param direction `"pos"`, `"neg"` or `"both"`.
#' @param query_size Genes in the query, or `NA`.
#' @return One-row data frame.
query_row <- function(key, fi, query_type, direction, query_size = NA_integer_) {
  data.frame(fit_key = key, factor_index = as.integer(fi), query_type = query_type,
             direction = direction, query_size = as.integer(query_size))
}

#' GSEA hits as `enrichment_results` rows
#'
#' @param res A trimmed fgsea() table.
#' @param key fit_key.
#' @param fi factor_index.
#' @param query_type Query type.
#' @param source MSigDB collection.
#' @param main_pathways Pathways collapsePathways() kept.
#' @return Data frame keyed by (fit_key, factor_index, query_type,
#'   direction = "both"), or `NULL`.
gsea_results <- function(res, key, fi, query_type, source, main_pathways) {
  if (is.null(res) || nrow(res) == 0) return(NULL)
  data.frame(fit_key = key, factor_index = as.integer(fi), query_type = query_type, direction = "both",
             term_id = res$pathway, source = source, term_name = res$pathway, p_value = res$padj,
             intersection_size = lengths(res$leadingEdge), term_size = res$size,
             nes = res$NES, es = res$ES, log2err = res$log2err,
             is_main_pathway = as.integer(res$pathway %in% (main_pathways %||% character(0))),
             genes = json_genes(res$leadingEdge))
}

#' ORA hits as `enrichment_results` rows
#'
#' @param res A trimmed fora() table.
#' @param key fit_key.
#' @param fi factor_index.
#' @param query_type Query type.
#' @param direction Gene direction tested.
#' @param source MSigDB collection.
#' @return Data frame, or `NULL`.
ora_results <- function(res, key, fi, query_type, direction, source) {
  if (is.null(res) || nrow(res) == 0) return(NULL)
  data.frame(fit_key = key, factor_index = as.integer(fi), query_type = query_type, direction = direction,
             term_id = res$pathway, source = source, term_name = res$pathway, p_value = res$padj,
             intersection_size = res$overlap, term_size = res$size,
             nes = NA_real_, es = NA_real_, log2err = NA_real_, is_main_pathway = NA_integer_,
             genes = json_genes(res$overlapGenes))
}

#' Rows from one trimmed run_fgsea_job() result
#'
#' @param x Trimmed job result.
#' @return `list(queries, results)`: Hallmark GSEA (`fgsea`), CoGAPS marker
#'   ORA (`cogaps_fora`), GSEA against the other collections (`gsea`) and
#'   top-gene ORA (`ora`).
fgsea_rows <- function(x) {
  key <- x$fit_key
  q <- list(); r <- list()
  for (g in x$gsea %||% list()) {
    q[[length(q) + 1]] <- query_row(key, g$factor_index, "fgsea", "both")
    r[[length(r) + 1]] <- gsea_results(g$result, key, g$factor_index, "fgsea", "HALLMARK", g$main_pathways)
  }
  for (fi_chr in names(x$fora)) {
    q[[length(q) + 1]] <- query_row(key, fi_chr, "cogaps_fora", "pos")
    r[[length(r) + 1]] <- ora_results(x$fora[[fi_chr]], key, fi_chr, "cogaps_fora", "pos", "HALLMARK")
  }
  for (g in x$local_gsea %||% list()) {
    q[[length(q) + 1]] <- query_row(key, g$factor_index, "gsea", "both")
    r[[length(r) + 1]] <- gsea_results(g$result, key, g$factor_index, "gsea", g$source, g$main_pathways)
  }
  for (o in x$local_ora %||% list()) {
    q[[length(q) + 1]] <- query_row(key, o$factor_index, "ora", o$direction, o$n_genes %||% NA)
    r[[length(r) + 1]] <- ora_results(o$result, key, o$factor_index, "ora", o$direction, o$source)
  }
  list(queries = do.call(rbind, q), results = do.call(rbind, r))
}

#' Rows from one trimmed run_wgcna_ora_job() result
#'
#' @param x Trimmed job result.
#' @return `list(queries, results)`, factor_index = module.
wgcna_ora_rows <- function(x) {
  key <- x$fit_key
  q <- list(); r <- list()
  for (o in x$ora %||% list()) {
    q[[length(q) + 1]] <- query_row(key, o$module, "ora", "pos", o$n_genes %||% NA)
    r[[length(r) + 1]] <- ora_results(o$result, key, o$module, "ora", "pos", o$source)
  }
  list(queries = do.call(rbind, q), results = do.call(rbind, r))
}

#' All of a dataset's enrichment rows (the `enrichment_<dataset>` target)
#'
#' @param fgsea_results Flat list of trimmed fgsea branch values.
#' @param wgcna_results Flat list of trimmed WGCNA ORA branch values.
#' @return `list(queries, results)`, queries unique per (fit, factor, type,
#'   direction) and results per (query, term), keeping the most significant
#'   hit when collections overlap.
enrichment_rows <- function(fgsea_results = list(), wgcna_results = list()) {
  sets <- c(lapply(Filter(Negate(is.null), fgsea_results), fgsea_rows),
            lapply(Filter(Negate(is.null), wgcna_results), wgcna_ora_rows))
  queries <- do.call(rbind, lapply(sets, `[[`, "queries"))
  results <- do.call(rbind, lapply(sets, `[[`, "results"))
  qcols <- c("fit_key", "factor_index", "query_type", "direction")
  if (!is.null(queries)) {
    queries <- queries[order(is.na(queries$query_size)), , drop = FALSE]
    queries <- queries[!duplicated(queries[, qcols]), , drop = FALSE]
  }
  if (!is.null(results)) {
    results <- results[order(results$p_value), , drop = FALSE]
    results <- results[!duplicated(results[, c(qcols, "term_id")]), , drop = FALSE]
  }
  list(queries = queries, results = results)
}

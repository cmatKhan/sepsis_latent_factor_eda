# Standalone, DBI-only batch functional-enrichment runner for
# R/ingest_results.R's optional `--run-enrichment` flag.
#
# Deliberately reimplements (rather than sources from) the app's
# equivalent logic in app/app.R / app/R/db_helpers.R -- app.R's own
# header states the app is "fully decoupled from ingest", and this keeps
# that boundary intact rather than having the ingest CLI reach into
# app/R/*. The duplication is intentional; keep the two in sync manually
# if the enrichment-caching scheme (enrichment_cache/enrichment_queried
# tables, factor_id keying) ever changes.
#
# Scope: every dataset x method x `status = 'ok'` fit whose family is
# `seed_sweep` or `param_grid` (i.e. an actual explorable fit, not a
# `maskcv` fit, which has no loadings) x every factor, plus every WGCNA
# module -- across the WHOLE db (every dataset previously ingested, not
# just the one just processed), per the "for all datasets and options"
# request this was built for.

library(DBI)
source(here::here("R/lib/ingest/symbol_mapping.R"))

#' get_fit()/get_factor_id()/enrichment_cached()/enrichment_store() below
#' mirror app/R/db_helpers.R's functions of the same name exactly (same
#' schema, same queries) -- see that file's comments for the schema
#' rationale.
ib_get_fit <- function(con, fit_id) {
  DBI::dbGetQuery(con, "SELECT * FROM fits WHERE fit_id = ?", params = list(fit_id))
}

ib_get_factor_id <- function(con, fit_id, factor_index) {
  DBI::dbGetQuery(con,
    "SELECT factor_id FROM factors WHERE fit_id = ? AND factor_index = ?",
    params = list(fit_id, factor_index))$factor_id
}

ib_load_loadings <- function(con, db_path, fit_id) {
  f <- ib_get_fit(con, fit_id)
  if (nrow(f) == 0 || is.na(f$loadings_file)) return(NULL)
  path <- resolve_artifact(f$loadings_file, db_path)
  if (!file.exists(path)) return(NULL)
  as.matrix(readRDS(path))
}

ib_enrichment_cached <- function(con, factor_id, query_type, direction) {
  hit <- DBI::dbGetQuery(con,
    "SELECT COUNT(*) AS n FROM enrichment_queried
     WHERE factor_id = ? AND query_type = ? AND direction = ?",
    params = list(factor_id, query_type, direction))$n > 0
  if (!hit) return(NULL)
  DBI::dbGetQuery(con,
    "SELECT source, term_id, term_name, p_value, intersection_size, term_size, query_size
     FROM enrichment_cache
     WHERE factor_id = ? AND query_type = ? AND direction = ?
     ORDER BY p_value",
    params = list(factor_id, query_type, direction))
}

ib_enrichment_store <- function(con, factor_id, query_type, direction, gost_result) {
  DBI::dbExecute(con,
    "INSERT OR REPLACE INTO enrichment_queried (factor_id, query_type, direction, queried_at)
     VALUES (?, ?, ?, datetime('now'))",
    params = list(factor_id, query_type, direction))
  if (!is.null(gost_result) && !is.null(gost_result$result) && nrow(gost_result$result) > 0) {
    r <- gost_result$result
    DBI::dbWriteTable(con, "enrichment_cache", data.frame(
      factor_id = factor_id, query_type = query_type, direction = direction,
      source = r$source, term_id = r$term_id, term_name = r$term_name,
      p_value = r$p_value, intersection_size = r$intersection_size,
      term_size = r$term_size,
      query_size = if (!is.null(r$query_size)) r$query_size else NA_integer_,
      genes = if (!is.null(r$intersection)) as.character(r$intersection) else NA_character_,
      queried_at = as.character(Sys.time())
    ), append = TRUE)
  }
  invisible(NULL)
}

#' Copies another factor's already-cached rows onto `factor_id` (marking
#' it queried) without calling g:Profiler -- see ib_find_or_reuse() below.
ib_copy_enrichment <- function(con, factor_id, query_type, direction, cached_rows) {
  DBI::dbExecute(con,
    "INSERT OR REPLACE INTO enrichment_queried (factor_id, query_type, direction, queried_at)
     VALUES (?, ?, ?, datetime('now'))",
    params = list(factor_id, query_type, direction))
  if (nrow(cached_rows) > 0) {
    DBI::dbWriteTable(con, "enrichment_cache", data.frame(
      factor_id = factor_id, query_type = query_type, direction = direction,
      source = cached_rows$source, term_id = cached_rows$term_id, term_name = cached_rows$term_name,
      p_value = cached_rows$p_value, intersection_size = cached_rows$intersection_size,
      term_size = cached_rows$term_size, genes = NA_character_,
      queried_at = as.character(Sys.time()), query_size = cached_rows$query_size
    ), append = TRUE)
  }
  invisible(NULL)
}

#' Cache lookup that also recognizes a numerically-identical loading
#' vector already queried under a DIFFERENT fit_id -- PCA gets a
#' distinct fit_id per rank, but prcomp(rank. = k) just truncates the
#' same underlying SVD (see R/methods/pca.R), so PC 3 at rank 7 and PC 3
#' at rank 10 are byte-identical. Restricted to PCA, the only method with
#' this guarantee (see app/app.R's identical find_or_reuse_enrichment()
#' for the same rationale).
ib_find_or_reuse <- function(con, db_path, fit_id, factor_index, qtype, direction) {
  factor_id <- ib_get_factor_id(con, fit_id, factor_index)
  cached <- ib_enrichment_cached(con, factor_id, qtype, direction)
  if (!is.null(cached)) return(cached)

  fit <- ib_get_fit(con, fit_id)
  if (nrow(fit) == 0 || fit$method != "pca") return(NULL)

  v <- ib_load_loadings(con, db_path, fit_id)
  if (is.null(v)) return(NULL)
  v <- v[, factor_index]

  siblings <- DBI::dbGetQuery(con,
    "SELECT fit_id FROM fits
     WHERE dataset_id = ? AND method = 'pca' AND family = 'seed_sweep'
       AND fit_id != ? AND status = 'ok'",
    params = list(fit$dataset_id, fit_id))$fit_id
  for (other_fit in siblings) {
    L2 <- ib_load_loadings(con, db_path, other_fit)
    if (is.null(L2)) next
    for (fi2 in seq_len(ncol(L2))) {
      v2 <- L2[, fi2]
      if (length(v) != length(v2) || !identical(names(v), names(v2))) next
      if (!isTRUE(all.equal(as.numeric(v), as.numeric(v2), tolerance = 1e-8))) next
      other_factor_id <- ib_get_factor_id(con, other_fit, fi2)
      hit <- ib_enrichment_cached(con, other_factor_id, qtype, direction)
      if (!is.null(hit)) {
        ib_copy_enrichment(con, factor_id, qtype, direction, hit)
        return(hit)
      }
    }
  }
  NULL
}

#' Runs (or reuses) enrichment for one factor's (qtype, direction)
#' combination and stores the result, unless already cached. `v` must
#' already be Ensembl-remapped by the caller (see run_all_enrichment()) --
#' THE canonical cross-dataset identifier for enrichment queries, same as
#' the slurm fgsea_grid/gprofiler_grid pipeline (R/ingest_jobs/
#' fgsea_job.R, gprofiler_job.R) -- not left as native platform ids, which
#' gprofiler2::gost() often can't recognize at all (e.g. raw microarray
#' probe ids).
ib_run_one <- function(con, db_path, fit_id, factor_index, v, qtype, direction, topn) {
  if (!is.null(ib_find_or_reuse(con, db_path, fit_id, factor_index, qtype, direction))) return(invisible(NULL))
  factor_id <- ib_get_factor_id(con, fit_id, factor_index)
  genes <- if (qtype == "ora") {
    if (direction == "neg") names(sort(v))[seq_len(min(topn, length(v)))]
    else names(sort(v, decreasing = TRUE))[seq_len(min(topn, length(v)))]
  } else {
    if (direction == "neg") names(sort(v)) else names(sort(v, decreasing = TRUE))
  }
  res <- tryCatch(
    gprofiler2::gost(
      query = genes, organism = "hsapiens", significant = TRUE,
      ordered_query = (qtype == "gsea"), correction_method = "fdr",
      sources = c("GO:BP", "GO:MF", "REAC", "KEGG", "WP")),
    error = function(e) {
      message("  g:Profiler query failed for fit ", fit_id, " factor ", factor_index,
              " (", qtype, "/", direction, "): ", conditionMessage(e))
      NULL
    })
  ib_enrichment_store(con, factor_id, qtype, direction, res)
  invisible(NULL)
}

#' Main entry point: batch-runs functional enrichment for every ok fit's
#' every factor across the WHOLE db (every dataset previously ingested),
#' skipping anything already cached (or reusable via the PCA cross-rank
#' shortcut above). `query_types` defaults to both ORA and GSEA, same
#' combinations a single factor's Level 3 Enrichment tab in the app can
#' otherwise trigger one at a time. `topn` is the ORA gene-list size
#' (irrelevant to GSEA, which always uses the full ranked list).
run_all_enrichment <- function(con, db_path, query_types = c("ora", "gsea"), topn = 100) {
  if (!requireNamespace("gprofiler2", quietly = TRUE)) {
    stop("Package 'gprofiler2' is required for --run-enrichment")
  }

  fits <- DBI::dbGetQuery(con,
    "SELECT fit_id, dataset_id, method, rank, n_factors FROM fits
     WHERE status = 'ok' AND family IN ('seed_sweep', 'param_grid')
       AND method != 'wgcna' AND n_factors IS NOT NULL AND n_factors > 0
     ORDER BY dataset_id, method, fit_id")

  # Ensembl maps -- THE canonical cross-dataset identifier for enrichment
  # queries (see R/lib/ingest/symbol_mapping.R's header), built directly
  # from config/*.yml here since this CLI path runs on the login node
  # (same assumption cache_dataset_matrix() makes -- see R/ingest_results.R's
  # header) and has real filesystem access to those configs, unlike the app
  # (see app/R/metadata_helpers.R::ensembl_map_for_dataset(), which instead
  # reads dataset_metadata_sources.ensembl_col to stay decoupled from ingest).
  all_ds_ids <- unique(c(fits$dataset_id, DBI::dbGetQuery(con,
    "SELECT DISTINCT dataset_id FROM fits WHERE status = 'ok' AND method = 'wgcna'")$dataset_id))
  available_cfg <- list.files("config", pattern = "_config\\.yml$")
  # Prefers the CACHED feature_metadata artifact (datasets.
  # feature_metadata_file, populated by cache_dataset_metadata()) over a
  # raw feature_metadata_path read -- same rationale as
  # R/create_ingest_slurm_bundle.R's ensembl_maps construction: this CLI's
  # own header says it should run wherever the raw data resolves, but
  # falling back gracefully here means it still works even when that
  # assumption doesn't hold (e.g. run from the cluster after the cache was
  # already synced over).
  ensembl_maps <- setNames(lapply(all_ds_ids, function(id) {
    match_idx <- match(tolower(paste0(id, "_config.yml")), tolower(available_cfg))
    if (is.na(match_idx)) return(NULL)
    fm <- NULL
    f <- DBI::dbGetQuery(con, "SELECT feature_metadata_file FROM datasets WHERE dataset_id = ?",
                          params = list(id))$feature_metadata_file
    if (length(f) == 1 && !is.na(f)) {
      path <- resolve_artifact(f, db_path)
      if (file.exists(path)) fm <- tryCatch(readRDS(path), error = function(e) NULL)
    }
    tryCatch(build_ensembl_map(yaml::read_yaml(file.path("config", available_cfg[match_idx])), fm = fm),
             error = function(e) NULL)
  }), all_ds_ids)

  n_fits <- nrow(fits)
  message("run_all_enrichment(): ", n_fits, " non-WGCNA ok fit(s) to consider")
  for (i in seq_len(n_fits)) {
    fit_id <- fits$fit_id[i]
    message("[", i, "/", n_fits, "] ", fits$dataset_id[i], " / ", fits$method[i],
            " / fit_id ", fit_id, " (", fits$n_factors[i], " factors)")
    L <- ib_load_loadings(con, db_path, fit_id)
    if (is.null(L)) {
      message("  no loadings artifact -- skipping")
      next
    }
    L_ens <- remap_to_ensembl(L, ensembl_maps[[fits$dataset_id[i]]])
    dirs <- if (fits$method[i] == "pca") c("pos", "neg") else "pos"
    for (fi in seq_len(ncol(L_ens))) {
      v <- L_ens[, fi]
      for (dir in dirs) {
        for (qt in query_types) {
          ib_run_one(con, db_path, fit_id, fi, v, qt, dir, topn)
        }
      }
    }
  }

  # WGCNA: ORA only, whole module membership (no topN/direction), same
  # treatment as the app's Level 3 module view.
  wgcna_fits <- DBI::dbGetQuery(con,
    "SELECT fit_id, dataset_id FROM fits WHERE status = 'ok' AND method = 'wgcna'
     ORDER BY dataset_id, fit_id")
  n_wg <- nrow(wgcna_fits)
  message("run_all_enrichment(): ", n_wg, " ok WGCNA fit(s) to consider")
  for (i in seq_len(n_wg)) {
    fit_id <- wgcna_fits$fit_id[i]
    ensembl_map <- ensembl_maps[[wgcna_fits$dataset_id[i]]]
    modules <- DBI::dbGetQuery(con,
      "SELECT DISTINCT module FROM wgcna_modules WHERE fit_id = ? ORDER BY module",
      params = list(fit_id))$module
    message("[", i, "/", n_wg, "] ", wgcna_fits$dataset_id[i], " / wgcna / fit_id ", fit_id,
            " (", length(modules), " modules)")
    for (mod in modules) {
      factor_id <- ib_get_factor_id(con, fit_id, mod)
      if (length(factor_id) != 1) next
      if (!is.null(ib_enrichment_cached(con, factor_id, "ora", "pos"))) next
      genes <- DBI::dbGetQuery(con, "SELECT gene FROM wgcna_modules WHERE fit_id = ? AND module = ?",
                                params = list(fit_id, mod))$gene
      genes_ens <- if (!is.null(ensembl_map)) {
        mapped <- ensembl_map[genes]
        unique(mapped[!is.na(mapped) & nzchar(mapped)])
      } else {
        genes
      }
      res <- tryCatch(
        gprofiler2::gost(
          query = genes_ens, organism = "hsapiens", significant = TRUE,
          ordered_query = FALSE, correction_method = "fdr",
          sources = c("GO:BP", "GO:MF", "REAC", "KEGG", "WP")),
        error = function(e) {
          message("  g:Profiler query failed for fit ", fit_id, " module ", mod, ": ", conditionMessage(e))
          NULL
        })
      ib_enrichment_store(con, factor_id, "ora", "pos", res)
    }
  }

  message("run_all_enrichment(): done")
  invisible(NULL)
}

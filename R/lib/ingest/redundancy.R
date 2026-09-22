# "Are we discovering too many patterns?" diagnostic, computed on ONE
# representative fit per rank (see representative_fit_ids() below), NEVER
# across seeds -- this is a WITHIN-fit question (how distinct are this
# fit's own K factors from each other), orthogonal to the existing
# seed-stability machinery in similarity.R/pairs.R (which compares factors
# ACROSS fits). Two complementary views, both off the plain loadings
# matrix (genes x factors) via prepare_loadings()/pair_similarities() from
# similarity.R:
#   1. within_fit_similarity(): the K x K signed-cosine matrix of a fit's
#      OWN factors against each other (off-diagonal = redundancy) --
#      meaningless for PCA/Tucker (orthogonal by construction, always ~0
#      off-diagonal), so run_all_redundancy() skips those.
#   2. pattern_markers(): CoGAPS delegates to the real
#      CoGAPS::patternMarkers() (needs the full CogapsResult object -- see
#      fits.raw_result_file, captured at ingest time specifically for
#      this). Every other method gets a generic analog: each gene
#      "claimed" by whichever factor it loads onto most strongly (row
#      z-scored across factors, to avoid non-negative NMF/CoGAPS loadings
#      being confounded by a gene's overall expression level). Shrinking/
#      zero marker counts for newly added factors as rank grows is the
#      "too many patterns" signal.
#
# As in R/lib/ingest/pairs.R, functions here come in a pure
# (`_from_universe`/`_for_fit`) form usable during a PARALLEL, no-DB-
# writer compute phase (see ingest_dataset.R::compute_ingest_bundle()),
# plus a thin DB-querying wrapper preserving the original API for the
# single-dataset direct path and merge-time use.

within_fit_similarity <- function(loadings) {
  prep <- prepare_loadings(loadings)
  pair_similarities(prep, prep)$cosine   # K x K, diagonal == 1
}

generic_pattern_markers <- function(loadings) {
  z <- t(scale(t(loadings)))             # row (gene) z-score across factors
  # A gene with (near-)identical loadings across every factor -- common for
  # sparse methods like sPCA, especially at low K, where most loadings are
  # exactly 0 -- gets zero column variance from scale(), producing an
  # all-NaN row here. which.max() on an all-NaN vector returns integer(0),
  # which breaks apply()'s result simplification (it silently returns a
  # list instead of an atomic vector, and everything downstream indexing
  # into it fails with "invalid subscript type 'list'"). Treat such a gene
  # as tied across all factors (z = 0 everywhere) rather than erroring --
  # which.max() then deterministically picks the first factor.
  z[!is.finite(z)] <- 0
  best <- apply(z, 1, which.max)
  data.frame(gene = rownames(loadings), factor_index = best,
             score = z[cbind(seq_len(nrow(z)), best)], stringsAsFactors = FALSE)
}

#' CoGAPS: the real patternMarkers(), flattened to the same
#' (gene, factor_index, score) shape generic_pattern_markers() returns.
#' `score` here is rank-within-pattern (patternMarkers() doesn't return a
#' comparable continuous statistic), not directly comparable across methods.
cogaps_pattern_markers <- function(cogaps_result) {
  pm <- CoGAPS::patternMarkers(cogaps_result)
  out <- do.call(rbind, lapply(seq_along(pm$PatternMarkers), function(k) {
    genes <- pm$PatternMarkers[[k]]
    if (length(genes) == 0) return(NULL)
    data.frame(gene = genes, factor_index = k, score = seq_along(genes), stringsAsFactors = FALSE)
  }))
  if (is.null(out)) out <- data.frame(gene = character(0), factor_index = integer(0), score = numeric(0))
  out
}

#' Per-fit scalar summary -- what the Level 1 redundancy-vs-rank plot uses.
redundancy_summary <- function(loadings, markers) {
  sim <- within_fit_similarity(loadings)
  diag(sim) <- NA
  list(
    max_offdiag = max(abs(sim), na.rm = TRUE),
    median_offdiag = median(abs(sim), na.rm = TRUE),
    n_factors_with_no_markers = length(setdiff(seq_len(ncol(loadings)), unique(markers$factor_index))),
    matrix = sim
  )
}

#' Compute one fit's redundancy diagnostic from already-resolved absolute
#' artifact paths (no DB access) -- `loadings_path` and `raw_result_path`
#' (CoGAPS only, NA otherwise) are absolute file paths. Returns
#' list(markers, summary) or NULL if the loadings can't be read.
run_redundancy_for_fit <- function(method, loadings_path, raw_result_path = NA_character_) {
  if (is.na(loadings_path) || !file.exists(loadings_path)) return(NULL)
  L <- as.matrix(readRDS(loadings_path))

  markers <- if (method == "cogaps" && !is.na(raw_result_path) && file.exists(raw_result_path)) {
    cogaps_raw <- tryCatch(readRDS(raw_result_path), error = function(e) NULL)
    if (is.null(cogaps_raw)) generic_pattern_markers(L) else cogaps_pattern_markers(cogaps_raw)
  } else {
    generic_pattern_markers(L)
  }
  list(markers = markers, summary = redundancy_summary(L, markers))
}

#' One method + dataset -> the fit_id(s) ingest should treat as
#' "representative" for every non-seed-sweep analysis (redundancy,
#' enrichment, projectr, pattern-drivers). PCA collapses to a SINGLE fit
#' (its max configured rank -- components are nested/prefix-consistent,
#' see R/methods/pca.R, so every smaller rank's PCs are already contained
#' in it). nmf/cogaps/ica pick the lowest-mse seed AT EACH rank (mirrors
#' the app's existing "best seed per rank" convention). sPCA has no seed
#' dimension either, but DOES have its own second swept parameter (`para`,
#' stored in fits.alpha -- see extract_result()'s spca branch) crossed
#' with rank/K -- picks the lowest-mse para AT EACH K, same idea as
#' nmf/cogaps/ica, rather than treating every (K, para) combination as
#' independently representative (confirmed directly: doing the latter made
#' sPCA 81.6% of all fgsea_grid rows across this project's real DB --
#' 1824 of 2236 -- and was the dominant contributor to fgsea_grid tasks
#' timing out/OOMing on the cluster). cp/tucker/wgcna have no second
#' parameter to collapse -- every ok fit is already "representative".
#'
#' Pure version: `fits` is a data.frame(id, method, family, rank, alpha,
#' mse, status) covering the FULL universe (already-merged + this
#' bundle's new fits) for one dataset -- same id-sign convention as
#' R/lib/ingest/pairs.R (positive = real fit_id, negative = local id).
select_representative_ids_from_universe <- function(fits, method) {
  f <- fits[fits$method == method & fits$status == "ok", , drop = FALSE]
  if (nrow(f) == 0) return(integer(0))

  if (method == "pca") {
    f <- f[order(-f$rank), , drop = FALSE]
    return(f$id[1])
  }
  if (method %in% c("nmf", "cogaps", "ica")) {
    f <- f[f$family == "seed_sweep", , drop = FALSE]
    if (nrow(f) == 0) return(integer(0))
    ranks <- sort(unique(f$rank))
    return(vapply(ranks, function(rk) {
      sub <- f[f$rank == rk, , drop = FALSE]
      sub$id[which.min(sub$mse)]
    }, integer(1)))
  }
  if (method == "spca") {
    ks <- sort(unique(f$rank))
    return(vapply(ks, function(k) {
      sub <- f[f$rank == k, , drop = FALSE]
      sub$id[which.min(sub$mse)]
    }, integer(1)))
  }
  f$id
}

#' Thin DB-querying wrapper -- unchanged public behavior.
representative_fit_ids <- function(con, dataset_id, method) {
  fits <- DBI::dbGetQuery(con,
    "SELECT fit_id AS id, method, family, rank, alpha, mse, status
     FROM fits WHERE dataset_id = ? AND method = ?", params = list(dataset_id, method))
  select_representative_ids_from_universe(fits, method)
}

#' Runs redundancy for every representative fit of every non-orthogonal
#' method (nmf/cogaps/spca/ica -- PCA/Tucker skipped, orthogonal by
#' construction; WGCNA has its own analogous cross-power module-overlap
#' view already, computed elsewhere). `force`/--recompute-redundancy
#' recomputes regardless of the cached artifact; otherwise a fit's
#' redundancy is only recomputed if this file (redundancy.R) has been
#' edited more recently than the cached artifact -- fit data itself can't
#' go stale (fit_ids are immutable; delete_family() already cascades), so
#' the only real staleness vector here is the ANALYSIS CODE changing.
#' @param project_root absolute project-root path used to locate THIS file
#'   (redundancy.R itself) for the staleness check below. Deliberately NOT
#'   here::here() -- that requires the `here` package (not installed in
#'   every container this can run inside, e.g. ingest_core) AND a
#'   resolvable project anchor from the CURRENT working directory, which
#'   inside a containerized slurm job is rslurm's own bundle directory, not
#'   the project root. Defaults to getwd(), correct when called from
#'   R/ingest_results.R's direct CLI use (run from the project root);
#'   run_ingest_core_compute_job() passes its bind-mounted PROJECT_ROOT
#'   explicitly via compute_ingest_bundle()'s own project_root argument.
run_all_redundancy <- function(con, db_path, dataset_id, force = FALSE, project_root = getwd()) {
  art_dir <- artifacts_dir(db_path, dataset_id)
  dir.create(art_dir, recursive = TRUE, showWarnings = FALSE)
  this_file <- file.path(project_root, "R/lib/ingest/redundancy.R")
  redundancy_src_mtime <- if (file.exists(this_file)) file.mtime(this_file) else Sys.time()

  for (method in c("nmf", "cogaps", "spca", "ica")) {
    fit_ids <- representative_fit_ids(con, dataset_id, method)
    for (fit_id in fit_ids) {
      existing <- DBI::dbGetQuery(con, "SELECT redundancy_file FROM fit_redundancy WHERE fit_id = ?",
                                   params = list(fit_id))
      have_cache <- nrow(existing) > 0 && !is.na(existing$redundancy_file[1])

      stale <- FALSE
      if (!force && have_cache) {
        art_path <- resolve_artifact(existing$redundancy_file[1], db_path)
        if (file.exists(art_path) && file.mtime(art_path) < redundancy_src_mtime) {
          message("  fit_redundancy for fit_id ", fit_id, " is STALE (redundancy.R edited since) -- recomputing")
          stale <- TRUE
        }
      }
      if (!force && !stale && have_cache) next

      f <- DBI::dbGetQuery(con, "SELECT * FROM fits WHERE fit_id = ?", params = list(fit_id))
      if (nrow(f) == 0 || is.na(f$loadings_file)) next
      out <- run_redundancy_for_fit(method, resolve_artifact(f$loadings_file, db_path),
                                     if (is.na(f$raw_result_file)) NA_character_ else resolve_artifact(f$raw_result_file, db_path))
      if (is.null(out)) next
      markers <- out$markers; summ <- out$summary

      fname <- sprintf("%s_fit%d_redundancy.rds", method, fit_id)
      saveRDS(summ$matrix, file.path(art_dir, fname))

      if (stale || have_cache) {
        DBI::dbExecute(con, "DELETE FROM pattern_markers WHERE fit_id = ?", params = list(fit_id))
        DBI::dbExecute(con, "DELETE FROM fit_redundancy WHERE fit_id = ?", params = list(fit_id))
      }
      DBI::dbWriteTable(con, "pattern_markers",
        cbind(fit_id = fit_id, markers[, c("factor_index", "gene", "score")]), append = TRUE)
      DBI::dbExecute(con,
        "INSERT INTO fit_redundancy (fit_id, max_offdiag_cosine, median_offdiag_cosine, n_factors_with_no_markers, redundancy_file)
         VALUES (?, ?, ?, ?, ?)",
        params = list(fit_id, summ$max_offdiag, summ$median_offdiag, summ$n_factors_with_no_markers,
                      file.path("stability_artifacts", dataset_id, fname)))
    }
  }
  invisible(NULL)
}

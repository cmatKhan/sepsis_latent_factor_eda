# Cross-dataset analyses as targets: projectR projections of each
# dataset's representative fits onto every other dataset (`projections`),
# and gene-space agreement between same-method bases of datasets in the
# same family (`gene_space_agreement`). Rows are keyed by fit_key and
# written by write_cross_dataset_db() after every dataset's own DB write.

#' Save a dataset's matrix next to the DB
#'
#' The `matrix_file_<dataset>` file target: projections onto this dataset read it.
#'
#' @param mat Feature x sample matrix.
#' @param dataset_id Dataset id.
#' @param db_path Path to the DB.
#' @return Path of `stability_artifacts/<dataset>/matrix.rds`.
dataset_matrix_file <- function(mat, dataset_id, db_path) {
  art_dir <- artifacts_dir(db_path, dataset_id)
  dir.create(art_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(art_dir, "matrix.rds")
  saveRDS(mat, path)
  path
}

#' Whether two datasets are in the same family
#'
#' @param families config/dataset_families.yml, parsed.
#' @param ds_a,ds_b Dataset ids.
#' @return `TRUE` if both belong to the same single family.
same_family <- function(families, ds_a, ds_b) {
  fam_a <- names(Filter(function(v) ds_a %in% v, families))
  fam_b <- names(Filter(function(v) ds_b %in% v, families))
  length(fam_a) == 1 && length(fam_b) == 1 && fam_a == fam_b
}

#' projectR jobs for one source dataset
#'
#' One row per (representative fit, other dataset), grouped (`tar_group`) by
#' target dataset so each branch projects every representative fit onto one
#' dataset. Within a family the intercept is left out.
#'
#' @param source_id Source dataset id.
#' @param parts Named list (by loadings method) of `list(ingest, representatives)`.
#' @param matrix_files Every dataset's matrix path (a file target, so unnamed:
#'   the dataset id is read from the path).
#' @param families config/dataset_families.yml, parsed.
#' @param db_path Path to the DB.
#' @return Job table, or a placeholder.
projectr_jobs <- function(source_id, parts, matrix_files, families, db_path) {
  # matrix_files is a file target, which drops names: the dataset id is the
  # artifact directory (stability_artifacts/<id>/matrix.rds).
  names(matrix_files) <- basename(dirname(matrix_files))
  rows <- list()
  for (method in names(parts)) {
    ing <- parts[[method]]$ingest
    keys <- parts[[method]]$representatives
    if (length(keys) == 0) next
    lf <- vapply(keys, fit_artifact, character(1), ingest = ing, kind = "loadings", db_path = db_path)
    sf <- vapply(keys, fit_artifact, character(1), ingest = ing, kind = "scores", db_path = db_path)
    ok <- !is.na(lf)
    if (!any(ok)) next
    for (target_id in setdiff(names(matrix_files), source_id)) {
      within <- same_family(families, source_id, target_id)
      rows[[length(rows) + 1]] <- data.frame(
        source_fit_key = keys[ok], source_method = method,
        loadings_file = unname(lf[ok]), scores_file = unname(sf[ok]),
        target_dataset_id = target_id, target_matrix_file = unname(matrix_files[[target_id]]),
        projection_type = if (within) "within_dataset" else "cross_dataset",
        include_intercept = !within)
    }
  }
  jobs <- placeholder_jobs(do.call(rbind, rows), "source_fit_key")
  jobs$tar_group <- as.integer(factor(jobs$target_dataset_id %||% NA_character_, exclude = NULL))
  jobs
}

#' Run one group of projections (one `projectr_<dataset>` branch)
#'
#' @param jobs The group's job rows.
#' @param source_id Source dataset id.
#' @param ensembl_maps Every dataset's Ensembl map.
#' @param source_matrix_file The source dataset's matrix path.
#' @param db_path Path to the DB.
#' @return `projections` rows (projection_row()), or `NULL`.
projection_rows <- function(jobs, source_id, ensembl_maps, source_matrix_file, db_path) {
  rows <- lapply(seq_len(nrow(jobs)), function(i) {
    projection_row(jobs[i, , drop = FALSE], source_id, ensembl_maps, source_matrix_file, db_path)
  })
  do.call(rbind, rows)
}

#' Run one projection
#'
#' Saves the full projectR result as an artifact.
#'
#' @param job One job row.
#' @param source_id Source dataset id.
#' @param ensembl_maps Every dataset's Ensembl map.
#' @param source_matrix_file The source dataset's matrix path.
#' @param db_path Path to the DB.
#' @return One `projections` row keyed by source fit_key, or `NULL` for a
#'   placeholder or a failed projection.
projection_row <- function(job, source_id, ensembl_maps, source_matrix_file, db_path) {
  if (!has_job(job$source_fit_key)) return(NULL)
  x <- run_projectr_job(job$source_fit_key, job$source_method, source_id, job$loadings_file,
                        job$target_dataset_id, job$target_matrix_file,
                        job$projection_type, job$include_intercept,
                        ensembl_maps[[source_id]], ensembl_maps[[job$target_dataset_id]],
                        source_mat_file = source_matrix_file,
                        source_scores_file = if (is.na(job$scores_file)) NULL else job$scores_file)
  if (is.null(x$result) || is.null(x$result$projection)) return(NULL)
  res <- x$result
  fname <- sprintf("projection_%s_tgt_%s.rds", job$source_fit_key, job$target_dataset_id)
  saveRDS(res, file.path(artifacts_dir(db_path, source_id), fname))
  summ <- function(v, f) if (is.null(v)) NA_real_ else f(v, na.rm = TRUE)
  data.frame(
    source_fit_key = job$source_fit_key, target_dataset_id = job$target_dataset_id,
    projection_type = job$projection_type, include_intercept = as.integer(isTRUE(job$include_intercept)),
    n_genes_matched = x$n_genes_matched, n_samples = ncol(res$projection),
    mean_r_squared = summ(res$r_squared, mean), median_r_squared = summ(res$r_squared, stats::median),
    mean_pval = summ(res$pval, mean), median_pval = summ(res$pval, stats::median),
    mean_pvar = summ(res$pvar, mean), median_pvar = summ(res$pvar, stats::median),
    path = artifact_rel(source_id, fname))
}

#' A method's "optimal" fit for gene-space comparison
#'
#' PCA's fit; sPCA's selected fit with the most variance explained; for seed
#' sweeps, the lowest-mse seed at the rank with the lowest mean mse (bootstrap
#' resamples excluded).
#'
#' @param ingest An ingest_method_fits() result.
#' @param representatives Representative fit_keys.
#' @return A fit_key, or `NULL`.
optimal_fit_key <- function(ingest, representatives) {
  f <- ingest$fits[ingest$fits$status == "ok" & ingest$fits$bootstrap == 0, , drop = FALSE]
  if (nrow(f) == 0) return(NULL)
  if (ingest$method == "pca") return(f$fit_key[1])
  if (ingest$method == "spca") {
    m <- ingest$metrics[ingest$metrics$name == "pev_sparse" & ingest$metrics$fit_key %in% representatives, ]
    return(if (nrow(m)) m$fit_key[which.max(m$value)] else NULL)
  }
  f <- f[!is.na(f$mse), , drop = FALSE]
  if (nrow(f) == 0) return(NULL)
  mean_mse <- tapply(f$mse, f$rank, mean)
  best_rank <- as.numeric(names(mean_mse)[which.min(mean_mse)])
  at <- f[f$rank == best_rank, , drop = FALSE]
  at$fit_key[which.min(at$mse)]
}

#' Gene-space agreement with the other datasets in a family
#'
#' For every other dataset in `parts` and every loadings method both have,
#' projects the source's optimal fit onto the target's (projectR on shared
#' Ensembl genes), Hungarian-matches by |value| and reports the mean matched
#' |projection|.
#'
#' @param source_id Source dataset id.
#' @param parts Named list (by dataset) of named lists (by method) of
#'   `list(ingest, representatives)`.
#' @param ensembl_maps Every dataset's Ensembl map.
#' @param db_path Path to the DB.
#' @return Data frame (source_fit_key, target_fit_key, n_genes_matched,
#'   mean_abs_diagonal), or `NULL`.
gene_space_rows <- function(source_id, parts, ensembl_maps, db_path) {
  rows <- list()
  load <- function(ds, ing, key) {
    L <- as.matrix(readRDS(fit_artifact(ing, key, "loadings", db_path)))
    if (length(ensembl_maps[[ds]])) remap_to_ensembl(L, ensembl_maps[[ds]]) else L
  }
  for (target_id in setdiff(names(parts), source_id)) {
    for (m in intersect(names(parts[[source_id]]), names(parts[[target_id]]))) {
      a <- parts[[source_id]][[m]]; b <- parts[[target_id]][[m]]
      ka <- optimal_fit_key(a$ingest, a$representatives)
      kb <- optimal_fit_key(b$ingest, b$representatives)
      if (is.null(ka) || is.null(kb)) next
      La <- load(source_id, a$ingest, ka)
      Lb <- load(target_id, b$ingest, kb)
      shared <- intersect(rownames(La), rownames(Lb))
      if (length(shared) < 3) next
      proj <- tryCatch(projectR::projectR(data = La[shared, , drop = FALSE], loadings = Lb[shared, , drop = FALSE],
                                          full = TRUE), error = function(e) NULL)
      if (is.null(proj)) next
      P <- proj$projection
      rows[[length(rows) + 1]] <- data.frame(
        source_fit_key = ka, target_fit_key = kb, n_genes_matched = length(shared),
        mean_abs_diagonal = mean(abs(P[hungarian_match_abs(P)]), na.rm = TRUE))
    }
  }
  do.call(rbind, rows)
}

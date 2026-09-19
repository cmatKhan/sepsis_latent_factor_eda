# Phase-2 slurm SINGLE job (NOT an array, NOT parallelized) -- g:Profiler's
# API is rate-limited, so every (fit, factor, direction) combination is
# looped over sequentially in one process. `targets` (one row per
# fit x factor x direction) and `ensembl_maps` are baked in as global
# objects at grid-build time; loadings are read directly from their
# artifact files (no DB access from the compute node). Remapping to
# Ensembl (not gene symbol) is the canonical cross-dataset identifier
# space here -- see R/lib/ingest/symbol_mapping.R's header for why.
# gprofiler2::gost() natively accepts Ensembl gene ids as `query` (no
# extra argument needed -- it auto-detects id type).
#
# `gprofiler_targets`/`ensembl_maps` are read as FREE VARIABLES, NOT
# function parameters -- see run_ingest_core_job()'s header (R/ingest_jobs/
# ingest_core_job.R) for why: this job is staged via submit_job_family()
# with jobs_df = NULL (slurm_call(), no `params`), so rslurm calls this
# function with ZERO arguments. Declaring them as formals (as an earlier
# version of this function did, named `targets`) shadows the global lookup
# -- a required-with-no-default formal fails outright with "argument ...
# is missing, with no default", and one WITH a default (like the earlier
# `symbol_maps = list()`) would silently fall back to that default instead
# of the real baked-in map, never even erroring. Named `gprofiler_targets`
# specifically (not the more generic `targets`) to avoid colliding with
# --stage core's own unrelated `targets` global of the same script run.
run_gprofiler_job <- function() {
  library(gprofiler2)
  n <- nrow(gprofiler_targets)
  results <- vector("list", n)
  loaded <- new.env(parent = emptyenv())   # cache loadings per fit_id within this one process

  # Per-combo diagnostics -- added because gprofiler2::gost()'s own
  # "No results to show / Please make sure that the organism is correct
  # or set significant = FALSE" message() (printed straight to the log,
  # not something we control) is ambiguous on its own: it fires both for
  # the EXPECTED case (top-100 loadings genes for this factor just don't
  # enrich significantly for anything) and for real problems (empty/
  # malformed gene id vector after remap_to_ensembl(), wrong organism,
  # transient API failure). Logging dataset/fit/factor/direction + gene
  # count + a few sample ids BEFORE each call, and the actual hit count
  # (or the real error, previously swallowed to NULL with no trace) AFTER
  # each call, makes it possible to tell from the log alone which of
  # those cases each "No results to show" line belongs to, and whether
  # it's happening for every submission or only some.
  for (i in seq_len(n)) {
    fid <- gprofiler_targets$fit_id[i]
    dsid <- gprofiler_targets$dataset_id[i]
    fi <- gprofiler_targets$factor_index[i]
    dir_i <- gprofiler_targets$direction[i]
    key <- as.character(fid)
    if (!exists(key, envir = loaded)) {
      L <- as.matrix(readRDS(gprofiler_targets$loadings_file[i]))
      assign(key, remap_to_ensembl(L, ensembl_maps[[dsid]]), envir = loaded)
    }
    L_ens <- get(key, envir = loaded)
    v <- L_ens[, fi]
    genes <- if (dir_i == "neg") names(sort(v))[seq_len(min(100, length(v)))]
             else names(sort(v, decreasing = TRUE))[seq_len(min(100, length(v)))]
    genes_ordered <- names(sort(v, decreasing = (dir_i != "neg")))

    message(sprintf(
      "[gprofiler %d/%d] dataset=%s fit=%d factor=%d direction=%s n_loadings=%d n_genes(top100)=%d n_valid_ensembl=%d sample_ids=%s",
      i, n, dsid, fid, fi, dir_i, length(v), length(genes),
      sum(grepl("^ENSG", genes)), paste(utils::head(genes, 3), collapse = ", ")
    ))
    if (length(genes) == 0) {
      message("  -- skipping: zero genes after remap_to_ensembl(); check ensembl_maps[[", dsid, "]]")
    }

    res_ora <- tryCatch(
      gprofiler2::gost(query = genes, organism = "hsapiens", significant = TRUE,
                        correction_method = "fdr", sources = c("GO:BP", "GO:MF", "REAC", "KEGG", "WP")),
      error = function(e) {
        message("  ORA gost() failed: ", conditionMessage(e))
        NULL
      })
    message("  ORA: ", if (is.null(res_ora)) "call failed (see error above)"
            else if (is.null(res_ora$result)) "no significant terms"
            else paste0(nrow(res_ora$result), " significant terms"))

    res_gsea <- tryCatch(
      gprofiler2::gost(query = genes_ordered, organism = "hsapiens", significant = TRUE, ordered_query = TRUE,
                        correction_method = "fdr", sources = c("GO:BP", "GO:MF", "REAC", "KEGG", "WP")),
      error = function(e) {
        message("  GSEA gost() failed: ", conditionMessage(e))
        NULL
      })
    message("  GSEA: ", if (is.null(res_gsea)) "call failed (see error above)"
            else if (is.null(res_gsea$result)) "no significant terms"
            else paste0(nrow(res_gsea$result), " significant terms"))

    results[[i]] <- list(dataset_id = dsid, method = gprofiler_targets$method[i],
                          fit_id = fid, factor_index = fi, direction = dir_i,
                          ora = res_ora, gsea = res_gsea)
    Sys.sleep(1)   # polite pacing -- API-bound, not compute-bound
  }
  results
}

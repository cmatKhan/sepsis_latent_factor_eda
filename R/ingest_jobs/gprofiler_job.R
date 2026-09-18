# Phase-2 slurm SINGLE job (NOT an array, NOT parallelized) -- g:Profiler's
# API is rate-limited, so every (fit, factor, direction) combination is
# looped over sequentially in one process. `targets` (one row per
# fit x factor x direction) and `symbol_maps` are baked in as global
# objects at grid-build time; loadings are read directly from their
# artifact files (no DB access from the compute node).
#
# `gprofiler_targets`/`symbol_maps` are read as FREE VARIABLES, NOT
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
  results <- vector("list", nrow(gprofiler_targets))
  loaded <- new.env(parent = emptyenv())   # cache loadings per fit_id within this one process

  for (i in seq_len(nrow(gprofiler_targets))) {
    fid <- gprofiler_targets$fit_id[i]
    key <- as.character(fid)
    if (!exists(key, envir = loaded)) {
      L <- as.matrix(readRDS(gprofiler_targets$loadings_file[i]))
      assign(key, remap_to_symbol(L, symbol_maps[[gprofiler_targets$dataset_id[i]]]), envir = loaded)
    }
    L_sym <- get(key, envir = loaded)
    fi <- gprofiler_targets$factor_index[i]
    v <- L_sym[, fi]
    genes <- if (gprofiler_targets$direction[i] == "neg") names(sort(v))[seq_len(min(100, length(v)))]
             else names(sort(v, decreasing = TRUE))[seq_len(min(100, length(v)))]

    res_ora <- tryCatch(
      gprofiler2::gost(query = genes, organism = "hsapiens", significant = TRUE,
                        correction_method = "fdr", sources = c("GO:BP", "GO:MF", "REAC", "KEGG", "WP")),
      error = function(e) NULL)
    res_gsea <- tryCatch(
      gprofiler2::gost(query = names(sort(v, decreasing = (gprofiler_targets$direction[i] != "neg"))),
                        organism = "hsapiens", significant = TRUE, ordered_query = TRUE,
                        correction_method = "fdr", sources = c("GO:BP", "GO:MF", "REAC", "KEGG", "WP")),
      error = function(e) NULL)

    results[[i]] <- list(dataset_id = gprofiler_targets$dataset_id[i], method = gprofiler_targets$method[i],
                          fit_id = fid, factor_index = fi, direction = gprofiler_targets$direction[i],
                          ora = res_ora, gsea = res_gsea)
    Sys.sleep(1)   # polite pacing -- API-bound, not compute-bound
  }
  results
}

# Phase-2 slurm SINGLE job (NOT an array, NOT parallelized) -- g:Profiler's
# API is rate-limited, so every (fit, factor, direction) combination is
# looped over sequentially in one process. `targets` (one row per
# fit x factor x direction) and `symbol_maps` are baked in as global
# objects at grid-build time; loadings are read directly from their
# artifact files (no DB access from the compute node).

run_gprofiler_job <- function(targets, symbol_maps = list()) {
  library(gprofiler2)
  results <- vector("list", nrow(targets))
  loaded <- new.env(parent = emptyenv())   # cache loadings per fit_id within this one process

  for (i in seq_len(nrow(targets))) {
    fid <- targets$fit_id[i]
    key <- as.character(fid)
    if (!exists(key, envir = loaded)) {
      L <- as.matrix(readRDS(targets$loadings_file[i]))
      assign(key, remap_to_symbol(L, symbol_maps[[targets$dataset_id[i]]]), envir = loaded)
    }
    L_sym <- get(key, envir = loaded)
    fi <- targets$factor_index[i]
    v <- L_sym[, fi]
    genes <- if (targets$direction[i] == "neg") names(sort(v))[seq_len(min(100, length(v)))]
             else names(sort(v, decreasing = TRUE))[seq_len(min(100, length(v)))]

    res_ora <- tryCatch(
      gprofiler2::gost(query = genes, organism = "hsapiens", significant = TRUE,
                        correction_method = "fdr", sources = c("GO:BP", "GO:MF", "REAC", "KEGG", "WP")),
      error = function(e) NULL)
    res_gsea <- tryCatch(
      gprofiler2::gost(query = names(sort(v, decreasing = (targets$direction[i] != "neg"))),
                        organism = "hsapiens", significant = TRUE, ordered_query = TRUE,
                        correction_method = "fdr", sources = c("GO:BP", "GO:MF", "REAC", "KEGG", "WP")),
      error = function(e) NULL)

    results[[i]] <- list(dataset_id = targets$dataset_id[i], method = targets$method[i],
                          fit_id = fid, factor_index = fi, direction = targets$direction[i],
                          ora = res_ora, gsea = res_gsea)
    Sys.sleep(1)   # polite pacing -- API-bound, not compute-bound
  }
  results
}

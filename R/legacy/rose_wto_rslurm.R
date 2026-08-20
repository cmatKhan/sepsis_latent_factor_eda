# wTO rslurm job submission for analyzing_rose_decomp_results.qmd (Q5).
#
# Reads the gene sets / expression inputs written by the qmd's `wto-gene-sets`
# and `wto-rep-measure-setup` chunks and submits the two rslurm bundles
# (wTO.Complete per gene set x day, and the repeated-measures wTO per gene
# set) out-of-band. Run this script directly on the cluster (e.g.
# `Rscript R/rose_wto_rslurm.R`); the qmd's collect chunks then read the
# results back from data/wto_complete/ and data/wto_rep_measure/.

library(here)
library(wTO)
library(rslurm)

wto_inputs      <- readRDS(here("data/wto_inputs.rds"))
overlap_sets    <- wto_inputs$overlap_sets
wto_k           <- wto_inputs$wto_k

wto_rep_inputs   <- readRDS(here("data/wto_rep_measure_inputs.rds"))
rep_measure_sets <- wto_rep_inputs$rep_measure_sets
vst_paired       <- wto_rep_inputs$vst_paired
id_paired        <- wto_rep_inputs$id_paired
wto_rep_k        <- wto_rep_inputs$wto_rep_k

decomp   <- readRDS(here("results/rose_decomp.rds"))
vst_day0 <- decomp$expr$vst_day0
vst_day2 <- decomp$expr$vst_day2

## ---- wTO.Complete: one job per (gene set, day) --------------------------

# Each worker returns a named list so the collect chunk can reconstruct the
# geneset/day mapping without relying on file-listing order.
run_wto_complete_job <- function(geneset_idx, day) {
  library(wTO)
  genes <- overlap_sets[[geneset_idx]]
  data  <- if (day == "day0") vst_day0 else vst_day2
  genes <- intersect(genes, rownames(data))
  res <- wTO.Complete(
    k = wto_k, n = 1000,
    Data              = as.data.frame(data),
    Overlap           = genes,
    method            = "p",
    method_resampling = "Bootstrap",
    pvalmethod        = "BH",
    plot              = FALSE
  )
  list(result = res, geneset = names(overlap_sets)[geneset_idx], day = day)
}

wto_complete_jobs <- expand.grid(
  geneset_idx      = seq_along(overlap_sets),
  day              = c("day0", "day2"),
  stringsAsFactors = FALSE
)

sjob_wto_complete <- slurm_apply(
  f              = run_wto_complete_job,
  params         = wto_complete_jobs,
  jobname        = "wto_complete",
  nodes          = nrow(wto_complete_jobs),
  cpus_per_node  = 1,
  global_objects = c("overlap_sets", "vst_day0", "vst_day2", "wto_k"),
  pkgs           = c("wTO"),
  libPaths       = "/ref/mblab/software/chasem/R/4.6.0",
  rscript_path   = "Rscript",
  sh_template    = "~/.rslurm/submit_sh.txt",
  slurm_options  = list(mem = "16G",
                        "cpus-per-task" = wto_k,
                        time = "08:00:00",
                        container = "docker://rocker/tidyverse:latest"),
  submit         = FALSE
)

## ---- wTO.rep_measure: one job per gene set -------------------------------

# wTO:::rmcor computes a repeated-measures (within-subject ANCOVA) correlation
# for one gene pair; this builds the full gene x gene matrix from it.
pairwise_rmcor <- function(mat, id_vec) {
  rmcor_fn <- getFromNamespace("rmcor", "wTO")
  g   <- ncol(mat)
  nms <- colnames(mat)
  Cor <- matrix(0, g, g, dimnames = list(nms, nms))
  if (g < 2) return(as.data.frame(Cor))
  idx  <- utils::combn(g, 2)
  vals <- vapply(seq_len(ncol(idx)), function(p) {
    v <- suppressWarnings(rmcor_fn(id_vec, mat[, idx[1, p]], mat[, idx[2, p]]))
    if (is.na(v)) v <- 0
    v
  }, numeric(1))
  Cor[cbind(idx[1, ], idx[2, ])] <- vals
  Cor[cbind(idx[2, ], idx[1, ])] <- vals
  as.data.frame(Cor)
}

# Reimplements wTO::wTO.rep_measure's algorithm exactly (same rmcor / wTO() /
# wTO.in.line() calls), but distributes the n bootstrap replicates across a
# PSOCK cluster instead of running them in a single-threaded for-loop.
#
# Two things forced a reimplementation rather than calling wTO.rep_measure
# directly: (1) the installed version crashes (native SIGABRT inside
# Rfast::Crossprod) whenever Overlap is a strict subset of Data's genes,
# because wTO::wTO() is then called on a non-square matrix; passing Data
# already restricted to the gene set of interest (Overlap == all of Data)
# avoids this, and is what this function assumes. (2) it has no built-in
# parallelism and its rmcor step is O(genes^2) per replicate. A fork-based
# parallelization (parallel::mclapply) was tried first and reliably hung —
# Rfast's C++ correlation code is not fork-safe once its internal thread
# pool has started — so this uses parallel::makeCluster(type = "PSOCK")
# instead, which starts independent worker processes rather than forking.
wTO_rep_measure_parallel <- function(Data, ID, sign = "sign", delta = 0.2, n = 1000,
                                      k = 10) {
  ID    <- as.factor(ID)
  Datat <- t(Data)

  Cor0     <- pairwise_rmcor(Datat, ID)
  wtomelt0 <- wTO::wTO(Cor0, sign)

  cl <- parallel::makeCluster(k, type = "PSOCK")
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterEvalQ(cl, library(wTO))
  parallel::clusterExport(cl, c("pairwise_rmcor", "Datat", "ID", "wtomelt0", "sign", "delta"),
                           envir = environment())
  parallel::clusterSetRNGStream(cl)

  U_list <- parallel::parLapply(cl, seq_len(n), function(b) {
    bootID    <- sample(levels(ID), replace = TRUE)
    Data_boot <- do.call(rbind, lapply(bootID, function(bi) subset(Datat, ID == bi)))
    Cor_b     <- pairwise_rmcor(Data_boot, ID)
    res       <- wTO::wTO(Cor_b, sign)
    (res < wtomelt0 - delta) + (res > wtomelt0 + delta)
  })

  U <- Reduce(`+`, U_list)

  wtomelt0_line <- wTO::wTO.in.line(wtomelt0)
  U_line        <- wTO::wTO.in.line(U)
  data.table::data.table(wtomelt0_line, pval = U_line$wTO / n)
}

# wTO_rep_measure_parallel already spans a k-worker PSOCK cluster internally
# for the n bootstrap replicates, so each rslurm job here requests k CPUs on
# a single node rather than being split further.
run_wto_rep_measure_job <- function(geneset_idx) {
  library(wTO)
  genes <- rep_measure_sets[[geneset_idx]]
  genes <- intersect(genes, rownames(vst_paired))
  dat   <- as.data.frame(vst_paired[genes, , drop = FALSE])
  res   <- wTO_rep_measure_parallel(Data = dat, ID = id_paired, sign = "sign",
                                     delta = 0.2, n = 1000, k = wto_rep_k)
  list(result = res, geneset = names(rep_measure_sets)[geneset_idx])
}

wto_rep_jobs <- data.frame(geneset_idx = seq_along(rep_measure_sets))

sjob_wto_rep_measure <- slurm_apply(
  f              = run_wto_rep_measure_job,
  params         = wto_rep_jobs,
  jobname        = "wto_rep_measure",
  nodes          = nrow(wto_rep_jobs),
  cpus_per_node  = 1,
  global_objects = c("rep_measure_sets", "vst_paired", "id_paired", "wto_rep_k",
                     "wTO_rep_measure_parallel", "pairwise_rmcor"),
  pkgs           = c("wTO", "data.table"),
  libPaths       = "/ref/mblab/software/chasem/R/4.6.0",
  rscript_path   = "Rscript",
  sh_template    = "~/.rslurm/submit_sh.txt",
  slurm_options  = list(mem = "16G",
                        "cpus-per-task" = wto_rep_k,
                        time = "12:00:00",
                        container = "docker://rocker/tidyverse:latest"),
  submit         = FALSE
)

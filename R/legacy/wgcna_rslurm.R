# WGCNA rslurm job submission (Langfelder & Horvath 2008,
# literature/wgcna.pdf). Genome-wide weighted co-expression network,
# fit before/independently of the `wTO` candidate-gene-set analysis in
# analyzing_rose_decomp_results.qmd §Q5 — see CLAUDE.md "Methods compared".
#
# Reads vst_day0/vst_day2 directly from results/rose_decomp.rds (already
# exported by rose.qmd); no new saveRDS was needed anywhere else.
#
# Three job families, all single-wave except module preservation:
#   1. Per-day network + module detection (blockwiseModules), one job per
#      day — the WGCNA analogue of wTO.Complete.
#   2. Consensus module detection across Day0+Day2 (blockwiseConsensusModules)
#      — the WGCNA analogue of wTO.rep_measure.
#   3. Module preservation (Day0 modules tested in Day2 data and vice versa)
#      — needs job family 1's per-day module colors already on disk, so it's
#      gated on data/wgcna_network/runs having both days' results. Run this
#      script once to submit families 1 and 2; once family 1 finishes on the
#      cluster, run it again to pick up family 3.
#
# Run out-of-band: `Rscript R/wgcna_rslurm.R`.

library(here)
library(WGCNA)
library(rslurm)

options(stringsAsFactors = FALSE)

decomp <- readRDS(here("results/rose_decomp.rds"))
vst_day0 <- decomp$expr$vst_day0
vst_day2 <- decomp$expr$vst_day2

wgcna_cpu_per_task <- 8L
wgcna_n_permutations <- 100L
wgcna_power_grid <- c(1:10, seq(12, 20, by = 2))

# Standard WGCNA tutorial convention: smallest power reaching the target
# scale-free-topology fit, falling back to the best fit tested if none does.
# This is a free-parameter choice with the same structure as PCA's rank
# elbow — see CLAUDE.md Q3.
choose_power <- function(fit_indices, target_r2 = 0.8) {
  ok <- fit_indices$SFT.R.sq >= target_r2
  if (any(ok)) fit_indices$Power[which(ok)[1]] else fit_indices$Power[which.max(fit_indices$SFT.R.sq)]
}

# Read a directory of per-job RDS result files (rslurm convention: each file
# holds a list whose first element is the job result). Returns list().
read_job_dir <- function(dir) {
  if (!dir.exists(dir)) {
    return(list())
  }
  files <- list.files(dir, pattern = "\\.[rR][dD][sS]$", full.names = TRUE)
  lapply(files, function(f) {
    x <- readRDS(f)
    if (is.list(x) && length(x) == 1) x[[1]] else x
  })
}

## ---- Family 1: per-day network + module detection ------------------------

run_wgcna_network_job <- function(day) {
  library(WGCNA)
  enableWGCNAThreads(nThreads = wgcna_cpu_per_task)
  options(stringsAsFactors = FALSE)

  mat <- if (day == "day0") vst_day0 else vst_day2
  datExpr <- t(mat)

  gsg <- goodSamplesGenes(datExpr, verbose = 0)
  if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]

  sft <- pickSoftThreshold(datExpr,
    powerVector = wgcna_power_grid,
    networkType = "signed", verbose = 0
  )
  power <- choose_power(sft$fitIndices)

  net <- blockwiseModules(
    datExpr,
    power = power, networkType = "signed", TOMType = "signed",
    minModuleSize = 30, mergeCutHeight = 0.25, numericLabels = TRUE,
    pamRespectsDendro = FALSE, saveTOMs = FALSE,
    nThreads = wgcna_cpu_per_task, verbose = 0
  )

  list(
    day = day, power = power, sft = sft$fitIndices, net = net,
    genes = colnames(datExpr)
  )
}

wgcna_network_jobs <- data.frame(day = c("day0", "day2"), stringsAsFactors = FALSE)

sjob_wgcna_network <- slurm_apply(
  f = run_wgcna_network_job,
  params = wgcna_network_jobs,
  jobname = "wgcna_network",
  nodes = nrow(wgcna_network_jobs),
  cpus_per_node = 1,
  global_objects = c(
    "vst_day0", "vst_day2", "wgcna_power_grid",
    "choose_power", "wgcna_cpu_per_task"
  ),
  pkgs = c("WGCNA"),
  libPaths = "/ref/mblab/software/chasem/R/4.6.0",
  rscript_path = "Rscript",
  sh_template = "~/.rslurm/submit_sh.txt",
  slurm_options = list(
    mem = "16G",
    "cpus-per-task" = wgcna_cpu_per_task,
    time = "04:00:00",
    container = "oras://community.wave.seqera.io/library/r-wgcna:1.74--5149c638df2976dd"
  ),
  submit = FALSE
)

## ---- Family 2: consensus module detection across Day0 + Day2 -------------

# A single job with no parameter grid — slurm_call (rather than slurm_apply)
# is the rslurm-idiomatic way to submit exactly one job.
run_wgcna_consensus_job <- function() {
  library(WGCNA)
  enableWGCNAThreads(nThreads = wgcna_cpu_per_task)
  options(stringsAsFactors = FALSE)

  multiExpr <- list(
    day0 = list(data = t(vst_day0)),
    day2 = list(data = t(vst_day2))
  )

  powers <- vapply(multiExpr, function(set) {
    sft <- pickSoftThreshold(set$data,
      powerVector = wgcna_power_grid,
      networkType = "signed", verbose = 0
    )
    choose_power(sft$fitIndices)
  }, numeric(1))

  net <- blockwiseConsensusModules(
    multiExpr,
    power = powers, networkType = "signed", TOMType = "signed",
    minModuleSize = 30, mergeCutHeight = 0.25, numericLabels = TRUE,
    nThreads = wgcna_cpu_per_task, verbose = 0
  )

  list(powers = powers, net = net, genes = colnames(multiExpr$day0$data))
}

sjob_wgcna_consensus <- slurm_call(
  f = run_wgcna_consensus_job,
  params = list(),
  jobname = "wgcna_consensus",
  global_objects = c(
    "vst_day0", "vst_day2", "wgcna_power_grid",
    "choose_power", "wgcna_cpu_per_task"
  ),
  pkgs = c("WGCNA"),
  libPaths = "/ref/mblab/software/chasem/R/4.6.0",
  rscript_path = "Rscript",
  sh_template = "~/.rslurm/submit_sh.txt",
  slurm_options = list(
    mem = "24G",
    "cpus-per-task" = wgcna_cpu_per_task,
    time = "06:00:00",
    container = "docker://rocker/tidyverse:latest"
  ),
  submit = FALSE
)

## ---- Family 3: module preservation (Day0 <-> Day2), gated on family 1 ----

wgcna_network_dir <- here("data/wgcna_network/runs")
network_runs <- read_job_dir(wgcna_network_dir)
net_day0 <- Find(function(x) x$day == "day0", network_runs)
net_day2 <- Find(function(x) x$day == "day2", network_runs)

if (is.null(net_day0) || is.null(net_day2)) {
  message(
    "wgcna_rslurm: ", wgcna_network_dir,
    " doesn't have both day0 and day2 results yet — submit/wait on ",
    "the network job first, then rerun this script for module preservation."
  )
} else {
  genes_day0 <- net_day0$genes
  genes_day2 <- net_day2$genes
  colors_day0 <- net_day0$net$colors
  colors_day2 <- net_day2$net$colors

  # Bidirectional: are Day0 modules preserved in Day2 data, and vice versa?
  # Zsummary >= 10 = strong evidence of preservation, 2-10 = weak/moderate,
  # <2 = no evidence (Langfelder & Horvath 2011 convention).
  run_wgcna_modpres_job <- function() {
    library(WGCNA)
    enableWGCNAThreads(nThreads = wgcna_cpu_per_task)
    options(stringsAsFactors = FALSE)

    multiExpr <- list(
      day0 = list(data = t(vst_day0)[, genes_day0, drop = FALSE]),
      day2 = list(data = t(vst_day2)[, genes_day2, drop = FALSE])
    )

    mp_day0_ref <- modulePreservation(
      multiExpr,
      multiColor = list(day0 = colors_day0),
      referenceNetworks = 1, nPermutations = wgcna_n_permutations,
      randomSeed = 42, networkType = "signed", quickCor = 0, verbose = 0
    )
    mp_day2_ref <- modulePreservation(
      multiExpr,
      multiColor = list(day2 = colors_day2),
      referenceNetworks = 2, nPermutations = wgcna_n_permutations,
      randomSeed = 42, networkType = "signed", quickCor = 0, verbose = 0
    )

    list(day0_modules_in_day2 = mp_day0_ref, day2_modules_in_day0 = mp_day2_ref)
  }

  sjob_wgcna_modpres <- slurm_call(
    f = run_wgcna_modpres_job,
    params = list(),
    jobname = "wgcna_modpres",
    global_objects = c(
      "vst_day0", "vst_day2", "genes_day0", "genes_day2",
      "colors_day0", "colors_day2", "wgcna_cpu_per_task",
      "wgcna_n_permutations"
    ),
    pkgs = c("WGCNA"),
    libPaths = "/ref/mblab/software/chasem/R/4.6.0",
    rscript_path = "Rscript",
    sh_template = "~/.rslurm/submit_sh.txt",
    slurm_options = list(
      mem = "16G",
      "cpus-per-task" = wgcna_cpu_per_task,
      time = "06:00:00",
      container = "docker://rocker/tidyverse:latest"
    ),
    submit = FALSE
  )
}

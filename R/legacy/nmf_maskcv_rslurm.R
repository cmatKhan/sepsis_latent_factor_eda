# NMF masking-CV rank selection for nmf_characterization.qmd (Day 0, Day 2,
# and combined Day0+Day2 bases).
#
# Companion to R/nmf_rslurm.R's full grid: that script fits at full quality
# with 10 seeds/rank, but scores reconstruction on the SAME entries the
# model was fit on (in-sample MSE), which is monotonically decreasing in
# rank and cannot itself select a rank without overfitting. However, it
# does provide a measure of the stability of the patterns over seeds.
# This script specifically tests for rank optimality by holds out of a
# random 1% of each basis's matrix entries (set to NA, which NNLM::nnmf's
# native missing-data handling and excludes those samples from the fit
# entirely -- confirmed via nnmf's `check.k` doc.
# nmf_characterization.qmd sources its elbow-rank auto-pick from THIS
# script's output, not the full grid's in-sample curve.

library(here)
library(rslurm)
library(NNLM)
library(DESeq2)
library(SummarizedExperiment)

rank_range <- 5:20
mask_frac <- 0.01
mask_seed <- 2024L

dds <- readRDS(here("results/dds.rds"))
vst_mat <- assay(varianceStabilizingTransformation(dds, blind = TRUE))

day0_samp <- colData(dds)$timepoint == "Day 0"
day2_samp <- colData(dds)$timepoint == "Day 2"

vst_day0 <- vst_mat[, day0_samp]
vst_day2 <- vst_mat[, day2_samp]
vst_combined <- cbind(vst_day0, vst_day2)

shift_nonneg <- function(mat) mat - matrixStats::rowMins(mat)

mat_nn_by_basis <- list(
  day0 = shift_nonneg(vst_day0),
  day2 = shift_nonneg(vst_day2),
  combined = shift_nonneg(vst_combined)
)

# One fixed held-out mask per basis, drawn once so every rank is scored
# against identical held-out cells.
set.seed(mask_seed)
mask_idx_by_basis <- lapply(mat_nn_by_basis, function(mat) {
  sample(length(mat), round(mask_frac * length(mat)))
})

nmf_maskcv_jobs <- expand.grid(
  basis = names(mat_nn_by_basis),
  rank = rank_range,
  stringsAsFactors = FALSE
)

run_nmf_maskcv_job <- function(basis, rank) {
  library(NNLM)
  mat_nn_basis <- mat_nn_by_basis[[basis]]
  mask_idx <- mask_idx_by_basis[[basis]]

  mat_masked <- mat_nn_basis
  mat_masked[mask_idx] <- NA

  fit <- nnmf(mat_masked, k = rank, max.iter = 10000, verbose = 0L)
  recon <- fit$W %*% fit$H
  mse <- mean((mat_nn_basis[mask_idx] - recon[mask_idx])^2)

  list(basis = basis, rank = rank, mse = mse)
}

sjob_nmf_maskcv <- slurm_apply(
  f = run_nmf_maskcv_job,
  params = nmf_maskcv_jobs,
  jobname = "nmf_maskcv",
  nodes = nrow(nmf_maskcv_jobs),
  cpus_per_node = 1,
  global_objects = c("mat_nn_by_basis", "mask_idx_by_basis"),
  pkgs = c("NNLM"),
  libPaths = "/ref/mblab/software/chasem/R/4.6.0",
  rscript_path = "Rscript",
  sh_template = "~/.rslurm/submit_sh.txt",
  slurm_options = list(
    mem = "8G",
    "cpus-per-task" = 1,
    time = "01:00:00",
    container = "docker://rocker/tidyverse:latest"
  ),
  submit = FALSE
)

# NMF pattern grid search for nmf_characterization.qmd (Day 0, Day 2, and
# combined Day0+Day2 bases).
#
# Standalone: reads results/dds.rds directly (same source R/pca.R uses) and
# recomputes VST itself, rather than depending on any rose.qmd-produced
# input file. rose.qmd's own NMF fitting (masking-CV rank selection, then
# separate multi-restart/subsample waves) is being retired/refactored into
# standalone components; this script and nmf_characterization.qmd are the
# NMF replacement, following R/pca.R's pattern.
#
# Grid: fit NNLM::nnmf at every rank in rank_range, 10 seeds per rank, for
# all three bases (day0, day2, combined) in a single combined job array —
# one sbatch submission instead of three, following the same "one job
# function parameterized by basis, all per-basis data passed as global
# objects" pattern R/rose_wto_rslurm.R already uses for its geneset x day
# grid. Each job returns its reconstruction error (mean squared full
# in-sample reconstruction error, no masking) alongside its fitted W/H, so
# nmf_characterization.qmd can (a) plot rank vs. reconstruction error with
# seed-based error bars for rank selection, and (b) pull the exact W
# matrix for whichever (basis, rank, seed) is chosen after inspecting that
# plot — no separate refit step.
#
# rank_range was originally narrowed to 5:10 (rather than a wide exploratory
# sweep) because the masking-CV rank search previously done in rose.qmd
# identified ~7 as the reconstruction-optimal NMF rank. It's now widened to
# 5:20 to match R/nmf_maskcv_rslurm.R's range: that script's masking-CV
# curve (data/rslurm_output/nmf/maskcv) picked an elbow at rank 17 for the
# day0 basis, outside this grid's old 5:10 ceiling, and the combined basis
# behaves differently from day0/day2 (its masking-CV MSE keeps improving
# past rank 10) — so a real rank decision needs full-quality seed-stability
# fits available at those higher ranks too, not just the masking-CV single
# fit. Note the masking-CV curve gets numerically unstable/noisy above
# ~rank 12-13 for day0/day2 (MSE swings between single digits and the
# hundreds run to run) — treat any elbow picked above that range with
# proportionate skepticism once these results are in, rather than at face
# value.
#
# Run out-of-band: `Rscript R/nmf_rslurm.R`.

library(here)
library(rslurm)
library(NNLM)
library(DESeq2)
library(SummarizedExperiment)

rank_range <- 5:20
nmf_seeds <- c(42L, 123L, 456L, 7L, 99L, 2024L, 8675309L, 271828L, 31415L, 90210L)

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

nmf_grid_jobs <- expand.grid(
  basis = names(mat_nn_by_basis),
  rank = rank_range,
  seed = nmf_seeds,
  stringsAsFactors = FALSE
)

run_nmf_grid_job <- function(basis, rank, seed) {
  library(NNLM)
  mat_nn_basis <- mat_nn_by_basis[[basis]]
  set.seed(seed)
  fit <- nnmf(mat_nn_basis, k = rank, max.iter = 10000, verbose = 0L)
  rownames(fit$W) <- rownames(mat_nn_basis)
  colnames(fit$W) <- paste0("Pattern_", seq_len(ncol(fit$W)))
  colnames(fit$H) <- colnames(mat_nn_basis)
  rownames(fit$H) <- colnames(fit$W)

  recon <- fit$W %*% fit$H
  mse <- mean((mat_nn_basis - recon)^2)
  list(basis = basis, rank = rank, seed = seed, mse = mse, W = fit$W, H = fit$H)
}

sjob_nmf_grid <- slurm_apply(
  f = run_nmf_grid_job,
  params = nmf_grid_jobs,
  jobname = "nmf_grid",
  nodes = nrow(nmf_grid_jobs),
  cpus_per_node = 1,
  global_objects = c("mat_nn_by_basis"),
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

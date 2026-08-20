# CoGAPS pattern grid search for cogaps_characterization.qmd (Day 0, Day 2,
# and combined Day0+Day2 bases).
#
# Standalone: reads results/dds.rds directly, same as R/nmf_rslurm.R and
# R/pca.R. Does NOT depend on NMF's rank or any rose.qmd-produced input
# file — rose.qmd's own CoGAPS fitting (reusing the NMF rank, plus an
# alphaA/alphaP seed-stability grid) is being retired/refactored into
# standalone components; this script and cogaps_characterization.qmd are
# the CoGAPS replacement.
#
# Grid: fit CoGAPS at every rank in rank_range, 10 seeds per rank, for all
# three bases (day0, day2, combined) in a single combined job array — one
# sbatch submission instead of three, following the same "one job function
# parameterized by basis, all per-basis data passed as global objects"
# pattern R/rose_wto_rslurm.R already uses for its geneset x day grid.
# Full nIterations = 15000L (these fits are used directly for
# characterization, not just for exploratory rank selection, so they're
# run at full quality rather than a cheaper reduced-iteration sweep).
#
# rank_range is narrowed to 5:10 (rather than a wide exploratory sweep):
# the masking-CV rank search previously done in rose.qmd identified ~7 as
# the reconstruction-optimal rank (for NMF; CoGAPS is a different —
# Bayesian — model, so this is a shared prior rather than a guarantee it
# lands the same, but it's the reasonable starting neighborhood). This
# grid exists to get seed-stability error bars in that neighborhood, not
# to re-discover the rank from scratch. Widen it again if that prior needs
# revisiting.
#
# Note: R/nmf_rslurm.R/R/nmf_maskcv_rslurm.R were later widened to 5:20
# after NMF's masking-CV elbow landed outside its old 5:10 ceiling (see
# their header comments). CoGAPS's grids were deliberately NOT widened to
# match — CoGAPS fits are far more expensive (15000 iterations,
# containerized, up to 6h/fit), and CoGAPS's own masking-CV curve (see
# R/cogaps_maskcv_rslurm.R) hasn't shown the same boundary-elbow problem
# that motivated NMF's change. So the two methods' rank_range can
# legitimately diverge here; don't "fix" this back to matching without
# first checking whether CoGAPS's masking-CV elbow is actually landing on
# a boundary.
#
# alphaA/alphaP are left at CoGAPS's own defaults (0.01 each) — the
# seed x alpha stability grid that used to tune these has been removed for
# now. To add it back: after a rank is picked from this grid search's
# reconstruction-error plot, add a second job family here that fixes
# nPatterns at the chosen rank and sweeps a small alphaA/alphaP grid across
# a few seeds, scoring each alpha by mean pairwise cosine similarity of
# featureLoadings across seeds (see git history / the pre-rewrite version
# of this file for the exact original `run_cogaps_seed_dX`-style
# implementation and its `align_A`-based stability scoring).
#
# Each job returns its reconstruction error (mean squared full in-sample
# error, featureLoadings %*% t(sampleFactors) vs. the input matrix)
# alongside the fitted CogapsResult, so cogaps_characterization.qmd can
# (a) plot rank vs. reconstruction error with seed-based error bars, and
# (b) pull the exact result for whichever (basis, rank, seed) is chosen
# after inspecting that plot — no separate refit step.
#
# Run out-of-band: `Rscript R/cogaps_rslurm.R`.

library(here)
library(rslurm)
library(CoGAPS)
library(DESeq2)
library(SummarizedExperiment)

rank_range          <- 5:10
cogaps_seeds        <- c(42L, 123L, 456L, 7L, 99L, 2024L, 8675309L, 271828L, 31415L, 90210L)
cogaps_cpu_per_task <- 10L

dds     <- readRDS(here("results/dds.rds"))
vst_mat <- assay(varianceStabilizingTransformation(dds, blind = TRUE))

day0_samp <- colData(dds)$timepoint == "Day 0"
day2_samp <- colData(dds)$timepoint == "Day 2"

vst_day0     <- vst_mat[, day0_samp]
vst_day2     <- vst_mat[, day2_samp]
vst_combined <- cbind(vst_day0, vst_day2)

shift_nonneg <- function(mat) mat - matrixStats::rowMins(mat)

mat_nn_by_basis <- list(day0 = shift_nonneg(vst_day0),
                        day2 = shift_nonneg(vst_day2),
                        combined = shift_nonneg(vst_combined))

cogaps_grid_jobs <- expand.grid(
  basis = names(mat_nn_by_basis),
  rank  = rank_range,
  seed  = cogaps_seeds,
  stringsAsFactors = FALSE
)

run_cogaps_grid_job <- function(basis, rank, seed, nsets) {
  library(CoGAPS)
  mat_nn_basis <- mat_nn_by_basis[[basis]]
  params <- CogapsParams(
    nPatterns   = rank,
    nIterations = 15000L,
    seed        = seed,
    distributed = "genome-wide"
    # alphaA/alphaP intentionally left at CoGAPS defaults — see the
    # header comment above for where/how to reintroduce alpha tuning.
  )
  params <- setDistributedParams(params, nSets = nsets)
  result <- tryCatch(
    CoGAPS(mat_nn_basis, params = params, nThreads = 1),
    error = function(e) NULL
  )
  if (is.null(result)) return(list(basis = basis, rank = rank, seed = seed, mse = NA_real_, result = NULL))

  recon <- result@featureLoadings %*% t(result@sampleFactors)
  mse   <- mean((mat_nn_basis - recon)^2)
  list(basis = basis, rank = rank, seed = seed, mse = mse, result = result)
}

sjob_cogaps_grid <- slurm_apply(
  f              = run_cogaps_grid_job,
  params         = data.frame(cogaps_grid_jobs, nsets = cogaps_cpu_per_task),
  jobname        = "cogaps_grid",
  nodes          = nrow(cogaps_grid_jobs),
  cpus_per_node  = 1,
  global_objects = c("mat_nn_by_basis"),
  pkgs           = c("CoGAPS"),
  libPaths       = "/ref/mblab/software/chasem/R/4.6.0",
  rscript_path   = "Rscript",
  sh_template    = "~/.rslurm/submit_sh.txt",
  slurm_options  = list(
    mem             = "10G",
    "cpus-per-task" = cogaps_cpu_per_task,
    time            = "06:00:00",
    container       = "docker://ghcr.io/fertiglab/cogaps:sha-3b3e002"
  ),
  submit = FALSE
)

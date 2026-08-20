# CoGAPS masking-CV rank/alpha selection for cogaps_characterization.qmd
# (Day 0, Day 2, and combined Day0+Day2 bases).
#
# Companion to R/cogaps_rslurm.R's full grid, same relationship as
# R/nmf_maskcv_rslurm.R has to R/nmf_rslurm.R: the full grid scores
# reconstruction in-sample (monotonically decreasing in rank, can't select a
# rank on its own) and leaves alphaA/alphaP at CoGAPS's defaults (0.01,
# 0.01). This script instead holds out a random 1% of each basis's matrix
# entries and scores reconstruction MSE only on those held-out cells, across
# a rank x alpha grid, so BOTH the rank and the alphaA/alphaP sparsity
# parameter get an out-of-sample-informed choice rather than a default/prior
# assumption. cogaps_characterization.qmd sources its elbow-rank auto-pick
# from THIS script's output, not the full grid's in-sample curve.
#
# Masking mechanism: unlike NNLM::nnmf, CoGAPS has no documented native
# NA/missing-data handling, so held-out cells are excluded via the
# `uncertainty` matrix instead -- a standard technique for CV with
# uncertainty-weighted Bayesian NMF (the model's likelihood weights each
# cell by 1/uncertainty^2, so an enormous uncertainty on a cell makes it
# contribute ~nothing to the fit, equivalent to exclusion). Kept cells get
# CoGAPS's textbook default convention of "10% of signal, floored" as their
# uncertainty (`pmax(0.1 * data, 0.1 * mean(data))`) -- CoGAPS's own
# internal default when `uncertainty = NULL` isn't documented precisely, but
# this only needs to be a *consistent* uncertainty across every fit in this
# grid for the rank/alpha comparison to be fair, not an exact match to
# CoGAPS's internal default.
#
# alpha_range sweeps alphaA and alphaP TOGETHER (one value applied to both,
# same convention the removed pre-refactor alpha-tuning grid used -- see
# R/cogaps_rslurm.R's header comment) across CoGAPS's default (0.01) and six
# other values spanning zero sparsity to considerably more sparsity than
# default.
#
# No multi-seed sweep here (unlike R/cogaps_rslurm.R's 10 seeds/rank) -- a
# single fit per (rank, alpha) is enough to pick a rank and alpha;
# seed-stability is a separate question already answered by the full grid
# at whichever rank/alpha this script selects. nIterations is kept at full
# quality (15000) rather than a cheaper reduced sweep, so a convergence
# shortfall doesn't get mistaken for a genuine rank/alpha effect on
# held-out error -- this does mean the grid (3 bases x 6 ranks x 7 alphas =
# 126 jobs) is not cheap; that's the tradeoff for a rank/alpha choice that's
# actually informed by the same iteration budget the final characterization
# fits use.
#
# rank_range matches R/cogaps_rslurm.R's 5:10 for consistency with the rest
# of the pipeline. If the elbow found here lands on a boundary (5 or 10),
# widen rank_range and rerun rather than trusting a boundary elbow.
#
# Run out-of-band: `Rscript R/cogaps_maskcv_rslurm.R`.

library(here)
library(rslurm)
library(CoGAPS)
library(DESeq2)
library(SummarizedExperiment)

rank_range          <- 5:10
alpha_range         <- c(0.001, 0, 0.01, 0.05, 0.1, 0.15, 0.2)
mask_frac           <- 0.01
mask_seed           <- 2024L
maskcv_seed         <- 42L   # single fixed CoGAPS seed -- this grid picks rank/alpha, not seed stability
cogaps_cpu_per_task <- 10L

dds     <- readRDS(here("results/dds.rds"))
vst_mat <- assay(varianceStabilizingTransformation(dds, blind = TRUE))

day0_samp <- colData(dds)$timepoint == "Day 0"
day2_samp <- colData(dds)$timepoint == "Day 2"

vst_day0     <- vst_mat[, day0_samp]
vst_day2     <- vst_mat[, day2_samp]
vst_combined <- cbind(vst_day0, vst_day2)

shift_nonneg <- function(mat) mat - matrixStats::rowMins(mat)

mat_nn_by_basis <- list(
  day0 = shift_nonneg(vst_day0),
  day2 = shift_nonneg(vst_day2),
  combined = shift_nonneg(vst_combined)
)

# One fixed held-out mask per basis, drawn once so every (rank, alpha)
# combination is scored against identical held-out cells.
set.seed(mask_seed)
mask_idx_by_basis <- lapply(mat_nn_by_basis, function(mat) {
  sample(length(mat), round(mask_frac * length(mat)))
})

# Default-uncertainty ("10% of signal, floored") matrix per basis, with
# held-out cells overwritten to an enormous value below in the job function
# so they're effectively excluded from that job's fit.
uncertainty_by_basis <- lapply(mat_nn_by_basis, function(mat) {
  pmax(0.1 * mat, 0.1 * mean(mat))
})

cogaps_maskcv_jobs <- expand.grid(
  basis = names(mat_nn_by_basis),
  rank  = rank_range,
  alpha = alpha_range,
  stringsAsFactors = FALSE
)

run_cogaps_maskcv_job <- function(basis, rank, alpha, nsets) {
  library(CoGAPS)
  mat_nn_basis <- mat_nn_by_basis[[basis]]
  mask_idx     <- mask_idx_by_basis[[basis]]

  unc <- uncertainty_by_basis[[basis]]
  unc[mask_idx] <- 1e4 * max(mat_nn_basis)

  params <- CogapsParams(
    nPatterns   = rank,
    nIterations = 15000L,
    seed        = maskcv_seed,
    alphaA      = alpha,
    alphaP      = alpha,
    distributed = "genome-wide"
  )
  params <- setDistributedParams(params, nSets = nsets)
  result <- tryCatch(
    CoGAPS(mat_nn_basis, params = params, uncertainty = unc, nThreads = 1),
    error = function(e) NULL
  )
  if (is.null(result)) return(list(basis = basis, rank = rank, alpha = alpha, mse = NA_real_))

  recon <- result@featureLoadings %*% t(result@sampleFactors)
  mse   <- mean((mat_nn_basis[mask_idx] - recon[mask_idx])^2)
  list(basis = basis, rank = rank, alpha = alpha, mse = mse)
}

sjob_cogaps_maskcv <- slurm_apply(
  f              = run_cogaps_maskcv_job,
  params         = data.frame(cogaps_maskcv_jobs, nsets = cogaps_cpu_per_task),
  jobname        = "cogaps_maskcv",
  nodes          = nrow(cogaps_maskcv_jobs),
  cpus_per_node  = 1,
  global_objects = c("mat_nn_by_basis", "mask_idx_by_basis", "uncertainty_by_basis", "maskcv_seed"),
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

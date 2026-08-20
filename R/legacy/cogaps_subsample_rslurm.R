# CoGAPS subject-subsample stability for cogaps_characterization.qmd (Day 0,
# Day 2, and combined Day0+Day2 bases).
#
# Companion to R/nmf_subsample_rslurm.R -- same rationale, same stratified
# subject-sampling design (see that script's header comment for the full
# reasoning: paired-subject-only pool, exact class composition via
# per-class-arm sampling, approximate gender/surv_90 via a secondary
# within-class stratified draw). Duplicated here rather than shared, per
# this repo's existing convention of standalone rslurm scripts (see e.g.
# `shift_nonneg`/`read_job_dir` duplicated across R/nmf_rslurm.R,
# R/cogaps_rslurm.R, and their masking-CV/subsample counterparts).
#
# Unlike R/cogaps_maskcv_rslurm.R, this grid does NOT sweep alphaA/alphaP --
# alphaA/alphaP are left at CoGAPS's defaults (0.01 each), since the
# question here is compositional stability (does a skewed class mix change
# the patterns?), not sparsity tuning. No seed sweep either (single fit per
# basis x rank x scenario; seed stability is answered separately by the
# full grid). rank_range matches R/cogaps_rslurm.R's 5:10 (54 fits total: 3
# bases x 6 ranks x 3 scenarios), at full nIterations = 15000 so a
# convergence shortfall doesn't get mistaken for a genuine
# composition-driven difference in the patterns.
#
# Run out-of-band: `Rscript R/cogaps_subsample_rslurm.R`.

library(here)
library(rslurm)
library(CoGAPS)
library(DESeq2)
library(SummarizedExperiment)
library(dplyr)

rank_range           <- 5:10
class_scenarios      <- c(hyper10 = 0.10, hyper30 = 0.30, hyper50 = 0.50)
n_subsample_subjects <- 60L
gender_frac_target   <- 0.30
surv_frac_target      <- 0.60
subsample_seed_by_scenario <- c(hyper10 = 1010L, hyper30 = 1030L, hyper50 = 1050L)
cogaps_seed          <- 42L   # single fixed CoGAPS seed -- this grid tests composition, not seed stability
cogaps_cpu_per_task  <- 10L

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

# --- subject-level metadata, restricted to the 93 Day0+Day2-paired subjects ---
subject_of   <- colData(dds)$subject_number
sample_of    <- colnames(dds)
timepoint_of <- colData(dds)$timepoint

subj_timepoints <- tibble::tibble(subject_number = subject_of, timepoint = timepoint_of) |>
  dplyr::distinct() |>
  dplyr::count(subject_number, name = "n_timepoints")
paired_subjects <- subj_timepoints$subject_number[subj_timepoints$n_timepoints == 2]

subj_meta_paired <- tibble::tibble(
  subject_number = subject_of, class = colData(dds)$class,
  gender = colData(dds)$gender, surv_90 = colData(dds)$surv_90
) |>
  dplyr::distinct() |>
  dplyr::filter(subject_number %in% paired_subjects)

stopifnot(nrow(subj_meta_paired) == length(paired_subjects))

# --- hierarchical stratified subject sampler (class exact, gender/surv_90 approximate) ---
lr_round <- function(raw, total) {
  fl  <- floor(raw)
  rem <- raw - fl
  n_add <- total - sum(fl)
  if (n_add > 0) {
    idx <- order(rem, decreasing = TRUE)[seq_len(n_add)]
    fl[idx] <- fl[idx] + 1
  }
  fl
}

sample_within_class <- function(pool_meta, n_target, gender_frac, surv_frac) {
  if (n_target <= 0) return(character(0))
  cell_defs <- expand.grid(gender = c("Female", "Male"), surv_90 = c(FALSE, TRUE),
                            stringsAsFactors = FALSE)
  cell_defs$p_cell <- ifelse(cell_defs$gender == "Female", gender_frac, 1 - gender_frac) *
    ifelse(cell_defs$surv_90, surv_frac, 1 - surv_frac)
  cell_defs$target <- lr_round(cell_defs$p_cell * n_target, n_target)
  cell_defs$key    <- paste(cell_defs$gender, cell_defs$surv_90)
  pool_meta$key    <- paste(pool_meta$gender, pool_meta$surv_90)

  sampled <- character(0)
  for (i in seq_len(nrow(cell_defs))) {
    pool_i <- pool_meta$subject_number[pool_meta$key == cell_defs$key[i]]
    take   <- min(length(pool_i), cell_defs$target[i])
    if (take > 0) sampled <- c(sampled, sample(pool_i, take))
  }
  shortfall <- n_target - length(sampled)
  if (shortfall > 0) {
    remaining <- setdiff(pool_meta$subject_number, sampled)
    take <- min(shortfall, length(remaining))
    if (take > 0) sampled <- c(sampled, sample(remaining, take))
  }
  sampled
}

stratified_subject_sample <- function(subj_meta, class_frac, n_subjects,
                                       gender_frac, surv_frac, seed) {
  set.seed(seed)
  n_hyper <- round(class_frac * n_subjects)
  n_hypo  <- n_subjects - n_hyper
  hyper_pool <- subj_meta |> dplyr::filter(class == "Hyperinflammatory")
  hypo_pool  <- subj_meta |> dplyr::filter(class == "Hypoinflammatory")
  stopifnot(n_hyper <= nrow(hyper_pool), n_hypo <= nrow(hypo_pool))
  c(sample_within_class(hyper_pool, n_hyper, gender_frac, surv_frac),
    sample_within_class(hypo_pool,  n_hypo,  gender_frac, surv_frac))
}

subject_subsamples <- lapply(names(class_scenarios), function(nm) {
  stratified_subject_sample(subj_meta_paired, class_scenarios[[nm]], n_subsample_subjects,
                             gender_frac_target, surv_frac_target,
                             subsample_seed_by_scenario[[nm]])
})
names(subject_subsamples) <- names(class_scenarios)

achieved_composition <- function(subjects) {
  m <- subj_meta_paired |> dplyr::filter(subject_number %in% subjects)
  list(class_frac_achieved  = mean(m$class == "Hyperinflammatory"),
       gender_frac_achieved = mean(m$gender == "Female"),
       surv_frac_achieved   = mean(m$surv_90 == TRUE))
}

for (nm in names(subject_subsamples)) {
  ac <- achieved_composition(subject_subsamples[[nm]])
  cat(sprintf("%s: target=%.2f | achieved hyper=%.2f female=%.2f survTRUE=%.2f | n=%d\n",
              nm, class_scenarios[[nm]], ac$class_frac_achieved, ac$gender_frac_achieved,
              ac$surv_frac_achieved, length(subject_subsamples[[nm]])))
}

# --- per-(basis, scenario) column sets ---
cols_for_subjects <- function(subjects, timepoints) {
  sample_of[subject_of %in% subjects & timepoint_of %in% timepoints]
}

basis_timepoints <- list(day0 = "Day 0", day2 = "Day 2", combined = c("Day 0", "Day 2"))

cols_by_basis_scenario <- lapply(names(mat_nn_by_basis), function(b) {
  setNames(lapply(names(subject_subsamples), function(nm) {
    cols_for_subjects(subject_subsamples[[nm]], basis_timepoints[[b]])
  }), names(subject_subsamples))
})
names(cols_by_basis_scenario) <- names(mat_nn_by_basis)

cogaps_subsample_jobs <- expand.grid(
  basis    = names(mat_nn_by_basis),
  rank     = rank_range,
  scenario = names(class_scenarios),
  stringsAsFactors = FALSE
)

run_cogaps_subsample_job <- function(basis, rank, scenario, nsets) {
  library(CoGAPS)
  cols   <- cols_by_basis_scenario[[basis]][[scenario]]
  mat_nn <- mat_nn_by_basis[[basis]][, cols]

  params <- CogapsParams(
    nPatterns   = rank,
    nIterations = 15000L,
    seed        = cogaps_seed,
    distributed = "genome-wide"
    # alphaA/alphaP intentionally left at CoGAPS defaults -- this grid tests
    # class-composition stability, not sparsity tuning (see
    # R/cogaps_maskcv_rslurm.R for the alpha sweep).
  )
  params <- setDistributedParams(params, nSets = nsets)
  result <- tryCatch(
    CoGAPS(mat_nn, params = params, nThreads = 1),
    error = function(e) NULL
  )

  ac <- achieved_composition(subject_subsamples[[scenario]])
  if (is.null(result)) {
    return(list(basis = basis, rank = rank, scenario = scenario,
                class_frac_target = class_scenarios[[scenario]],
                class_frac_achieved = ac$class_frac_achieved,
                gender_frac_achieved = ac$gender_frac_achieved,
                surv_frac_achieved = ac$surv_frac_achieved,
                subjects = subject_subsamples[[scenario]], result = NULL))
  }
  list(basis = basis, rank = rank, scenario = scenario,
       class_frac_target = class_scenarios[[scenario]],
       class_frac_achieved = ac$class_frac_achieved,
       gender_frac_achieved = ac$gender_frac_achieved,
       surv_frac_achieved = ac$surv_frac_achieved,
       subjects = subject_subsamples[[scenario]], result = result)
}

sjob_cogaps_subsample <- slurm_apply(
  f              = run_cogaps_subsample_job,
  params         = data.frame(cogaps_subsample_jobs, nsets = cogaps_cpu_per_task),
  jobname        = "cogaps_subsample",
  nodes          = nrow(cogaps_subsample_jobs),
  cpus_per_node  = 1,
  global_objects = c("mat_nn_by_basis", "cols_by_basis_scenario", "subject_subsamples",
                      "subj_meta_paired", "class_scenarios", "achieved_composition",
                      "cogaps_seed"),
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

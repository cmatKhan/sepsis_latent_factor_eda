# NMF subject-subsample stability for nmf_characterization.qmd (Day 0, Day 2,
# and combined Day0+Day2 bases).
#
# A third stability check alongside R/nmf_rslurm.R's full grid (seed
# stability) and R/nmf_maskcv_rslurm.R's masking-CV grid (rank selection):
# does the fitted pattern set hold up when the fitting cohort's *class*
# composition (Hyperinflammatory vs. Hypoinflammatory, `class` in
# data/rose_se.rds -- the primary grouping variable throughout this repo)
# is skewed away from the cohort's roughly 50/50 split? Three deliberately
# biased subject subsamples are drawn -- 10%, 30%, 50% Hyperinflammatory --
# while holding gender and surv_90 close to their overall cohort proportions
# (~30% Female / ~40% surv_90 == FALSE) so a skewed class mix isn't
# confounded with a skewed gender or survival mix.
#
# Subsampling is by SUBJECT, not by sample column, and the sampling pool is
# restricted to the 93 subjects with BOTH a Day 0 and a Day 2 sample (of 128
# total) -- so "sample by subject_number" always yields both timepoints for
# a selected subject, per the class/gender/surv_90 composition request.
# This pool is shared across all three bases: the day0/day2 bases just pull
# that subject's single relevant column, while the combined basis pulls
# both. One subsample is drawn per class-mix scenario (not per basis, not
# per rank) and reused everywhere that scenario appears, so a day0 vs. day2
# vs. combined comparison at a fixed scenario is looking at literally the
# same subjects.
#
# Quota construction: class composition is hit EXACTLY by drawing
# round(class_frac * n) Hyperinflammatory and the rest Hypoinflammatory
# directly from each class's subject pool (38 Hyper / 55 Hypo available
# among the 93 paired subjects -- comfortably above what's needed at
# n_subsample_subjects = 60 for any of the three scenarios, including the
# tight case class_frac = 0.10 -> 54 of 55 available Hypo subjects). Gender
# and surv_90 are then approximated via a secondary 2x2 (gender x surv_90)
# stratified draw WITHIN each class arm (target cell counts from the
# independence-product formula, largest-remainder rounded, with same-class-
# arm-only top-up on shortfall) -- this keeps the class marginal exact even
# when a gender/surv_90 sub-cell is short (e.g. only 3 Hypoinflammatory
# Female surv_90=FALSE subjects exist), at the cost of gender/surv_90 being
# "roughly" rather than exactly on target. Validated against the real data:
# all three scenarios land within a few points of the 30%/60% gender/surv_90
# targets at n_subsample_subjects = 60.
#
# No alpha applies to NMF; no seed sweep (single fit per basis x rank x
# scenario -- seed stability is a separate question already answered by the
# full grid); rank_range matches R/nmf_rslurm.R's 5:10 for consistency with
# the rest of the pipeline (54 fits total: 3 bases x 6 ranks x 3 scenarios).
#
# Run out-of-band: `Rscript R/nmf_subsample_rslurm.R`.

library(here)
library(rslurm)
library(NNLM)
library(DESeq2)
library(SummarizedExperiment)
library(dplyr)

rank_range <- 5:10
class_scenarios <- c(hyper10 = 0.10, hyper30 = 0.30, hyper50 = 0.50)
n_subsample_subjects <- 60L
gender_frac_target <- 0.30 # cohort-wide ~30% Female
surv_frac_target <- 0.60 # cohort-wide ~60% surv_90 == TRUE
subsample_seed_by_scenario <- c(hyper10 = 1010L, hyper30 = 1030L, hyper50 = 1050L)

dds <- readRDS(here("results/dds.rds"))
vst_mat <- assay(varianceStabilizingTransformation(dds, blind = TRUE))

day0_samp <- colData(dds)$timepoint == "Day 0"
day2_samp <- colData(dds)$timepoint == "Day 2"

vst_day0 <- vst_mat[, day0_samp]
vst_day2 <- vst_mat[, day2_samp]
vst_combined <- cbind(vst_day0, vst_day2)

shift_nonneg <- function(mat) mat - matrixStats::rowMins(mat)

# Nonneg-shifted per-basis matrices, shift computed on the FULL basis (not
# per-subsample) so the same shift constant applies across every scenario --
# subsample column-subsetting happens after this, per job.
mat_nn_by_basis <- list(
  day0 = shift_nonneg(vst_day0),
  day2 = shift_nonneg(vst_day2),
  combined = shift_nonneg(vst_combined)
)

# --- subject-level metadata, restricted to the 93 Day0+Day2-paired subjects ---
subject_of <- colData(dds)$subject_number
sample_of <- colnames(dds)
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
  fl <- floor(raw)
  rem <- raw - fl
  n_add <- total - sum(fl)
  if (n_add > 0) {
    idx <- order(rem, decreasing = TRUE)[seq_len(n_add)]
    fl[idx] <- fl[idx] + 1
  }
  fl
}

sample_within_class <- function(pool_meta, n_target, gender_frac, surv_frac) {
  if (n_target <= 0) {
    return(character(0))
  }
  cell_defs <- expand.grid(
    gender = c("Female", "Male"), surv_90 = c(FALSE, TRUE),
    stringsAsFactors = FALSE
  )
  cell_defs$p_cell <- ifelse(cell_defs$gender == "Female", gender_frac, 1 - gender_frac) *
    ifelse(cell_defs$surv_90, surv_frac, 1 - surv_frac)
  cell_defs$target <- lr_round(cell_defs$p_cell * n_target, n_target)
  cell_defs$key <- paste(cell_defs$gender, cell_defs$surv_90)
  pool_meta$key <- paste(pool_meta$gender, pool_meta$surv_90)

  sampled <- character(0)
  for (i in seq_len(nrow(cell_defs))) {
    pool_i <- pool_meta$subject_number[pool_meta$key == cell_defs$key[i]]
    take <- min(length(pool_i), cell_defs$target[i])
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
  n_hypo <- n_subjects - n_hyper
  hyper_pool <- subj_meta |> dplyr::filter(class == "Hyperinflammatory")
  hypo_pool <- subj_meta |> dplyr::filter(class == "Hypoinflammatory")
  stopifnot(n_hyper <= nrow(hyper_pool), n_hypo <= nrow(hypo_pool))
  c(
    sample_within_class(hyper_pool, n_hyper, gender_frac, surv_frac),
    sample_within_class(hypo_pool, n_hypo, gender_frac, surv_frac)
  )
}

subject_subsamples <- lapply(names(class_scenarios), function(nm) {
  stratified_subject_sample(
    subj_meta_paired, class_scenarios[[nm]], n_subsample_subjects,
    gender_frac_target, surv_frac_target,
    subsample_seed_by_scenario[[nm]]
  )
})
names(subject_subsamples) <- names(class_scenarios)

achieved_composition <- function(subjects) {
  m <- subj_meta_paired |> dplyr::filter(subject_number %in% subjects)
  list(
    class_frac_achieved = mean(m$class == "Hyperinflammatory"),
    gender_frac_achieved = mean(m$gender == "Female"),
    surv_frac_achieved = mean(m$surv_90 == TRUE)
  )
}

for (nm in names(subject_subsamples)) {
  ac <- achieved_composition(subject_subsamples[[nm]])
  cat(sprintf(
    "%s: target=%.2f | achieved hyper=%.2f female=%.2f survTRUE=%.2f | n=%d\n",
    nm, class_scenarios[[nm]], ac$class_frac_achieved, ac$gender_frac_achieved,
    ac$surv_frac_achieved, length(subject_subsamples[[nm]])
  ))
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

nmf_subsample_jobs <- expand.grid(
  basis = names(mat_nn_by_basis),
  rank = rank_range,
  scenario = names(class_scenarios),
  stringsAsFactors = FALSE
)

# scenario_metadata <- lapply(
#   names(subject_subsamples),
#   function(nm) {
#     subj_meta_paired |>
#       dplyr::filter(subject_number %in% subject_subsamples[[nm]]) |>
#       dplyr::mutate(scenario = nm)
#   }
# )
# names(scenario_metadata) <- names(subject_subsamples)

run_nmf_subsample_job <- function(basis, rank, scenario) {
  library(NNLM)
  cols <- cols_by_basis_scenario[[basis]][[scenario]]
  mat_nn <- mat_nn_by_basis[[basis]][, cols]

  fit <- nnmf(mat_nn, k = rank, max.iter = 10000, verbose = 0L)
  rownames(fit$W) <- rownames(mat_nn)
  colnames(fit$W) <- paste0("Pattern_", seq_len(ncol(fit$W)))
  colnames(fit$H) <- colnames(mat_nn)
  rownames(fit$H) <- colnames(fit$W)

  ac <- achieved_composition(subject_subsamples[[scenario]])
  list(
    basis = basis, rank = rank, scenario = scenario,
    class_frac_target = class_scenarios[[scenario]],
    class_frac_achieved = ac$class_frac_achieved,
    gender_frac_achieved = ac$gender_frac_achieved,
    surv_frac_achieved = ac$surv_frac_achieved,
    subjects = subject_subsamples[[scenario]],
    W = fit$W, H = fit$H
  )
}

sjob_nmf_subsample <- slurm_apply(
  f = run_nmf_subsample_job,
  params = nmf_subsample_jobs,
  jobname = "nmf_subsample",
  nodes = nrow(nmf_subsample_jobs),
  cpus_per_node = 1,
  global_objects = c(
    "mat_nn_by_basis", "cols_by_basis_scenario", "subject_subsamples",
    "subj_meta_paired", "class_scenarios", "achieved_composition"
  ),
  pkgs = c("NNLM"),
  libPaths = "/ref/mblab/software/chasem/R/4.6.0",
  rscript_path = "Rscript",
  sh_template = "~/.rslurm/submit_sh.txt",
  slurm_options = list(
    mem = "8G",
    "cpus-per-task" = 1,
    time = "01:00:00",
    out = "logs/%a.out",
    container = "docker://rocker/tidyverse:latest"
  ),
  submit = FALSE
)

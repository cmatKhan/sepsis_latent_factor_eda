# Phase 0 POC, step 4 (see ../README.md). Only run this after
# ../01_raw_crew_poc.R has passed. Drives R/methods/pca.R's REAL,
# unmodified pca_registry through a targets + crew.cluster pipeline, to
# confirm the full stack (not just raw crew) works together on the real
# cluster. No Apptainer -- workers get R via `spack load r@4.6.1`, same as
# ../01_raw_crew_poc.R (see its header comment for why).
#
# Expected to be run interactively from an RStudio Server session on the
# cluster, with this project open. Because here::here() (via
# sepsis_timecourse_factors_eda.Rproj) resolves to the repo root
# regardless of RStudio's starting directory, `targets::tar_make()` can be
# run from the R console with the working directory set to EITHER this
# directory (02_targets_poc/, so _targets.R is found automatically) or the
# repo root (passing script/store paths explicitly) -- e.g. from the R
# console:
#   setwd(here::here("poc/phase0_slurm_apptainer/02_targets_poc"))
#   targets::tar_make()

# ---- EDIT THESE (kept in sync with ../01_raw_crew_poc.R) ----
# Must be set before any library() call, AND before tar_option_set(),
# since tar_make()'s default callr subprocess re-sources this whole file
# from scratch.
user_lib_path    <- "/ref/mblab/software/chasem/R-rstudio/4.6"
r_spec           <- "r@4.6.1"
partition        <- NULL
script_directory <- here::here("poc/phase0_slurm_apptainer/02_targets_poc/job_scripts")
# ---------------------------------------------------------------

.libPaths(c(user_lib_path, .libPaths()))

library(targets)
library(crew)
library(crew.cluster)
library(here)

if (!nzchar(Sys.which("sbatch"))) {
  stop("`sbatch` not found on PATH -- see the sanity check in ../01_raw_crew_poc.R; ",
       "run that step first if you haven't already.")
}

# Reuses the real project code unmodified -- the whole point of this POC
# is proving targets/crew.cluster can drive the EXISTING method-registry
# convention (R/lib/method_registry.R) with zero changes to it.
source(here::here("R/methods/pca.R"))

# See ../01_raw_crew_poc.R's header comment for why script_lines is built
# this way (spack load + explicit R_LIBS_USER export, matching
# crew_options_slurm()'s own documented `module load R` pattern).
build_spack_script_lines <- function(r_spec, user_lib_path) {
  c(
    sprintf('eval $(spack load --sh %s)', r_spec),
    sprintf('export R_LIBS_USER="%s"', user_lib_path)
  )
}

controller <- crew_controller_slurm(
  name = "pca",
  workers = 2,
  options_cluster = crew_options_slurm(
    script_lines             = build_spack_script_lines(r_spec, user_lib_path),
    memory_gigabytes_per_cpu = 4,
    cpus_per_task            = 1,
    time_minutes             = 30,
    partition                = partition,
    script_directory         = script_directory,
    verbose                  = TRUE
  ),
  seconds_launch = 120
)

tar_option_set(controller = crew_controller_group(controller))

list(
  # deployment = "main": runs in tar_make()'s own process, NOT dispatched
  # to a crew worker -- confirmed necessary (see plan file): once a global
  # controller is set, "worker" is the default deployment for EVERY
  # target, so lightweight/structural targets need this set explicitly or
  # they run pointlessly on a compute node too.
  tar_target(pca_grid, pca_registry$build_grid(pca_registry$defaults), deployment = "main"),

  # Dynamic branching over grid rows -- replaces rslurm's array-task/
  # max_array_size batching; crew's auto-scaling worker pool absorbs the
  # granularity problem instead.
  tar_target(
    pca_fit,
    {
      mat <- matrix(rnorm(2000), nrow = 100)  # POC stand-in for a real loaded matrix
      pca_registry$fn(rank = pca_grid$rank)
    },
    pattern = map(pca_grid),
    resources = tar_resources(crew = tar_resources_crew(controller = "pca"))
  ),

  tar_target(pca_summary, {
    data.frame(rank = vapply(pca_fit, function(x) x$rank, numeric(1)),
               mse  = vapply(pca_fit, function(x) x$mse, numeric(1)))
  }, deployment = "main")
)

# Phase 0 POC, step 3 (see README.md). Submits ONE real SLURM job via
# crew_controller_slurm() -- no Apptainer, no `targets` yet, just raw
# `crew`. Workers get R via `spack load r@4.6.1` (the same `script_lines`
# mechanism crew.cluster's own docs show for `module load R`), since this
# cluster's SLURM client tools and R are both available directly via spack
# once a worker's job script runs -- no container boundary to cross at
# all. This is simpler than, and was adopted in place of, an earlier
# Apptainer-wrapping version of this POC (see git history / the plan
# file): it's a closer match to how `targets`+`crew.cluster` are actually
# meant to be used, and sidesteps an entire class of container-vs-host
# library-compatibility issues that approach kept surfacing.
#
# Expected to be run interactively from an RStudio Server session on the
# cluster, with this project open (so here::here() resolves to the repo
# root via sepsis_timecourse_factors_eda.Rproj regardless of whatever
# directory RStudio happened to start in) -- source this file, or run it
# line by line, from the R console. Run 00_preflight_check.sh from
# RStudio's Terminal tab first.
#
# Success looks like a "POC PASSED" message within ~1-2 minutes. Report
# back the full console output either way -- see README.md.

# ---- EDIT THESE ----
# Where packages (crew/crew.cluster/targets/etc.) live -- used for BOTH
# this session (prepended to .libPaths() below) and the worker (exported
# as R_LIBS_USER in script_lines further down, since a fresh SLURM job's
# shell starts with none of this session's state and needs it set
# explicitly).
user_lib_path    <- "/ref/mblab/software/chasem/R-rstudio/4.6"
r_spec           <- "r@4.6.1"   # the spack spec `spack load --sh` resolves
partition        <- NULL         # today's rslurm templates don't set one either
script_directory <- here::here("poc/phase0_slurm_apptainer/job_scripts")
# ---------------------

.libPaths(c(user_lib_path, .libPaths()))

library(crew)
library(crew.cluster)
library(here)

# ---- sanity checks (fail fast, before spending 120s waiting) ----
if (!nzchar(Sys.which("sbatch"))) {
  stop("`sbatch` not found on PATH in this session -- crew_controller_slurm() submits\n",
       "jobs by calling sbatch directly from wherever this R process runs, so this\n",
       "session's node needs it. Load whatever spack/environment setup gives you sbatch\n",
       "(the same one that already makes `sbatch` available in your reconfigured\n",
       "RStudio launch) before sourcing this file.")
}
detected_host <- utils::head(nanonext::ip_addr(), n = 1L)
message("sbatch found at: ", Sys.which("sbatch"))
message("crew will advertise this session as reachable at: ", detected_host, " (auto-detected via nanonext::ip_addr()).")
message("If this node is multi-homed (e.g. a separate public-facing and internal\n",
        "  cluster-fabric interface), and compute nodes can't reach ", detected_host, ",\n",
        "  pass `host = \"<internal-ip>\"` explicitly to crew_controller_slurm() below\n",
        "  instead of relying on auto-detection.")

#' crew.cluster's generated job script is the `#SBATCH` header lines,
#' followed by `script_lines` verbatim, followed by one final line:
#' `Rscript -e '<crew::crew_worker(...) call>'` (confirmed by reading
#' crew.cluster:::crew_class_launcher_cluster$public_methods$launch_worker's
#' source). So `script_lines` just needs to put the right `Rscript` on
#' PATH and point it at the right library before that line runs -- exactly
#' the same shape as crew_options_slurm()'s own documented
#' `script_lines = "module load R"` example, just via spack instead of
#' environment modules, and with R_LIBS_USER exported explicitly rather
#' than assumed, since a fresh job's shell won't have this session's
#' .libPaths() state.
build_spack_script_lines <- function(r_spec, user_lib_path) {
  c(
    sprintf('eval $(spack load --sh %s)', r_spec),
    sprintf('export R_LIBS_USER="%s"', user_lib_path)
  )
}

options_cluster <- crew_options_slurm(
  script_lines             = build_spack_script_lines(r_spec, user_lib_path),
  memory_gigabytes_per_cpu = 4,     # matches EARLI_config.yml's slurm.pca.mem (4G)
  cpus_per_task            = 1,     # matches slurm.pca.cpus_per_task
  time_minutes             = 30,    # matches slurm.pca.time ("00:30:00")
  partition                = partition,
  script_directory         = script_directory,
  verbose                  = TRUE
)

controller <- crew_controller_slurm(
  name           = "phase0_poc",
  workers        = 1,
  options_cluster = options_cluster,
  seconds_launch = 120   # generous -- a real queued job may take longer to start than a local fake one
)

controller$start()

message("Submitted. Pushing a task that only succeeds if it's really running on a real")
message("SLURM-allocated compute node, with the right R and the right libPaths...")
controller$push(
  command = {
    fit <- prcomp(t(mat), rank. = 3, center = TRUE, scale. = FALSE)
    list(
      node      = Sys.info()[["nodename"]],
      r_version = R.version.string,
      lib_paths = .libPaths(),
      rank      = ncol(fit$x),
      mse       = mean((t(mat) - fit$x %*% t(fit$rotation))^2)
    )
  },
  data = list(mat = matrix(rnorm(500), nrow = 50))
)

deadline <- Sys.time() + 180
result <- NULL
while (Sys.time() < deadline) {
  controller$wait(seconds_timeout = 5)
  out <- controller$pop()
  if (!is.null(out)) { result <- out; break }
}

if (is.null(result)) {
  cat("\nPOC INCONCLUSIVE: no result within 180s.\n")
  cat("Most likely cause: the compute node can't connect back to this process over\n")
  cat("the network, or the SLURM job never started (check `squeue`, and the generated\n")
  cat("script/log files under: ", script_directory, "\n")
  controller$terminate()
  quit(status = 1)
}

cat("\n=== raw result row ===\n")
print(result)

submitting_node <- Sys.info()[["nodename"]]
cat("\nSubmitting (RStudio Server) node: ", submitting_node, "\n")

ok <- isTRUE(result$error %in% c(NA_character_)) || is.na(result$error)
if (!ok) {
  cat("\nPOC FAILED: task errored -- see 'error'/'trace' columns above.\n")
  controller$terminate()
  quit(status = 1)
}

payload <- result$result[[1]]
cat("\nWorker-reported compute node: ", payload$node, "\n")
cat("Worker-reported R version: ", payload$r_version, "\n")
cat("Worker-reported .libPaths(): ", paste(payload$lib_paths, collapse = ", "), "\n")

if (identical(payload$node, submitting_node)) {
  cat("\nWARNING: worker ran on the SAME node as this session -- that's unexpected for a\n")
  cat("real SLURM allocation (job may not have actually been scheduled elsewhere; check\n")
  cat("`sacct`/`squeue`). Not necessarily a failure if this cluster's queue legitimately\n")
  cat("placed the job back on this same node, but worth a second look.\n")
}
if (!(user_lib_path %in% payload$lib_paths)) {
  cat("\nWARNING: user_lib_path did not show up in the worker's .libPaths() -- the\n")
  cat("R_LIBS_USER export in script_lines may not have taken effect.\n")
}

cat("\nPOC PASSED: crew.cluster + spack-based script_lines works on this cluster.\n")
controller$terminate()

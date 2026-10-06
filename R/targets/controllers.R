# crew controllers for the targets pipeline: one crew_controller_slurm()
# per entry in config/pipeline.yml's `controllers:`, combined with
# crew_controller_group() and selected per target via
# tar_resources_crew(controller = "<name>") -- the "heterogeneous workers"
# pattern from the targets manual's crew chapter.

#' Read config/pipeline.yml
#'
#' @param path Path to the pipeline config.
#' @return The parsed config (datasets, data root, DB path, cluster setup,
#'   controllers).
read_pipeline_config <- function(path = "config/pipeline.yml") {
  pipeline <- yaml::read_yaml(path)
  stopifnot(
    "config/pipeline.yml needs a `datasets` list" = length(pipeline$datasets) > 0,
    "config/pipeline.yml needs a `controllers` block" = length(pipeline$controllers) > 0
  )
  pipeline
}

#' Shell lines a worker's SLURM script runs before starting R
#'
#' crew.cluster's `script_lines` hook (its docs' `module load R` example): the
#' spack activation `rstudio_spack.sbatch` also uses -- or the config's
#' `cluster.setup_lines` on a cluster without spack -- then `cd` to the project
#' so R starts renv from `.Rprofile`.
#'
#' @param cluster The pipeline config's `cluster:` block.
#' @param project_root Project directory.
#' @return Character vector of shell lines.
worker_script_lines <- function(cluster, project_root = getwd()) {
  setup <- if (!is.null(cluster$setup_lines)) {
    unlist(cluster$setup_lines)
  } else {
    view <- file.path(cluster$spack_env, "view")
    c(
      sprintf('eval "$(spack env activate --sh %s)"', cluster$spack_env),
      # rstudio_spack.sbatch sets this too: spack env activate doesn't, and
      # packages compiled against the view's libraries need it to load.
      sprintf('export LD_LIBRARY_PATH="%s/lib:%s/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"', view, view),
      "unset R_LIBS R_LIBS_USER R_LIBS_SITE"
    )
  }
  c(setup, sprintf('cd "%s"', project_root))
}

#' A controller's worker limit when config/pipeline.yml sets none: no cap
#' in practice (above SLURM's MaxArraySize, 20001), so SLURM's scheduler
#' and per-user limits decide what runs when.
UNCAPPED_WORKERS <- 50000L

#' Build the pipeline's crew controllers
#'
#' One `crew_controller_slurm()` per `controllers:` entry (plus
#' `<name>_backup` where a `backup:` block is given), combined with
#' `crew_controller_group()`; targets pick one with
#' `tar_resources_crew(controller = "<name>")`. A controller without `workers:`
#' gets UNCAPPED_WORKERS (docs: Running the workflow, "SLURM resources").
#'
#' @param pipeline The parsed pipeline config.
#' @param project_root Project directory (worker logs go to `logs/crew/`).
#' @return A crew controller group.
build_controller_group <- function(pipeline, project_root = getwd()) {
  log_dir <- file.path(project_root, "logs", "crew")
  dir.create(file.path(log_dir, "scripts"), recursive = TRUE, showWarnings = FALSE)
  script_lines <- worker_script_lines(pipeline$cluster, project_root)

  slurm_controller <- function(name, spec, backup = NULL) {
    crew.cluster::crew_controller_slurm(
      name = name,
      workers = spec$workers %||% UNCAPPED_WORKERS,
      seconds_idle = spec$seconds_idle %||% 30,
      # Soft wall time: after this the worker finishes its current task and
      # exits, and crew launches a fresh one if tasks remain -- so a worker
      # isn't killed by SLURM's --time mid-fit. Needs
      # time_minutes * 60 >= seconds_wall + the longest single fit.
      seconds_wall = spec$seconds_wall %||% (spec$time_minutes * 60 / 2),
      # tasks_max = 1: the worker runs one task and exits -- for methods
      # whose single fits can take hours (spca), so a long fit never starts
      # late in a worker's life.
      tasks_max = spec$tasks_max %||% Inf,
      crashes_max = spec$crashes_max %||% 5L,
      backup = backup,
      options_cluster = crew.cluster::crew_options_slurm(
        script_lines = script_lines,
        script_directory = file.path(log_dir, "scripts"),
        log_output = file.path(log_dir, "%x_%j.out"),
        log_error = NULL,  # stderr goes to log_output too
        memory_gigabytes_required = spec$memory_gigabytes_required,
        cpus_per_task = spec$cpus_per_task,
        time_minutes = spec$time_minutes,
        partition = pipeline$cluster$partition
      )
    )
  }

  controllers <- list()
  for (name in names(pipeline$controllers)) {
    spec <- pipeline$controllers[[name]]
    # A `backup:` block (crew's backup controllers, see the "Backup
    # controllers" section of crew's groups vignette): a task whose worker
    # crashes `crashes_max` times in a row on this controller -- e.g. SLURM
    # killing it for exceeding --mem -- is retried on "<name>_backup",
    # whose spec is this one with the `backup:` fields overridden. targets
    # (>= 1.10.0.9002) does the retry automatically.
    backup <- NULL
    if (!is.null(spec$backup)) {
      backup_spec <- utils::modifyList(spec[setdiff(names(spec), c("backup", "crashes_max"))], spec$backup)
      backup <- slurm_controller(paste0(name, "_backup"), backup_spec)
    }
    controllers[[name]] <- slurm_controller(name, spec, backup)
    if (!is.null(backup)) controllers[[paste0(name, "_backup")]] <- backup
  }
  do.call(crew::crew_controller_group, unname(controllers))
}

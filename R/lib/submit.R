# Generic slurm_apply()/slurm_call() submission wrapper. Every method's
# orchestrator script builds a job grid + job function, then calls this
# instead of hand-assembling slurm_options each time.

library(rslurm)

if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

# write_submit_all_script() -- see R/lib/submit_all_script.R.
source(here::here("R/lib/submit_all_script.R"))

# Path each job's bundle is bind-mounted to inside its container, once
# apptainer performs the bind (see `build_apptainer_rscript_path()` and
# R/README.md's squashfs section) -- arbitrary but must agree with nothing
# else, since nothing outside this function needs to know it.
.rslurm_squashfs_mount <- "/mnt/rslurm_bundle"

#' Build the `rscript_path` passed to slurm_apply()/slurm_call(): a full
#' `apptainer run ...` invocation (not just "Rscript") that runs
#' `rscript_cmd` inside `container`, with its working directory set to
#' this job family's directory inside the squashfs archive bind-mounted
#' via `$SQUASHFS_PATH` (a plain shell-variable reference, left
#' unexpanded here -- it's resolved when the generated submit.sh actually
#' runs, using whatever `submit_all_<dataset_id>.sh` exported at
#' submission time; see R/lib/submit_all_script.R). `apptainer` itself
#' must already be on PATH by the time this line runs -- submit_sh.txt/
#' submit_single_sh.txt take care of that via `spack load`.
#'
#' `lib_paths` (typically `cluster_cfg$libPaths`) is ALSO explicitly
#' bind-mounted here (source == destination, so paths match exactly what
#' rslurm's `libPaths` argument tells `slurm_run.R` to add to
#' `.libPaths()`). This is necessary, not just belt-and-suspenders:
#' apptainer does not bind arbitrary host paths like `/ref/...` by
#' default, and `.libPaths()` silently *drops* any path that isn't a real,
#' readable directory from inside the container -- confirmed on HTCF via
#' `apptainer exec <image> Rscript -e '.libPaths(c(.libPaths(), "/ref/..."));
#' print(.libPaths())'`: without an explicit bind, the appended path
#' vanishes and installed packages there (e.g. NNLM for NMF) are never
#' found, even though rslurm dutifully generated the `.libPaths()` line.
#'
#' No `--cleanenv`: the container must inherit the job's full environment
#' -- Slurm sets `SLURM_ARRAY_TASK_ID` (which `slurm_run.R` reads via
#' `Sys.getenv()` to know which rows of `params.RDS` this array task
#' handles) on the job's own process, not via anything we explicitly
#' `--env`-forward. An earlier version of this function added `--cleanenv`
#' to investigate an unrelated WGCNA container failure; that turned out
#' not to be the cause (confirmed by direct testing: identical failure
#' with or without it) and it silently broke every array job instead
#' (`SLURM_ARRAY_TASK_ID` vanishing inside `--cleanenv` produced
#' `Sys.getenv() == ""` -> `NA` -> `.rslurm_istart:.rslurm_iend : NA/NaN
#' argument`). Don't reintroduce `--cleanenv` without also explicitly
#' `--env`-forwarding every `SLURM_*` variable any job function might
#' need, not just the ones we happen to already know about.
build_apptainer_rscript_path <- function(container, jobname, rscript_cmd, lib_paths = character(0)) {
  binds <- c(
    sprintf('-B "$SQUASHFS_PATH:%s:image-src=/"', .rslurm_squashfs_mount),
    if (length(lib_paths) > 0) sprintf('-B "%s:%s"', lib_paths, lib_paths)
  )
  paste(c(
    "apptainer run",
    sprintf('--pwd "%s/_rslurm_%s"', .rslurm_squashfs_mount, jobname),
    binds,
    container,
    rscript_cmd
  ), collapse = " ")
}

#' Submit one job family.
#'
#' @param f job function (as passed to slurm_apply/slurm_call)
#' @param jobs_df data.frame of per-job parameters (params for slurm_apply);
#'   pass NULL (or a zero-row/zero-col data.frame) for a single job with no
#'   grid, which dispatches to slurm_call instead.
#' @param jobname passed through to rslurm
#' @param global_objects character vector of object names f depends on
#' @param pkgs character vector of packages f depends on
#' @param cluster_cfg one method's (or network backend's) entry from
#'   cluster_config.yml -- must have mem/cpus_per_task/time/container/
#'   libPaths/sh_template/rscript_path, and may optionally have
#'   sh_template_single (used instead of sh_template when this family
#'   dispatches to slurm_call -- see param `jobs_df` -- falls back to
#'   sh_template if not given). `container` is not passed through as a
#'   native `--container=` SBATCH option (removed -- see
#'   config/rslurm_templates/submit_sh.txt); instead it's baked into the
#'   `apptainer run ...` command this function constructs as the actual
#'   `rscript_path` (see `build_apptainer_rscript_path()` below).
#' @param slurm_options_extra named list merged into slurm_options, for any
#'   job-family-specific overrides (e.g. a container override)
#' @param submit passed through to rslurm; defaults to FALSE (build the
#'   sbatch materials without actually submitting)
#' @param output_dir directory rslurm's `_rslurm_<jobname>` bundle is
#'   created in -- defaults to the current working directory (rslurm's own
#'   default). Created if it doesn't exist. Restored to the prior working
#'   directory afterward regardless of success/failure.
submit_job_family <- function(f, jobs_df, jobname, global_objects = character(0),
                               pkgs = character(0), cluster_cfg,
                               slurm_options_extra = list(), submit = FALSE,
                               output_dir = getwd()) {
  # No `container` here -- it's baked into rscript_path (below) as part of
  # an explicit `apptainer run ...` invocation instead of a native
  # `--container=` SBATCH option.
  slurm_options <- utils::modifyList(
    list(
      mem             = cluster_cfg$mem,
      "cpus-per-task" = cluster_cfg$cpus_per_task,
      time            = cluster_cfg$time
    ),
    slurm_options_extra
  )

  is_single_job <- is.null(jobs_df) || nrow(jobs_df) == 0
  sh_template <- if (is_single_job) {
    cluster_cfg$sh_template_single %||% cluster_cfg$sh_template
  } else {
    cluster_cfg$sh_template
  }
  # Resolve to absolute paths *before* the setwd() below -- relative
  # template paths (e.g. "config/rslurm_templates/submit_sh.txt") are
  # relative to the project root this script was run from, not output_dir.
  sh_template <- normalizePath(sh_template, mustWork = TRUE)

  # Default r_template points at this project's local copies of rslurm's
  # own templates, whose only change is writing results to
  # Sys.getenv("RSLURM_OUTPUT_DIR") instead of a relative path -- see
  # config/rslurm_templates/slurm_run_R.txt / slurm_run_single_R.txt and
  # R/README.md's squashfs packaging section for why. Falls back to
  # rslurm's own bundled templates untouched if this project's copies
  # aren't present for some reason.
  r_template <- if (is_single_job) {
    here::here("config/rslurm_templates/slurm_run_single_R.txt")
  } else {
    here::here("config/rslurm_templates/slurm_run_R.txt")
  }
  if (file.exists(r_template)) r_template <- normalizePath(r_template) else r_template <- NULL

  common_args <- list(
    f              = f,
    jobname        = jobname,
    global_objects = global_objects,
    pkgs           = pkgs,
    libPaths       = cluster_cfg$libPaths,
    rscript_path   = build_apptainer_rscript_path(cluster_cfg$container, jobname, cluster_cfg$rscript_path,
                                                   lib_paths = cluster_cfg$libPaths),
    sh_template    = sh_template,
    r_template     = r_template,
    slurm_options  = slurm_options,
    submit         = submit
  )

  # rslurm has no output-directory argument of its own -- it always writes
  # `_rslurm_<jobname>` under the current working directory -- so this
  # temporarily switches into `output_dir` for the duration of the call.
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  old_wd <- setwd(output_dir)
  on.exit(setwd(old_wd), add = TRUE)

  if (is_single_job) {
    do.call(slurm_call, c(common_args, list(params = list())))
  } else {
    do.call(slurm_apply, c(common_args, list(
      params        = jobs_df,
      nodes         = nrow(jobs_df),
      cpus_per_node = 1
    )))
  }
}

#' Save a named list of sjob objects for this dataset/method, so a later
#' step can find each job's output directory without re-running setup.
save_sjobs <- function(sjobs, dataset_id, method, results_root = here::here("results/rslurm")) {
  dir.create(file.path(results_root, dataset_id), recursive = TRUE, showWarnings = FALSE)
  out_path <- file.path(results_root, dataset_id, paste0(method, "_sjobs.rds"))
  saveRDS(sjobs, out_path)
  message("Saved sjob metadata to ", out_path)
  invisible(out_path)
}

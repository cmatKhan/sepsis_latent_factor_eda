# Auto-discovers method registries from R/methods/*.R, so adding a new
# method requires nothing more than adding a file there -- see R/README.md's
# "Adding a new method" section.
#
# Convention: R/methods/<name>.R must define exactly one object named
# `<name>_registry`, a flat list with:
#   - `global_object`: the matrix/tensor variable name the job function
#     expects (e.g. "mat", "mat_nn", "tnsr").
#   - `jobname`: passed to submit_job_family()/rslurm (e.g. "pca_grid").
#   - `fn`: the job function (as passed to slurm_apply/slurm_call).
#   - `build_grid`: function(params) -> data.frame ready for
#     submit_job_family()'s `jobs_df` -- the method script's own sweep/
#     cross logic, however it needs to build one (typically just
#     `expand.grid(...)` over the specific arguments the job function
#     needs -- see R/methods/pca.R for the simple case, R/methods/cogaps.R
#     for one with a nested sub-block). `params` is `defaults` merged with
#     the dataset config's overrides (see R/create_slurm_bundle.R).
#   - `defaults` (optional, default `list()`): named list of default
#     parameter values, overridable per-key by the dataset config.
#   - `pkgs` (optional): packages the job function needs.
#   - `resource_defaults` (optional): function(params, slurm_cfg) -> params,
#     for defaults derived from cluster config (e.g. CoGAPS's nSets/
#     nThreads).
#
# Two further fields are optional and default to FALSE when absent:
#   - `network`: TRUE if this method's config/slurm entries nest under
#     methods$network$<name>/slurm$network$<name> instead of the flat
#     methods$<name>/slurm$<name> (WGCNA today).
#   - `requires_subject_timepoint`: TRUE if this method needs
#     dataset.subject_id_col/timepoint_col/sample_metadata_path (CP/Tucker
#     today) -- checked generically in R/lib/metadata.R instead of a
#     hardcoded method-name list.
#
# There is no per-method "family" concept anymore -- each method describes
# exactly one kind of run (a stability-design second family was explored
# and deliberately removed; see git history if reviving that).
#
# Each file is source()'d into globalenv (same as the explicit source()
# calls this replaces) -- rslurm's slurm_apply()/slurm_call() only ship a
# job function's *own* global objects across (see R/methods/ica.R's header
# for why every job function is self-contained), but the job function
# itself, and the registry list referencing it, must still be findable in
# globalenv by name.

#' Discover every method registry under `dir`, sourcing each file into
#' globalenv and validating the resulting `<name>_registry` object.
#'
#' @param dir directory to scan for `*.R` files (default R/methods/)
#' @return a named list, keyed by method name, of
#'   `list(registry = <name>_registry, network = TRUE/FALSE)` -- exactly the
#'   shape R/create_slurm_bundle.R's per-method loop expects.
discover_method_registries <- function(dir = here::here("R/methods")) {
  files <- sort(list.files(dir, pattern = "\\.R$", full.names = TRUE))
  if (length(files) == 0) stop("No method files found under ", dir)

  manifest <- list()
  for (f in files) {
    key <- tools::file_path_sans_ext(basename(f))
    registry_var <- paste0(key, "_registry")

    # Clear any stale object of this name before sourcing -- otherwise, in
    # an interactive session that already ran discover_method_registries()
    # once, a leftover `<name>_registry` from a PRIOR successful source()
    # would silently mask this file forgetting to (re)define it.
    if (exists(registry_var, envir = globalenv(), inherits = FALSE)) {
      rm(list = registry_var, envir = globalenv())
    }
    source(f, local = FALSE)

    if (!exists(registry_var, envir = globalenv(), inherits = FALSE)) {
      stop(
        "R/methods/", basename(f), " must define a `", registry_var, "` list ",
        "(a file named `<name>.R` must define `<name>_registry`) -- see ",
        "R/README.md's \"Adding a new method\" section."
      )
    }

    registry <- get(registry_var, envir = globalenv())
    validate_method_registry(registry, key)
    manifest[[key]] <- list(registry = registry, network = isTRUE(registry$network))
  }
  manifest
}

#' Fail fast (before any slurm job is built) on a malformed registry, rather
#' than surfacing a cryptic error deep inside R/create_slurm_bundle.R's loop
#' -- same eager-validation philosophy as R/lib/metadata.R's dataset config
#' checks.
validate_method_registry <- function(registry, key) {
  required <- c("global_object", "jobname", "fn", "build_grid")
  missing <- setdiff(required, names(registry))
  if (length(missing) > 0) {
    stop("`", key, "_registry` is missing required field(s): ", paste(missing, collapse = ", "))
  }
  if (!is.character(registry$global_object) || length(registry$global_object) != 1) {
    stop("`", key, "_registry$global_object` must be a single string")
  }
  if (!is.function(registry$fn) || !is.function(registry$build_grid)) {
    stop("`", key, "_registry`'s `fn` and `build_grid` must both be functions")
  }
  if (!is.null(registry$defaults) &&
      (!is.list(registry$defaults) || (length(registry$defaults) > 0 && is.null(names(registry$defaults))))) {
    stop("`", key, "_registry$defaults` must be a named list")
  }
  invisible(registry)
}

# Method discovery. Every R/methods/<name>.R defines:
#   - `<name>_registry`, the fitting contract (fit function, parameter grid,
#     defaults) -- the fit targets depend on it, so editing it refits;
#   - `<name>_ingest`, the ingestion contract (method properties and
#     extract()) -- only ingestion depends on it.
# discover_method_registries() sources every file and validates both. Field
# lists: docs, Methods, "Adding a method".

#' Run one fit, recording failures and wall time
#'
#' Wraps a fit function's body (every method except sPCA) so a failing fit
#' returns its error message instead of stopping the rest of its batch;
#' ingestion records it as failed. A fit-target dependency: editing this
#' function refits every method that uses it.
#'
#' @param expr The fit's body, evaluated lazily.
#' @param base Fields known before fitting (e.g. rank, seed), returned with
#'   `error` on failure.
#' @return The fit's result list plus `elapsed` (seconds), or `base` plus
#'   `error` and `elapsed` when the fit errored.
capture_fit <- function(expr, base = list()) {
  t0 <- proc.time()[["elapsed"]]
  out <- tryCatch(expr, error = function(e) c(base, list(error = conditionMessage(e))))
  out$elapsed <- proc.time()[["elapsed"]] - t0
  out
}

#' Discover the method files
#'
#' Sources every `R/methods/*.R` into the global environment and validates the
#' `<name>_registry` and `<name>_ingest` each must define (the file name and
#' the variable names must match).
#'
#' @param dir Directory of method files.
#' @return Named list (by method) of `list(registry, ingest, network)`, the
#'   shape the target factories (R/targets/factories.R) take.
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
    for (v in paste0(key, c("_registry", "_ingest"))) {
      if (exists(v, envir = globalenv(), inherits = FALSE)) rm(list = v, envir = globalenv())
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
    ingest_var <- paste0(key, "_ingest")
    if (!exists(ingest_var, envir = globalenv(), inherits = FALSE)) {
      stop("R/methods/", basename(f), " must define `", ingest_var, "` (the ingest contract) -- ",
           "see R/README.md's \"Adding a new method\" section.")
    }
    spec <- get(ingest_var, envir = globalenv())
    validate_method_ingest(spec, key)
    manifest[[key]] <- list(registry = registry, ingest = spec, network = isTRUE(registry$network))
  }
  manifest
}

#' Validate a fitting contract
#'
#' Fails when _targets.R is read rather than mid-run: required fields
#' (`global_object`, `jobname`, `fn`, `build_grid`) and their types.
#'
#' @param registry A `<name>_registry` list.
#' @param key The method name.
#' @return `registry`, invisibly.
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

#' Validate an ingestion contract
#'
#' Required fields: `family` (`"seed_sweep"` or `"param_grid"`),
#' `sign_ambiguous`, `has_loadings`, `resource_params` and `extract`
#' (docs: Methods, "Adding a method").
#'
#' @param spec A `<name>_ingest` list.
#' @param key The method name.
#' @return `spec`, invisibly.
validate_method_ingest <- function(spec, key) {
  required <- c("family", "sign_ambiguous", "has_loadings", "resource_params", "extract")
  missing <- setdiff(required, names(spec))
  if (length(missing) > 0) {
    stop("`", key, "_ingest` is missing required field(s): ", paste(missing, collapse = ", "))
  }
  if (!spec$family %in% c("seed_sweep", "param_grid")) {
    stop("`", key, "_ingest$family` must be \"seed_sweep\" or \"param_grid\"")
  }
  if (!is.function(spec$extract)) stop("`", key, "_ingest$extract` must be a function")
  invisible(spec)
}

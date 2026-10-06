# Glue between the method registries (R/methods/*.R) and the fit targets:
# building each method's grid, running a batch of fits, and reading fits
# back. The fit targets depend on these functions -- editing their code
# refits every method (comments don't count).

#' A method's resolved parameters for one dataset
#'
#' The registry's `defaults` overlaid with the dataset config's method block
#' (merge_named_list()), then `resource_defaults()` (e.g. CoGAPS `nThreads` =
#' the controller's `cpus_per_task`).
#'
#' @param registry A `<name>_registry`.
#' @param meta The dataset config.
#' @param method Method name.
#' @param cpus_per_task The method controller's CPUs per worker.
#' @return Named list of parameters for `build_grid()`.
resolve_method_params <- function(registry, meta, method, cpus_per_task) {
  method_meta <- if (isTRUE(registry$network)) meta$methods$network[[method]] else meta$methods[[method]]
  resolved <- merge_named_list(registry$defaults, method_meta)
  if (!is.null(registry$resource_defaults)) {
    resolved <- registry$resource_defaults(resolved, list(cpus_per_task = cpus_per_task))
  }
  resolved
}

#' Group grid rows into batches
#'
#' Marks `batch_size` consecutive rows with one `tar_group` so the fit target
#' branches over batches (`iteration = "group"`): each branch costs the main
#' process about a second, more than a PCA fit takes.
#'
#' @param grid Grid data frame, one row per fit.
#' @param batch_size Rows per batch.
#' @return `grid` with a `tar_group` column.
add_batches <- function(grid, batch_size) {
  grid$tar_group <- as.integer(ceiling(seq_len(nrow(grid)) / batch_size))
  grid
}

#' A method's parameter grid for one dataset (the `grid_<m>_<dataset>` target)
#'
#' @param registry A `<name>_registry`.
#' @param meta The dataset config.
#' @param method Method name.
#' @param cpus_per_task The method controller's CPUs per worker.
#' @param batch_size Rows per batch.
#' @return Grid data frame with a `tar_group` column, one row per fit.
build_method_grid <- function(registry, meta, method, cpus_per_task, batch_size = 1L) {
  resolved <- resolve_method_params(registry, meta, method, cpus_per_task)
  add_batches(registry$build_grid(resolved), batch_size)
}

#' Second-stage grid from the first stage's fits
#'
#' For a registry with `refine_grid` (sPCA: points around each K's
#' Index-of-Sparseness peak).
#'
#' @param registry A `<name>_registry`.
#' @param meta The dataset config.
#' @param method Method name.
#' @param cpus_per_task The method controller's CPUs per worker.
#' @param batch_size Rows per batch.
#' @param fits The first stage's fit target value (a list of batches).
#' @param pca_fits The dataset's PCA fit target value.
#' @return Grid data frame with a `tar_group` column.
build_refine_grid <- function(registry, meta, method, cpus_per_task, batch_size, fits, pca_fits) {
  resolved <- resolve_method_params(registry, meta, method, cpus_per_task)
  grid <- registry$refine_grid(flatten_fits(fits), flatten_fits(pca_fits)[[1]]$result, resolved)
  if (is.null(grid) || nrow(grid) == 0) stop(method, ": refine_grid() returned no rows")
  add_batches(grid, batch_size)
}

#' Selection table over both stages' fits
#'
#' @param registry A `<name>_registry` with `select_fits`.
#' @param fits,refine_fits The two stages' fit target values.
#' @param pca_fits The dataset's PCA fit target value.
#' @return The registry's selection table (sPCA: spca_select_fits()).
select_method_fits <- function(registry, fits, refine_fits, pca_fits) {
  registry$select_fits(c(flatten_fits(fits), flatten_fits(refine_fits)),
                       flatten_fits(pca_fits)[[1]]$result)
}

#' Flatten a fit target's batches
#'
#' @param x A fit target value: a list of batches, each a list of
#'   `list(params, result)`.
#' @return One `list(params, result)` per grid row.
flatten_fits <- function(x) unlist(x, recursive = FALSE)

#' Columns build_method_grid() adds for batching -- not job-function args.
GRID_BATCH_COLUMNS <- "tar_group"

#' Fit every grid row in one batch (one dynamic branch)
#'
#' Fit functions read their matrix from a free variable named
#' `registry$global_object` (`mat` or `mat_nn`); it's bound in a new enclosing
#' environment for the call rather than assigned globally.
#'
#' @param registry A `<name>_registry`.
#' @param batch The batch's grid rows.
#' @param mat The input matrix.
#' @return One `list(params, result)` per row.
run_method_fits <- function(registry, batch, mat) {
  fn <- registry$fn
  environment(fn) <- list2env(
    stats::setNames(list(mat), registry$global_object),
    parent = environment(fn)
  )
  params <- batch[, setdiff(names(batch), GRID_BATCH_COLUMNS), drop = FALSE]
  lapply(seq_len(nrow(params)), function(i) {
    row <- params[i, , drop = FALSE]
    # list-columns (CoGAPS's params/distributed_params/run) arrive as a
    # length-1 list per row; unwrap them.
    args <- lapply(as.list(row), function(x) if (is.list(x)) x[[1]] else x)
    list(params = row, result = do.call(fn, args))
  })
}

#' A method's fits for one dataset, from the targets store
#'
#' @param dataset_id Dataset id.
#' @param method Method name.
#' @param store targets store directory.
#' @return One `list(params, result)` per grid row, in grid order, including a
#'   refine stage's fits.
read_method_fits <- function(dataset_id, method, store = targets::tar_config_get("store")) {
  stage_names(dataset_id, method, "fit", store) |>
    lapply(targets::tar_read_raw, store = store) |>
    lapply(flatten_fits) |>
    unlist(recursive = FALSE)
}

#' A method's grid for one dataset, from the targets store
#'
#' @param dataset_id Dataset id.
#' @param method Method name.
#' @param store targets store directory.
#' @return Grid data frame without the batching column, both stages,
#'   matching read_method_fits()'s order.
read_method_grid <- function(dataset_id, method, store = targets::tar_config_get("store")) {
  grids <- lapply(stage_names(dataset_id, method, "grid", store), targets::tar_read_raw, store = store)
  grid <- do.call(rbind, grids)
  grid[, setdiff(names(grid), GRID_BATCH_COLUMNS), drop = FALSE]
}

#' Built target names for a method's stages
#'
#' @param dataset_id Dataset id.
#' @param method Method name.
#' @param prefix `"fit"` or `"grid"`.
#' @param store targets store directory.
#' @return The built targets among `<prefix>_<method>_<key>` and
#'   `<prefix>_<method>_refine_<key>`.
stage_names <- function(dataset_id, method, prefix, store) {
  key <- target_key(dataset_id)
  names <- c(paste(prefix, method, key, sep = "_"), paste(prefix, method, "refine", key, sep = "_"))
  built <- targets::tar_meta(store = store, fields = "name")$name
  intersect(names, built)
}

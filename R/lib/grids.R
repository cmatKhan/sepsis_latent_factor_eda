# Generic (method-agnostic) job-grid builders. Each returns a data.frame
# suitable for slurm_apply()'s `params` argument -- one row per job. There's
# no "basis" dimension here: one config = one matrix = one grid.

#' param x seed grid -- the "seed-sweep" stability design (does the
#' discovered factorization change across random seeds, at a given
#' parameter value?).
build_grid_seed_sweep <- function(param_values, seeds, param_name = "rank") {
  grid <- expand.grid(param = param_values, seed = seeds, stringsAsFactors = FALSE)
  names(grid)[names(grid) == "param"] <- param_name
  grid
}

#' param [x extra_grid dims] grid -- the "masking-CV" stability design
#' (which parameter value best reconstructs held-out entries?). `extra` is
#' an optional named list of additional vectors to cross in (e.g. CoGAPS's
#' alpha_range).
build_grid_masking_cv <- function(param_values, param_name = "rank", extra = NULL) {
  dims <- c(list(param = param_values), extra)
  grid <- do.call(expand.grid, c(dims, stringsAsFactors = FALSE))
  names(grid)[names(grid) == "param"] <- param_name
  grid
}

#' Arbitrary parameter grid -- used by network methods (WGCNA's power_grid,
#' wTO's n/delta grid), which don't have a masking-CV analogue. `param_grid`
#' is a named list of vectors to cross.
build_grid_param <- function(param_grid) {
  do.call(expand.grid, c(param_grid, stringsAsFactors = FALSE))
}

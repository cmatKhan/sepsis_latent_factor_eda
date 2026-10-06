# Merging a method's config overrides onto its registry `defaults`
# (docs: Methods, "Adding a method").

#' Null default
#'
#' `a` unless it is `NULL`, then `b` (base R's `%||%`, defined here for older R).
#'
#' @param a A value, possibly `NULL`.
#' @param b The fallback.
#' @return `a`, or `b` when `a` is `NULL`.
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Merge config overrides onto a method's defaults
#'
#' Layers `override`'s entries onto `base` by name (plain `c()` would keep both
#' copies of a duplicated name). Recurses one level into entries that are named
#' lists in both, such as CoGAPS's `params` sub-block, so overriding one nested
#' key (`params.nPatterns`) keeps the other nested defaults.
#'
#' @param base Named list of defaults (a registry's `defaults`), or `NULL`.
#' @param override Named list of overrides (the dataset config's method block), or `NULL`.
#' @return The merged named list.
merge_named_list <- function(base, override) {
  base <- base %||% list()
  override <- override %||% list()
  for (nm in names(override)) {
    base_val <- base[[nm]]
    override_val <- override[[nm]]
    if (is.list(base_val) && !is.null(names(base_val)) &&
        is.list(override_val) && !is.null(names(override_val))) {
      base[[nm]] <- merge_named_list(base_val, override_val)
    } else {
      base[[nm]] <- override_val
    }
  }
  base
}

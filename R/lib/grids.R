# Small shared helper used by R/create_slurm_bundle.R to resolve a
# method's config overrides against its script-defined `defaults`. See
# R/README.md's "Adding a new method" for the registry schema
# (`defaults` + `build_grid`).

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Layer `override`'s keys on top of `base` by NAME (unlike plain `c()`,
#' which just concatenates and would leave duplicate-named entries with
#' the base's value found first) -- used to merge a method's `defaults`
#' with the dataset config's overrides.
#'
#' Recurses one level into any key whose value is itself a NAMED list in
#' both `base` and `override` (e.g. CoGAPS's `params` sub-block) --
#' otherwise, overriding just one nested key (say `params.nPatterns`)
#' would wipe out every other default nested under that same key (e.g.
#' `params.seed`/`params.nIterations`) instead of layering on top of them.
#' Every other method here has no nested list values in `defaults` at
#' all, so this recursion is a no-op for them -- plain top-level
#' replacement, same as before.
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

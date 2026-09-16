# Shared metadata/config loading + validation.
#
# `dataset_metadata.yml` describes ONE dataset x ONE matrix, resolved one of
# three ways (see R/lib/matrices.R::load_input_matrix() for the precedence):
#   1. `dataset.matrix_path` -- an already-built matrix, used as-is.
#   2. `dataset.preprocessing_script` -- a script YOU write that reads the
#      raw parquet trio, does arbitrary manipulation, and writes the result
#      out (see R/README.md's "Preprocessing script contract"). Any
#      sample/feature filtering belongs here, not as a post-hoc id list --
#      if you need more than a bare long-to-wide pivot, write a script.
#   3. neither given -- the framework reads `dataset.expression_path`
#      (long format) and pivots it wide with no other transformation.
#
# There is no "bases" concept here: this framework always runs on exactly
# one matrix per config. If you want to run, say, Day 0 and Day 2 samples
# separately, build/point at two matrices and write two dataset_metadata.yml
# files (one per run) -- see R/README.md.
#
# `methods:` and `slurm:` (SLURM resource settings, formerly a separate
# cluster_config.yml -- merged in so different datasets can use different
# resources) both live in the SAME file. A method is enabled purely by
# being PRESENT under `methods:` -- there is no top-level `enabled:` flag.
#
# Each method's block IS its overrides, directly -- e.g. `methods.pca.rank`
# -- layered on top of that method's own script-defined `defaults` (see
# R/methods/<name>.R). There's no wrapper key: a method now describes
# exactly one kind of run (a stability-design second family, and the
# `full`/`enabled`/`params` wrapper that came with it, were explored and
# deliberately removed; see git history if reviving either). A key only
# takes effect if that method's `build_grid()` actually references it --
# see each R/methods/<name>.R file for exactly which keys it wires
# through, and its own `defaults` for what happens if you omit one.
#
# Because the schema is now effectively opaque per-method (this file has
# no way to know which keys a given method's build_grid() cares about),
# validation here is limited to the boolean-token-key footgun below, plus
# checking that a configured method actually exists -- both eager
# (before any slurm_apply()/slurm_call() is constructed) so a typo'd path
# fails fast on the login node instead of surfacing as a cryptic error
# deep in a batch array.

library(yaml)

if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

#' Read and validate a dataset_metadata.yml file (methods + slurm sections
#' included).
#'
#' @param manifest optional method manifest, as returned by
#'   `discover_method_registries()` (see R/lib/method_registry.R) --
#'   used both to check that every configured method was actually
#'   discovered, and to generically check per-method requirements (e.g.
#'   cp/tucker's subject/timepoint columns) via each registry's
#'   `requires_subject_timepoint` flag, instead of hardcoding method names
#'   here. Defaults to NULL for callers (e.g. R/legacy/*.R) that never
#'   source R/methods/*.R at all -- both checks are skipped entirely when
#'   no manifest is supplied (see validate_dataset_metadata() below); the
#'   cp/tucker check falls back to a hardcoded name check instead so those
#'   (already-stale) legacy callers keep their original guardrail.
read_dataset_metadata <- function(path, manifest = NULL) {
  meta <- yaml::read_yaml(path)
  validate_dataset_metadata(meta, manifest)
  meta
}

validate_dataset_metadata <- function(meta, manifest = NULL) {
  stopifnot(
    "dataset_metadata.yml must have a top-level `dataset` block" = !is.null(meta$dataset),
    "dataset block must have an `id`" = !is.null(meta$dataset$id)
  )

  matrix_path <- meta$dataset$matrix_path
  script_path <- meta$dataset$preprocessing_script

  if (!is.null(matrix_path)) {
    if (!file.exists(matrix_path)) {
      stop("dataset$matrix_path does not exist: ", matrix_path)
    }
  } else if (!is.null(script_path)) {
    if (!file.exists(script_path)) {
      stop("dataset$preprocessing_script does not exist: ", script_path)
    }
  } else {
    if (is.null(meta$dataset$expression_path) || !file.exists(meta$dataset$expression_path)) {
      stop("None of dataset$matrix_path, dataset$preprocessing_script are set, so ",
           "dataset$expression_path must point to an existing long-format ",
           "expression parquet file to pivot by default. See R/README.md's ",
           "\"Preparing your input matrix\".")
    }
  }

  if (is.null(meta$methods) || length(meta$methods) == 0) {
    stop("dataset_metadata.yml must have at least one method under `methods:` ",
         "(presence = enabled -- there is no `enabled:` flag)")
  }

  # Catches a real R `yaml` package footgun (confirmed, 2.3.12): an
  # UNQUOTED key that happens to be a YAML 1.1 legacy boolean token (y/Y/
  # yes/Yes/YES/n/N/no/No/NO/true/True/TRUE/false/False/FALSE/on/On/ON/
  # off/Off/OFF) parses as a literal TRUE/FALSE, not the string you wrote --
  # e.g. a real argument name `n` would silently become the key `FALSE`,
  # which downstream becomes the column name "FALSE." (R's make.names()
  # escaping the reserved word) and breaks the tool call with "unused
  # argument". Quote any such key in YAML (`"n": [...]`) -- see
  # R/README.md's config reference. Applied recursively (e.g. CoGAPS's
  # nested `params`/`distributed_params`/`run`), since there's no longer a
  # fixed wrapper depth to stop at.
  check_bool_key_typo_recursive <- function(block, label) {
    for (nm in names(block)) {
      if (nm %in% c("TRUE", "FALSE")) {
        stop(label, " has a key literally named `", nm, "` -- almost certainly an unquoted ",
             "YAML boolean-token key (y/n/yes/no/on/off/true/false, in any case) that got ",
             "coerced instead of staying a literal string. Quote it in the YAML, e.g. `\"n\": ",
             "[...]` instead of `n: [...]`. See R/README.md's config reference.")
      }
      val <- block[[nm]]
      if (is.list(val) && !is.null(names(val))) check_bool_key_typo_recursive(val, paste0(label, ".", nm))
    }
  }

  # Only checks "does this method exist" when a manifest is supplied (i.e.
  # from R/create_slurm_bundle.R, which sources R/methods/*.R before
  # calling read_dataset_metadata()) -- meaningless without it. Legacy
  # callers (manifest = NULL) skip that half; the boolean-token-key check
  # always runs regardless.
  for (nm in names(meta$methods)) {
    if (nm == "network") {
      for (backend in names(meta$methods$network)) {
        if (!is.null(manifest) && is.null(manifest[[backend]]$registry)) {
          stop("methods.network.", backend, " is configured but no such method was discovered under R/methods/")
        }
        check_bool_key_typo_recursive(meta$methods$network[[backend]], paste0("methods.network.", backend))
      }
    } else {
      if (!is.null(manifest) && is.null(manifest[[nm]]$registry)) {
        stop("methods.", nm, " is configured but no such method was discovered under R/methods/")
      }
      check_bool_key_typo_recursive(meta$methods[[nm]], paste0("methods.", nm))
    }
  }

  # Tensor methods (cp/tucker today) need dataset:-level subject/timepoint
  # columns to reshape the standard matrix into a genes x subjects x
  # timepoints array -- see R/lib/tensors.R::build_tensor(). Presence-only
  # check here (fails fast); the columns' actual existence in
  # sample_metadata_path is checked at build_tensor() runtime.
  #
  # Driven by each method's registry (`requires_subject_timepoint`, see
  # R/lib/method_registry.R) rather than a hardcoded cp/tucker name list --
  # this is the ONLY reason `manifest` is threaded through to this
  # function, and only fires when a manifest is actually supplied (i.e.
  # from R/create_slurm_bundle.R, which sources R/methods/*.R before
  # calling read_dataset_metadata()). Callers that never discover method
  # registries at all (R/legacy/*.R) fall back to the original cp/tucker-
  # only check so they keep the same guardrail they've always had.
  if (!is.null(manifest)) {
    check_subject_timepoint <- function(nm) {
      reg <- manifest[[nm]]$registry
      if (is.null(reg) || !isTRUE(reg$requires_subject_timepoint)) return(invisible())
      if (is.null(meta$dataset$subject_id_col) || is.null(meta$dataset$timepoint_col)) {
        stop("methods.", nm, " is configured but dataset.subject_id_col and ",
             "dataset.timepoint_col are not both set -- required to build the genes x ",
             "subjects x timepoints tensor. See R/README.md's config reference.")
      }
      if (is.null(meta$dataset$sample_metadata_path)) {
        stop("methods.", nm, " is configured but dataset.sample_metadata_path ",
             "is not set -- required to look up each sample's subject/timepoint")
      }
    }
    for (nm in setdiff(names(meta$methods), "network")) check_subject_timepoint(nm)
  } else if (!is.null(meta$methods$cp) || !is.null(meta$methods$tucker)) {
    if (is.null(meta$dataset$subject_id_col) || is.null(meta$dataset$timepoint_col)) {
      stop("methods.cp / methods.tucker are configured but dataset.subject_id_col and ",
           "dataset.timepoint_col are not both set -- required to build the genes x ",
           "subjects x timepoints tensor. See R/README.md's config reference.")
    }
    if (is.null(meta$dataset$sample_metadata_path)) {
      stop("methods.cp / methods.tucker are configured but dataset.sample_metadata_path ",
           "is not set -- required to look up each sample's subject/timepoint")
    }
  }

  validate_slurm_config(meta$slurm, names(meta$methods))
  invisible(meta)
}

#' Every method actually present in `methods:` must have a matching
#' `slurm:` entry with the fields submit_job_family() needs.
#'
#' Also guards against a real R `yaml` package (2.3.12, confirmed) footgun:
#' `<<: *anchor` merge keys DO NOT let sibling keys in the same mapping
#' override the merged-in values -- every key after `<<:` is silently
#' dropped instead. A config that used anchors this way had EVERY method's
#' slurm entry silently collapse to whichever one anchor originally
#' defined, container included. It's normal/expected for SOME methods to
#' intentionally share a container (e.g. pca/nmf both just need
#' tidyverse) -- but if there are 3+ methods configured and EVERY single
#' one ends up with the identical container, that's the exact symptom of
#' the merge-key bug, not a real intentional setup, so it's flagged.
validate_slurm_config <- function(slurm, method_names) {
  if (is.null(slurm)) stop("dataset_metadata.yml must have a top-level `slurm` block ",
                            "(SLURM resource settings per method -- formerly cluster_config.yml)")
  required <- c("mem", "cpus_per_task", "time", "container", "libPaths", "sh_template", "rscript_path")
  check_entry <- function(entry, label) {
    missing <- setdiff(required, names(entry))
    if (length(missing) > 0) {
      stop("slurm entry '", label, "' is missing required field(s): ", paste(missing, collapse = ", "))
    }
  }
  containers <- list()
  for (nm in method_names) {
    if (nm == "network") {
      for (backend in names(slurm$network)) {
        check_entry(slurm$network[[backend]], paste0("network.", backend))
        containers[[paste0("network.", backend)]] <- slurm$network[[backend]]$container
      }
    } else {
      if (is.null(slurm[[nm]])) stop("slurm.", nm, " is missing (methods.", nm, " is configured)")
      check_entry(slurm[[nm]], nm)
      containers[[nm]] <- slurm[[nm]]$container
    }
  }
  container_vals <- unlist(containers)
  if (length(container_vals) >= 3 && length(unique(container_vals)) == 1) {
    warning("Every configured method's slurm entry ([", paste(names(containers), collapse = ", "),
            "]) resolves to the SAME container ('", container_vals[[1]], "') -- this is the exact ",
            "symptom of R yaml's `<<: *anchor` merge-key bug, where sibling keys after `<<:` are ",
            "silently dropped (see this function's comment / R/README.md's config reference). ",
            "Write each slurm entry out in full rather than using anchors, unless every method ",
            "genuinely does share one container.", call. = FALSE)
  }
  invisible(slurm)
}

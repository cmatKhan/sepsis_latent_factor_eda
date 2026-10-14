# Study design for variancePartition and dream (R/targets/variance.R): the
# optional `design:` block of a dataset config, turned into the model data,
# formulas and contrasts.
#
# One cell-means design for every dataset: `cell` is group x time (time alone
# or group alone when the other is absent), so each contrast of interest is a
# difference of cell means, and ragged group x time layouts (a control arm
# sampled at fewer timepoints) need no special case.

DESIGN_FIELDS <- c("platform", "subject", "time", "time_order", "group", "group_at_first_time",
                   "covariates", "vp_covariates", "numeric", "drop", "keep", "set_na", "recode",
                   "contrasts", "de_covariates", "control_probes")

#' Validate a config's `design:` block
#'
#' Checks field names and types only; that the columns exist is checked
#' against the sample metadata by prepare_design().
#'
#' @param design The parsed `design:` block.
#' @return `design`, invisibly (errors on the first problem).
validate_design <- function(design) {
  bad <- setdiff(names(design), DESIGN_FIELDS)
  if (length(bad) > 0) stop("design: unknown field(s) ", paste(bad, collapse = ", "))
  if (!isTRUE(design$platform %in% c("array", "rnaseq"))) {
    stop("design.platform must be 'array' or 'rnaseq'")
  }
  if (is.null(design[["time"]]) && is.null(design[["group"]])) {
    stop("design needs a `time` or a `group` column (cell = group x time)")
  }
  if (!is.null(design[["time"]]) && is.null(design$time_order)) {
    stop("design.time is set, so design.time_order (chronological levels) is required")
  }
  if (isTRUE(design$group_at_first_time) && (is.null(design[["group"]]) || is.null(design[["subject"]]))) {
    stop("design.group_at_first_time needs design.group and design.subject")
  }
  if (!is.null(design[["subject"]]) && is.null(design[["time"]])) {
    stop("design.subject is set without design.time; repeated measures need a time column")
  }
  for (f in c("subject", "time", "group")) {
    if (!is.null(design[[f]]) && !(is.character(design[[f]]) && length(design[[f]]) == 1)) {
      stop("design.", f, " must be one column name")
    }
  }
  if (length(setdiff(design$vp_covariates, design$numeric)) > 0) {
    stop("design.vp_covariates must also be listed in design.numeric (they are fixed effects)")
  }
  cp <- design$control_probes
  if (!is.null(cp) && !(is.character(cp$id_pattern) || (is.character(cp$column) && length(cp$values) > 0))) {
    stop("design.control_probes needs `id_pattern`, or `column` and `values`")
  }
  if (isTRUE(cp$covariate) && !identical(design$platform, "array")) {
    stop("design.control_probes.covariate applies to arrays only")
  }
  not_cov <- setdiff(design$numeric, c(design$covariates, design$vp_covariates, design$de_covariates))
  if (length(not_cov) > 0) {
    stop("design.numeric lists column(s) that aren't covariates: ", paste(not_cov, collapse = ", "))
  }
  for (f in c("drop", "keep", "set_na", "recode", "contrasts", "control_probes")) {
    if (!is.null(design[[f]]) && (!is.list(design[[f]]) || is.null(names(design[[f]])))) {
      stop("design.", f, " must be a mapping")
    }
  }
  vars <- design_variables(design)
  if (any(make.names(vars) != vars)) {
    stop("design column names must be syntactic R names (they go in model formulas): ",
         paste(vars[make.names(vars) != vars], collapse = ", "))
  }
  invisible(design)
}

#' Columns a design uses in its models
#'
#' @param design A `design:` block.
#' @return Character vector: subject, time, group, covariates, variance-model-only
#'   covariates, DE covariates.
design_variables <- function(design) {
  unique(c(design[["subject"]], design[["time"]], design[["group"]], design$covariates, design$vp_covariates,
           design$de_covariates))
}

#' Missing-value recoding shared by every model column
#'
#' @param x A metadata column.
#' @return `x` as character with `"NA"` and `""` set to `NA`.
as_model_character <- function(x) {
  x <- as.character(x)
  x[x %in% c("NA", "")] <- NA_character_
  x
}

#' Build the model data, formulas and contrasts for one dataset
#'
#' Applies the design's `drop`/`keep` filters and `set_na`/`recode`
#' recodings, optionally fixes the group at each subject's first timepoint
#' (`group_at_first_time`), coerces and standardizes `numeric`
#' covariates, makes every other model column a factor (time in `time_order`),
#' drops samples with a missing model variable (counted per column), and builds
#' the cell factor with make.names()-safe levels.
#'
#' Formulas follow the variancePartition vignette: categorical variables as
#' random effects, numeric ones fixed, constituents included with the
#' group:time interaction; dream gets `~ 0 + cell + (1|subject)` plus only the
#' config's `de_covariates`.
#'
#' @param meta The resolved dataset config (with a `design:` block).
#' @param sample_meta Sample metadata data frame (read_sample_metadata()).
#' @param control_factor compute_control_factor() output, or `NULL`. When
#'   given, its scores become the numeric covariate `control_factor` in both
#'   the variance and the dream formulas, and the share of the control
#'   probes' variance it captures is kept as `control_factor_fraction`.
#' @return List: `spec` (the design block), `platform`, `data` (model data frame,
#'   rownames the sample ids), `vp_formula`, `dream_formula`, `cell_levels`
#'   (readable cell level -> safe level), `contrasts` (build_contrasts()),
#'   `cancor` (cancor_pairs()), `dropped` (data frame of reason and count),
#'   `n_samples_metadata`, `n_group_recoded` (samples whose group
#'   `group_at_first_time` changed), `numeric_centers` and `numeric_scales`
#'   (each numeric covariate's mean and SD, used to standardize it),
#'   `control_factor_fraction`.
prepare_design <- function(meta, sample_meta, control_factor = NULL) {
  spec <- meta$design
  validate_design(spec)
  sid <- meta$dataset$sample_id_col %||% "sample_id"
  df <- as.data.frame(sample_meta)
  if (!is.null(control_factor)) {
    df$control_factor <- unname(control_factor$factor[as.character(df[[sid]])])
    spec$vp_covariates <- c(unlist(spec$vp_covariates), "control_factor")
    spec$de_covariates <- c(unlist(spec$de_covariates), "control_factor")
    spec$numeric <- c(unlist(spec$numeric), "control_factor")
  }
  needed <- c(sid, design_variables(spec), names(spec$drop), names(spec$keep), names(spec$set_na),
              names(spec$recode))
  missing_cols <- setdiff(needed, names(df))
  if (length(missing_cols) > 0) {
    stop(meta$dataset$id, ": design column(s) not in the sample metadata: ",
         paste(missing_cols, collapse = ", "))
  }
  n0 <- nrow(df)
  dropped <- list()
  note <- function(reason, n) if (n > 0) dropped[[length(dropped) + 1]] <<- data.frame(reason = reason, n = n)

  for (col in names(spec$drop)) {
    hit <- as.character(df[[col]]) %in% as.character(unlist(spec$drop[[col]]))
    note(paste0("drop ", col, " in [", paste(unlist(spec$drop[[col]]), collapse = ", "), "]"), sum(hit))
    df <- df[!hit, , drop = FALSE]
  }
  for (col in names(spec$keep)) {
    hit <- as.character(df[[col]]) %in% as.character(unlist(spec$keep[[col]]))
    note(paste0("keep ", col, " in [", paste(unlist(spec$keep[[col]]), collapse = ", "), "]"), sum(!hit))
    df <- df[hit, , drop = FALSE]
  }

  vars <- design_variables(spec)
  numeric_vars <- unlist(spec$numeric)
  for (v in vars) {
    x <- as_model_character(df[[v]])
    if (!is.null(spec$set_na[[v]])) x[x %in% as.character(unlist(spec$set_na[[v]]))] <- NA_character_
    for (from in names(spec$recode[[v]])) x[x %in% from] <- as.character(spec$recode[[v]][[from]])
    if (v %in% numeric_vars) {
      num <- suppressWarnings(as.numeric(x))
      bad <- !is.na(x) & is.na(num)
      if (any(bad)) {
        stop(meta$dataset$id, ": design.numeric column '", v, "' has non-numeric values: ",
             paste(utils::head(unique(x[bad]), 5), collapse = ", "))
      }
      df[[v]] <- num
    } else {
      df[[v]] <- x
    }
  }
  for (v in vars) {
    miss <- is.na(df[[v]])
    note(paste0("missing ", v), sum(miss))
    df <- df[!miss, , drop = FALSE]
  }
  # Standardized: centring keeps a numeric covariate's coefficient from
  # correlating with the intercept (colinearityScore() was ~0.98 for uncentred
  # age, near the vignette's 0.99 limit), and scaling keeps lme4 from refusing
  # predictors on very different scales (ROSE sequencing depth ~1e7 next to
  # percentages). A variance fraction, var(X * beta), is unchanged by either.
  numeric_centers <- vapply(numeric_vars, function(v) mean(df[[v]]), numeric(1))
  numeric_scales <- vapply(numeric_vars, function(v) stats::sd(df[[v]]), numeric(1))
  for (v in numeric_vars) df[[v]] <- (df[[v]] - numeric_centers[[v]]) / numeric_scales[[v]]

  if (!is.null(spec[["time"]])) {
    order <- as.character(unlist(spec$time_order))
    unknown <- setdiff(unique(df[[spec[["time"]]]]), order)
    if (length(unknown) > 0) {
      stop(meta$dataset$id, ": time values not in design.time_order: ", paste(unknown, collapse = ", "))
    }
    df[[spec[["time"]]]] <- droplevels(factor(df[[spec[["time"]]]], levels = order))
  }
  # A group that is a time-varying state (e.g. "SIRS resolved" at day 3) made
  # into a fixed arm: every sample takes its subject's first-timepoint value.
  n_group_recoded <- 0L
  if (isTRUE(spec$group_at_first_time)) {
    first <- tapply(seq_len(nrow(df)), df[[spec[["subject"]]]],
                    function(i) i[which.min(as.integer(df[[spec[["time"]]]][i]))])
    arm <- df[[spec[["group"]]]][unlist(first)[as.character(df[[spec[["subject"]]]])]]
    n_group_recoded <- sum(arm != df[[spec[["group"]]]])
    df[[spec[["group"]]]] <- arm
  }
  for (v in setdiff(vars, c(spec[["time"]], numeric_vars))) df[[v]] <- droplevels(factor(df[[v]]))

  ids <- as.character(df[[sid]])
  if (anyDuplicated(ids)) stop(meta$dataset$id, ": duplicated sample ids in the sample metadata")
  data <- df[, vars, drop = FALSE]
  rownames(data) <- ids

  readable <- if (!is.null(spec[["group"]]) && !is.null(spec[["time"]])) {
    droplevels(interaction(data[[spec[["group"]]]], data[[spec[["time"]]]], sep = ".", lex.order = TRUE))
  } else if (!is.null(spec[["time"]])) data[[spec[["time"]]]] else data[[spec[["group"]]]]
  cell_levels <- stats::setNames(make.names(levels(readable), unique = TRUE), levels(readable))
  data$cell <- factor(unname(cell_levels[as.character(readable)]), levels = unname(cell_levels))

  out <- list(
    spec = spec,
    platform = spec$platform,
    data = data,
    vp_formula = vp_formula(spec),
    dream_formula = dream_formula(spec),
    cell_levels = cell_levels,
    dropped = if (length(dropped) > 0) do.call(rbind, dropped) else data.frame(reason = character(0), n = integer(0)),
    n_samples_metadata = n0,
    n_group_recoded = n_group_recoded,
    numeric_centers = numeric_centers,
    numeric_scales = numeric_scales,
    control_factor_fraction = control_factor$fraction %||% NA_real_
  )
  out$contrasts <- build_contrasts(out)
  out$cancor <- cancor_pairs(spec, data)
  out
}

#' variancePartition formula
#'
#' @param spec A `design:` block.
#' @return `~ (1|subject) + (1|time) + (1|group) + (1|group:time) + (1|cov) + numeric`
#'   (`vp_covariates` among the numeric terms), absent terms left out.
vp_formula <- function(spec) terms_formula(vp_terms(spec))

#' The variance model's terms
#'
#' @param spec A `design:` block.
#' @return List: `random` (variables, `group:time` for the interaction) and
#'   `fixed` (numeric covariates, including `vp_covariates`).
vp_terms <- function(spec) {
  numeric_vars <- unlist(spec$numeric)
  list(random = c(spec[["subject"]], spec[["time"]], spec[["group"]], interaction_term(spec),
                  setdiff(unlist(spec$covariates), numeric_vars)),
       fixed = unique(c(intersect(unlist(spec$covariates), numeric_vars), unlist(spec$vp_covariates))))
}

#' The group x time interaction term's name
#'
#' @param spec A `design:` block.
#' @return `"group:time"`, or `NULL` without both.
interaction_term <- function(spec) {
  if (!is.null(spec[["group"]]) && !is.null(spec[["time"]])) paste0(spec[["group"]], ":", spec[["time"]])
}

#' A formula from vp_terms()-shaped terms
#'
#' @param terms List of `random` and `fixed` term names.
#' @return `~ (1|random) + ... + fixed`, or `NULL` when there are no terms.
terms_formula <- function(terms) {
  rhs <- c(if (length(terms$random)) paste0("(1|", terms$random, ")"), terms$fixed)
  if (!length(rhs)) return(NULL)
  stats::as.formula(paste("~", paste(rhs, collapse = " + ")), env = globalenv())
}

#' A formula as one line of text
#'
#' @param f A formula, or `NULL`.
#' @return Character, `NA` for `NULL`.
formula_text <- function(f) if (is.null(f)) NA_character_ else paste(deparse(f), collapse = "")

#' dream formula
#'
#' @param spec A `design:` block.
#' @return `~ 0 + cell + (1|subject)` (no subject term when cross-sectional),
#'   plus any `de_covariates` as fixed effects.
dream_formula <- function(spec) {
  terms <- c("0", "cell", if (!is.null(spec[["subject"]])) paste0("(1|", spec[["subject"]], ")"),
             unlist(spec$de_covariates))
  stats::as.formula(paste("~", paste(terms, collapse = " + ")), env = globalenv())
}

#' Contrasts for dream
#'
#' Time course: for each group (or all samples when there is no group) with at
#' least `min_pairs` subjects sampled at 2+ of its timepoints, each later
#' timepoint against the group's first in `time_order`, plus an omnibus spec
#' (the group's contrast names, for a joint F-test at ingest). Groups below
#' `min_pairs` stay in the fit without contrasts. Cross-sectional: the config's
#' `contrasts:` (`{name: "levelA - levelB"}`, readable group levels), else
#' every pairwise group comparison (later level minus earlier).
#'
#' @param design_obj prepare_design()'s list (needs `spec`, `data`, `cell_levels`,
#'   `dream_formula`).
#' @param min_pairs Minimum subjects with repeated samples per tested group.
#' @return List: `table` (name, kind, expression, readable), `omnibus` (name ->
#'   contrast names), `L` (makeContrastsDream() matrix, `NULL` when there are
#'   no contrasts), `skipped` (groups left without contrasts and why).
build_contrasts <- function(design_obj, min_pairs = 3) {
  spec <- design_obj$spec
  data <- design_obj$data
  cl <- design_obj$cell_levels
  coef <- function(readable) paste0("cell", cl[[readable]])
  safe <- function(x) make.names(x)
  exprs <- character(0); readable <- character(0)
  omnibus <- list(); skipped <- character(0)

  if (!is.null(spec[["time"]]) && !is.null(spec[["subject"]])) {
    groups <- if (is.null(spec[["group"]])) "all" else levels(data[[spec[["group"]]]])
    for (g in groups) {
      in_g <- if (is.null(spec[["group"]])) rep(TRUE, nrow(data)) else data[[spec[["group"]]]] == g
      time_g <- droplevels(data[[spec[["time"]]]][in_g])
      if (nlevels(time_g) < 2) { skipped[g] <- "one timepoint"; next }
      subj_g <- droplevels(data[[spec[["subject"]]]][in_g])
      n_paired <- sum(tapply(as.character(time_g), subj_g, function(x) length(unique(x))) >= 2)
      if (n_paired < min_pairs) { skipped[g] <- paste0(n_paired, " subjects with >= 2 timepoints"); next }
      t0 <- levels(time_g)[1]
      prefix <- if (is.null(spec[["group"]])) "" else paste0(g, ".")
      names_g <- character(0)
      for (t in levels(time_g)[-1]) {
        nm <- if (is.null(spec[["group"]])) paste0(safe(t), "vs", safe(t0))
              else paste0(safe(g), "_", safe(t), "vs", safe(t0))
        exprs[nm] <- paste(coef(paste0(prefix, t)), "-", coef(paste0(prefix, t0)))
        readable[nm] <- paste0(prefix, t, " - ", prefix, t0)
        names_g <- c(names_g, nm)
      }
      omnibus[[if (is.null(spec[["group"]])) "time_omnibus" else paste0(safe(g), "_time_omnibus")]] <- names_g
    }
  } else if (!is.null(spec[["group"]])) {
    pairs <- if (!is.null(spec$contrasts)) {
      lapply(spec$contrasts, function(e) {
        sides <- trimws(strsplit(e, " - ", fixed = TRUE)[[1]])
        if (length(sides) != 2 || !all(sides %in% names(cl))) {
          stop("design.contrasts '", e, "' must be 'levelA - levelB' with levels from: ",
               paste(names(cl), collapse = ", "))
        }
        sides
      })
    } else {
      combos <- utils::combn(names(cl), 2, simplify = FALSE)
      stats::setNames(lapply(combos, rev),
                      vapply(combos, function(p) paste0(safe(p[2]), "vs", safe(p[1])), character(1)))
    }
    for (nm in names(pairs)) {
      exprs[nm] <- paste(coef(pairs[[nm]][1]), "-", coef(pairs[[nm]][2]))
      readable[nm] <- paste(pairs[[nm]], collapse = " - ")
    }
  }

  table <- data.frame(name = names(exprs), kind = "pairwise", expression = unname(exprs),
                      readable = unname(readable[names(exprs)]))
  L <- if (length(exprs) > 0) {
    variancePartition::makeContrastsDream(design_obj$dream_formula, data, contrasts = exprs)
  }
  list(table = table, omnibus = omnibus, L = L, skipped = skipped)
}

#' Canonical correlations between the design's variables
#'
#' variancePartition::canCorPairs() on the model variables (the interaction
#' left out). Warnings are not stored in the value: subject-level variables
#' nested in subject already show as correlation 1 in the matrix, and a
#' stored once-per-session notice (lme4's "'subbars' has moved to
#' reformulas") made the design value depend on which dataset's design was
#' built first in the R session. That notice is muffled; other warnings
#' reach the target's warnings.
#'
#' @param spec A `design:` block.
#' @param data The model data frame.
#' @return List: `matrix` (or `NULL` on error), `error`.
cancor_pairs <- function(spec, data) {
  vars <- setdiff(design_variables(spec), unlist(spec$de_covariates))
  if (length(vars) < 2) return(list(matrix = NULL, error = NA_character_))
  res <- tryCatch(
    withCallingHandlers(
      variancePartition::canCorPairs(stats::as.formula(paste("~", paste(vars, collapse = " + "))), data),
      warning = function(w) if (grepl("subbars", conditionMessage(w))) invokeRestart("muffleWarning")
    ),
    error = function(e) structure(conditionMessage(e), class = "cancor_error")
  )
  if (inherits(res, "cancor_error")) return(list(matrix = NULL, error = unclass(res)))
  list(matrix = res, error = NA_character_)
}

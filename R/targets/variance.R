# Variance partition and differential expression per dataset
# (variancePartition's fitVarPartModel() and dream()), for every pipeline
# dataset whose config has a `design:` block. Follows the package vignettes:
# variance fractions with ML, dream with REML and Kenward-Roger /
# Satterthwaite degrees of freedom (`ddf = "adaptive"`), categorical variables
# as random effects. The full results are saved as target values, so what the
# DB ingests from them can change without refitting.

#' What is kept from each gene's variancePartition fit
#'
#' fitVarPartModel()'s `fxn`: the per-gene lmer fits are dropped (about 10 GB
#' for 20k genes x 1000 samples, per the vignette). Every call is
#' namespace-qualified because BiocParallel's Snow workers don't attach
#' packages.
#'
#' @param fit One gene's fitted model (merMod, or lm with no random effect).
#' @return List: `variances` (calcVarPart() absolute; the fractions are these
#'   over their sum), `sigma`, `nobs`, `singular`, `colinearity`.
vp_record <- function(fit) {
  list(
    variances = variancePartition::calcVarPart(fit, returnFractions = FALSE),
    sigma = stats::sigma(fit),
    nobs = stats::nobs(fit),
    singular = if (inherits(fit, "merMod")) lme4::isSingular(fit) else FALSE,
    colinearity = tryCatch(as.numeric(variancePartition::colinearityScore(fit)),
                           error = function(e) NA_real_)
  )
}

#' Parallel backend for one task
#'
#' @param cpus Worker processes (the controller's `cpus_per_task`).
#' @return A BiocParallel SnowParam, or SerialParam for one CPU.
variance_bpparam <- function(cpus) {
  if (cpus > 1) BiocParallel::SnowParam(workers = cpus, type = "SOCK", progressbar = FALSE)
  else BiocParallel::SerialParam()
}

#' Expression ready for the models (the `norm_<dataset>` target)
#'
#' RNA-seq: voomWithDreamWeights() with the dream formula (an EList whose
#' precision weights both fitVarPartModel() and dream() use). Arrays: the
#' collapsed log-scale matrix unchanged.
#'
#' @param expr prepare_expression() output.
#' @param design_obj prepare_design() output.
#' @param cpus Worker processes.
#' @return List: `values` (EList or matrix), `seconds`, `resources`
#'   (job_resources()).
normalize_for_models <- function(expr, design_obj, cpus) {
  t0 <- proc.time()[["elapsed"]]
  values <- if (expr$platform == "rnaseq") {
    variancePartition::voomWithDreamWeights(expr$values, design_obj$dream_formula, design_obj$data,
                                            BPPARAM = variance_bpparam(cpus))
  } else {
    expr$values
  }
  list(values = values, seconds = proc.time()[["elapsed"]] - t0, resources = job_resources(cpus))
}

#' Variance partition (the `vp_<dataset>` target)
#'
#' fitVarPartModel() with ML (`REML = FALSE`, the vignette's choice for
#' variance fractions), keeping vp_record() per gene.
#'
#' @param norm normalize_for_models() output.
#' @param design_obj prepare_design() output.
#' @param genes prepare_expression()'s `genes` table.
#' @param cpus Worker processes.
#' @return List: `formula`, `genes`, `varPart` (gene x term fractions),
#'   `variances` (gene x term absolute), `flags` (per gene sigma, nobs,
#'   singular, colinearity), `n_genes_input` (genes modelled, failed ones
#'   included), `errors`, `reml`, `weights_used`, `seconds`, `resources`
#'   (job_resources()), `versions`.
run_variance_partition <- function(norm, design_obj, genes, cpus) {
  t0 <- proc.time()[["elapsed"]]
  records <- variancePartition::fitVarPartModel(
    norm$values, design_obj$vp_formula, design_obj$data, REML = FALSE, fxn = vp_record,
    showWarnings = FALSE, BPPARAM = variance_bpparam(cpus))
  errors <- attr(records, "errors")
  variances <- do.call(rbind, lapply(records, `[[`, "variances"))
  rownames(variances) <- names(records)
  flags <- data.frame(
    gene = names(records),
    sigma = vapply(records, `[[`, numeric(1), "sigma"),
    nobs = vapply(records, function(r) as.numeric(r$nobs), numeric(1)),
    singular = vapply(records, `[[`, logical(1), "singular"),
    colinearity = vapply(records, function(r) as.numeric(r$colinearity)[1], numeric(1)),
    row.names = NULL
  )
  list(
    formula = design_obj$vp_formula,
    genes = genes,
    varPart = as.data.frame(variances / rowSums(variances)),
    variances = as.data.frame(variances),
    flags = flags,
    n_genes_input = nrow(model_expression(norm)),
    errors = errors,
    reml = FALSE,
    weights_used = inherits(norm$values, "EList") && !is.null(norm$values$weights),
    seconds = proc.time()[["elapsed"]] - t0,
    resources = job_resources(cpus),
    versions = variance_versions()
  )
}

#' Differential expression (the `dream_<dataset>` target)
#'
#' dream() on the cell-means formula with the design's contrasts, then
#' eBayes(). The full fit is kept so any coefficient or joint test can be
#' derived later; dream_tables() builds the design's contrasts' tables.
#'
#' @param norm normalize_for_models() output.
#' @param design_obj prepare_design() output.
#' @param cpus Worker processes.
#' @return List: `fit_eb` (eBayes-moderated MArrayLM2; holds every field of
#'   the dream fit), `L`, `contrast_table`, `omnibus`, `formula`, `ddf`,
#'   `errors`, `seconds`, `resources` (job_resources()), `versions`. The
#'   values and data it was fit on are the `norm_` and `design_` targets.
run_dream <- function(norm, design_obj, cpus) {
  t0 <- proc.time()[["elapsed"]]
  L <- design_obj$contrasts$L
  fit <- variancePartition::dream(norm$values, design_obj$dream_formula, design_obj$data, L = L,
                                  ddf = "adaptive", BPPARAM = variance_bpparam(cpus))
  errors <- attr(fit, "errors")
  fit_eb <- variancePartition::eBayes(fit)
  list(
    fit_eb = fit_eb,
    L = L,
    contrast_table = design_obj$contrasts$table,
    omnibus = design_obj$contrasts$omnibus,
    formula = design_obj$dream_formula,
    ddf = "adaptive",
    errors = errors,
    seconds = proc.time()[["elapsed"]] - t0,
    resources = job_resources(cpus),
    versions = variance_versions()
  )
}

#' Every contrast's and joint test's full table (the `de_tables_<dataset>` target)
#'
#' The one place dream's fit is read: topTable() for each contrast and, for
#' each group with more than one contrast, their joint F-test.
#'
#' @param dream run_dream() output.
#' @return List: `tables` (named topTable() data frames, genes in fit order),
#'   `contrasts` (contrast, kind, readable, n_genes, n_fdr05), and the fit's
#'   `seconds`, `errors`, `n_errors`, `resources`.
dream_tables <- function(dream) {
  omni <- Filter(function(x) length(x) > 1, dream$omnibus)
  coefs <- c(stats::setNames(as.list(dream$contrast_table$name), dream$contrast_table$name), omni)
  tables <- lapply(coefs, function(cf) {
    variancePartition::topTable(dream$fit_eb, coef = cf, number = Inf, sort.by = "none")
  })
  # A design can have no contrasts (no group with enough repeated subjects).
  contrasts <- data.frame(
    contrast = names(coefs) %||% character(0),
    kind = rep(c("pairwise", "omnibus"), c(nrow(dream$contrast_table), length(omni))),
    readable = c(dream$contrast_table$readable, vapply(omni, paste, character(1), collapse = ", ")),
    n_genes = rep(nrow(dream$fit_eb$coefficients), length(coefs)),
    n_fdr05 = vapply(tables, function(t) sum(t$adj.P.Val < 0.05, na.rm = TRUE), numeric(1)),
    row.names = NULL)
  list(tables = tables, contrasts = contrasts, seconds = dream$seconds, errors = dream$errors,
       n_errors = length(dream$errors), resources = dream$resources)
}

#' Resources this SLURM job has used so far
#'
#' Read from the job's cgroup (v2), which, unlike `seff`/`sacct`, counts the
#' BiocParallel Snow workers: they detach from the job's process tree, so
#' SLURM's accounting misses them, but they stay in its cgroup. With
#' `tasks_max: 1` a worker job runs one task, so this is that task's use.
#'
#' @param cpus CPUs the task was given (for the efficiency ratio).
#' @return List: `job_id`, `cpus`, `memory_peak_gb`, `cpu_seconds` (`NA`
#'   outside SLURM or without cgroup v2).
job_resources <- function(cpus) {
  dir <- file.path("/sys/fs/cgroup/system.slice/slurmstepd.scope",
                   paste0("job_", Sys.getenv("SLURM_JOB_ID")))
  read_num <- function(f, key = NULL) tryCatch({
    x <- readLines(file.path(dir, f))
    if (!is.null(key)) x <- sub(paste0("^", key, " "), "", grep(paste0("^", key, " "), x, value = TRUE))
    as.numeric(x[1])
  }, error = function(e) NA_real_, warning = function(w) NA_real_)
  list(job_id = Sys.getenv("SLURM_JOB_ID"), cpus = cpus,
       memory_peak_gb = read_num("memory.peak") / 1024^3,
       cpu_seconds = read_num("cpu.stat", "usage_usec") / 1e6)
}

#' Package versions recorded with each result
#'
#' @return Named character vector.
variance_versions <- function() {
  pkgs <- c("variancePartition", "lme4", "limma", "edgeR", "BiocParallel")
  vapply(pkgs, function(p) as.character(utils::packageVersion(p)), character(1))
}

#' Which design variables are constant within subject
#'
#' The vignette's nesting note: subject-level variables (sex, group, age) take
#' share from the subject term, so the between-subject total is their sum.
#'
#' @param design_obj prepare_design() output.
#' @return Variable names (excluding subject itself); empty when cross-sectional.
subject_level_variables <- function(design_obj) {
  spec <- design_obj$spec
  if (is.null(spec[["subject"]])) return(character(0))
  data <- design_obj$data
  vars <- setdiff(design_variables(spec), c(spec[["subject"]], spec[["time"]], unlist(spec$de_covariates)))
  Filter(function(v) all(tapply(as.character(data[[v]]), data[[spec[["subject"]]]],
                                function(x) length(unique(x)) == 1)), vars)
}

#' One dataset's summary (the `vp_summary_<dataset>` target)
#'
#' Per term: distribution of the variance fraction and the share of genes above
#' 5% and 10%; the derived per-gene totals (category_totals()); fit flags; DE
#' genes per contrast and joint test at FDR 5%.
#'
#' @param dataset_id Dataset id.
#' @param design_obj prepare_design() output.
#' @param expr prepare_expression() output.
#' @param norm normalize_for_models() output.
#' @param vp run_variance_partition() output.
#' @param de dream_tables() output.
#' @param dupcor run_duplicate_correlation() output, or `NULL`.
#' @param composition analyze_composition() output, or `NULL`.
#' @return List: `dataset` (one row), `terms`, `contrasts`.
summarize_variance <- function(dataset_id, design_obj, expr, norm, vp, de, dupcor = NULL,
                               composition = NULL) {
  spec <- design_obj$spec
  vpm <- vp$varPart
  derived <- category_totals(design_obj, vpm)
  all_terms <- c(as.list(vpm), derived)
  kind <- function(t) {
    if (t == "Residuals") "residual" else if (t %in% names(derived)) "derived"
    else if (t %in% unlist(spec$numeric)) "fixed" else "random"
  }
  terms <- do.call(rbind, lapply(names(all_terms), function(t) {
    x <- all_terms[[t]]
    q <- stats::quantile(x, c(0.1, 0.5, 0.9), na.rm = TRUE, names = FALSE)
    data.frame(dataset_id = dataset_id, term = t, kind = kind(t), n_genes = sum(!is.na(x)),
               mean = mean(x, na.rm = TRUE), q10 = q[1], median = q[2], q90 = q[3],
               share_gt05 = mean(x > 0.05, na.rm = TRUE), share_gt10 = mean(x > 0.10, na.rm = TRUE))
  }))
  contrasts <- if (nrow(de$contrasts)) cbind(dataset_id = dataset_id, de$contrasts) else
    cbind(dataset_id = character(0), de$contrasts)

  technical <- if (design_obj$platform == "array") "control_mean" else "log2_lib_size"
  dataset <- data.frame(
    dataset_id = dataset_id, platform = design_obj$platform,
    vp_formula = formula_text(design_obj$vp_formula),
    dream_formula = formula_text(design_obj$dream_formula),
    n_samples = nrow(design_obj$data),
    n_subjects = if (is.null(spec[["subject"]])) NA_integer_ else nlevels(design_obj$data[[spec[["subject"]]]]),
    n_dropped = sum(design_obj$dropped$n), n_features_raw = expr$n_features_raw, n_genes = expr$n_genes,
    n_genes_fit = nrow(vpm), n_genes_failed = vp$n_genes_input - nrow(vpm),
    share_singular = mean(vp$flags$singular),
    max_colinearity = max(vp$flags$colinearity, na.rm = TRUE),
    share_colinearity_gt99 = mean(vp$flags$colinearity > 0.99, na.rm = TRUE),
    n_dream_errors = de$n_errors,
    dupcor_consensus = if (is.null(dupcor)) NA_real_ else dupcor$consensus,
    median_subject_fraction = if (is.null(spec[["subject"]])) NA_real_ else stats::median(vpm[[spec[["subject"]]]]),
    rho_PC1_nl = composition_rho(composition, "PC1", "nl_score"),
    rho_sd_nl = composition_rho(composition, "expr_sd", "nl_score"),
    rho_mean_nl = composition_rho(composition, "expr_mean", "nl_score"),
    rho_sd_technical = composition_rho(composition, "expr_sd", technical),
    rho_mean_technical = composition_rho(composition, "expr_mean", technical),
    n_control_probes = if (is.null(composition)) NA_integer_ else composition$control$n_probes,
    sex_separation = if (is.null(composition)) NA_real_ else composition$sex_check$separation,
    sex_agreement = if (is.null(composition)) NA_real_ else composition$sex_check$agreement,
    n_subjects_sex_discordant = if (is.null(composition)) NA_integer_ else composition$sex_check$n_subjects_discordant,
    nl_subject_fraction = trait_fraction(composition, "nl_score", spec[["subject"]]),
    nl_time_fraction = trait_fraction(composition, "nl_score", spec[["time"]]),
    nl_group_fraction = trait_fraction(composition, "nl_score", spec[["group"]]),
    subject_level_vars = paste(subject_level_variables(design_obj), collapse = ","),
    norm_seconds = norm$seconds, vp_seconds = vp$seconds, dream_seconds = de$seconds,
    norm_mem_gb = norm$resources$memory_peak_gb %||% NA_real_,
    vp_mem_gb = vp$resources$memory_peak_gb %||% NA_real_,
    dream_mem_gb = de$resources$memory_peak_gb %||% NA_real_,
    vp_cpu_efficiency = cpu_efficiency(vp),
    dream_cpu_efficiency = cpu_efficiency(de)
  )
  list(dataset = dataset, terms = terms, contrasts = contrasts)
}

#' A Spearman correlation from analyze_composition()
#'
#' @param composition analyze_composition() output, or `NULL`.
#' @param a,b Trait names.
#' @return rho, or `NA` when either trait is absent.
composition_rho <- function(composition, a, b) {
  if (is.null(composition)) return(NA_real_)
  cc <- composition$correlations
  hit <- (cc$x == a & cc$y == b) | (cc$x == b & cc$y == a)
  if (any(hit)) cc$rho[hit][1] else NA_real_
}

#' A sample-level trait's variance fraction for one term
#'
#' @param composition analyze_composition() output, or `NULL`.
#' @param trait Trait name (e.g. `"nl_score"`).
#' @param term Term name, or `NULL`.
#' @return The fraction, or `NA`.
trait_fraction <- function(composition, trait, term) {
  tv <- composition$trait_varpart
  if (is.null(tv) || is.null(term) || !term %in% names(tv) || !trait %in% tv$trait) return(NA_real_)
  tv[[term]][tv$trait == trait]
}

#' CPU efficiency of a variance task
#'
#' @param x A run_variance_partition() or dream_tables() result.
#' @return CPU seconds over (elapsed seconds x CPUs), or `NA`.
cpu_efficiency <- function(x) {
  r <- x$resources
  if (is.null(r)) return(NA_real_)
  r$cpu_seconds / (x$seconds * r$cpus)
}

#' The dataset whose design applies to a dataset's fits
#'
#' A dataset with its own `design:` block is its own parent; a per-timepoint
#' subset uses its family's parent (config/dataset_families.yml) when that
#' parent is in the pipeline and has a design.
#'
#' @param dataset_id Dataset id.
#' @param with_design variance_datasets() output.
#' @param families config/dataset_families.yml, parsed.
#' @return The parent dataset id, or `NULL` when no design applies.
variance_parent <- function(dataset_id, with_design, families) {
  if (dataset_id %in% with_design) return(dataset_id)
  parent <- family_of(families, dataset_id)
  if (length(parent) == 1 && parent %in% with_design) parent else NULL
}

#' Pipeline datasets with a `design:` block
#'
#' @param pipeline The parsed pipeline config.
#' @return Dataset ids.
variance_datasets <- function(pipeline) {
  Filter(function(id) !is.null(yaml::read_yaml(dataset_config_path(id))$design), pipeline$datasets)
}

#' Variance targets for one dataset
#'
#' `design_`, `expr_`, `norm_`, `vp_`, `dream_`, `de_tables_`, `composition_`
#' (R/targets/composition.R), `vp_explore_`, `vp_summary_` and `vp_qc_<dataset>`
#' (R/targets/variance_qc.R), plus `dupcor_<dataset>` for repeated-measures
#' datasets, `controls_<dataset>` for arrays with `control_probes`, and
#' `control_factor_<dataset>` when those are a covariate.
#' Reuses tar_dataset_fits()'s `meta_` and `sample_meta_` targets; nothing here
#' depends on a fit, so the factor fits are never invalidated by this module.
#'
#' @param dataset_id Dataset id.
#' @param pipeline The parsed pipeline config.
#' @return List of target objects.
tar_dataset_variance <- function(dataset_id, pipeline) {
  key <- target_key(dataset_id)
  name <- function(...) paste(c(..., key), collapse = "_")
  sym <- function(...) as.symbol(name(...))
  on <- function(controller) targets::tar_resources(crew = targets::tar_resources_crew(controller = controller))
  cpus <- as.integer(pipeline$controllers$vp$cpus_per_task %||% 1L)
  # The duplicateCorrelation comparison needs repeated measures.
  spec <- yaml::read_yaml(dataset_config_path(dataset_id))$design
  repeated <- !is.null(spec[["subject"]])
  rnaseq <- identical(spec$platform, "rnaseq")
  dupcor <- if (repeated) sym("dupcor") else NULL
  has_controls <- !is.null(spec$control_probes)
  control_cov <- isTRUE(spec$control_probes$covariate)
  out <- list(
    targets::tar_target_raw(
      name("design"),
      substitute(prepare_design(META, SM, CF),
                 list(META = sym("meta"), SM = sym("sample_meta"),
                      CF = if (control_cov) sym("control_factor") else NULL)),
      deployment = "main"
    ),
    targets::tar_target_raw(
      name("expr"),
      substitute(prepare_expression(META, DESIGN), list(META = sym("meta"), DESIGN = sym("design"))),
      resources = on("prep")
    ),
    # Arrays pass through on the main process; RNA-seq runs voom on a worker.
    targets::tar_target_raw(
      name("norm"),
      substitute(normalize_for_models(EXPR, DESIGN, CPUS),
                 list(EXPR = sym("expr"), DESIGN = sym("design"), CPUS = cpus)),
      deployment = if (rnaseq) "worker" else "main",
      resources = if (rnaseq) on("vp") else targets::tar_resources()
    ),
    targets::tar_target_raw(
      name("vp"),
      substitute(run_variance_partition(NORM, DESIGN, EXPR$genes, CPUS),
                 list(NORM = sym("norm"), DESIGN = sym("design"), EXPR = sym("expr"), CPUS = cpus)),
      resources = on("vp")
    ),
    targets::tar_target_raw(
      name("dream"),
      substitute(run_dream(NORM, DESIGN, CPUS), list(NORM = sym("norm"), DESIGN = sym("design"), CPUS = cpus)),
      resources = on("vp")
    ),
    targets::tar_target_raw(
      name("de_tables"),
      substitute(dream_tables(DREAM), list(DREAM = sym("dream"))),
      # The only reader of the dream fit. variancePartition is loaded before
      # the fit is read: a dim() call on an MArrayLM2 before the namespace is
      # loaded returns NULL, and stays NULL for the session, which breaks
      # topTable()'s joint F-test.
      packages = "variancePartition",
      resources = on("ingest")
    ),
    targets::tar_target_raw(
      name("composition"),
      substitute(analyze_composition(META, DESIGN, EXPR, NORM, read_composition_markers(MARKERS), CONTROLS),
                 list(META = sym("meta"), DESIGN = sym("design"), EXPR = sym("expr"), NORM = sym("norm"),
                      MARKERS = quote(composition_markers),
                      CONTROLS = if (has_controls) sym("controls") else NULL)),
      resources = on("ingest")
    ),
    targets::tar_target_raw(
      name("vp_explore"),
      substitute(explore_variance(ID, META, DESIGN, NORM, VP, COMP, MAT),
                 list(ID = dataset_id, META = sym("meta"), DESIGN = sym("design"), NORM = sym("norm"),
                      VP = sym("vp"), COMP = sym("composition"), MAT = sym("mat"))),
      resources = on("ingest")
    ),
    targets::tar_target_raw(
      name("vp_summary"),
      substitute(summarize_variance(ID, DESIGN, EXPR, NORM, VP, DE, DUPCOR, COMP),
                 list(ID = dataset_id, DESIGN = sym("design"), EXPR = sym("expr"), NORM = sym("norm"),
                      VP = sym("vp"), DE = sym("de_tables"), DUPCOR = dupcor, COMP = sym("composition"))),
      resources = on("ingest")
    ),
    targets::tar_target_raw(
      name("vp_qc"),
      substitute(write_variance_qc(ID, DESIGN, NORM, VP, DE, SUMMARY, DUPCOR, COMP),
                 list(ID = dataset_id, DESIGN = sym("design"), NORM = sym("norm"), VP = sym("vp"),
                      DE = sym("de_tables"), SUMMARY = sym("vp_summary"), DUPCOR = dupcor,
                      COMP = sym("composition"))),
      format = "file",
      resources = on("ingest")
    )
  )
  if (has_controls) {
    out <- c(out, list(targets::tar_target_raw(
      name("controls"),
      substitute(control_probe_matrix(META), list(META = sym("meta"))),
      resources = on("prep")
    )))
  }
  if (control_cov) {
    out <- c(out, list(targets::tar_target_raw(
      name("control_factor"),
      substitute(compute_control_factor(CONTROLS), list(CONTROLS = sym("controls"))),
      deployment = "main"
    )))
  }
  if (repeated) {
    out <- c(out, list(targets::tar_target_raw(
      name("dupcor"),
      substitute(run_duplicate_correlation(EXPR, NORM, DESIGN),
                 list(EXPR = sym("expr"), NORM = sym("norm"), DESIGN = sym("design"))),
      resources = on("ingest")
    )))
  }
  out
}

#' The cross-dataset variance targets
#'
#' `composition_markers` (the marker-gene file) and `variance_explore`
#' (combine_variance_explore() over every dataset with a design).
#'
#' @param dataset_ids Datasets with variance targets.
#' @return List of target objects.
tar_variance_cross_dataset <- function(dataset_ids) {
  named <- function(prefix) as.call(c(quote(list), stats::setNames(
    lapply(dataset_ids, function(id) as.symbol(paste0(prefix, target_key(id)))), dataset_ids)))
  list(
    targets::tar_target_raw("composition_markers", "config/composition_markers.yml", format = "file",
                            deployment = "main"),
    targets::tar_target_raw(
      "variance_explore",
      substitute(combine_variance_explore(EXPLORE, SUMMARIES),
                 list(EXPLORE = named("vp_explore_"), SUMMARIES = named("vp_summary_"))),
      deployment = "main"
    )
  )
}

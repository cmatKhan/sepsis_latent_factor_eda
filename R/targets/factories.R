# Target factories (targets manual, "Target factories"): every target for one
# dataset, built from its config and the method registries, plus the
# cross-dataset targets. tar_target_raw() rather than tar_map() because each
# method's fit target needs its own packages and crew controller. The full
# list of targets per dataset is in the docs (Running the workflow, "Targets
# per dataset").

#' A dataset id made safe for target names
#'
#' @param dataset_id Dataset id.
#' @return The id with characters other than letters, digits and `_` replaced by `_`.
target_key <- function(dataset_id) gsub("[^A-Za-z0-9_]", "_", dataset_id)

#' Find a dataset's config file
#'
#' Config file names don't always match ids (`ANEMONES_day1_config.yml`
#' declares `ANEMONES_DAY1`), so this matches on the declared `dataset.id`.
#'
#' @param dataset_id Dataset id.
#' @return Path to the one `config/*_config.yml` declaring it (errors if not exactly one).
dataset_config_path <- function(dataset_id) {
  files <- Sys.glob("config/*_config.yml")
  ids <- vapply(files, function(f) yaml::read_yaml(f)$dataset$id %||% NA_character_, character(1))
  hit <- files[!is.na(ids) & ids == dataset_id]
  if (length(hit) != 1) stop("expected exactly one config/*_config.yml with dataset.id '", dataset_id,
                             "', found ", length(hit))
  hit
}

#' Methods configured for a dataset
#'
#' A method runs when it appears under `methods:` (network methods under
#' `methods.network`).
#'
#' @param meta The dataset config.
#' @param manifest discover_method_registries() output.
#' @return Method names, in manifest order.
configured_methods <- function(meta, manifest) {
  flat <- setdiff(names(meta$methods), "network")
  nested <- names(meta$methods$network)
  intersect(names(manifest), c(flat, nested))
}

#' The files a dataset's matrix is built from
#'
#' The config, preprocessing script, prebuilt matrix and raw data files; the
#' `inputs_<dataset>` file target, so editing any of them rebuilds that dataset.
#'
#' @param config_path Path to the dataset config.
#' @param meta The resolved dataset config.
#' @return Character vector of existing paths (errors on a missing one).
dataset_input_files <- function(config_path, meta) {
  ds <- meta$dataset
  files <- c(config_path, ds$matrix_path, ds$preprocessing_script,
             ds$expression_path, ds$sample_metadata_path, ds$feature_metadata_path)
  missing <- files[!file.exists(files)]
  if (length(missing) > 0) {
    stop("Dataset '", ds$id, "' references files that don't exist on this machine: ",
         paste(missing, collapse = ", "),
         " -- check config/pipeline.yml's data_root (or SEPSIS_DATA_ROOT).")
  }
  unname(files)
}

#' Every target for one dataset
#'
#' Inputs, matrix, and per method the grid and fit targets (plus sPCA's refine
#' stage), then tar_dataset_ingest()'s targets.
#'
#' @param dataset_id Dataset id.
#' @param manifest discover_method_registries() output.
#' @param pipeline The parsed pipeline config.
#' @param data_root Data root directory.
#' @return List of target objects.
tar_dataset_fits <- function(dataset_id, manifest, pipeline, data_root) {
  config_path <- dataset_config_path(dataset_id)
  # Read once at definition time to decide which targets exist; the
  # meta_<key> target re-reads it at run time, so edits are tracked.
  meta <- read_dataset_metadata(config_path, manifest = manifest,
                                data_root = data_root)
  key <- target_key(dataset_id)
  name <- function(...) paste(c(..., key), collapse = "_")
  methods <- setdiff(configured_methods(meta, manifest), unlist(pipeline$skip_methods))

  unknown <- setdiff(methods, names(pipeline$controllers))
  if (length(unknown) > 0) {
    stop("No controller in config/pipeline.yml for method(s): ", paste(unknown, collapse = ", "))
  }

  targets <- list(
    targets::tar_target_raw(
      name("inputs"),
      as.call(c(quote(c), as.list(dataset_input_files(config_path, meta)))),
      format = "file",
      deployment = "main"
    ),
    targets::tar_target_raw(
      name("meta"),
      substitute(
        {
          INPUTS
          read_dataset_metadata(CONFIG, data_root = DATA_ROOT)
        },
        list(INPUTS = as.symbol(name("inputs")), CONFIG = config_path, DATA_ROOT = data_root)
      ),
      deployment = "main"
    ),
    targets::tar_target_raw(
      name("mat"),
      substitute(load_input_matrix(META), list(META = as.symbol(name("meta")))),
      resources = targets::tar_resources(crew = targets::tar_resources_crew(controller = "prep"))
    )
  )

  needs_nonneg <- any(vapply(methods, function(m) isTRUE(manifest[[m]]$registry$needs_nonneg), logical(1)))
  if (needs_nonneg) {
    targets <- c(targets, list(targets::tar_target_raw(
      name("mat_nn"),
      substitute(shift_nonneg(MAT), list(MAT = as.symbol(name("mat")))),
      packages = "matrixStats",
      deployment = "main"
    )))
  }

  for (method in methods) {
    registry <- manifest[[method]]$registry
    registry_sym <- as.symbol(paste0(method, "_registry"))
    mat_sym <- as.symbol(if (isTRUE(registry$needs_nonneg)) name("mat_nn") else name("mat"))
    cpus <- pipeline$controllers[[method]]$cpus_per_task
    batch <- as.integer(pipeline$controllers[[method]]$batch_size %||% 1L)

    # One dynamic branch per batch of grid rows, on the method's controller.
    fit_target <- function(fit_name, grid_name) {
      targets::tar_target_raw(
        fit_name,
        substitute(
          run_method_fits(REGISTRY, GRID, MAT),
          list(REGISTRY = registry_sym, GRID = as.symbol(grid_name), MAT = mat_sym)
        ),
        pattern = substitute(map(GRID), list(GRID = as.symbol(grid_name))),
        iteration = "list",
        packages = registry$pkgs %||% character(0),
        resources = targets::tar_resources(crew = targets::tar_resources_crew(controller = method))
      )
    }

    targets <- c(targets, list(
      targets::tar_target_raw(
        name("grid", method),
        substitute(
          build_method_grid(REGISTRY, META, METHOD, CPUS, BATCH),
          list(REGISTRY = registry_sym, META = as.symbol(name("meta")), METHOD = method,
               CPUS = cpus, BATCH = batch)
        ),
        iteration = "group",
        deployment = "main"
      ),
      fit_target(name("fit", method), name("grid", method))
    ))

    # Optional data-driven second stage (registry `refine_grid`, e.g. spca's
    # Index-of-Sparseness search): its grid depends on the first stage's
    # fits and the dataset's PCA fit; `select_fits` then tabulates both.
    if (!is.null(registry$refine_grid)) {
      if (!"pca" %in% methods) {
        stop(dataset_id, ": method '", method, "' has a refine stage that needs the dataset's ",
             "PCA fit, but methods.pca isn't configured.")
      }
      targets <- c(targets, list(
        targets::tar_target_raw(
          name("grid", method, "refine"),
          substitute(
            build_refine_grid(REGISTRY, META, METHOD, CPUS, BATCH, FITS, PCA),
            list(REGISTRY = registry_sym, META = as.symbol(name("meta")), METHOD = method,
                 CPUS = cpus, BATCH = batch, FITS = as.symbol(name("fit", method)),
                 PCA = as.symbol(name("fit", "pca")))
          ),
          iteration = "group",
          deployment = "main"
        ),
        fit_target(name("fit", method, "refine"), name("grid", method, "refine")),
        targets::tar_target_raw(
          name("selection", method),
          substitute(
            select_method_fits(REGISTRY, FITS, REFINE, PCA),
            list(REGISTRY = registry_sym, FITS = as.symbol(name("fit", method)),
                 REFINE = as.symbol(name("fit", method, "refine")),
                 PCA = as.symbol(name("fit", "pca")))
          ),
          deployment = "main"
        )
      ))
    }
  }
  c(targets, tar_dataset_ingest(dataset_id, methods, manifest, pipeline, name))
}

#' Ingestion targets for one dataset
#'
#' Per method: ingest, pairs, representatives, redundancy, ICASSO and
#' enrichment targets; WGCNA's ORA, per-fit kME and dataset diagnostics; per
#' dataset: metadata, matrix file, diagnostics, enrichment rows, drivers and
#' `db_<dataset>`, which writes the dataset's DB rows.
#'
#' @param dataset_id Dataset id.
#' @param methods The dataset's configured methods.
#' @param manifest discover_method_registries() output.
#' @param pipeline The parsed pipeline config.
#' @param name Function building this dataset's target names.
#' @return List of target objects.
tar_dataset_ingest <- function(dataset_id, methods, manifest, pipeline, name) {
  db_path <- pipeline$db_path %||% "results/targets/stability.sqlite"
  on_ingest <- targets::tar_resources(crew = targets::tar_resources_crew(controller = "ingest"))
  on_fgsea <- targets::tar_resources(crew = targets::tar_resources_crew(controller = "fgsea"))
  fgsea_cpus <- pipeline$controllers$fgsea$cpus_per_task %||% 1L
  sym <- function(...) as.symbol(name(...))
  spec_sym <- function(m) as.symbol(paste0(m, "_ingest"))
  targets <- list()
  parts <- list()
  fgsea_parts <- list()
  wgcna_parts <- list()

  for (method in methods) {
    registry <- manifest[[method]]$registry
    spec <- manifest[[method]]$ingest
    has_refine <- !is.null(registry$refine_grid)
    fits_expr <- if (has_refine) {
      substitute(c(flatten_fits(A), flatten_fits(B)),
                 list(A = sym("fit", method), B = sym("fit", method, "refine")))
    } else {
      substitute(flatten_fits(A), list(A = sym("fit", method)))
    }
    targets <- c(targets, list(
      targets::tar_target_raw(
        name("ingest", method),
        substitute(ingest_method_fits(FITS, SPEC, METHOD, ID, DB),
                   list(FITS = fits_expr, SPEC = spec_sym(method), METHOD = method, ID = dataset_id,
                        DB = db_path)),
        resources = on_ingest
      ),
      targets::tar_target_raw(
        name("pairs", method),
        substitute(compute_method_pairs(INGEST, SPEC, DB),
                   list(INGEST = sym("ingest", method), SPEC = spec_sym(method), DB = db_path)),
        packages = c("clue", "mclust"),
        resources = on_ingest
      ),
      targets::tar_target_raw(
        name("representatives", method),
        substitute(representative_keys(INGEST, SPEC, SELECTION),
                   list(INGEST = sym("ingest", method), SPEC = spec_sym(method),
                        SELECTION = if (has_refine) sym("selection", method) else NULL)),
        deployment = "main"
      )
    ))
    part <- list(ingest = sym("ingest", method), pairs = sym("pairs", method),
                 representatives = sym("representatives", method),
                 selection = if (has_refine) sym("selection", method) else NULL)

    if (isTRUE(spec$has_loadings)) {
      if (method != "pca") {   # PCA's components are orthogonal by construction
        targets <- c(targets, list(targets::tar_target_raw(
          name("redundancy", method),
          substitute(compute_method_redundancy(INGEST, REPS, DB),
                     list(INGEST = sym("ingest", method), REPS = sym("representatives", method), DB = db_path)),
          packages = if (method == "cogaps") "CoGAPS" else character(0),
          resources = on_ingest
        )))
        part$redundancy <- sym("redundancy", method)
      }
      if (method == "ica") {
        targets <- c(targets, list(targets::tar_target_raw(
          name("icasso", method),
          substitute(compute_method_icasso(INGEST, DB), list(INGEST = sym("ingest", method), DB = db_path)),
          resources = on_ingest
        )))
        part$icasso <- sym("icasso", method)
      }
      # Enrichment: one branch per representative fit, trimmed to the
      # significant terms (R/targets/enrichment.R).
      targets <- c(targets, list(
        targets::tar_target_raw(
          name("fgsea_jobs", method),
          substitute(fgsea_jobs(INGEST, REPS, RED, DB),
                     list(INGEST = sym("ingest", method), REPS = sym("representatives", method),
                          RED = part$redundancy, DB = db_path)),
          deployment = "main"
        ),
        targets::tar_target_raw(
          name("fgsea", method),
          substitute(
            if (has_job(JOBS$fit_key)) {
              trim_enrichment_result(run_fgsea_job(ID, METHOD, JOBS$fit_key, JOBS$loadings_file, msigdb_pathways,
                                                   ENSEMBL, JOBS$cogaps_marker_genes[[1]], CPUS,
                                                   both_directions = SIGNED))
            },
            list(ID = dataset_id, METHOD = method, JOBS = sym("fgsea_jobs", method),
                 ENSEMBL = sym("ensembl_map"), CPUS = fgsea_cpus, SIGNED = isTRUE(spec$sign_ambiguous))),
          pattern = substitute(map(J), list(J = sym("fgsea_jobs", method))),
          iteration = "list",
          packages = c("fgsea", "BiocParallel"),
          resources = on_fgsea
        )
      ))
      fgsea_parts <- c(fgsea_parts, sym("fgsea", method))
    } else {
      # WGCNA: module ORA per fit, kME per fit, dataset-level diagnostics.
      targets <- c(targets, list(
        targets::tar_target_raw(
          name("wgcna_ora_jobs"),
          substitute(wgcna_ora_jobs(INGEST), list(INGEST = sym("ingest", method))),
          deployment = "main"
        ),
        targets::tar_target_raw(
          name("wgcna_ora"),
          substitute(
            if (has_job(JOBS$fit_key)) {
              trim_enrichment_result(run_wgcna_ora_job(ID, JOBS$fit_key, JOBS$module_genes[[1]],
                                                       JOBS$universe_genes[[1]], msigdb_pathways, ENSEMBL, CPUS))
            },
            list(ID = dataset_id, JOBS = sym("wgcna_ora_jobs"), ENSEMBL = sym("ensembl_map"),
                 CPUS = fgsea_cpus)),
          pattern = substitute(map(J), list(J = sym("wgcna_ora_jobs"))),
          iteration = "list",
          packages = c("fgsea", "BiocParallel"),
          resources = on_fgsea
        ),
        targets::tar_target_raw(
          name("wgcna_fit_jobs"),
          substitute(unique(INGEST$modules$fit_key) %||% NA_character_, list(INGEST = sym("ingest", method))),
          deployment = "main"
        ),
        targets::tar_target_raw(
          name("wgcna_fit_diag"),
          substitute(if (has_job(KEY)) wgcna_fit_diagnostics(KEY, INGEST, MAT, DB),
                     list(KEY = sym("wgcna_fit_jobs"), INGEST = sym("ingest", method), MAT = sym("mat"),
                          DB = db_path)),
          pattern = substitute(map(J), list(J = sym("wgcna_fit_jobs"))),
          iteration = "list",
          resources = on_ingest
        ),
        targets::tar_target_raw(
          name("wgcna_dataset_diag"),
          substitute(wgcna_dataset_diagnostics(MAT, META, SM),
                     list(MAT = sym("mat"), META = sym("meta"), SM = sym("sample_meta"))),
          packages = "WGCNA",
          resources = on_ingest
        )
      ))
      wgcna_parts <- c(wgcna_parts, sym("wgcna_ora"))
      part$wgcna_fits <- sym("wgcna_fit_diag")
    }
    parts[[method]] <- as.call(c(quote(list), Filter(Negate(is.null), part)))
  }

  loadings_methods <- Filter(function(m) isTRUE(manifest[[m]]$ingest$has_loadings), methods)
  rep_list <- as.call(c(quote(list), lapply(setNames(nm = loadings_methods), function(m) sym("representatives", m))))
  ing_list <- as.call(c(quote(list), lapply(setNames(nm = loadings_methods), function(m) sym("ingest", m))))
  c(targets, list(
    targets::tar_target_raw(
      name("sample_meta"),
      substitute(read_sample_metadata(META), list(META = sym("meta"))),
      deployment = "main"
    ),
    targets::tar_target_raw(
      name("metadata_files"),
      substitute(write_metadata_files(META, SM, DB), list(META = sym("meta"), SM = sym("sample_meta"), DB = db_path)),
      format = "file",
      deployment = "main"
    ),
    # A file target: projections onto this dataset read it.
    targets::tar_target_raw(
      name("matrix_file"),
      substitute(dataset_matrix_file(MAT, ID, DB), list(MAT = sym("mat"), ID = dataset_id, DB = db_path)),
      format = "file",
      deployment = "main"
    ),
    targets::tar_target_raw(
      name("matrix_diagnostics"),
      substitute(matrix_diagnostics_row(MAT), list(MAT = sym("mat"))),
      resources = on_ingest
    ),
    targets::tar_target_raw(
      name("ensembl_map"),
      substitute(build_ensembl_map(META), list(META = sym("meta"))),
      deployment = "main"
    ),
    targets::tar_target_raw(
      name("enrichment"),
      # Each part is a pattern's list of branch values: flatten one level.
      substitute(enrichment_rows(unlist(F, recursive = FALSE), unlist(W, recursive = FALSE)),
                 list(F = as.call(c(quote(list), fgsea_parts)), W = as.call(c(quote(list), wgcna_parts)))),
      deployment = "main"
    ),
    targets::tar_target_raw(
      name("drivers"),
      substitute(pattern_driver_rows(INGESTS, REPS, MAT, SM, META$dataset$sample_id_col %||% "sample_id", ID, DB),
                 list(INGESTS = ing_list, REPS = rep_list, MAT = sym("mat"), SM = sym("sample_meta"),
                      META = sym("meta"), ID = dataset_id, DB = db_path)),
      packages = "projectR",
      resources = on_ingest
    ),
    targets::tar_target_raw(
      name("db"),
      substitute(
        write_dataset_db(DB, db_schema, META, c(MF, MD), METHODS, METHOD_TABLE, ENRICH,
                         MATDIAG, WGCNA, DRIVERS),
        list(DB = db_path, META = sym("meta"), MF = sym("matrix_file"), MD = sym("metadata_files"),
             METHODS = as.call(c(quote(list), parts)), ENRICH = sym("enrichment"),
             MATDIAG = sym("matrix_diagnostics"),
             WGCNA = if ("wgcna" %in% methods) sym("wgcna_dataset_diag") else NULL,
             DRIVERS = sym("drivers"))),
      packages = c("DBI", "RSQLite"),
      deployment = "main"
    )
  ))
}

#' The `methods` table rows
#'
#' @param manifest discover_method_registries() output.
#' @return Data frame (method, family, sign_ambiguous, has_loadings), from
#'   each method's ingestion contract.
method_table <- function(manifest) {
  do.call(rbind, lapply(names(manifest), function(m) {
    s <- manifest[[m]]$ingest
    data.frame(method = m, family = s$family, sign_ambiguous = isTRUE(s$sign_ambiguous),
               has_loadings = isTRUE(s$has_loadings))
  }))
}

#' Cross-dataset targets
#'
#' The DB schema file target, gathered Ensembl maps and matrix files, dataset
#' families, per dataset its projectR jobs and gene-space agreement, and
#' `db_cross_dataset`, which writes them after every dataset's DB write.
#'
#' @param pipeline The parsed pipeline config.
#' @param manifest discover_method_registries() output.
#' @param data_root Data root directory.
#' @return List of target objects.
tar_cross_dataset <- function(pipeline, manifest, data_root) {
  db_path <- pipeline$db_path %||% "results/targets/stability.sqlite"
  ids <- unlist(pipeline$datasets)
  keys <- setNames(target_key(ids), ids)
  families <- yaml::read_yaml("config/dataset_families.yml")
  loadings_methods_of <- lapply(setNames(nm = ids), function(id) {
    meta <- read_dataset_metadata(dataset_config_path(id), manifest = manifest, data_root = data_root)
    ms <- setdiff(configured_methods(meta, manifest), unlist(pipeline$skip_methods))
    Filter(function(m) isTRUE(manifest[[m]]$ingest$has_loadings), ms)
  })
  named_list <- function(prefix, which = ids) {
    as.call(c(quote(list), lapply(which, function(id) as.symbol(paste(prefix, keys[[id]], sep = "_")))
              |> setNames(which)))
  }
  parts_of <- function(id) {
    as.call(c(quote(list), lapply(setNames(nm = loadings_methods_of[[id]]), function(m) {
      substitute(list(ingest = I, representatives = R),
                 list(I = as.symbol(paste("ingest", m, keys[[id]], sep = "_")),
                      R = as.symbol(paste("representatives", m, keys[[id]], sep = "_"))))
    })))
  }

  targets <- list(
    targets::tar_target_raw("db_schema", quote(c("R/db/schema.sql", "R/db/views.sql")), format = "file",
                            deployment = "main"),
    targets::tar_target_raw("ensembl_maps", named_list("ensembl_map"), deployment = "main"),
    targets::tar_target_raw("matrix_files", substitute(unlist(X), list(X = named_list("matrix_file"))),
                            format = "file", deployment = "main"),
    targets::tar_target_raw("dataset_families_file", "config/dataset_families.yml", format = "file",
                            deployment = "main"),
    targets::tar_target_raw("dataset_families", quote(yaml::read_yaml(dataset_families_file)),
                            deployment = "main")
  )
  on_projectr <- targets::tar_resources(crew = targets::tar_resources_crew(controller = "projectr"))
  projection_syms <- list(); gene_space_syms <- list()
  for (id in ids) {
    key <- keys[[id]]
    nm <- function(prefix) paste(prefix, key, sep = "_")
    if (length(ids) > 1) {
      targets <- c(targets, list(
        targets::tar_target_raw(
          nm("projectr_jobs"),
          substitute(projectr_jobs(ID, PARTS, matrix_files, dataset_families, DB),
                     list(ID = id, PARTS = parts_of(id), DB = db_path)),
          iteration = "group",
          deployment = "main"
        ),
        targets::tar_target_raw(
          nm("projectr"),
          substitute(projection_rows(JOBS, ID, ensembl_maps, MATRIX_FILE, DB),
                     list(JOBS = as.symbol(nm("projectr_jobs")), ID = id,
                          MATRIX_FILE = as.symbol(nm("matrix_file")), DB = db_path)),
          pattern = substitute(map(J), list(J = as.symbol(nm("projectr_jobs")))),
          iteration = "list",
          packages = "projectR",
          resources = on_projectr
        )
      ))
      projection_syms <- c(projection_syms, as.symbol(nm("projectr")))
    }
    kin <- Filter(function(b) b != id && same_family(families, id, b), ids)
    if (length(kin)) {
      targets <- c(targets, list(targets::tar_target_raw(
        nm("gene_space"),
        substitute(gene_space_rows(ID, PARTS, ensembl_maps, DB),
                   list(ID = id, DB = db_path,
                        PARTS = as.call(c(quote(list), lapply(setNames(nm = c(id, kin)), parts_of))))),
        packages = "projectR",
        resources = targets::tar_resources(crew = targets::tar_resources_crew(controller = "ingest"))
      )))
      gene_space_syms <- c(gene_space_syms, as.symbol(nm("gene_space")))
    }
  }
  c(targets, list(targets::tar_target_raw(
    "db_cross_dataset",
    substitute(write_cross_dataset_db(DB, WRITTEN, unlist(P, recursive = FALSE), G),
               list(DB = db_path, WRITTEN = named_list("db"),
                    P = as.call(c(quote(list), projection_syms)),
                    G = as.call(c(quote(list), gene_space_syms)))),
    packages = c("DBI", "RSQLite"),
    deployment = "main"
  )))
}

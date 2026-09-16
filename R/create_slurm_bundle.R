# Create a directory that may be executed, via the submit_all_<dataset>.sh
# script inside of that directory, via slurm
#
# Usage (interactive):
#   output_dir <- here("slurm_bundles")
#   dataset_metadata_path <- "config/GSE110487_config.yml"
#   # optional
#   only_methods = c("pca", "cogaps")
#   source("R/create_slurm_bundle.R")

library(here)
library(optparse)
source(here("R/lib/metadata.R"))
source(here("R/lib/matrices.R"))
source(here("R/lib/tensors.R"))
source(here("R/lib/grids.R"))
source(here("R/lib/submit.R"))
source(here("R/lib/method_registry.R"))

option_list <- list(
    make_option(c("--config"),
        type = "character",
        default = NULL,
        help = "path to the dataset configuration",
        metavar = "path",
        required = TRUE
    ),
    make_option(c("--only"),
        type = "character", default = NULL,
        help = "Comma-separated list of methods to run (e.g., pca,cogaps)",
        metavar = "methods"
    ),
    make_option(c("--output"),
        type = "character",
        default = here("slurm_bundles"),
        help = "path to the output directory. output will be to <output>/<config basename>",
        metavar = "path"
    )
)

parser <- OptionParser(
    usage = "Usage: Rscript %prog [options]",
    option_list = option_list,
    description = paste0(
        "Given a dataset config, use rslurm to create ",
        "a directory in ./slurm_bundles/<yaml_basename> ",
        "with slurm bundles for each of the configured methods ",
        "that may be run on a slurm cluster"
    )
)

if (!interactive()) {
    parsed_args <- parse_args(parser)
    if (!file.exists(parsed_args$config)) {
        stop(sprintf("The config path does not exist: %s", parsed_args$config))
    }
    dataset_metadata_path <- parsed_args$config
    if (!dir.exists(parsed_args$output)) {
        message("creating output directory at: ", parsed_args$output)
        dir.create(parsed_args$output)
    }

    output_dir <- parsed_args$output

    only_methods <- if (!is.null(parsed_args$only)) {
        strsplit(parsed_args$only, ",")[[1]]
    } else {
        NULL
    }
}

if (interactive()) {
    if (!file.exists(dataset_metadata_path)) {
        stop("`dataset_metadata_path` DNE: ", dataset_metadata_path)
    }

    if (!exists("output_dir")) {
        message("`output_dir` set to 'slurm_bundles'")
        output_dir <- "slurm_bundles"
    }
    only_methods <- if (exists("only_methods")) {
        strsplit(only_methods, ",")[[1]]
    } else {
        NULL
    }
}

# this stores a list where hte top level names are
# methods discovered automatically from R/methods.
# See R/README.md#adding-a-new-method for more details
# about the information that is stored for each method
METHOD_MANIFEST <- discover_method_registries()

dataset_meta <- read_dataset_metadata(dataset_metadata_path, manifest = METHOD_MANIFEST)
dataset_id <- dataset_meta$dataset$id
job_output_dir <- file.path(output_dir, dataset_id)

# detect which methods are active in the config and according to `only_methods`.
# This is a subset of METHOD_MANIFEST
active_specs <- list()
for (name in names(METHOD_MANIFEST)) {
    if (!is.null(only_methods) && !(name %in% only_methods)) next
    spec <- METHOD_MANIFEST[[name]]
    # This extracts the method specification from the dataset config
    method_meta <- if (spec$network) dataset_meta$methods$network[[name]] else dataset_meta$methods[[name]]
    # if method_meta is null, then this particular method is not enabled for
    # this dataset and it is skipped
    if (is.null(method_meta)) next
    slurm_cfg <- if (spec$network) dataset_meta$slurm$network[[name]] else dataset_meta$slurm[[name]]
    active_specs[[name]] <- list(spec = spec, method_meta = method_meta, slurm_cfg = slurm_cfg)
}

# The raw matrix is always needed; the non-negative-shifted variant and/or
# the tensor reshape are each built at most once, and only if some active
# method's registry actually asks for them (registry$needs_nonneg /
# registry$needs_tensor) -- no lazy cache required now that "active" is
# known up front.
raw_mat <- load_input_matrix(dataset_meta)
# if any active method needs the non-negative matrix, set this flag to TRUE.
# else FALSE
needs_nonneg <- any(vapply(active_specs, function(s) isTRUE(s$spec$registry$needs_nonneg), logical(1)))
# if any active method needs a tensor matrix (ie it has multiple timepionts),
# set this flag to TRUE. else FALSE
needs_tensor <- any(vapply(active_specs, function(s) isTRUE(s$spec$registry$needs_tensor), logical(1)))
nn_mat <- if (needs_nonneg) shift_nonneg(raw_mat) else NULL
tnsr_arr <- if (needs_tensor) build_tensor(raw_mat, dataset_meta) else NULL

sjobs <- list()

for (name in names(active_specs)) {
    active <- active_specs[[name]]
    spec <- active$spec
    method_meta <- active$method_meta
    slurm_cfg <- active$slurm_cfg
    registry <- spec$registry

    mat_name <- registry$global_object
    mat_value <- if (isTRUE(registry$needs_tensor)) {
        tnsr_arr
    } else if (isTRUE(registry$needs_nonneg)) {
        nn_mat
    } else {
        raw_mat
    }
    assign(mat_name, mat_value, envir = .GlobalEnv)

    message("== ", toupper(name), " ==")

    # method_meta is the dataset config's ENTIRE override for this method
    # (e.g. `dataset_meta$methods$pca`) -- no per-method wrapper key to
    # unwrap anymore, since each method now describes exactly one kind of
    # run. merge_named_list() recurses one level into nested named-list
    # values (e.g. CoGAPS's `params` sub-block), so overriding one nested
    # key doesn't wipe out its siblings' defaults.
    resolved <- merge_named_list(registry$defaults, method_meta)
    if (!is.null(registry$resource_defaults)) resolved <- registry$resource_defaults(resolved, slurm_cfg)

    grid <- registry$build_grid(resolved)

    sjob <- submit_job_family(
        f              = registry$fn,
        jobs_df        = grid,
        jobname        = registry$jobname,
        global_objects = mat_name,
        pkgs           = registry$pkgs %||% character(0),
        cluster_cfg    = slurm_cfg,
        output_dir     = job_output_dir
    )
    message("  [", registry$jobname, "] ", nrow(grid), " job(s)")

    sjobs[[name]] <- sjob
    save_sjobs(sjob, dataset_id = dataset_id, method = name)
}

# Regenerates slurm_bundles/<dataset_id>/submit_all_<dataset_id>.sh from
# every job family recorded so far (this run's methods plus any others
# already run for this dataset) -- safe/expected after every setup run.
write_submit_all_script(dataset_id, job_output_dir)

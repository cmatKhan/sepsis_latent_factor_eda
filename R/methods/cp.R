# CP / CANDECOMP-PARAFAC (rTensor::cp()) method registry -- a genuinely
# different METHOD from Tucker for every purpose in this framework (own
# config block, own ingest branch/artifacts) despite sharing
# R/lib/tensors.R's build_tensor() and the rTensor package.
#
# Deterministic given `num_components` (no seed/init argument).
#
# Produces THREE factor matrices, tied together by one shared rank (unlike
# Tucker, whose modes can have different ranks): gene loadings, SUBJECT-mode
# "scores" (rownames = subject ids from dataset.subject_id_col -- NOT
# sample ids), and time-mode loadings (rownames = timepoint levels). `mse`
# is `1 - norm_percent/100` so it's on the same "lower is better" scale as
# every other method.
#
# Requires `dataset.subject_id_col`/`dataset.timepoint_col` (see
# R/lib/tensors.R::build_tensor()) -- validated in R/lib/metadata.R via
# `requires_subject_timepoint` below.
#
# Operates on the global object `tnsr` (an rTensor Tensor, built once via
# build_tensor() and shared with Tucker -- see R/create_slurm_bundle.R).

run_cp_job <- function(num_components) {
  library(rTensor)
  result <- tryCatch(cp(tnsr, num_components = num_components), error = function(e) NULL)
  if (is.null(result)) {
    return(list(rank = num_components, mse = NA_real_, converged = FALSE, result = NULL))
  }
  gene_loadings <- result$U[[1]]; rownames(gene_loadings) <- dimnames(tnsr@data)$gene
  subject_loadings <- result$U[[2]]; rownames(subject_loadings) <- dimnames(tnsr@data)$subject
  time_loadings <- result$U[[3]]; rownames(time_loadings) <- dimnames(tnsr@data)$time
  colnames(gene_loadings) <- colnames(subject_loadings) <- colnames(time_loadings) <-
    paste0("Component_", seq_len(num_components))

  list(rank = num_components, mse = 1 - result$norm_percent / 100, converged = result$conv,
       loadings = gene_loadings, scores = subject_loadings, time_loadings = time_loadings,
       # lambdas -- the per-component scaling; U's columns are unit-norm, so
       # without this the fitted tensor can't be reconstructed and the
       # loadings have no real relative magnitude
       lambdas = result$lambdas)
}

cp_registry <- list(
  needs_tensor = TRUE,
  requires_subject_timepoint = TRUE,
  global_object = "tnsr",
  jobname = "cp_grid",
  fn = run_cp_job,
  pkgs = "rTensor",
  defaults = list(num_components = c(3, 4, 5, 6, 8, 10)),
  build_grid = function(p) expand.grid(num_components = p$num_components, stringsAsFactors = FALSE)
)

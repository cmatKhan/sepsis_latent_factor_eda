# Tucker (rTensor::tucker()) method registry -- a genuinely different
# METHOD from CP for every purpose in this framework, despite sharing
# R/lib/tensors.R's build_tensor() and the rTensor package -- see
# R/methods/cp.R's header for why.
#
# Deterministic given its ranks (no seed/init argument).
#
# Unlike CP, Tucker allows a DIFFERENT rank per mode -- rTensor's real
# argument is `ranks = c(r1, r2, r3)`, a length-3 vector assembled here
# from three independent config keys (`rank_genes`/`rank_subjects`/
# `rank_time`, each independently sweepable, crossed via `expand.grid()`
# like any other set of arguments -- these three don't need to be coupled
# to each other, unlike sPCA's K/para).
#
# Produces THREE factor matrices with their OWN ranks: gene loadings
# (rank_genes columns), SUBJECT-mode "scores" (rank_subjects columns,
# rownames = subject ids -- NOT sample ids, same caveat as CP), and
# time-mode loadings (rank_time columns, rownames = timepoint levels),
# plus the (small) core tensor -- kept as an artifact, not surfaced
# anywhere yet. `mse` is `1 - norm_percent/100`.
#
# Requires `dataset.subject_id_col`/`dataset.timepoint_col` -- validated
# in R/lib/metadata.R via `requires_subject_timepoint` below. Operates on
# the global object `tnsr` (shared with CP -- see R/create_slurm_bundle.R).

run_tucker_job <- function(rank_genes, rank_subjects, rank_time) {
  library(rTensor)
  ranks <- c(rank_genes, rank_subjects, rank_time)
  result <- tryCatch(tucker(tnsr, ranks = ranks), error = function(e) NULL)
  if (is.null(result)) {
    return(list(rank_genes = rank_genes, rank_subjects = rank_subjects, rank_time = rank_time,
                mse = NA_real_, converged = FALSE, result = NULL))
  }
  gene_loadings <- result$U[[1]]; rownames(gene_loadings) <- dimnames(tnsr@data)$gene
  subject_loadings <- result$U[[2]]; rownames(subject_loadings) <- dimnames(tnsr@data)$subject
  time_loadings <- result$U[[3]]; rownames(time_loadings) <- dimnames(tnsr@data)$time
  colnames(gene_loadings) <- paste0("Component_", seq_len(rank_genes))
  colnames(subject_loadings) <- paste0("Component_", seq_len(rank_subjects))
  colnames(time_loadings) <- paste0("Component_", seq_len(rank_time))

  list(rank_genes = rank_genes, rank_subjects = rank_subjects, rank_time = rank_time,
       mse = 1 - result$norm_percent / 100, converged = result$conv,
       loadings = gene_loadings, scores = subject_loadings, time_loadings = time_loadings,
       # core -- the one Tucker-specific object with no CP analogue.
       # all_resids -- same iteration-by-iteration residual-norm trace as
       # CP's own (?tucker's own recommended convergence check).
       core = result$Z, all_resids = result$all_resids)
}

tucker_registry <- list(
  needs_tensor = TRUE,
  requires_subject_timepoint = TRUE,
  global_object = "tnsr",
  jobname = "tucker_grid",
  fn = run_tucker_job,
  pkgs = "rTensor",
  defaults = list(rank_genes = c(10, 20, 30), rank_subjects = c(3, 4), rank_time = 3),
  build_grid = function(p) {
    expand.grid(rank_genes = p$rank_genes, rank_subjects = p$rank_subjects,
                rank_time = p$rank_time, stringsAsFactors = FALSE)
  }
)

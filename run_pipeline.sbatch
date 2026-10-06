#!/bin/bash
#SBATCH --job-name=targets_main
#SBATCH --output=logs/targets_main_%j.out
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=7-00:00:00

# Run the targets pipeline's main process (`tar_make()`) as a SLURM job,
# for runs that outlive an interactive RStudio Server session (e.g. CoGAPS
# fits, 12 h each). The main process only coordinates: it submits the
# worker jobs through crew.cluster and stays alive until they finish, so
# it needs little CPU/memory but a long wall time.
#
# usage, from the project root:
#   sbatch run_pipeline.sbatch                         # tar_make() everything outdated
#   sbatch run_pipeline.sbatch 'starts_with("fit_pca_")'   # a subset (tidyselect)
#
# Follow progress with `tail -f logs/targets_main_<jobid>.out`, or from R
# with targets::tar_progress_summary() / tar_poll().

set -euo pipefail

# ------ Get R 4.6 on PATH ------
# HTCF: the same spack environment rstudio_spack.sbatch launches RStudio from.
SPACK_ENV=/ref/mblab/software/chasem/spack_envs/rstudio-4.6.1
eval "$(spack env activate --sh "$SPACK_ENV")"
export LD_LIBRARY_PATH="$SPACK_ENV/view/lib:$SPACK_ENV/view/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
unset R_LIBS R_LIBS_USER R_LIBS_SITE
# Without spack, replace the four lines above with whatever provides R 4.6,
# e.g.:
#   module load R/4.6.1
# and set `cluster.setup_lines` in config/pipeline.yml to the same thing so
# the workers match. Run renv::restore() once on the new cluster first.

# The project root: renv activates from its .Rprofile, and _targets.R and
# _targets/ live here.
cd "${SLURM_SUBMIT_DIR}"

if [ $# -ge 1 ]; then
  Rscript -e "targets::tar_make(names = $1)"
else
  Rscript -e 'targets::tar_make()'
fi

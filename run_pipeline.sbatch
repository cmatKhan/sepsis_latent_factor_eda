#!/bin/bash
#SBATCH --job-name=targets_main
#SBATCH --output=logs/targets_main_%j.out
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=7-00:00:00

# Run the targets pipeline's main process (`tar_make()`) as a SLURM job,
# for runs that outlive an interactive RStudio Server session (e.g. CoGAPS
# fits, 12 h each). The main process only coordinates: it submits the
# worker jobs through crew.cluster and stays alive until they finish, so
# it needs little CPU but a long wall time. Its memory grows with the
# number of tasks pushed to crew and not yet collected (~10.6 MB each; a
# full run's opening burst of ~1,400 tasks reached 14.75 GB and crashed a
# 16G request, 2026-10-14); a background loop records the job's memory per
# process and the number of worker jobs, every 30 s, to
# logs/targets_main_<jobid>_mem.tsv.
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

# ------ Memory log ------
# One row per process in this job's cgroup every 30 s: time, cgroup total
# (GB), worker jobs of this user in SLURM (pending + running), pid, RSS
# (MB), command.
MEM_LOG="logs/targets_main_${SLURM_JOB_ID}_mem.tsv"
CGROUP="/sys/fs/cgroup$(awk -F: '$1 == "0" {print $3}' /proc/self/cgroup)"
echo -e "time\tcgroup_gb\tworker_jobs\tpid\trss_mb\tcommand" > "$MEM_LOG"
(
  set +eo pipefail  # a process exiting between listing and ps must not stop the log
  while true; do
    now=$(date +%FT%T)
    total=$(awk '{printf "%.2f", $1 / 1024^3}' "$CGROUP/memory.current" 2>/dev/null || echo NA)
    workers=$(squeue -u "$USER" -h -o %j 2>/dev/null | grep -c '^crew-worker' || true)
    for pid in $(cat "$CGROUP/cgroup.procs" 2>/dev/null); do
      ps -o rss=,args= -p "$pid" 2>/dev/null |
        awk -v t="$now" -v g="$total" -v w="$workers" -v p="$pid" \
          '{rss = $1; $1 = ""; printf "%s\t%s\t%s\t%s\t%.0f\t%s\n", t, g, w, p, rss / 1024, substr($0, 2, 120)}'
    done >> "$MEM_LOG"
    sleep 30
  done
) &
MEM_LOGGER=$!
trap 'kill $MEM_LOGGER 2>/dev/null' EXIT

if [ $# -ge 1 ]; then
  Rscript -e "targets::tar_make(names = $1)"
else
  Rscript -e 'targets::tar_make()'
fi

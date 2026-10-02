#!/bin/sh
# Cheap, no-SLURM-job-submitted-yet sanity check, split into two parts:
# (1) does THIS session (the one launching the POC) have what it needs,
# and (2) does the exact `spack load` + `R_LIBS_USER` sequence the
# generated SLURM job script will run actually resolve the needed
# packages, in a FRESH shell (so a stale assumption from this session's
# own environment doesn't hide a problem that would only show up once a
# real job runs).
#
# Expected to be run from RStudio Server's Terminal tab (View > Terminal,
# or the "Terminal" tab next to "Console") in the same session you'll run
# 01_raw_crew_poc.R from -- that way "this session" below really is the
# one crew_controller_slurm() will call sbatch from.
#
# See README.md for what to do based on the results.

set -eu

# ---- EDIT THIS (kept in sync with 01_raw_crew_poc.R) ----
USER_LIB_PATH="/ref/mblab/software/chasem/R-rstudio/4.6"
R_SPEC="r@4.6.1"
# --------------------

echo "=== 1. This session (where you'll launch the POC from) ==="
echo
echo -n "sbatch on PATH: "
if command -v sbatch >/dev/null 2>&1; then
  command -v sbatch
else
  echo "NOT FOUND -- crew_controller_slurm() calls sbatch directly from this"
  echo "  session, so the POC can't submit anything until this resolves."
fi
echo
echo "R packages needed in THIS session (crew/crew.cluster for step 3; also"
echo "targets/igraph for step 4), with $USER_LIB_PATH prepended to .libPaths():"
Rscript -e '
  .libPaths(c("'"$USER_LIB_PATH"'", .libPaths()))
  for (pkg in c("crew", "crew.cluster", "mirai", "nanonext", "targets", "igraph", "here")) {
    ok <- requireNamespace(pkg, quietly = TRUE)
    cat(sprintf("  %-14s %s\n", pkg, ok))
  }
'

echo
echo "=== 2. What a real SLURM job's script_lines will actually run ==="
echo
echo "Simulating, in a FRESH shell (new process, none of this session's state),"
echo "exactly what crew.cluster's generated job script will run before its own"
echo "trailing 'Rscript -e ...' line:"
echo
echo "  eval \$(spack load --sh $R_SPEC)"
echo "  export R_LIBS_USER=\"$USER_LIB_PATH\""
echo

sh -c "
  eval \$(spack load --sh $R_SPEC)
  export R_LIBS_USER=\"$USER_LIB_PATH\"
  echo 'Rscript resolved to:' \$(command -v Rscript || echo NOT FOUND)
  Rscript -e '
    cat(\"R version: \", R.version.string, \"\n\", sep = \"\")
    cat(\".libPaths(): \", paste(.libPaths(), collapse = \", \"), \"\n\", sep = \"\")
    for (pkg in c(\"crew\", \"crew.cluster\", \"mirai\", \"nanonext\", \"targets\", \"igraph\")) {
      ok <- requireNamespace(pkg, quietly = TRUE)
      cat(sprintf(\"  %-14s %s\n\", pkg, ok))
    }
  '
"

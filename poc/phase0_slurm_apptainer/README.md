# Phase 0 POC: does `targets` + `crew.cluster` work on the real cluster?

This directory is a throwaway proof-of-concept, not part of the pipeline. It exists to
answer one question before any real migration work starts: **can `crew.cluster` launch
and manage SLURM worker jobs for a `targets` pipeline, end to end, on the real cluster?**

No Apptainer/containers are involved -- an earlier version of this POC wrapped workers in
our per-method `.sif` images, but after reconfiguring the RStudio launch (via spack) to
have `sbatch` directly on `PATH`, there's no container boundary to cross at all: workers
just get R via `spack load r@4.6.1` (the exact same `script_lines` mechanism
`crew_options_slurm()`'s own docs show for `module load R`), the same way the rest of this
cluster's spack-based R setup already works. This is simpler and a closer match to how
`targets`+`crew.cluster` are actually meant to be used.

## Expected workflow: run this from your (spack-configured) RStudio Server session on the cluster

`crew_controller_slurm()` calls `sbatch` directly, from wherever the R process that
creates it runs -- there's no separate "submission step" to run elsewhere. That means:

- Whatever session you run these scripts from needs `sbatch` on its `PATH` -- true now
  that your RStudio launch is reconfigured via spack.
- That session's R process has to stay alive for as long as the run takes, since compute
  nodes connect back to it over the network to report results (this is `mirai`'s
  host/worker model). For this short POC (a couple minutes) that's a non-issue either
  way; it only matters once Phase 1 runs the real pipeline for hours.
- All three files here use `here::here(...)` (via
  `sepsis_timecourse_factors_eda.Rproj`) for every path, so it doesn't matter what
  directory RStudio happens to open in -- just have this project open.
- Workers run on bare compute nodes (no container), getting R via
  `eval $(spack load --sh r@4.6.1)` and an explicit `R_LIBS_USER` export in
  `script_lines` -- a fresh SLURM job's shell starts with none of your interactive
  session's state, so this is set explicitly rather than assumed.

## What you'll need to edit before running

Each of `00_preflight_check.sh`, `01_raw_crew_poc.R`, and `02_targets_poc/_targets.R` has
an `# ---- EDIT THESE ----` block near the top (keep them in sync if you change one):

- `user_lib_path` / `USER_LIB_PATH` -- defaults to `/ref/mblab/software/chasem/R-rstudio/4.6`,
  where `crew`/`crew.cluster`/`targets`/etc. are installed. Used both for this session's
  own `.libPaths()` and exported as `R_LIBS_USER` for workers.
- `r_spec` / `R_SPEC` -- defaults to `r@4.6.1`, the spack spec `spack load --sh` resolves
  on a worker's job script.
- `partition` -- left `NULL` (today's rslurm templates don't set one either). Set it if
  your cluster needs one.

## Steps

1. **Preflight check** (no SLURM job submitted yet). From RStudio's **Terminal** tab (not
   the Console -- this is a shell script):

   ```sh
   bash poc/phase0_slurm_apptainer/00_preflight_check.sh
   ```

   Part 1 (this session) should show `sbatch` resolving to a real path and all of
   `crew`/`crew.cluster`/`mirai`/`nanonext`/`targets`/`igraph`/`here` as `TRUE`.

   Part 2 runs the *exact* `eval $(spack load --sh r@4.6.1)` / `export R_LIBS_USER=...`
   sequence that will appear in a real generated job script, in a fresh shell (not this
   session's own environment), and checks the same packages there. This is the one that
   actually predicts whether step 3 will work -- if it shows packages missing while part 1
   shows them present, something about how spack/R_LIBS_USER resolve differs between an
   interactive session and a fresh job shell, and is worth tracking down before going
   further.

2. **Raw crew + SLURM round trip** (the actual test -- submits one real SLURM job via
   `crew_controller_slurm()`). From the RStudio **Console**, with this project open:

   ```r
   source("poc/phase0_slurm_apptainer/01_raw_crew_poc.R")
   ```

   It prints `sbatch`/network sanity info up front, then should print a `POC PASSED`
   line within ~1-2 minutes, confirming:
   - a real SLURM job was submitted and started on a compute node (checked via
     `Sys.info()[["nodename"]]`, which should differ from the submitting node's own),
   - `spack load r@4.6.1` + the `R_LIBS_USER` export correctly resolved `crew`/etc. on
     that worker,
   - a real `prcomp()` fit executed and its result came back over the network.

   If it prints `POC INCONCLUSIVE` (no result within 180s) instead of hanging forever,
   that most likely means the worker process on the compute node can't connect back to
   this session over the network, or the job never started -- check `squeue`, and the
   generated job scripts/logs under `poc/phase0_slurm_apptainer/job_scripts/`. A hang or
   inconclusive result is itself a useful (if less convenient) finding -- report it.

3. **Only if step 2 passes**, try the fuller `targets`-based version, from the Console:

   ```r
   setwd(here::here("poc/phase0_slurm_apptainer/02_targets_poc"))
   targets::tar_make()
   ```

   This drives `R/methods/pca.R`'s *real, unmodified* `pca_registry` through dynamic
   branching (`pattern = map(...)`), the same construct that would replace today's
   `expand.grid()`-based grid building. `targets::tar_read(pca_summary)` afterward will
   print results if it succeeds.

## What to report back

Whichever step you get to, paste back:
- which step you reached,
- pass / fail / inconclusive,
- the console output (especially any error message).

That's enough to know whether to proceed to Phase 1 as planned, adjust the approach, or
dig into something specific.

# TODO

1. **Remove the hard-coded requirement to use spack, and the machine-specific paths in the
   code.** Nothing here should assume HTCF: the user gets R 4.6.1 running themselves
   (`docs/setup.qmd` already starts from that), and the code finds R, the data and the DB
   without a path baked in. Where they are today:

   - **spack**
     - `R/targets/controllers.R`: `worker_script_lines()` builds `spack env activate` and the
       `LD_LIBRARY_PATH` export from `cluster$spack_env` whenever `cluster.setup_lines` is
       null, so spack is the default; the roxygen comment names `rstudio_spack.sbatch`.
     - `config/pipeline.yml`: `cluster.spack_env: /ref/mblab/software/chasem/spack_envs/rstudio-4.6.1`
       and the comments around it.
     - `run_pipeline.sbatch`: lines 29-34 set `SPACK_ENV` and activate it before `tar_make()`.
     - `README.md` (lines 16-29, the spack environment and its activation commands) and
       `docs/workflow.qmd` ("Running on a scheduler other than SLURM", item 3) still describe it.
     - `CLAUDE.md`: the HTCF `srun`/spack command and the note about the two spack versions.
       That file is instructions for the AI assistant; decide whether the HTCF commands stay.
   - **Paths**
     - `config/pipeline.yml`: `data_root: /scratch/mblab/chasem/hf_sepsis_collection`. The
       committed config should not name a machine; `SEPSIS_DATA_ROOT` (see
       `default_data_root()`, `R/lib/data_paths.R`) or a local value should supply it, with a
       clear error when neither is set. `docs/setup.qmd` ("Raw data") quotes the HTCF path.
     - `R/preprocessing/*.R` (33 scripts): the interactive fallback inside each
       `if (!framework_run)` block sets `expression_path`, `sample_metadata_path` and
       `feature_metadata_path` to `~/projects/hf_sepsis_collection/<dataset>/...`. Use
       `default_data_root()` or stop with a message.
     - Re-run `grep -rn "/home/\|/scratch/\|/ref/\|~/projects" R config app run_pipeline.sbatch`
       when done; the list above is from that search.

   End state: R is assumed to be on `PATH`; `cluster.setup_lines` is optional (null means
   nothing runs before R starts except the `cd` into the project); there is no `spack_env`
   key; `run_pipeline.sbatch` assumes R on `PATH`; the paths come from the environment or a
   local config.

   Caution when making the change: the preprocessing scripts are on the matrix path, which
   invalidates fits (CLAUDE.md, "Don't refit by accident"), and `worker_script_lines()` sits
   beside the controller code. Run `targets::tar_outdated()` afterwards and confirm no
   `fit_` or `grid_` target is outdated before launching anything.

2. **Adjust controller resources in `config/pipeline.yml` (and `run_pipeline.sbatch`) to what
   the 2026-10-14/15 full runs used.** Measurements and their scaling with dataset size are in
   `docs/workflow.qmd` ("Expected resource use"); the scripts are
   `logs/dev/resource_analysis.R` and `logs/dev/resource_predict.R`. Target: ~75% of the
   request at the largest observed peak.

   | Controller | Now | Observed peak | Change to | What changing it reruns |
   |---|---|---|---|---|
   | `fgsea` | 32 GB, 8 CPUs | 3.9 GB; ~2.5 cores used (31% of 8) | 6 GB, 3 CPUs | memory: nothing. CPUs: every `fgsea_` target and what follows (enrichment, `db_`), because `cpus_per_task` is passed into the fgsea commands (`fgsea_cpus`, `R/targets/factories.R`); no fits |
   | `cogaps` | 18 GB | 10.3 GB | 14 GB | nothing (do **not** change `cpus_per_task`: it is CoGAPS's nSets and changes the grid, refitting CoGAPS) |
   | `wgcna` | 6 GB | 2.8 GB (bounded: the gene × gene matrices at the 8,000-gene cap) | 4 GB | nothing |
   | `ingest` | 8 GB | 2.3 GB | 3 GB | nothing |
   | `projectr` | 4 GB, `time_minutes: 120` | 2.1 GB; a branch took 112 min, and two workers were retried on the backup (a late-starting long branch overran the 2 h limit) | 3 GB, `time_minutes: 240` | nothing |
   | `ica` | 0.5 GB | 0.9 GB (12 retries on `ica_backup`) | 1 GB | nothing |
   | `pca` | 0.5 GB | 0.55 GB (1 retry on `pca_backup`) | 1 GB | nothing |
   | `spca`, `nmf`, `vp`, `prep` | 4 / 0.5 / 12 / 6 GB | 3.0 / 0.3 / 8.6 / 3.4 GB | keep | |
   | main process (`run_pipeline.sbatch --mem`) | 64 GB | 17.7 GB (2026-10-14 opening burst) | 24 GB | nothing (not a target) |

   Memory, `time_minutes` and `workers` are not hashed by targets; `batch_size` and
   `cpus_per_task` are (they reach the grids and some commands), so leave those alone unless
   the rerun is intended. Run `targets::tar_outdated()` after the edit and confirm only what the
   table says (nothing, or the fgsea chain) is outdated. Update the measurement comments in
   `pipeline.yml` (several still quote EARLI-only numbers from 2026-10-02/04) at the same time.

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

# CLAUDE.md

Guidance for working in this repository. Project documentation is the Quarto site in `docs/`
(setup, data, methods, database, workflow, app, function reference); read the relevant page
rather than relying on memory, and keep it current when behavior changes.

## Working style

- **Ask before hacking anything.** Prefer the idioms of the targets, crew and crew.cluster
  vignettes; if something needs a workaround, say so and ask first.
- **No backward compatibility.** Keep the codebase focused: delete legacy code rather than keep
  it runnable (git has the history). Correctness is shown by parity checks, not by preserving
  old paths.
- **Don't commit unless asked.** Work is on branch `target_rewrite`, not merged to `main` until
  signed off.
- **Explain unfamiliar methods descriptively** (what the method does, what the choice changes),
  and give a recommendation rather than a survey of options.
- **Never report work as complete while its jobs are still running.** Give status (what's done,
  running, queued), not a summary.
- **Report outcomes faithfully**, including when an earlier statement turned out wrong.

## Running things on HTCF

- R comes from the spack environment. For an ad-hoc compute-node command:
  ```sh
  srun --mem=8G --time=60 bash -c 'eval "$(/ref/mblab/software/spack/bin/spack env activate --sh /ref/mblab/software/chasem/spack_envs/rstudio-4.6.1)"; unset R_LIBS R_LIBS_USER R_LIBS_SITE; Rscript ...'
  ```
  Login-node `Rscript` also works for light reads (renv activates from `.Rprofile`).
- **Long runs go through `sbatch run_pipeline.sbatch`**, never a background tool command (those
  stop at 2 hours). The main `tar_make()` process must outlive the run: workers connect back to
  it, and killing it kills them.
- Background watchers also stop at 2 hours; re-arm them rather than assuming a run finished.
- **Scripts that compute nodes read must live on shared storage** (`logs/dev/`), not `/tmp`.
  Use the session scratchpad for local-only temp files.
- **Whole-DB scans** (row counts, integrity or foreign-key checks on the 3.6 GB DB) are too slow
  from the login node over the network filesystem: run them in an sbatch job that copies the DB
  to node-local disk first (see `logs/dev/verify_db.sbatch`).

## Jobs and resources

- **Log every job id** you submit (session scratchpad `jobs/main_jobs.txt`), with a note of what
  it is.
- Check usage with `seff <jobid>`, or for arrays
  `script -qc "stty cols 100; python /scratch/mblab/chasem/seff-array/seff-array.py <id>" /dev/null`
  (seff-array needs a TTY, and misparses ReqMem values like `512M`).
- **Targets:** about 75% average memory and CPU use, never below 50%, erring high. Size memory
  for the largest dataset (about 8,000 genes x 518 samples), not EARLI.
- **No worker caps just to limit submissions.** SLURM's scheduler and per-user limits (200
  running jobs) handle queuing. A `workers:` cap needs another reason, stated in a comment
  (pca/nmf/ica: tasks shorter than a worker's ~30 s start-up).
- Out-of-memory crashes are handled by each controller's `backup:` (crew backup controllers).

## Don't refit by accident

The fits take days (sPCA alone about a day across datasets). They are invalidated by any code
change to what the `fit_*` and `grid_*` targets depend on:

- each `<name>_registry` and its fit function (`R/methods/*.R`);
- `run_method_fits()`, `flatten_fits()`, `build_method_grid()`, `resolve_method_params()`,
  `add_batches()`, `GRID_BATCH_COLUMNS` (`R/targets/run_method_fit.R`), `merge_named_list()`
  (`R/lib/grids.R`), `capture_fit()` (`R/lib/method_registry.R`);
- the matrix path: `load_input_matrix()`, `run_preprocessing_script()` (including error-message
  strings inside them): every `mat_` target reruns, and the fits too unless each rebuilt matrix
  is bit-identical;
- a controller's `batch_size` or `cpus_per_task` (they change the grid).

Editing `read_dataset_metadata()` or a config's comments only reruns the cheap `meta_` targets:
their value is unchanged, so targets stops there. Deferred changes and what each would rerun are
in `TODOS.md`.

Rules:
- Ingestion-side changes go in each method's `<name>_ingest` or in `R/targets/`, never the
  registry.
- Comments and whitespace don't invalidate (targets hashes functions without them); roxygen
  edits are always safe.
- After any edit, check `targets::tar_outdated()` shows no `fit_`/`grid_` targets before
  launching anything.
- **Store more at fit time.** When a fit function does change, return the whole model object
  where cheap, so ingestion can change later without refitting.

## The database

- A derived output of the targets store, with no migrations. Changing `R/db/schema.sql` or
  `views.sql` means bumping `SCHEMA_VERSION` (`R/db/connect.R`), deleting the DB and its
  `stability_artifacts/`, and rerunning `tar_make()` (only ingestion reruns).
- Those two SQL files are a file target, so even a comment edit reruns every dataset's DB write
  (about an hour). Their header comments still point to the old `R/db/SCHEMA.md`; fix that with
  the next real schema change.
- Read-only references: `results/targets_v0/` (the targets DB before the redesign) and
  `/scratch/mblab/chasem/r_projects/results/stability.sqlite` (the rslurm-era DB).

## targets lessons from the rewrite

- **File targets drop names** (`format = "file"` returns bare paths); recover ids from the path
  (see `projectr_jobs()`).
- **targets can't branch over an empty table** -- it stops the whole `tar_make()`, not just one
  target. Job tables return `placeholder_jobs()` and branches check `has_job()`.
- **`c(list(a, b))` doesn't flatten** a list of patterns' branch lists; use
  `unlist(x, recursive = FALSE)`.
- **Each branch costs the main process about a second**, so batch tiny tasks (`batch_size`;
  projectR grouped per target dataset).
- Light targets use `deployment = "main"`; storage and retrieval are on workers.
- **Dry-run ingestion changes on one dataset** before a full run (`logs/dev/test_ingest_earli.R`
  calls the ingestion functions directly into a scratch DB). Include the combining steps
  (enrichment rows, projections): the EARLI dry run skipped them and missed a bug that wrote
  zero enrichment rows.

## Scientific decisions to keep

- **sPCA:** K 2-12; `para_coarse` 0.02-0.5 plus IS refinement (Guerra-Urzola et al. 2021);
  `lambda` a single value (1e-6); `max_hours` 8, after which a fit is recorded as
  `killed_not_converged`.
- **WGCNA:** one block over all genes (`maxBlockSize = 50000`); power swept, not auto-picked.
- **PCA:** one fit per dataset at the maximum rank; first-n slicing happens downstream.
- **ICA:** seed sweep plus bootstrap resamples, for ICASSO.
- **Stability storage:** matched pairs plus context (runner-up, margin, unmatched-cell null
  stats, histograms), not every factor x factor cell. WGCNA kME: own module, runner-up, analytic
  and permutation-null p-values.
- **Keep everything computed** in the DB, even what the app doesn't read yet.

## Open items

See `TODOS.md` for deferred changes and their rerun costs. In short:

- **App:** update `app/` to the new DB schema (its queries target the old layout; its comments
  still reference removed files).
- **gains:** its preprocessing script is a stub; pooling the four accessions
  (E-MTAB-4421/4451/5273/5274) needs the user's decision.
- **ROSE ×3:** data isn't on the cluster. **MARS, SHIP-TREND:** preprocessing scripts don't
  exist.
- `dev/phase1_parity.R` and `dev/phase2_parity.R` compare against the old schema and are now
  obsolete.
- Nothing from the targets rewrite is committed yet.

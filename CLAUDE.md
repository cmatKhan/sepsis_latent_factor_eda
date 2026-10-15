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
- **Docs point to code, never the reverse.** No doc paths (`docs/*.qmd`, page names) in code
  strings or comments, in `R/` or `config/`: each docs section names the scripts and functions
  behind it instead. Code strings (error messages) are hashed by targets, so a doc pointer there
  forces reruns when the docs move; comments aren't hashed, but follow the same rule.

## Running things on HTCF

- R comes from the spack environment. For an ad-hoc compute-node command:
  ```sh
  srun --mem=8G --time=60 bash -c 'eval "$(spack env activate --sh /ref/mblab/software/chasem/spack_envs/rstudio-4.6.1)"; unset R_LIBS R_LIBS_USER R_LIBS_SITE; Rscript ...'
  ```
  (`spack` on PATH is spack 1.1; `/ref/mblab/software/spack/bin/spack` is 0.17 and fails.)
  Login-node `Rscript` also works for light reads (renv activates from `.Rprofile`).
- **Long runs go through `sbatch run_pipeline.sbatch`**, never a background tool command (those
  stop at 2 hours). The main `tar_make()` process must outlive the run: workers connect back to
  it, and killing it kills them.
- Background watchers also stop at 2 hours; re-arm them rather than assuming a run finished.
- **Scripts that compute nodes read must live on shared storage** (`logs/dev/`), not `/tmp`.
  Use the session scratchpad for local-only temp files.
- **Whole-DB scans** (row counts, integrity or foreign-key checks on the multi-GB DB) are too
  slow from the login node over the network filesystem: run them in an sbatch job that copies
  the DB to node-local disk first.
- **Run a subset by name** (`sbatch run_pipeline.sbatch '<tidyselect>'`): targets builds only
  those and their dependencies. Check `tar_outdated(names = ...)` with the same selection first.

## Jobs and resources

- **Log every job id** you submit (`logs/jobs/main_jobs.txt`, shared storage and git-ignored; not
  the session scratchpad, which lives under `/tmp` and was wiped mid-session once), with a note of what
  it is.
- Check usage with `seff <jobid>`, or for arrays
  `script -qc "stty cols 100; python /scratch/mblab/chasem/seff-array/seff-array.py <id>" /dev/null`
  (seff-array needs a TTY, and misparses ReqMem values like `512M`).
- **Targets:** about 75% average memory and CPU use, never below 50%, erring high. Size memory
  for the largest dataset (SHIP-TREND, 8,000 genes x 991 samples), not EARLI.
- **The `vp` controller's tasks record their job's cgroup peak memory and CPU seconds** in
  `$resources` (`job_resources()`). `seff` missed the earlier BiocParallel `SnowParam` workers
  (they left the process tree); the workers are now forked (`MulticoreParam`, because
  SnowParam's ports collided between tasks starting together on a node), so check whether
  `seff` agrees with `$resources`.
- **No worker caps just to limit submissions.** SLURM's scheduler and per-user limits (200
  running jobs) handle queuing. A `workers:` cap needs another reason, stated in a comment
  (pca/nmf/ica: tasks shorter than a worker's ~30 s start-up).
- Out-of-memory crashes are handled by each controller's `backup:` (crew backup controllers).

## Don't refit by accident

The fits take days (sPCA alone about a day across datasets). They are invalidated by any code
change to what the `fit_*` and `grid_*` targets depend on:

- each `<name>_registry` and its fit function (`R/methods/*.R`);
- `run_method_fits()`, `fit_subjects()`, `flatten_fits()`, `build_method_grid()`,
  `build_refine_grid()`, `select_method_fits()`, `resolve_method_params()`, `add_batches()`,
  `GRID_BATCH_COLUMNS` (`R/targets/run_method_fit.R`), `merge_named_list()` (`R/lib/grids.R`),
  `capture_fit()` (`R/lib/method_registry.R`);
- the `subjects_<dataset>` and `sample_meta_<dataset>` targets (`read_sample_metadata()`,
  `read_parquet_plain()`): every fit of that dataset reads them;
- the matrix path: `load_input_matrix()`, `run_preprocessing_script()` (including error-message
  strings inside them): every `mat_` target reruns, and the fits too unless each rebuilt matrix
  is bit-identical;
- a controller's `batch_size` or `cpus_per_task` (they change the grid).

Editing `read_dataset_metadata()` or a config's comments only reruns the cheap `meta_` targets:
their value is unchanged, so targets stops there. The repeated-measures module (`R/lib/design.R`,
`R/lib/expression.R`, `R/targets/variance*.R`, `composition.R`) never touches a fit.

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
  (about an hour).
- Schema 2 (2026-10-13): repeated-measures tables, enrichment gene lists and pattern markers as
  artifacts. The earlier (schema 1) results are on the user's local machine, not the cluster.
- Schema 3 (2026-10-14): `projections` gains centred scores (`transfer_r2`, `native_r2`,
  `relative_r2`, `null_r2`, `n_components`; `projection_scores()`, `projection_scores_<dataset>`
  targets). projectR's own R² (inflated by the mean profile) is no longer stored; the app reads
  `relative_r2`.

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
- **Dry-run ingestion changes on one dataset** before a full run (`logs/dev/test_ingest.R`
  calls the ingestion functions directly into a scratch DB). Include the combining steps
  (enrichment rows, projections): an earlier EARLI dry run skipped them and missed a bug that
  wrote zero enrichment rows.
- **Load a package before reading its S4 objects.** `dim()` on a variancePartition dream fit
  (`MArrayLM2`) read before the package's namespace is loaded returns `NULL`, and stays `NULL`
  for the session; give targets that read such objects `packages = "variancePartition"`.
- **Skipping a method (`skip_methods:`) and running to the end deletes its artifact files**
  (`db_cross_dataset`'s orphan cleanup: the DB no longer references them), but its `ingest_`
  targets stay up to date, because ingestion writes artifacts as a side effect, not as file
  targets. Before adding the method back, `tar_invalidate()` its `ingest_`, `pairs_`,
  `representatives_`, `redundancy_`, `fgsea_jobs_`, `fgsea_` and `factor_vp_` targets for every
  dataset that had been ingested (2026-10-14: hfgp-500fg's sPCA). Its fits survive.
- **Don't store once-per-session messages in a target value.** A captured lme4 deprecation
  notice made design values depend on build order and reran everything downstream.

## Scientific decisions to keep

- **sPCA:** K 2-12; a scan at `K_scan` 5 over para 0.2-0.5, then every K at the scan's
  Index-of-Sparseness winner +/- 0.05 (Guerra-Urzola et al. 2021). On the earlier full grid every
  winner was at para 0.25-0.50 and moved by at most 0.10 across K within a dataset; the replay of
  this design was within 0.05 of every K's winner at >= 80% of Ks in 30/30 datasets for about
  half the compute. `lambda` a single value (1e-6); `max_hours` 8 (then `killed_not_converged`);
  convergence from spca()'s own trace (`iterations`, `converged`).
- **WGCNA:** one block over all genes (`maxBlockSize = 50000`); power swept, not auto-picked.
- **PCA:** one fit per dataset at the maximum rank; first-n slicing happens downstream.
- **ICA:** seed sweep plus bootstrap resamples **of subjects** (all of a drawn person's samples),
  for ICASSO.
- **WGCNA storage:** per-gene connectivity and the scale-free fit are stored per fit; the TOM is
  not (recomputable from the matrix and parameters).
- **Repeated measures** (2026-10-13 findings, `docs/methods.qmd`): time is 1-4% of variance in
  the top 8,000 genes, between-person 35-60%; leukocyte composition is a stable person trait;
  within-person change is mostly unstructured. So: fit the full datasets and all per-timepoint
  subsets as they are; resample by subject; treat within-family agreement as temporal
  reproducibility, not replication; describe every factor with a score-level variance
  partition (subsets use the parent's design without subject/time). No within-person-centred
  fits.
- **Variance models:** variancePartition with ML, dream with REML (`ddf = "adaptive"`);
  categorical variables random, numeric covariates fixed and standardized; arrays use the
  deposited log values (no variance filter); RNA-seq filterByExpr + TMM + voomWithDreamWeights.
  No between-array normalization by default: composition shifts are real biology. ANEMONES'
  spike-in controls enter as a covariate (their first PC), because a per-sample control offset
  over-corrected.
- **Data quality is flagged in the app, not fixed in the data:** GSE54514's annotations don't
  match its arrays (sex agreement 0.50) but it stays in; sex-discordant samples stay.
- **Stability storage:** matched pairs plus context (runner-up, margin, unmatched-cell null
  stats, histograms), not every factor x factor cell. WGCNA kME: own module, runner-up, analytic
  and permutation-null p-values.
- **Keep everything computed** in the DB, even what the app doesn't read yet.

## Open items

- **App:** not yet updated to schema 2 (repeated-measures tables; enrichment genes and pattern
  markers now artifacts; the stale comments in `app/app.R` and `app/R/db_helpers.R` naming
  removed files). Its Time structure page should switch to the stored variance results, and the
  flags for GSE54514 and sex-discordant samples belong there. The user develops the app locally.
- **gains:** its preprocessing script is a stub; pooling the four accessions
  (E-MTAB-4421/4451/5273/5274) needs the user's decision. Adding it reruns every projection
  (a new projectR target dataset).
- **Full re-run** pending (schema 2, new sPCA grid, subject bootstrap, MARS/SHIP-TREND/ROSE
  subsets); afterwards tune the controllers from it (sPCA on the new grid, every method on
  SHIP-TREND and MARS, `vp` from the cgroup records).
- Nothing from the targets rewrite is committed yet.

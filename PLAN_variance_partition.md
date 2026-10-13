# variancePartition + dream module: implementation plan

Written to be taken over on the cluster (SLURM), by you or a fresh Claude session with no memory of
the conversation that produced it. Everything the implementer needs is here: the decisions, the
facts established locally, the prototype's measured numbers, the files to write, the targets, the
output to save, the tests, and what to bring back. Nothing in this plan has been implemented in the
repository yet. This file is untracked: commit or copy it to the cluster (the repo rules say not to
commit unless asked, so that is your call).

## 1. Goal and decisions

**Goal.** Replace `R/de/` (limma + `duplicateCorrelation()`, hand-run, outside the pipeline) with a
targets module built on **variancePartition** and its **dream** function, following the package's
vignettes. The module computes the full result once and saves it; the ingest stage turns it into DB
tables and artifacts without rerunning any model. Results guide the rest of the work (how to handle
timepoints, what the factor analysis should do), so getting data out first is the priority.

**Why.**
- The Time structure page (app) showed most time-course variance is between subjects. What matters
  for trajectories is how variance splits across subject, time, group and group x time.
  variancePartition estimates that per gene with a mixed model; our own method-of-moments estimator
  agrees with it in simulation but is not a published method.
- `duplicateCorrelation()` assumes one subject variance for all genes. The between-subject share
  varies a lot across genes (ANEMONES: the middle 80% of genes spans about 0.34-0.70). dream
  estimates it per gene and shares the design with the variance partition.

**Decisions made with the user.**
- Cover the 12 datasets the old scripts ran (table in section 4). ROSE and MARS are in the pipeline
  only after the separate "add datasets" work (section 11); their `design:` blocks are written now.
- Drop the three follow-on analyses that lived in `R/de/` (fgsea of contrasts, clustering of
  time-omnibus genes, DE-vs-factor correlation). The saved dream output keeps every statistic, so
  they can be rebuilt later at ingest.
- First-pass VP formulas: subject, time, group, group:time, plus sex and age where present. Outcome,
  treatment, severity and cell-proportion columns are left out on purpose (they can mediate the time
  effect). The dream (DE) model gets no covariates in this first pass (as before).
- Do exactly what is best for variancePartition, not what is convenient for the factor analysis:
  no variance filter on arrays, own normalization path per platform, genes may differ from the
  factor analysis' top-8,000.
- No backward compatibility: delete `R/de/` and its documentation; document the current state only.
- Don't commit unless asked.

## 2. Prerequisites (do these first)

**A. The pipeline currently cannot load.** `_targets.R` line 33 sources `app/R/metadata_helpers.R`,
which was deleted during the app port. `generic_association_scan()` (used by
`wgcna_dataset_diagnostics()`, `R/targets/diagnostics.R` line 115) lived there. Recover it:
`git show HEAD:app/R/metadata_helpers.R` (HEAD still has the file), put `generic_association_scan()`
(and only what it needs) in `R/lib/association.R`, and replace the `app/R/metadata_helpers.R` entry in
`tar_source()` with `R/lib/association.R`. It needs `WGCNA::corPvalueStudent`. `targets` was not
installed on the machine where this was written, so how `tar_source` fails on a missing path was not
observed; confirm the fix with `targets::tar_manifest()`.

**B. Packages.** variancePartition 1.42.0 is installed in the project library on the machine where
this was written (it is not in `renv.lock`). Its dependencies present there: `lme4`, `lmerTest`,
`pbkrtest`, `remaCor`, `BiocParallel`, `edgeR`, `limma`. On the cluster:
`renv::restore()` then `renv::install("bioc::variancePartition")`, then record only the new packages
(`renv::snapshot(packages = c("variancePartition", "lme4", "lmerTest", "pbkrtest", "remaCor", "edgeR", "limma"), update = FALSE)`
or equivalent; check the lockfile diff). Do not take a blanket snapshot blindly: renv 1.3.1 was
loaded against a lockfile configured for 1.3.0 here, and a snapshot would change that line.

**C. Where the vignettes are** (read them; the plan follows them): in the installed package,
`system.file("doc", package = "variancePartition")`: `variancePartition.Rmd` (main tutorial, incl.
interpretation, "Which variables should be included?", RNA-seq processing, interactions),
`dream.Rmd` (DE, contrasts, `z.std`, KR vs Satterthwaite, comparison with `duplicateCorrelation`),
`rnd_effects.Rmd` (why ML for fractions, REML for dream, why categoricals are random),
`errors.Rmd` (every warning/error and its fix), `FAQ.Rmd`, `additional_visualization.Rmd`.

## 3. Rules taken from the vignettes, and how each is applied

| Vignette guidance | Applied as |
|---|---|
| RNA-seq: `DGEList`, expression filter, TMM, `voom`; for dream, `voomWithDreamWeights(dge, form, metadata)`; variancePartition uses voom precision weights automatically | counts -> `edgeR::DGEList` -> `filterByExpr()` -> TMM (`normLibSizes`, formerly `calcNormFactors`) -> `voomWithDreamWeights(dge, dream_formula, data)`; that one voom object feeds both `dream()` and `fitVarPartModel()` |
| Arrays: any standard normalization | the raw log-scale matrix (`pivot_expression_long()`), one probe per gene symbol (highest mean probe), no variance filter |
| Continuous variables fixed, categorical random; intercept required; fixed-effect categorical with many levels is unstable | numeric covariates (`age`) fixed; everything else `(1\|x)` |
| Interactions: include the constituents | `(1\|group:time)` with `(1\|group)`, `(1\|time)` |
| Variables constant within a person (sex, group) lower the apparent share of Individual; the person-level total is their sum | terms reported separately; "between-subject total" derived at ingest, not modelled |
| Fractions with ML (`REML = FALSE`); dream uses REML, Kenward-Roger for <= 20 samples else Satterthwaite (`ddf = "adaptive"`) | package defaults for each |
| `fitVarPartModel()` keeps all lmer fits (about 10 GB for 20k genes x 1000 samples); use `fxn` | `fxn = vp_record` (section 6); no lmer object stored |
| `canCorPairs()`, `colinearityScore()`; warnings in `errors.Rmd` | stored per dataset / per gene; gene-level errors saved |
| Do not read meaning into the residual fraction | not used as a metric anywhere |
| dream: tested variables fixed; contrasts via `makeContrastsDream()`; `dream()` -> `eBayes()` -> `topTable()`; report `z.std` / `F.std`; BiocParallel `SnowParam` | `~ 0 + cell + (1\|subject)`; `topTable(fit, coef = ...)` at ingest; `SnowParam(cpus, "SOCK")` |
| A VP run says which variables to add to the dream model | dream formula gets covariates via config only after the first VP results |

## 4. Design spec in each dataset config (`design:` block)

Add an optional `design:` block to `config/<ID>_config.yml` (the pipeline reads it with
`read_dataset_metadata()`; validate in `R/lib/metadata.R`). Fields: `platform` (`array`/`rnaseq`);
`subject` (omit for cross-sectional); `time`; `time_order` (explicit, chronological; required whenever
alphabetical is not chronological); `group`; `covariates` (list; non-numeric ones are random effects);
`numeric` (covariates fitted as numbers, coerced with `as.numeric`); `drop`/`keep`
(`{column: [values]}` sample filters); `set_na` (`{column: [values]}` recoded to NA); `contrasts`
(optional override `{name: "levelA - levelB"}`; default is generated, section 6). "NA", "" and NA are
treated as missing; samples with missing model variables are dropped once, consistently, and counted.

Metadata checked locally (types, levels, missing). Write these blocks:

```yaml
# ANEMONES (array)   disease: Control(30, Day 1 only)/Sepsis/SIRS; sex F/M; time levels sort chronologically
design: {platform: array, subject: patient_id, time: timepoint, time_order: ["Day: 1","Day: 2","Day: 5","Day: Discharge"], group: disease, covariates: [sex]}
# CORTICUS (array)   "Post(24h)" sorts before "Pre": time_order is REQUIRED; age is text -> numeric; treatment/responder/survive_28_days exist but are left out
design: {platform: array, subject: patient_id, time: timepoint, time_order: [Pre, "Post(24h)"], covariates: [sex, age], numeric: [age]}
# GSE13904 (array)   clinical_status: Control, Sepsis, Septic Shock, SIRS, SIRS resolved; no sex/age columns
design: {platform: array, subject: patient_id, time: timepoint, time_order: [day1, day3], group: clinical_status}
# GSE54514 (array)   day is integer 1-5 (treat as categorical); disease_status healthy / sepsis survivor / sepsis nonsurvivor; group_id is HC_1... (a person id); age numeric
design: {platform: array, subject: group_id, time: day, time_order: [1,2,3,4,5], group: disease_status, covariates: [sex, age], numeric: [age]}
# GSE95233 (array)   22 samples have timepoint "NA" (drop); sex has a level H (4 samples): set to NA unless you decide otherwise; age is text
design: {platform: array, subject: patient_id, time: timepoint, time_order: [D01,D02,D03], covariates: [sex, age], numeric: [age], drop: {timepoint: ["NA"]}, set_na: {sex: ["H"]}}
# GSE110487 (rnaseq)  only response_to_treatment besides subject/time
design: {platform: rnaseq, subject: patient_id, time: timepoint, time_order: [T1, T2]}
# GSE273700 (rnaseq)
design: {platform: rnaseq, subject: patient_id, time: timepoint, time_order: [Day1, Day8]}
# ROSE (rnaseq; pending dataset pass)  class Hyper/Hypo-inflammatory; gender; scr_age integer
design: {platform: rnaseq, subject: subject_number, time: timepoint, time_order: ["Day 0","Day 2"], group: class, covariates: [gender, scr_age], numeric: [scr_age]}
# EARLI (rnaseq, cross-sectional)  lca_label Hyper/Hypo; sex factor; age numeric
design: {platform: rnaseq, group: lca_label, covariates: [sex, age], numeric: [age]}
# MARS (array, cross-sectional; pending dataset pass)  healthy_control is logical
design: {platform: array, group: healthy_control, covariates: [sex, age], numeric: [age]}
# dilgom (array, cross-sectional)  age_band 25-34 ... 65-74 (5 levels); the old analysis compared 65-74 with 25-34
design: {platform: array, group: age_band, covariates: [sex], contrasts: {Age65_74_vs_25_34: "65-74 - 25-34"}}
# hfgp-500fg (rnaseq, cross-sectional)  positive control: sex (top DE hits should be Y-chromosome genes)
design: {platform: rnaseq, group: sex, covariates: [age], numeric: [age]}
```

Facts to keep in mind:
- Replicates at one subject-timepoint are almost absent (ANEMONES 2 cells of 445, GSE95233 2 of 122,
  GSE110487 3 of 59), so technical noise cannot be separated from subject-specific change; do not add
  a `(1|subject:time)` term.
- Subject-level covariates (sex, age, group) are constant within subject, so `canCorPairs()` will report
  correlation 1 with the subject term and warn; that is the nesting the vignette describes, not an
  error. Expect `isSingular` for many genes (random effects at the boundary); record it, do not "fix" it.
- The Time structure page's earlier numbers (time-only model, variance components by method of moments):
  between-subject share ANEMONES 0.55, CORTICUS 0.55, GSE13904 0.57, GSE110487 0.43, GSE273700 0.42,
  GSE95233 0.36, GSE54514 0.02 (unexplained: group_id is a clean person id, healthy subjects have two
  samples; check this dataset first). Use these as a sanity reference for the new numbers.

## 5. Code to write

Reuse before writing: `pivot_expression_long()` (`R/lib/matrices.R`), `read_sample_metadata()`
(`R/targets/diagnostics.R`), `build_symbol_map()` (`R/lib/ingest/symbol_mapping.R`), `with_transaction`,
`append_rows`, `artifacts_dir()`, `artifact_rel()` (`R/targets/db_write.R`). **Port from `R/de/de_helpers.R`
before deleting it:** `collapse_to_symbol()`, `build_cell_means_design()` (the `make.names` level map)
and the contrast logic of `run_de_timecourse()` / `run_de_case_control()`.

- `R/lib/design.R`
  - `validate_design(design)` (called from `validate_dataset_metadata()`).
  - `prepare_design(meta, sample_meta)` -> `list(spec, data, vp_formula, dream_formula, cell_levels, contrasts, cancor, dropped)`:
    set rownames of `data` to the sample ids (variancePartition warns otherwise); apply `drop`/`keep`/`set_na`;
    coerce `numeric`; convert the other model columns to factors (time with `time_order`, logical -> factor);
    drop samples with missing model variables and count them; `cell` = `interaction(group, time, sep = ".", lex.order = TRUE)`
    (time alone / group alone when the other is absent), with `make.names` safe levels and a readable map.
  - Formulas: `vp_formula = ~ (1|subject) + (1|time) + (1|group) + (1|group:time) + (1|cov) ... + numeric`
    (omit absent terms; cross-sectional has no subject/time); `dream_formula = ~ 0 + cell + (1|subject)`
    (`~ 0 + cell` when cross-sectional; covariates only if the config adds `de_covariates`).
  - `build_contrasts(design_obj)`: time-course: for each group level (or "all") with >= 3 subjects having
    >= 2 distinct timepoints in that group, one contrast per later timepoint vs the group's first timepoint
    in `time_order` (`cell<group.t> - cell<group.t0>`), plus an `omnibus` spec (name -> its contrast names,
    for the joint F at ingest); case/control: all pairwise group contrasts or the config override. Build
    the matrix with `variancePartition::makeContrastsDream(dream_formula, data, contrasts = ...)`.
  - `cancor_pairs(vp_formula, data)`: `variancePartition::canCorPairs()` wrapped so its warning is captured, not raised.
- `R/lib/expression.R`: `collapse_to_symbol()` (ported), `prepare_expression(meta, design_obj)`: arrays ->
  collapsed log matrix (no variance filter); RNA-seq -> integer counts -> `DGEList` -> `filterByExpr(dge,
  model.matrix(~ 0 + cell, data))` -> TMM; samples aligned to `design_obj$data`; keep the Ensembl id as the
  gene key and carry a symbol column (from `build_symbol_map()`) for display.
- `R/lib/association.R`: `generic_association_scan()` (prerequisite A).
- `R/targets/variance.R`: the factory (section 7) and the compute functions:
  - `normalize_for_models(expr, design_obj, cpus)`: RNA-seq -> `variancePartition::voomWithDreamWeights(dge, dream_formula,
    data, BPPARAM = SnowParam(cpus, "SOCK"))`; arrays -> the matrix. For cross-sectional RNA-seq confirm
    `voomWithDreamWeights` works with a fixed-effect-only formula (it should fall back to a fixed model); if it
    does not, use `limma::voom(dge, model.matrix(dream_formula, data))`.
  - `run_variance_partition(norm, design_obj, cpus)`: `fitVarPartModel(norm, vp_formula, data, REML = FALSE,
    fxn = vp_record, showWarnings = FALSE, BPPARAM = SnowParam(cpus, "SOCK"))`.
    **Every call inside `vp_record` must be namespace-qualified** (the Snow workers do not attach
    packages; the prototype failed with "could not find function calcVarPart" until it was):
    ```r
    vp_record <- function(fit) list(
      fractions   = variancePartition::calcVarPart(fit),
      variances   = variancePartition::calcVarPart(fit, returnFractions = FALSE),
      sigma       = stats::sigma(fit), nobs = stats::nobs(fit),
      singular    = lme4::isSingular(fit),
      colinearity = tryCatch(variancePartition::colinearityScore(fit), error = function(e) NA_real_))
    ```
    Fixed-effects-only models (no random term) are `lm` fits: guard `isSingular` with `inherits(fit, "merMod")`.
  - `run_dream(norm, design_obj, cpus)`: `dream(norm, dream_formula, data, L, BPPARAM = ...)` then `eBayes()`;
    default `ddf = "adaptive"`. `dream` returns an `MArrayLM2` (coefficients, `rdf`, `sigma`, `stdev.unscaled`,
    `residuals`, `hatvalues`, `cov.coefficients.list`, `t`, `p.value`, `formula`, `data`, ...); about 1.1 MB per 300
    genes x 94 samples, so roughly 60 MB for a 16k-gene array dataset and a few hundred MB for ANEMONES: acceptable,
    but keep an eye on the store.

## 6. What is saved (the "full" output)

Both targets are ordinary targets values (shared store), so ingest can change what it derives without a rerun.
- `design_<ds>`: the validated design object of `prepare_design()` (data frame used, formulas, contrasts spec and
  matrix, `canCorPairs` matrix, dropped-sample counts and reasons).
- `vp_<ds>`: `list(formula, data, genes (id + symbol), records (per-gene `vp_record`), varPart (fractions as a gene x term
  data frame, from `records`), errors (`attr(., "errors")` when present), reml, weights_used, versions)`.
  The lmer objects are deliberately not stored (the vignette's own warning about memory); `vp_record` holds what
  `extractVarPart()` uses plus the absolute variances.
- `dream_<ds>`: `list(fit, fit_eb, L, contrast_table (name, kind, expression), omnibus spec, formula, ddf, norm (voom
  EList or matrix), data, errors, versions)`. `topTable()` is not precomputed: ingest builds it for any coefficient or joint test.

## 7. Targets, controller and resources

New factory `tar_dataset_variance(dataset_id, pipeline, data_root)` in `R/targets/variance.R`, appended in
`_targets.R` beside `tar_dataset_fits`, for every pipeline dataset whose config has a `design:` block. It reuses the
existing `meta_<key>` and `sample_meta_<key>` targets. Target names (key = `target_key(dataset_id)`):

| Target | Controller | Depends on | Notes |
|---|---|---|---|
| `design_<key>` | main (`deployment = "main"`) | `meta`, `sample_meta` | cheap; stops the run early on a bad design |
| `expr_<key>` | `prep` | `meta`, `design` | raw counts/probes; memory like `mat_<key>` |
| `norm_<key>` | `vp` | `expr`, `design` | `voomWithDreamWeights` for RNA-seq (measured 453 s serial for 15,432 genes x 62 samples); pass-through for arrays |
| `vp_<key>` | `vp` | `norm`, `design` | |
| `dream_<key>` | `vp` | `norm`, `design` | |

These do not depend on any `fit_*`/`grid_*` target, so the (days-long) factor fits are not invalidated by anything in
this module; still run `targets::tar_outdated()` after the edits, per the repo rule.
`vp` controller (new entry in `config/pipeline.yml`, picked up by `build_controller_group()`): start with
`cpus_per_task: 8`, `memory_gigabytes_required: 24`, `time_minutes: 1440`, `tasks_max: 1` (one long task per worker),
`crashes_max: 1`, `backup: {memory_gigabytes_required: 48, crashes_max: 2}`. These are guesses to be replaced from
`seff` after the first run (repo target: about 75% average memory and CPU). `SnowParam(workers = cpus_per_task, "SOCK")`
inside the task, with `cpus` read from `pipeline$controllers$vp$cpus_per_task` as the fgsea targets do.

**Measured on the machine where this was written (serial, one core, real data):**

| Dataset | Step | Time |
|---|---|---|
| CORTICUS (arrays, 94 samples, 16,285 genes after collapse) | `fitVarPartModel` | 20 s per 300 genes (0.067 s/gene, so about 18 min for all) |
| | `dream` (one contrast) | 37.6 s per 300 genes (0.125 s/gene, about 34 min for all) |
| GSE110487 (RNA-seq, 62 samples, 39,376 -> 15,432 genes) | `voomWithDreamWeights`, all genes | 453 s |
| | `fitVarPartModel` / `dream`, 300 genes | 9.5 s / 16.8 s |

`SnowParam(4)` on 300 genes was slower (38 s) than serial (20 s): worker start-up is about 20 s, so it only pays off on
the full gene set. Larger datasets (ANEMONES 447 samples, ROSE 221, GSE13904 227) were not measured; per-gene cost
grows with samples and subjects, so measure ANEMONES first and size the controller from it.

## 8. Ingest and database (minimal, so the saved output is usable)

`vp_ingest_<key>` (main, light), supplied to the dataset's `db_<key>` write (`write_dataset_db()` gets a `variance`
argument; `tar_dataset_ingest()` passes it when the target exists). New tables in `R/db/schema.sql`:
- `variance_partition_models(dataset_id PK, formula, platform, n_samples, n_genes, n_genes_failed, share_singular, reml, weights, package_version, fitted_at)`
- `variance_partition_terms(dataset_id, term, kind ('fixed'/'random'/'residual'), n_genes, mean, median, q10, q90; PK (dataset_id, term))` (fractions)
- `de_contrasts(dataset_id, contrast, kind ('pairwise'/'omnibus'), expression, n_genes, n_fdr05, ddf; PK (dataset_id, contrast))`
Gene-level output as `dataset_artifacts` (extend the `kind` CHECK constraint): `variance_partition` (`variance_partition.rds`:
fractions, variances, flags per gene, `canCorPairs`, genes with symbols) and `de_results` (`de_results.rds`: all-gene
`topTable` per contrast with `z.std`, and the joint-F tables for the omnibus specs). Bump `SCHEMA_VERSION`
(`R/db/connect.R`) once for all schema changes in flight (this and any from the deferred TODOS pass), and update
`docs/database.qmd`. Derive the "between-subject total" (subject + sex + group + age, per the vignette's nesting note)
here, not in the model.

## 9. Removal

Delete every file in `R/de/` (21). Then fix: `docs/methods.qmd` (replace the "Differential expression (outside the
pipeline)" section with one describing this module: design block, what is computed, what is saved, how to rerun),
`docs/index.qmd` (the `R/de/` row and the methods row text), `docs/data.qmd` (the `design:` reference),
`docs/workflow.qmd` (new targets, the `vp` controller), `docs/database.qmd` (tables, artifacts), `config/dataset_metadata.example.yml`
(`design:` block), `CLAUDE.md` (Open items; scientific decisions: formulas, ML vs REML, no variance filter), and delete the
old DE results reference (`results/de/`) from any text. Document the current state only.

## 10. Testing, in order (cluster)

Log every submitted job id (CLAUDE.md: session scratchpad `jobs/main_jobs.txt`).
1. **Load check.** `Rscript -e 'targets::tar_manifest()'` builds (also proves prerequisite A). `tar_outdated()` shows
   no `fit_`/`grid_` targets.
2. **Config.** Every `design:` block through `read_dataset_metadata()` with the data root set; columns exist; counts of
   dropped samples printed per dataset (GSE95233 should drop the 22 NA-timepoint samples plus any with sex H set to NA).
3. **Dry run on one small dataset, interactive** (`srun --mem=16G --cpus-per-task=4 --time=120 ...`, command in
   CLAUDE.md), a script `logs/dev/test_variance_dataset.R` that sources the new functions and runs CORTICUS end to end
   without targets. Regression numbers from the prototype (first 300 genes after collapse; same formula
   `~ (1|patient_id) + (1|timepoint) + (1|sex) + age`, `time_order` Pre, Post(24h)): median fractions patient_id 0.461,
   sex 0.000, timepoint 0.000, age 0.010, Residuals 0.497; 82% of fits singular; dream contrast Post(24h) vs Pre top gene
   ADGRE1, logFC -0.907, t -10.57, z.std -7.56. For GSE110487 (RNA-seq): `filterByExpr` keeps 15,432 of 39,376;
   first 300 genes `~ (1|patient_id) + (1|timepoint)` medians patient_id 0.409, timepoint 0.000, Residuals 0.587; top
   dream gene 100128420 (logFC 0.553, P 8.1e-4). (Prototype used all samples with `pivot_expression_long()`; no
   filtering beyond that.) Confirm `fxn` works under `SnowParam` (qualified calls) and that results equal the serial run.
4. **Cross-checks that validate the choice** (do before trusting the numbers): variance fractions against the earlier
   estimator for ANEMONES with a subject + time model (between-subject about 0.55); dream against limma
   `duplicateCorrelation()` p-values (`variancePartition::plotCompareP`, as `dream.Rmd` does) and the consensus
   correlation against the median subject fraction, on CORTICUS, ANEMONES and GSE95233.
5. **Targets, one dataset.** `sbatch run_pipeline.sbatch 'c("design_CORTICUS", "expr_CORTICUS", "norm_CORTICUS", "vp_CORTICUS", "dream_CORTICUS")'`
   (a tidyselect subset; do not use a pattern that also matches the `fit_` targets), watch `logs/targets_main_<jobid>.out` and `logs/crew/`; then GSE110487 (RNA-seq) and EARLI
   (cross-sectional, no random effect), then ANEMONES (the largest array time-course; measure time and memory with
   `seff`), then everything.
6. **Ingest dry run.** For one dataset, `vp_ingest` into a scratch DB (like `logs/dev/test_ingest_earli.R`), tables and
   artifacts readable; check the artifact loads in a fresh R session.

## 11. After the run: bring back, and what is deferred

Bring back (or leave in the repo): for each dataset the median/10th/90th percentile fraction per term, the share of genes
with a time or group:time fraction above 5% and 10%, the between-subject total, singular share, genes failed, the number
of DE genes per contrast at FDR 5%, runtimes and `seff` for the `vp` controller. These decide how to treat timepoints and the factor analysis.

Deferred, unchanged by this plan (the earlier planning is in the conversation and `TODOS.md`): the TODOS pass (registry cleanup,
`fit_failed(power)`, sPCA trace-based convergence, sPCA K = 5 scan with config-overridable defaults, doc-pointer strings, DB schema
headers and an `enrichment_genes` artifact, metadata arrow fix, comment cleanup); adding ROSE x3 / MARS / SHIP-TREND (gains stays out);
the projection-score change (store the doubly-centered transfer score, the target's native rank-K PCA value, a shuffled null and their
ratio; `logs/dev/projection_metric*.R` has the working analysis); and the app change that replaces the Time structure page's live
computation with the stored results.

## 12. Risks and things not verified

- `dream()` object size and memory for the larger datasets (ANEMONES, GSE13904, ROSE); the controller sizes are guesses.
- `voomWithDreamWeights` is slow (7.5 min for 62 samples) and scales with samples and genes; ROSE (221 samples) and EARLI are the larger RNA-seq cases.
- How `targets`/`crew` behave with multi-hour tasks holding a Snow cluster inside a crew worker was not tested; if Snow sockets misbehave
  under the worker, `MulticoreParam` is the alternative (the vignette recommends Snow mainly for memory).
- GSE54514's near-zero between-subject share is unexplained; GSE95233 `sex = H`; CORTICUS `time_order`.
- Many singular fits (about 80% in the CORTICUS prototype) are expected with 2-level random effects and subject-nested covariates; the fractions are still
  usable, but report the share and consider whether a term should be fixed or dropped after the first results.
- Nothing here was run through `targets` (not installed where this was written); expect to debug the factory against the existing `tar_dataset_fits`.

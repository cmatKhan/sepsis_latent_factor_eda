# Differential expression (limma/voom)

Standalone, local (no cluster needed -- every dataset here is small enough
to run on a laptop) differential-expression analysis, one script per
configured BASE dataset (the whole dataset, never a per-timepoint slice
config like `ANEMONES_day1_config.yml`). The scientific target is genes
whose expression differentiates patients **over time** within a sepsis
time-course -- a repeated-measures question, not a one-shot case/control
test -- for the 8 datasets that have real repeated-measures structure. The
other 6 base datasets have no time axis at all; they get a simpler
case/control comparison instead (see the table below).

This lives entirely outside the main ingest/app pipeline: no `limma`,
`edgeR`, or `voom` usage existed anywhere in this repo before this
directory, and none of it is wired into `results/stability.sqlite` or the
Shiny app. Each script reads the sample/feature metadata parquet trio
directly (the same way `R/preprocessing/*.R` scripts do), not through the
app's DB.

## Why not just call `load_input_matrix()` / each dataset's `preprocessing_script`

Every dataset's `preprocessing_script` (used by the factor-analysis side
of this pipeline) applies a transform inappropriate for DE: RNA-seq
scripts run DESeq2's VST (voom needs raw counts, not variance-stabilized
data -- feeding VST output into `voom()` double-transforms and produces a
nonsensical mean-variance trend), and every script (RNA-seq or array)
pre-filters to a fixed top-N-by-variance gene set, which is circular for
DE (it pre-selects for the same kind of signal a DE test is trying to
detect, biasing power toward whichever genes already looked variable
before you controlled for the covariates you're now testing).

Instead, every `R/de/*.R` script calls
`R/lib/matrices.R::pivot_expression_long()` directly -- bypassing
`load_input_matrix()`'s preference for a configured `preprocessing_script`
entirely -- to get the RAW, untransformed feature x sample matrix. This is
what makes raw integer counts recoverable for RNA-seq datasets despite
every RNA-seq config pointing at a VST-applying `preprocessing_script`,
and what avoids the top-variance circularity for microarray datasets.
`R/de/de_helpers.R::filter_low_variance()` applies a much milder,
non-circular filter instead (drop the bottom quartile by TOTAL variance
across every sample, not conditioned on the group/time effect being
tested) -- purely for speed and multiple-testing power, the array analogue
of `edgeR::filterByExpr()` for RNA-seq.

## Platform: voom vs. plain limma

Checked directly against each dataset's raw `expression_path` (not just
config comments): RNA-seq datasets (EARLI, GSE110487, GSE273700, ROSE,
hfgp-500fg) have genuine raw integer counts with no negatives, confirmed
recoverable before each one's `preprocessing_script` applies DESeq2's VST
-- these use `edgeR::DGEList()` + `filterByExpr()` + `calcNormFactors()` +
`limma::voom()`. Every other dataset (ANEMONES, CORTICUS, GSE13904,
GSE54514, GSE95233, dilgom, gains, MARS, SHIP-TREND) is microarray,
already on a continuous/log2(-ish) scale (several already
negative/mean-centered, i.e. already normalized upstream of this repo) --
these use plain `limma::lmFit()`/`eBayes()`, no voom.

## Statistical design (the 8 time-course datasets)

Single shared "cell-means" recipe in `R/de/de_helpers.R::run_de_timecourse()`:

- **Design**: `~ 0 + cell`, where `cell` is `interaction(group, time)` (or
  bare `time` if the dataset has no group covariate) -- one coefficient
  per (group x time) combination actually present. This is what "include
  the condition in the model" (see below) becomes concretely: one uniform
  recipe that degenerates cleanly when there's no group, handles ragged
  group x time cells (e.g. a healthy-control arm only sampled at 2 of a
  sepsis arm's 5 timepoints) with no special-casing, and lets every
  contrast of interest be written directly as a difference of cell means
  rather than fighting interaction-term coefficient interpretation.
- **Blocking**: `limma::duplicateCorrelation()` on the subject id --
  handles a single repeated-measures random effect without spending a
  degree of freedom per subject (which a subject fixed effect would need,
  and which cannot handle an unbalanced/ragged design at all). RNA-seq
  datasets use the standard two-pass voom + duplicateCorrelation recipe
  (fit once for preliminary weights, re-estimate correlation, refit).
- **Contrasts reported**: for every group level with at least
  `min_pairs` (default 3) subjects sampled at >=2 of that group's
  timepoints, an overall F-test across all of that group's timepoint
  contrasts (`<group>_time_omnibus` -- the headline "does expression
  change over time within this group" result, an ANOVA-style test built
  from CONTRASTS between cell means, not a raw test of the means
  themselves) plus each pairwise timepoint contrast against that group's
  earliest timepoint. A group level below `min_pairs` still contributes
  its samples to the fit (nothing is dropped) but gets no reported
  contrast.

## Case/control design (the other 6 datasets)

`R/de/de_helpers.R::run_de_case_control()`: `~ 0 + group`, no blocking,
one contrast (or every pairwise comparison, if more than 2 levels) between
the most sepsis-relevant levels available for that dataset.

## Two decisions that shaped every script here

1. **Where a time-course dataset has a healthy-control/other-condition
   arm** (ANEMONES's `disease`, GSE13904's `clinical_status`, GSE54514's
   `disease_status`), that label is included in the model as `group_col`,
   not subsetted out -- a directive decision, not the default this
   analysis would have picked on its own (subsetting to sepsis-only was
   the initially-recommended option). Datasets with no such arm
   (CORTICUS, GSE110487, GSE273700) use no group covariate at all --
   there's no "other condition" to include.
2. **No covariate-maximizing** -- age, sex, treatment arm, severity
   scores, etc. are deliberately left OUT of every model even where
   available (CORTICUS's `treatment`, GSE54514's `apache_ii`, ROSE's
   cell-composition fractions, ...). This is a first pass, not the most
   complete model each dataset could support -- see the per-dataset
   scripts' own header comments for exactly which covariates exist but
   aren't used.

## Per-dataset summary

| Dataset | Platform | Design | Subject col | Time / group |
|---|---|---|---|---|
| ANEMONES | array | time-course | patient_id | timepoint x disease (Sepsis/SIRS/Control) |
| CORTICUS | array | time-course | patient_id | timepoint (Pre/Post24h); no group |
| GSE110487 | rnaseq | time-course | patient_id | timepoint (T1/T2); no group |
| GSE13904 | array | time-course | patient_id | timepoint x clinical_status (5 levels) |
| GSE273700 | rnaseq | time-course | patient_id | timepoint (Day1/Day8); no group |
| GSE54514 | array | time-course | group_id | day x disease_status (healthy/survivor/nonsurvivor) |
| GSE95233 | array | time-course | patient_id | timepoint; 22 no-timepoint controls dropped (structural, not a subsetting choice) |
| ROSE | rnaseq | time-course | subject_number | timepoint (Day 0/Day 2) x class (Hypo/Hyperinflammatory) |
| EARLI | rnaseq | case/control | -- | lca_label (Hypo/Hyper) |
| MARS | array | case/control | -- | healthy_control -> Sepsis vs Healthy |
| gains | array | case/control | -- | srs_group (SRS1/SRS2); bespoke matrix load (partitioned files, see script header) |
| dilgom | array | case/control | -- | age_band (oldest vs youngest) -- **no sepsis-relevant covariate exists**, illustrative only |
| hfgp-500fg | rnaseq | case/control | -- | sex -- **no sepsis-relevant covariate exists**, illustrative only (doubles as a positive-control QC check: top hits are Y-chromosome genes, as expected) |
| SHIP-TREND | -- | none | -- | **no sample_metadata.parquet exists at all** -- script documents this and does not run |

## Output

`results/de/<DATASET>/`:
- `topTable_<contrast>.csv` -- one per contrast (limma `topTable()`
  output; `symbol` column added for RNA-seq datasets whose native feature
  id is Entrez/Ensembl, via `R/lib/ingest/symbol_mapping.R::build_symbol_map()`
  -- the same shared helper `de_factor_helpers.R`'s gene-id reconciliation
  uses, not a separate DE-specific copy).
- `fit.rds` -- the full post-`eBayes()` `MArrayLM` object.
- `de_bundle.rds` -- `list(fit, expr, time, subject, group)`: the
  PRE-contrast cell-means fit (`fit$coefficients` is a genes x
  (group.time) matrix of per-cell mean expression), the expression matrix
  the fit ran on, and the sample-level factors -- everything
  `run_de_clustering.R`/`run_de_factor_correlation.R` need, with no
  recomputation. `time`/`subject` are `NULL` for the 6 case/control
  datasets (no time axis).
- `diag_pca.png` -- sample PCA colored by the main grouping factor.
- `diag_pval_hist.png` -- one p-value histogram per contrast (should be
  roughly uniform with a spike near 0; a spike near 1 instead usually
  means a missing covariate or unaccounted-for correlation).
- `diag_voom_mean_var.png` -- RNA-seq only, the post-`eBayes` mean-variance
  trend (`limma::plotSA()`), should be flat.

## Follow-on analyses

Three further analyses build on the `topTable`/`de_bundle` outputs above --
each is its own driver script, run once, reading already-computed DE
results rather than refitting anything.

### 1. fgsea enrichment (`R/de/run_de_fgsea.R`, `R/de/de_fgsea_helpers.R`)

Mirrors `R/ingest_jobs/fgsea_job.R`'s conventions for the factor-analysis
pipeline exactly -- same 6 msigdbr collections (Hallmark, GO:BP, GO:MF,
KEGG, REACTOME, WikiPathways), same `fgsea(minSize=10, maxSize=500)`
parameters, same canonical Ensembl gene-id space. Only **pairwise**
contrasts get fgsea (files with "omnibus" in the name are skipped) --
GSEA needs a signed ranking statistic, and an F-test has no sign; each
pairwise contrast's `t` column is the rank vector. Since DE's own gene-id
convention isn't always Ensembl (array datasets collapse to gene SYMBOL;
GSE110487 uses Entrez), `build_de_ensembl_map()` remaps per-dataset first
-- for array datasets this builds a NEW symbol-to-Ensembl lookup directly
from `feature_metadata.parquet` (the existing `build_ensembl_map()` is
probe-keyed, not symbol-keyed). Confirmed healthy mapping rates across
datasets (58-100%, never near-zero, which would indicate a broken join).
Output: `fgsea_<contrast>.csv` (all 6 sources stacked, `padj < 0.05` rows
only, same filter convention `enrichment_cache` uses).

### 2. Clustering + time recovery (`R/de/run_de_clustering.R`)

For each of the 8 time-course datasets, per group (or the whole dataset
where there's no group covariate):
- **Gene clustering by temporal shape**: significant genes (`adj.P.Val <
  0.05` in that group's `*_time_omnibus` contrast) pulled straight from
  `de_bundle.rds`'s `fit$coefficients`, z-scored, `hclust()`'d
  (correlation distance) into 6 clusters -- `gene_clusters_<contrast>.csv`
  + `gene_cluster_heatmap_<contrast>.png`. Real output shows clean,
  visually distinct trajectory shapes (monotonic increase, early spike
  then decline, late-onset change, etc.), not noise.
- **Does DE-gene expression actually recover time structure?** Samples
  clustered on the same significant genes' RAW per-sample expression,
  cut to the number of true timepoint levels, scored against the true
  time labels via `mclust::adjustedRandIndex()`. Real result across all
  8 datasets: mostly near-zero ARI, meaningfully positive only for
  GSE54514's smaller subgroups (healthy: 0.23, sepsis nonsurvivor: 0.32).
  This is an honest, expected finding, not a bug: `duplicateCorrelation()`
  deliberately isolates the (real, statistically significant) WITHIN-
  subject time effect from the much larger BETWEEN-subject variation --
  naive unsupervised clustering of raw expression has no such adjustment,
  so it's dominated by between-subject differences except where the time
  trajectory itself is unusually strong (e.g. a rapidly-deteriorating
  nonsurvivor subgroup). Summary across all datasets/groups:
  `results/de/sample_time_clustering_summary.csv`.

### 3. DE-vs-factor correlation (`R/de/run_de_factor_correlation.R`, `R/de/de_factor_helpers.R`)

For every dataset with both DE results and at least one ingested
factor-analysis method (11/13 -- MARS and gains were never ingested into
`results/stability.sqlite` at all, since their preprocessing pipelines
are incomplete/broken; both are skipped with a clear message, not
silently), correlates that dataset's DE contrasts against PCA/NMF/
CoGAPS/sPCA/ICA using ONE fixed rank per method (10, or K=10 for sPCA --
deliberately sidesteps model-selection for this first pass; see
`pick_fit_at_rank()`). Two tests, both per (factor, contrast):
- **Overlap** (`overlap_test()`): Fisher's exact test between a factor's
  gene set -- non-zero genes for sPCA, top-100 by \|weight\| for every
  other (dense) method -- and a contrast's significant genes, over the
  shared gene universe. Run for every contrast including the omnibus
  ones (a binary overlap test doesn't need a sign). Output:
  `factor_overlap.csv`.
- **GSEA-style** (`factor_gsea_test()`): rank every gene by a contrast's
  signed `t` statistic, test whether each factor's gene set is enriched
  at either end via `fgsea()` -- pairwise contrasts only (needs a sign).
  Output: `factor_gsea.csv`.

Gene-id reconciliation here is dataset-native (no Ensembl remap needed,
unlike Part 1) -- RNA-seq datasets' loadings already match DE's
identifier directly; array datasets' loadings are remapped from probe id
to gene symbol via `R/lib/ingest/symbol_mapping.R::remap_to_symbol()`
(an existing display-only helper, reused here for a real computational
join). Real results confirmed strong and biologically coherent, e.g.
ANEMONES: several PCA/ICA components enriched at `padj` as low as 1e-49
against the Day-5/Discharge-vs-Day-1 Sepsis contrasts, and one PCA
component's top-100 genes 100% overlapping a Sepsis time-omnibus hit list
(Fisher's exact `p = 6e-34`).

## Re-running

```bash
Rscript R/de/ANEMONES_de.R          # per-dataset DE (14 scripts)
Rscript R/de/run_de_fgsea.R              # Part 1, all datasets at once
Rscript R/de/run_de_clustering.R         # Part 2, the 8 time-course datasets
Rscript R/de/run_de_factor_correlation.R # Part 3, all datasets at once
```

Each per-dataset DE script is self-contained (sources `R/lib/matrices.R`
and `R/de/de_helpers.R` itself, no shared setup step needed).
`duplicateCorrelation()` on the larger array datasets (ANEMONES: 447
samples; GSE13904: 227; ROSE: 221) is the slow step -- 10-30 minutes each
on a laptop. Every other per-dataset script finishes in well under a
minute. The three follow-on drivers each run in a few minutes (fgsea's
`msigdbr` fetch is the slowest single step, done once per run).

-- Stability DB schema.
-- Created by open_stability_db() (R/db/connect.R) on a new database, which
-- then stamps PRAGMA user_version with SCHEMA_VERSION. There are no
-- migrations: the DB is rebuilt from the targets store.
--
-- Statements are separated by a line holding only ";;" (connect.R splits
-- on it), since triggers or CHECKs may contain semicolons.

-- ---- dataset level ---------------------------------------------------------

CREATE TABLE datasets (
  dataset_id              TEXT PRIMARY KEY,
  description             TEXT,
  n_samples               INTEGER NOT NULL,
  n_genes                 INTEGER NOT NULL,
  sample_id_col           TEXT,
  feature_id_col          TEXT,
  ensembl_col             TEXT,
  symbol_col              TEXT,
  subject_id_col          TEXT,
  sample_metadata_source  TEXT,
  feature_metadata_source TEXT,
  written_at              TEXT NOT NULL
)
;;
CREATE TABLE dataset_artifacts (
  dataset_id TEXT NOT NULL REFERENCES datasets(dataset_id) ON DELETE CASCADE,
  kind       TEXT NOT NULL CHECK (kind IN ('matrix', 'sample_metadata', 'feature_metadata',
                                           'variance_partition', 'de_results', 'composition',
                                           'enrichment_genes', 'pattern_markers')),
  path       TEXT NOT NULL,
  PRIMARY KEY (dataset_id, kind)
)
;;
CREATE TABLE matrix_diagnostics (
  dataset_id         TEXT PRIMARY KEY REFERENCES datasets(dataset_id) ON DELETE CASCADE,
  val_min            REAL,
  val_max            REAL,
  val_mean           REAL,
  val_sd             REAL,
  skew_p50           REAL,
  skew_p95           REAL,
  kurtosis_p50       REAL,
  kurtosis_p95       REAL,
  min_singular_value REAL,
  max_singular_value REAL,
  condition_number   REAL,
  max_sample_cor     REAL,
  computed_at        TEXT NOT NULL
)
;;

-- ---- methods and fits ------------------------------------------------------

CREATE TABLE methods (
  method         TEXT PRIMARY KEY,
  family         TEXT NOT NULL CHECK (family IN ('seed_sweep', 'param_grid')),
  sign_ambiguous INTEGER NOT NULL CHECK (sign_ambiguous IN (0, 1)),
  has_loadings   INTEGER NOT NULL CHECK (has_loadings IN (0, 1))
)
;;
CREATE TABLE fits (
  fit_id         INTEGER PRIMARY KEY,
  fit_key        TEXT NOT NULL UNIQUE,
  dataset_id     TEXT NOT NULL REFERENCES datasets(dataset_id) ON DELETE CASCADE,
  method         TEXT NOT NULL REFERENCES methods(method),
  rank           INTEGER,
  seed           INTEGER,
  bootstrap      INTEGER NOT NULL DEFAULT 0 CHECK (bootstrap IN (0, 1)),
  representative INTEGER NOT NULL DEFAULT 0 CHECK (representative IN (0, 1)),
  mse            REAL,
  n_factors      INTEGER,
  status         TEXT NOT NULL CHECK (status IN ('ok', 'failed')),
  error          TEXT,
  written_at     TEXT NOT NULL
)
;;
CREATE INDEX idx_fits_dataset_method_rank ON fits(dataset_id, method, rank)
;;
CREATE TABLE fit_params (
  fit_id     INTEGER NOT NULL REFERENCES fits(fit_id) ON DELETE CASCADE,
  name       TEXT NOT NULL,
  num_value  REAL,
  text_value TEXT,
  PRIMARY KEY (fit_id, name)
)
;;
CREATE TABLE fit_metrics (
  fit_id INTEGER NOT NULL REFERENCES fits(fit_id) ON DELETE CASCADE,
  name   TEXT NOT NULL,
  value  REAL,
  PRIMARY KEY (fit_id, name)
)
;;
CREATE TABLE fit_artifacts (
  fit_id INTEGER NOT NULL REFERENCES fits(fit_id) ON DELETE CASCADE,
  kind   TEXT NOT NULL CHECK (kind IN ('loadings', 'scores', 'diag', 'raw', 'dendrogram',
                                       'redundancy', 'kme')),
  path   TEXT NOT NULL,
  PRIMARY KEY (fit_id, kind)
)
;;
CREATE TABLE factors (
  factor_id          INTEGER PRIMARY KEY,
  fit_id             INTEGER NOT NULL REFERENCES fits(fit_id) ON DELETE CASCADE,
  factor_index       INTEGER NOT NULL,
  label              TEXT,
  n_genes            INTEGER,
  stability_cosine   REAL,
  stability_pearson  REAL,
  stability_spearman REAL,
  kurtosis           REAL,
  excess_kurtosis    REAL,
  UNIQUE (fit_id, factor_index)
)
;;

-- ---- pair stability --------------------------------------------------------

CREATE TABLE fit_pairs (
  fit_pair_id INTEGER PRIMARY KEY,
  fit_a       INTEGER NOT NULL REFERENCES fits(fit_id) ON DELETE CASCADE,
  fit_b       INTEGER NOT NULL REFERENCES fits(fit_id) ON DELETE CASCADE,
  same_rank   INTEGER NOT NULL CHECK (same_rank IN (0, 1)),
  n_factors_a INTEGER NOT NULL,
  n_factors_b INTEGER NOT NULL,
  ari         REAL,
  UNIQUE (fit_a, fit_b),
  CHECK (fit_a < fit_b)
)
;;
CREATE INDEX idx_fit_pairs_b ON fit_pairs(fit_b)
;;
CREATE TABLE factor_matches (
  match_id    INTEGER PRIMARY KEY,
  fit_pair_id INTEGER NOT NULL REFERENCES fit_pairs(fit_pair_id) ON DELETE CASCADE,
  factor_a    INTEGER NOT NULL REFERENCES factors(factor_id) ON DELETE CASCADE,
  factor_b    INTEGER NOT NULL REFERENCES factors(factor_id) ON DELETE CASCADE,
  UNIQUE (fit_pair_id, factor_a)
)
;;
CREATE INDEX idx_factor_matches_a ON factor_matches(factor_a)
;;
CREATE INDEX idx_factor_matches_b ON factor_matches(factor_b)
;;
CREATE TABLE match_scores (
  match_id  INTEGER NOT NULL REFERENCES factor_matches(match_id) ON DELETE CASCADE,
  metric    TEXT NOT NULL CHECK (metric IN ('cosine', 'pearson', 'spearman', 'jaccard')),
  value     REAL NOT NULL,
  runner_up REAL,
  margin    REAL,
  PRIMARY KEY (match_id, metric)
)
;;
CREATE TABLE fit_pair_null (
  fit_pair_id INTEGER NOT NULL REFERENCES fit_pairs(fit_pair_id) ON DELETE CASCADE,
  metric      TEXT NOT NULL CHECK (metric IN ('cosine', 'pearson', 'spearman', 'jaccard')),
  n           INTEGER NOT NULL,
  mean        REAL,
  sd          REAL,
  median      REAL,
  p95         REAL,
  max         REAL,
  PRIMARY KEY (fit_pair_id, metric)
)
;;
CREATE TABLE similarity_histograms (
  dataset_id  TEXT NOT NULL REFERENCES datasets(dataset_id) ON DELETE CASCADE,
  method      TEXT NOT NULL REFERENCES methods(method),
  level_a     INTEGER NOT NULL,
  level_b     INTEGER NOT NULL,
  metric      TEXT NOT NULL CHECK (metric IN ('cosine', 'pearson', 'spearman', 'jaccard')),
  bin         INTEGER NOT NULL,
  bin_lo      REAL NOT NULL,
  bin_hi      REAL NOT NULL,
  n_matched   INTEGER NOT NULL,
  n_unmatched INTEGER NOT NULL,
  PRIMARY KEY (dataset_id, method, level_a, level_b, metric, bin)
)
;;

-- ---- WGCNA -----------------------------------------------------------------

CREATE TABLE gene_modules (
  fit_id      INTEGER NOT NULL REFERENCES fits(fit_id) ON DELETE CASCADE,
  gene        TEXT NOT NULL,
  module      INTEGER NOT NULL,
  kme_own     REAL,
  next_module INTEGER,
  kme_next    REAL,
  margin      REAL,
  p_analytic  REAL,
  p_perm      REAL,
  PRIMARY KEY (fit_id, gene)
)
;;
CREATE INDEX idx_gene_modules_module ON gene_modules(fit_id, module)
;;
CREATE TABLE module_kme_summary (
  fit_id INTEGER NOT NULL REFERENCES fits(fit_id) ON DELETE CASCADE,
  module INTEGER NOT NULL,
  grp    TEXT NOT NULL CHECK (grp IN ('member', 'nonmember', 'permutation')),
  n      INTEGER NOT NULL,
  mean   REAL,
  sd     REAL,
  median REAL,
  p95    REAL,
  p99    REAL,
  PRIMARY KEY (fit_id, module, grp)
)
;;
CREATE TABLE module_kme_histograms (
  fit_id INTEGER NOT NULL REFERENCES fits(fit_id) ON DELETE CASCADE,
  module INTEGER NOT NULL,
  grp    TEXT NOT NULL CHECK (grp IN ('member', 'nonmember', 'permutation')),
  bin    INTEGER NOT NULL,
  bin_lo REAL NOT NULL,
  bin_hi REAL NOT NULL,
  count  INTEGER NOT NULL,
  PRIMARY KEY (fit_id, module, grp, bin)
)
;;
CREATE TABLE wgcna_sft (
  dataset_id     TEXT NOT NULL REFERENCES datasets(dataset_id) ON DELETE CASCADE,
  power          INTEGER NOT NULL,
  sft_r_sq       REAL,
  slope          REAL,
  truncated_r_sq REAL,
  mean_k         REAL,
  median_k       REAL,
  max_k          REAL,
  PRIMARY KEY (dataset_id, power)
)
;;
CREATE TABLE gene_significance (
  dataset_id TEXT NOT NULL REFERENCES datasets(dataset_id) ON DELETE CASCADE,
  field      TEXT NOT NULL,
  gene       TEXT NOT NULL,
  test       TEXT NOT NULL CHECK (test IN ('spearman', 'kruskal')),
  statistic  REAL,
  p_value    REAL,
  PRIMARY KEY (dataset_id, field, gene)
)
;;

-- ---- ICA / ICASSO ----------------------------------------------------------

CREATE TABLE icasso_clusters (
  dataset_id           TEXT NOT NULL REFERENCES datasets(dataset_id) ON DELETE CASCADE,
  rank                 INTEGER NOT NULL,
  cluster_id           INTEGER NOT NULL,
  iq                   REAL,
  n_members            INTEGER NOT NULL,
  centrotype_factor_id INTEGER REFERENCES factors(factor_id) ON DELETE CASCADE,
  dendro_file          TEXT,
  PRIMARY KEY (dataset_id, rank, cluster_id)
)
;;
CREATE TABLE icasso_membership (
  factor_id  INTEGER PRIMARY KEY REFERENCES factors(factor_id) ON DELETE CASCADE,
  dataset_id TEXT NOT NULL,
  rank       INTEGER NOT NULL,
  cluster_id INTEGER NOT NULL,
  intra_sim  REAL,
  FOREIGN KEY (dataset_id, rank, cluster_id)
    REFERENCES icasso_clusters(dataset_id, rank, cluster_id) ON DELETE CASCADE
)
;;

-- ---- enrichment ------------------------------------------------------------

CREATE TABLE enrichment_queries (
  query_id   INTEGER PRIMARY KEY,
  factor_id  INTEGER NOT NULL REFERENCES factors(factor_id) ON DELETE CASCADE,
  query_type TEXT NOT NULL,
  direction  TEXT NOT NULL CHECK (direction IN ('pos', 'neg', 'both')),
  query_size INTEGER,
  queried_at TEXT NOT NULL,
  UNIQUE (factor_id, query_type, direction)
)
;;
CREATE TABLE enrichment_results (
  query_id          INTEGER NOT NULL REFERENCES enrichment_queries(query_id) ON DELETE CASCADE,
  term_id           TEXT NOT NULL,
  source            TEXT NOT NULL,
  term_name         TEXT,
  p_value           REAL NOT NULL,
  intersection_size INTEGER,
  term_size         INTEGER,
  nes               REAL,
  es                REAL,
  log2err           REAL,
  is_main_pathway   INTEGER CHECK (is_main_pathway IN (0, 1)),
  PRIMARY KEY (query_id, term_id)
)
;;
CREATE INDEX idx_enrichment_results_p ON enrichment_results(query_id, p_value)
;;

-- ---- cross-dataset ---------------------------------------------------------

CREATE TABLE projections (
  projection_id     INTEGER PRIMARY KEY,
  source_fit_id     INTEGER NOT NULL REFERENCES fits(fit_id) ON DELETE CASCADE,
  target_dataset_id TEXT NOT NULL REFERENCES datasets(dataset_id) ON DELETE CASCADE,
  projection_type   TEXT NOT NULL CHECK (projection_type IN ('within_dataset', 'cross_dataset')),
  include_intercept INTEGER NOT NULL CHECK (include_intercept IN (0, 1)),
  n_genes_matched   INTEGER,
  n_samples         INTEGER,
  mean_r_squared    REAL,
  median_r_squared  REAL,
  mean_pval         REAL,
  median_pval       REAL,
  mean_pvar         REAL,
  median_pvar       REAL,
  path              TEXT NOT NULL,
  UNIQUE (source_fit_id, target_dataset_id, include_intercept)
)
;;
CREATE INDEX idx_projections_target ON projections(target_dataset_id)
;;
CREATE TABLE gene_space_agreement (
  source_fit_id     INTEGER NOT NULL REFERENCES fits(fit_id) ON DELETE CASCADE,
  target_fit_id     INTEGER NOT NULL REFERENCES fits(fit_id) ON DELETE CASCADE,
  n_genes_matched   INTEGER NOT NULL,
  mean_abs_diagonal REAL NOT NULL,
  PRIMARY KEY (source_fit_id, target_fit_id)
)
;;

-- ---- representative-fit analyses -------------------------------------------

CREATE TABLE pattern_drivers (
  driver_id            INTEGER PRIMARY KEY,
  factor_id            INTEGER NOT NULL REFERENCES factors(factor_id) ON DELETE CASCADE,
  grouping_col         TEXT NOT NULL,
  group1_level         TEXT NOT NULL,
  group2_level         TEXT NOT NULL,
  mode                 TEXT NOT NULL CHECK (mode IN ('CI', 'PV')),
  n_genes_considered   INTEGER,
  n_significant_shared INTEGER,
  path                 TEXT NOT NULL,
  computed_at          TEXT NOT NULL,
  UNIQUE (factor_id, grouping_col, group1_level, group2_level, mode)
)
;;
CREATE TABLE fit_redundancy (
  fit_id                    INTEGER PRIMARY KEY REFERENCES fits(fit_id) ON DELETE CASCADE,
  max_offdiag_cosine        REAL,
  median_offdiag_cosine     REAL,
  n_factors_with_no_markers INTEGER
)
;;

-- ---- repeated measures: design, sample QC, variance partition ---------------
-- Written for datasets with a `design:` block and, from the parent's design,
-- for their per-timepoint subsets (R/targets/variance.R, variance_explore.R,
-- composition.R).

CREATE TABLE dataset_design (
  dataset_id           TEXT PRIMARY KEY REFERENCES datasets(dataset_id) ON DELETE CASCADE,
  parent_dataset_id    TEXT NOT NULL,
  is_subset            INTEGER NOT NULL CHECK (is_subset IN (0, 1)),
  platform             TEXT,
  subject_col          TEXT,
  time_col             TEXT,
  time_order           TEXT,
  group_col            TEXT,
  covariates           TEXT,
  technical_covariates TEXT,
  vp_formula           TEXT,
  dream_formula        TEXT,
  n_samples            INTEGER,
  n_subjects           INTEGER,
  n_dropped            INTEGER
)
;;
CREATE TABLE samples (
  dataset_id       TEXT NOT NULL REFERENCES datasets(dataset_id) ON DELETE CASCADE,
  sample_id        TEXT NOT NULL,
  subject          TEXT,
  timepoint        TEXT,
  grp              TEXT,
  cell             TEXT,
  expr_mean        REAL,
  expr_sd          REAL,
  expr_q10         REAL,
  expr_median      REAL,
  expr_q90         REAL,
  pc1              REAL,
  pc2              REAL,
  pc3              REAL,
  neutrophil_score REAL,
  lymphocyte_score REAL,
  nl_score         REAL,
  control_mean     REAL,
  control_factor   REAL,
  log2_lib_size    REAL,
  norm_factor      REAL,
  chrY_score       REAL,
  xist             REAL,
  inferred_sex     TEXT,
  recorded_sex     TEXT,
  PRIMARY KEY (dataset_id, sample_id)
)
;;
CREATE TABLE variance_summary (
  dataset_id                      TEXT PRIMARY KEY REFERENCES datasets(dataset_id) ON DELETE CASCADE,
  n_genes                         INTEGER,
  n_genes_failed                  INTEGER,
  n_factor_genes                  INTEGER,
  share_singular                  REAL,
  max_colinearity                 REAL,
  dupcor_consensus                REAL,
  between_median                  REAL,
  time_median                     REAL,
  residual_median                 REAL,
  technical_median                REAL,
  repeatability_median            REAL,
  time_within_median              REAL,
  share_time_gt05                 REAL,
  share_time_gt10                 REAL,
  share_time_gt20                 REAL,
  between_vw_all                  REAL,
  time_vw_all                     REAL,
  residual_vw_all                 REAL,
  between_vw_factor_genes         REAL,
  time_vw_factor_genes            REAL,
  residual_vw_factor_genes        REAL,
  between_vw_top                  REAL,
  time_vw_top                     REAL,
  residual_vw_top                 REAL,
  n_subjects_repeated             INTEGER,
  n_changes                       INTEGER,
  change_pc1_share                REAL,
  change_pc1_cor_nl               REAL,
  median_r2_nl_all                REAL,
  median_r2_nl_factor_genes       REAL,
  median_r2_technical_all         REAL,
  median_r2_technical_factor_genes REAL,
  rho_pc1_nl                      REAL,
  rho_sd_nl                       REAL,
  rho_mean_nl                     REAL,
  rho_sd_technical                REAL,
  rho_mean_technical              REAL,
  n_control_probes                INTEGER,
  control_probes                  TEXT,
  control_factor_fraction         REAL,
  sex_separation                  REAL,
  sex_agreement                   REAL,
  n_subjects_sex_discordant       INTEGER,
  nl_subject_fraction             REAL,
  nl_time_fraction                REAL,
  nl_group_fraction               REAL,
  vp_seconds                      REAL,
  dream_seconds                   REAL,
  fitted_at                       TEXT NOT NULL
)
;;
CREATE TABLE variance_terms (
  dataset_id TEXT NOT NULL REFERENCES datasets(dataset_id) ON DELETE CASCADE,
  term       TEXT NOT NULL,
  category   TEXT NOT NULL,
  n_genes    INTEGER,
  mean       REAL,
  q10        REAL,
  median     REAL,
  q90        REAL,
  share_gt05 REAL,
  share_gt10 REAL,
  PRIMARY KEY (dataset_id, term)
)
;;
CREATE TABLE variance_categories (
  dataset_id      TEXT NOT NULL REFERENCES datasets(dataset_id) ON DELETE CASCADE,
  category        TEXT NOT NULL,
  present         INTEGER NOT NULL CHECK (present IN (0, 1)),
  mean            REAL,
  median          REAL,
  vw_all          REAL,
  vw_factor_genes REAL,
  vw_top          REAL,
  PRIMARY KEY (dataset_id, category)
)
;;
CREATE TABLE variance_bins (
  dataset_id TEXT NOT NULL REFERENCES datasets(dataset_id) ON DELETE CASCADE,
  binning    TEXT NOT NULL CHECK (binning IN ('expression_quartile', 'variance_rank')),
  bin        TEXT NOT NULL,
  n_genes         INTEGER,
  between_subject REAL,
  time            REAL,
  residual        REAL,
  PRIMARY KEY (dataset_id, binning, bin)
)
;;
CREATE TABLE trait_variance (
  dataset_id TEXT NOT NULL REFERENCES datasets(dataset_id) ON DELETE CASCADE,
  trait      TEXT NOT NULL,
  term       TEXT NOT NULL,
  category   TEXT NOT NULL,
  fraction   REAL,
  PRIMARY KEY (dataset_id, trait, term)
)
;;
CREATE TABLE composition_correlations (
  dataset_id TEXT NOT NULL REFERENCES datasets(dataset_id) ON DELETE CASCADE,
  x          TEXT NOT NULL,
  y          TEXT NOT NULL,
  rho        REAL,
  p_value    REAL,
  PRIMARY KEY (dataset_id, x, y)
)
;;
CREATE TABLE de_contrasts (
  dataset_id TEXT NOT NULL REFERENCES datasets(dataset_id) ON DELETE CASCADE,
  contrast   TEXT NOT NULL,
  kind       TEXT NOT NULL CHECK (kind IN ('pairwise', 'omnibus')),
  readable   TEXT,
  n_genes    INTEGER,
  n_fdr05    INTEGER,
  PRIMARY KEY (dataset_id, contrast)
)
;;
CREATE TABLE factor_variance_models (
  dataset_id TEXT NOT NULL REFERENCES datasets(dataset_id) ON DELETE CASCADE,
  method     TEXT NOT NULL REFERENCES methods(method),
  formula    TEXT NOT NULL,
  n_factors  INTEGER,
  PRIMARY KEY (dataset_id, method)
)
;;
CREATE TABLE factor_variance (
  factor_id INTEGER NOT NULL REFERENCES factors(factor_id) ON DELETE CASCADE,
  term      TEXT NOT NULL,
  category  TEXT NOT NULL,
  fraction  REAL,
  PRIMARY KEY (factor_id, term)
)
;;
CREATE TABLE variance_time_genes (
  dataset_id TEXT NOT NULL REFERENCES datasets(dataset_id) ON DELETE CASCADE,
  symbol          TEXT NOT NULL,
  time_all        REAL,
  between_subject REAL,
  PRIMARY KEY (dataset_id, symbol)
)
;;
CREATE TABLE time_gene_recurrence (
  symbol     TEXT PRIMARY KEY,
  n_datasets INTEGER NOT NULL,
  datasets   TEXT NOT NULL
)

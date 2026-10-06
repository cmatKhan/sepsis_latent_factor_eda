-- Stability DB schema. Documented table by table in R/db/SCHEMA.md.
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
  kind       TEXT NOT NULL CHECK (kind IN ('matrix', 'sample_metadata', 'feature_metadata')),
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
  genes             TEXT,
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
CREATE TABLE pattern_markers (
  factor_id INTEGER NOT NULL REFERENCES factors(factor_id) ON DELETE CASCADE,
  gene      TEXT NOT NULL,
  score     REAL NOT NULL,
  PRIMARY KEY (factor_id, gene)
)
;;
CREATE TABLE fit_redundancy (
  fit_id                    INTEGER PRIMARY KEY REFERENCES fits(fit_id) ON DELETE CASCADE,
  max_offdiag_cosine        REAL,
  median_offdiag_cosine     REAL,
  n_factors_with_no_markers INTEGER
)

-- Convenience views over R/db/schema.sql (documented in R/db/SCHEMA.md).
-- Statements are separated by a line holding only ";;".

-- fits + method family, with the commonly filtered parameters pivoted out
-- of fit_params.
CREATE VIEW v_fits AS
SELECT f.*,
       m.family,
       (SELECT num_value FROM fit_params p WHERE p.fit_id = f.fit_id AND p.name = 'para')   AS para,
       (SELECT num_value FROM fit_params p WHERE p.fit_id = f.fit_id AND p.name = 'lambda') AS lambda,
       (SELECT num_value FROM fit_params p WHERE p.fit_id = f.fit_id AND p.name = 'power')  AS power,
       (SELECT num_value FROM fit_params p WHERE p.fit_id = f.fit_id AND p.name = 'alpha')  AS alpha
FROM fits f
JOIN methods m ON m.method = f.method
;;
-- Every match in both orientations (this fit/factor -> the other), with
-- each metric's score, so per-fit and per-factor queries need no UNION.
CREATE VIEW v_factor_matches AS
SELECT fp.fit_pair_id, fm.match_id, fa.dataset_id, fa.method, fp.same_rank,
       fp.fit_a AS fit_id, fa.rank AS rank, fa.seed AS seed, xa.factor_id AS factor_id,
       xa.factor_index AS factor_index,
       fp.fit_b AS other_fit_id, fb.rank AS other_rank, fb.seed AS other_seed,
       xb.factor_id AS other_factor_id, xb.factor_index AS other_factor_index,
       ms.metric, ms.value, ms.runner_up, ms.margin
FROM factor_matches fm
JOIN fit_pairs fp ON fp.fit_pair_id = fm.fit_pair_id
JOIN match_scores ms ON ms.match_id = fm.match_id
JOIN factors xa ON xa.factor_id = fm.factor_a
JOIN factors xb ON xb.factor_id = fm.factor_b
JOIN fits fa ON fa.fit_id = fp.fit_a
JOIN fits fb ON fb.fit_id = fp.fit_b
UNION ALL
SELECT fp.fit_pair_id, fm.match_id, fb.dataset_id, fb.method, fp.same_rank,
       fp.fit_b, fb.rank, fb.seed, xb.factor_id, xb.factor_index,
       fp.fit_a, fa.rank, fa.seed, xa.factor_id, xa.factor_index,
       ms.metric, ms.value, ms.runner_up, ms.margin
FROM factor_matches fm
JOIN fit_pairs fp ON fp.fit_pair_id = fm.fit_pair_id
JOIN match_scores ms ON ms.match_id = fm.match_id
JOIN factors xa ON xa.factor_id = fm.factor_a
JOIN factors xb ON xb.factor_id = fm.factor_b
JOIN fits fa ON fa.fit_id = fp.fit_a
JOIN fits fb ON fb.fit_id = fp.fit_b
;;
-- What the DB holds per (dataset, method); replaces the old `ingests` table.
CREATE VIEW v_ingest_summary AS
SELECT f.dataset_id, f.method, m.family,
       COUNT(*)                     AS n_fits,
       SUM(f.status = 'ok')         AS n_ok,
       SUM(f.status = 'failed')     AS n_failed,
       SUM(f.representative)        AS n_representative,
       MAX(f.written_at)            AS written_at
FROM fits f
JOIN methods m ON m.method = f.method
GROUP BY f.dataset_id, f.method

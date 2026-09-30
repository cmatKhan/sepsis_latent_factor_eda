# SHIP-TREND has NO sample_metadata_path at all (confirmed: no
# sample_metadata.parquet exists anywhere for this dataset -- see
# config/SHIP-TREND_config.yml's own blocking TODO). With zero clinical/
# demographic covariates, there is no group to contrast and no way to
# construct any differential-expression comparison, case/control or
# otherwise -- unlike dilgom/hfgp-500fg, which at least have age/sex to
# fall back on. This script exists (per the user's direction that every
# configured dataset gets a script) purely to document that limitation
# loudly rather than fabricate a fake contrast; it deliberately does not
# attempt to load a matrix or run any model.
#
# If a sample_metadata.parquet is ever added for this dataset (see the
# config's own TODO), this script should be replaced with a real
# case/control (or, if a real subject/timepoint structure turns up,
# time-course) analysis following the same pattern as R/de/MARS_de.R /
# R/de/dilgom_de.R.

stop(
  "SHIP-TREND has no sample_metadata.parquet at all (see config/SHIP-TREND_config.yml's ",
  "own blocking TODO) -- there is no clinical/demographic covariate of any kind to build a ",
  "differential-expression contrast from. This script intentionally does not run. See this ",
  "file's header comment."
)

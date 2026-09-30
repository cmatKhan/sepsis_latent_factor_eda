# Shared dataset manifest for the R/de/run_de_*.R drivers -- dataset_id ->
# (config path, platform). Platform determines which gene-id reconciliation
# path a driver takes (see R/de/de_fgsea_helpers.R / R/de/de_factor_helpers.R
# header comments). SHIP-TREND is deliberately excluded -- it has no DE
# results at all (see R/de/SHIP-TREND_de.R).
DE_DATASETS <- list(
  ANEMONES     = list(config = "config/ANEMONES_config.yml",     platform = "array"),
  CORTICUS     = list(config = "config/CORTICUS_config.yml",     platform = "array"),
  GSE110487    = list(config = "config/GSE110487_config.yml",    platform = "rnaseq"),
  GSE13904     = list(config = "config/GSE13904_config.yml",     platform = "array"),
  GSE273700    = list(config = "config/GSE273700_config.yml",    platform = "rnaseq"),
  GSE54514     = list(config = "config/GSE54514_config.yml",     platform = "array"),
  GSE95233     = list(config = "config/GSE95233_config.yml",     platform = "array"),
  ROSE         = list(config = "config/ROSE_config.yml",         platform = "rnaseq"),
  EARLI        = list(config = "config/EARLI_config.yml",        platform = "rnaseq"),
  MARS         = list(config = "config/MARS_config.yml",         platform = "array"),
  gains        = list(config = "config/gains_config.yml",        platform = "array"),
  dilgom       = list(config = "config/dilgom_config.yml",       platform = "array"),
  `hfgp-500fg` = list(config = "config/hfgp-500fg_config.yml",   platform = "rnaseq")
)

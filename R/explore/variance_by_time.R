# Between-person versus time: plots of the repeated-measures summaries the
# pipeline stores (variance_explore and the vp_explore_<dataset> targets,
# R/targets/variance_explore.R), for deciding how factor analysis treats
# repeated measures.
#
# Interactive: open in RStudio (project root as working directory) and run
# section by section. Writes figures and tables to
# results/variance/exploration/.
#
# How to read the numbers:
#   * Fractions are per gene and sum to 1 (variancePartition, ML fits).
#     Categories: subject (the person), subject-level (variables fixed within
#     a person: sex, age, disease arm), time, group x time, technical
#     (ANEMONES' control-probe factor, ROSE's sequencing QC), residual
#     (within-person variation the model doesn't explain: day-to-day
#     fluctuation, person-specific change, noise).
#   * Medians weight every gene equally; factor methods are driven by
#     high-variance genes, so the variance-weighted shares over the genes the
#     factor analysis uses ("factor genes") are closer to what a fit sees.
#   * Repeatability = between / (between + residual): the share of a
#     single-timepoint snapshot that is stable person-level difference.
#   * Datasets whose sample annotations don't match their arrays show low
#     sex_agreement in the headline (GSE54514: 0.50).

library(targets)
library(dplyr)
library(tidyr)
library(ggplot2)
source("R/targets/factories.R")          # target_key()
source("R/targets/variance_explore.R")   # VARIANCE_CATEGORIES

out_dir <- "results/variance/exploration"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
save_plot <- function(p, name, width = 9, height = 6) {
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width = width, height = height, dpi = 120)
  p
}
category_colours <- stats::setNames(c("#1b6ca8", "#6aa6d6", "#9a77cf", "#d9822b", "#f2b56b", "#7f7f7f",
                                      "#b5b5b5", "#e6e6e6"), VARIANCE_CATEGORIES)

# ---- 1. Load ---------------------------------------------------------------

ve <- tar_read(variance_explore)
datasets <- ve$headline$dataset_id
repeated <- ve$headline$dataset_id[ve$headline$repeated]
genes <- bind_rows(lapply(datasets, function(id) {
  g <- tar_read_raw(paste0("vp_explore_", target_key(id)))$genes
  cbind(dataset_id = id, g)
}))

# ---- 2. Headline ------------------------------------------------------------

compact <- ve$headline |>
  transmute(dataset_id, repeated,
            between_median = round(between_median, 2), time_median = round(time_median, 3),
            between_factor_genes = round(between_vw_factor_genes, 2),
            time_factor_genes = round(time_vw_factor_genes, 3),
            residual_factor_genes = round(residual_vw_factor_genes, 2),
            repeatability = round(repeatability_median, 2),
            genes_time_gt10 = scales::percent(share_time_gt10, 0.1),
            sex_agreement = round(sex_agreement, 2))
print(compact, width = Inf)
write.csv(ve$headline, file.path(out_dir, "headline.csv"), row.names = FALSE)

# ---- 3. Where the variance sits ---------------------------------------------

bars <- ve$categories |>
  select(dataset_id, category, `mean over genes` = mean, `variance-weighted` = vw_all,
         `variance-weighted, factor genes` = vw_factor_genes) |>
  pivot_longer(-c(dataset_id, category), names_to = "view", values_to = "share") |>
  group_by(dataset_id, view) |> mutate(share = share / sum(share, na.rm = TRUE)) |> ungroup() |>
  mutate(category = factor(category, levels = names(category_colours)),
         dataset_id = factor(dataset_id, levels = c(repeated, setdiff(datasets, repeated))))
save_plot(ggplot(bars, aes(dataset_id, share, fill = category)) +
  geom_col(width = 0.8) + facet_wrap(~ view, ncol = 1) +
  scale_fill_manual(values = category_colours, drop = FALSE) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "share of expression variance", fill = NULL,
       title = "Where expression variance sits (time courses first)") +
  theme_bw() + theme(axis.text.x = element_text(angle = 30, hjust = 1)), "variance_composition", height = 9)

# ---- 4. Per gene: between vs time ------------------------------------------

save_plot(ggplot(filter(genes, dataset_id %in% repeated), aes(between, time_all)) +
  geom_bin2d(bins = 50) + scale_fill_viridis_c(trans = "log10") +
  facet_wrap(~ dataset_id) + geom_abline(linetype = 2, colour = "grey40") +
  labs(x = "between-person fraction", y = "time + group x time fraction",
       title = "Between-person vs time, per gene") + theme_bw(), "between_vs_time_scatter", height = 7)

save_plot(ggplot(filter(genes, dataset_id %in% repeated), aes(repeatability, colour = factor_gene)) +
  geom_density() + facet_wrap(~ dataset_id, scales = "free_y") +
  scale_colour_manual(values = c(`TRUE` = "#c0392b", `FALSE` = "grey50"),
                      labels = c(`TRUE` = "factor-analysis genes", `FALSE` = "other genes")) +
  labs(x = "repeatability: between / (between + residual)", colour = NULL,
       title = "How much of a single-timepoint snapshot is stable between-person variation") +
  theme_bw(), "repeatability", height = 7)

# ---- 5. By gene variance rank and expression level --------------------------

ranks <- ve$by_variance_rank |> pivot_longer(c(between, time, residual), names_to = "component",
                                            values_to = "median_fraction")
save_plot(ggplot(filter(ranks, dataset_id %in% repeated),
                 aes(variance_rank, median_fraction, colour = component, group = component)) +
  geom_line() + geom_point() + facet_wrap(~ dataset_id) +
  scale_colour_manual(values = c(between = "#1b6ca8", time = "#d9822b", residual = "grey50")) +
  labs(x = "gene variance rank (1 = most variable)", y = "median fraction") + theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)), "fractions_by_variance_rank", height = 7)

expr <- ve$by_expression |> pivot_longer(c(between, time, residual), names_to = "component",
                                        values_to = "median_fraction")
save_plot(ggplot(filter(expr, dataset_id %in% repeated),
                 aes(expression_quartile, median_fraction, colour = component, group = component)) +
  geom_line() + geom_point() + facet_wrap(~ dataset_id) +
  scale_colour_manual(values = c(between = "#1b6ca8", time = "#d9822b", residual = "grey50")) +
  labs(x = "mean expression quartile (Q1 = lowest)", y = "median fraction",
       title = "Fractions by expression level (low-count noise shows as residual)") + theme_bw(),
  "fractions_by_expression", height = 7)

# ---- 6. Within-person change and time-driven genes --------------------------

print(ve$change, width = Inf)
cat("\nTime-driven genes (time + group x time >", ve$threshold, ") per dataset:\n")
print(table(ve$time_genes$dataset_id))
cat("\nTime-driven in the most datasets:\n")
print(head(ve$time_gene_recurrence[ve$time_gene_recurrence$n_datasets > 1, ], 40), row.names = FALSE)
write.csv(ve$time_gene_recurrence, file.path(out_dir, "time_gene_recurrence.csv"), row.names = FALSE)

cat("\nFigures and tables in", out_dir, "\n")

# Time structure: how much of each time-course dataset's variance is between
# subjects, an average timepoint effect, or other within-subject variation.
# This decides how much of what a factor finds can be about change over time.
# Variance components are estimated per gene (or per factor) with
# variance_components() (R/data_helpers.R); intervals resample subjects.

SHARE_COLORS <- c(`between subjects` = "#4C78A8", `average time effect` = "#F58518",
                  `other within-subject` = "#B8B8B8")

share_long <- function(d, key) {
  out <- rbind(
    data.frame(key = d[[key]], part = "between subjects", share = d$between),
    data.frame(key = d[[key]], part = "average time effect", share = d$time),
    data.frame(key = d[[key]], part = "other within-subject", share = d$resid))
  out$part <- factor(out$part, levels = names(SHARE_COLORS))
  out
}

#' Stacked horizontal or vertical bars of the three parts, with the 95%
#' interval of the between-subject part drawn at its top edge.
share_plot <- function(d, key, horizontal = TRUE) {
  d$key <- factor(d[[key]], levels = unique(d[[key]]))
  long <- share_long(d, key); long$key <- factor(long$key, levels = levels(d$key))
  stack <- position_stack(reverse = TRUE)
  p <- if (horizontal) ggplot(long, aes(share, key, fill = part)) + geom_col(width = 0.7, position = stack) +
         geom_errorbar(data = d, aes(y = key, xmin = between_lo, xmax = between_hi), inherit.aes = FALSE, orientation = "y", width = 0.3) +
         scale_x_continuous(labels = scales::percent, expand = c(0, 0)) + coord_cartesian(xlim = c(0, 1)) + labs(x = "share of total variance", y = NULL)
       else ggplot(long, aes(key, share, fill = part)) + geom_col(width = 0.8, position = stack) +
         geom_errorbar(data = d, aes(x = key, ymin = between_lo, ymax = between_hi), inherit.aes = FALSE, width = 0.4) +
         scale_y_continuous(labels = scales::percent, expand = c(0, 0)) + coord_cartesian(ylim = c(0, 1)) + labs(x = NULL, y = "share of the factor's score variance")
  p + scale_fill_manual(values = SHARE_COLORS, name = NULL) + theme_minimal(base_size = 12) +
    theme(legend.position = "top", axis.text.x = if (horizontal) element_text() else element_text(angle = 60, hjust = 1))
}

mod_time_ui <- function(id) {
  ns <- NS(id)
  tagList(
    p(class = "text-muted",
      "Each time-course dataset's variance is split into three parts. ", tags$strong("Between subjects:"),
      " stable differences between people. ", tags$strong("Average time effect:"),
      " the change every subject shares between timepoints. ", tags$strong("Other within-subject:"),
      " change that differs between subjects, plus measurement noise (the two cannot be separated without technical replicates).",
      " Factors can only describe trajectories to the extent that the last two parts are large.",
      " Parts are variance components estimated per gene and summed (so high-variance genes count more); the whisker is a 95% interval from resampling subjects (not shown for the time effect, whose interval is unreliable near zero)."),
    card(card_header("All time-course datasets (pooled across timepoints)"),
         plotOutput(ns("overview"), height = "340px"), DTOutput(ns("table"))),
    card(card_header(textOutput(ns("gene_title"))),
         p(class = "text-muted", "Each gene's own between-subject share. The line is the pooled value in the bars above.",
           " Right: the same genes weighted by their variance, which is what the pooled value reflects."),
         plotOutput(ns("genes"), height = "300px")),
    card(card_header(textOutput(ns("factor_title"))), uiOutput(ns("factor_note")), plotOutput(ns("factors"), height = "380px")))
}

mod_time_server <- function(id, state, focus, con) {
  moduleServer(id, function(input, output, session) {
    tbl <- reactive(withProgress(message = "Splitting variance", value = 0.5, time_structure_table(con)))
    output$overview <- renderPlot({
      d <- tbl(); validate(need(!is.null(d) && nrow(d) > 0, "No time-course datasets with subject ids were found."))
      share_plot(d[order(d$between), ], "dataset")
    })
    output$table <- renderDT({
      d <- tbl(); req(!is.null(d))
      ci <- function(e, lo, hi) sprintf("%.2f [%.2f, %.2f]", e, lo, hi)
      out <- data.frame(Dataset = d$dataset, Subjects = d$subjects, Samples = d$samples,
                        `Subjects with repeats` = d$repeated, Timepoints = d$timepoints,
                        `Between subjects` = ci(d$between, d$between_lo, d$between_hi),
                        `Average time effect` = sprintf("%.2f", d$time),
                        `Other within-subject` = ci(d$resid, d$resid_lo, d$resid_hi), check.names = FALSE)
      datatable(out[order(d$between), ], rownames = FALSE, selection = "none", options = list(dom = "t", paging = FALSE))
    })

    fam <- reactive(family_of(state$dataset()))
    pooled_id <- reactive(fam()$pooled %||% state$dataset())
    output$gene_title <- renderText(paste0("Across genes: ", pooled_id()))
    output$genes <- renderPlot({
      f <- fam(); validate(need(!is.null(f) && length(f$siblings) > 0, "This dataset has no timepoint siblings."))
      g <- gene_between_fractions(con, f$pooled); validate(need(!is.null(g), "No subject ids for this dataset."))
      d <- tbl(); pooled_value <- d$between[d$dataset == f$pooled]
      long <- rbind(data.frame(between = g$between, w = 1 / nrow(g), view = "each gene counts once"),
                    data.frame(between = g$between, w = g$weight, view = "weighted by gene variance"))
      ggplot(long, aes(between, weight = w)) + geom_histogram(breaks = seq(0, 1, 0.05), fill = "#4C78A8", colour = "white") +
        geom_vline(xintercept = pooled_value, linetype = 2) +
        facet_wrap(~view, nrow = 1) + scale_x_continuous(labels = scales::percent) +
        scale_y_continuous(labels = scales::percent) +
        labs(x = "between-subject share of the gene's variance", y = "share of genes (or of variance)") + theme_minimal(base_size = 12)
    })

    output$factor_title <- renderText(paste0("Factors of the selected fit: ", state$dataset(), " / ", focus_label(focus())))
    samples <- reactive(time_structure_samples(con, state$dataset()))
    output$factor_note <- renderUI({
      f <- fam()
      if (is.null(f) || length(f$siblings) == 0) return(p(class = "text-muted", "This dataset has no timepoint siblings, so it has no within-subject structure to split."))
      if (f$pooled != state$dataset())
        return(p(class = "text-muted", "This is a single-timepoint dataset, so every subject appears once. Select the pooled dataset ", tags$strong(f$pooled), " in the sidebar to split its factors."))
      NULL
    })
    factor_split <- reactive({
      s <- samples(); validate(need(!is.null(s), "Select a pooled time-course dataset in the sidebar to see how its factors split."))
      S <- focus_scores(con, focus()); validate(need(!is.null(S), "No sample scores for this fit."))
      s <- s[s$sample_id %in% rownames(S), , drop = FALSE]
      bs <- bootstrap_fractions(S[s$sample_id, , drop = FALSE], s$subject, s$time, pooled = FALSE)
      data.frame(factor = colnames(S), between = bs$est$between, time = bs$est$time, resid = bs$est$resid,
                 between_lo = bs$lo$between, between_hi = bs$hi$between)
    }) |> bindCache(state$dataset(), focus()$fit_id, focus()$n_factors)
    output$factors <- renderPlot({
      d <- factor_split()
      share_plot(d, "factor", horizontal = FALSE) +
        labs(subtitle = "Factors that are mostly between-subject describe stable differences between people, not change over time.")
    })
  })
}

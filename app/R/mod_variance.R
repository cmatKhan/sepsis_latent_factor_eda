# Question 1: how much of each dataset's variance is change over time, and how much is
# between patients? The answer comes from the stored variancePartition results
# (R/variance_helpers.R); the tabs below it are the evidence and the sample checks that
# decide how far to trust it. A per-timepoint subset shows its parent's results.

# A card that sizes to its content (no fill), expandable to full screen.
vcard <- function(title, ..., min_height = NULL) {
  card(min_height = min_height, full_screen = TRUE, fill = FALSE, card_header(title),
       card_body(fillable = FALSE, fill = FALSE, ...))
}

# The all-datasets overview: not driven by the sidebar.
mod_variance_overview_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h4("How much of the variance is change over time?"),
    p(class = "text-muted",
      "Each bar splits a dataset's gene expression variance into what differs between patients,",
      "what changes with time, technical covariates, and the rest (change within a person that is not",
      "shared by time, plus noise)."),
    uiOutput(ns("notices")),
    card(fill = FALSE, full_screen = TRUE,
      card_header(
        div(class = "d-flex justify-content-between align-items-center",
            span("Variance by source, all time-course datasets"),
            div(style = "width: 340px;",
                selectInput(ns("geneset"), NULL, choices = VARIANCE_GENESETS, width = "100%")))),
      card_body(fillable = FALSE, fill = FALSE, plotOutput(ns("headline"), height = "480px", fill = FALSE)),
      uiOutput(ns("headline_note"))),
    vcard("Summary", DT::DTOutput(ns("headline_table"))),
    vcard("PC1 and the leukocyte marker genes",
          uiOutput(ns("pc1_text")),
          DT::DTOutput(ns("pc1_table"))),
    vcard("Time-driven genes seen in several datasets", DT::DTOutput(ns("recurrence")))
  )
}

mod_variance_overview_server <- function(id, con) {
  moduleServer(id, function(input, output, session) {
    ready <- variance_ready(con)
    design <- if (ready) design_table(con) else NULL
    summ <- if (ready) variance_summary_table(con) else NULL
    ok <- !is.null(design) && !is.null(summ) && nrow(summ) > 0
    output$notices <- renderUI({
      if (!ready || !ok)
        missing_notice("Variance results.", "The database has no repeated-measures tables (schema 2); rerun the pipeline.")
    })
    headline <- reactive({
      req(ok)
      tc <- intersect(time_course_parents(design), summ$dataset_id)
      validate(need(length(tc) > 0, "No time-course datasets have variance results."))
      if (input$geneset == "median") {
        s <- summ[match(tc, summ$dataset_id), ]
        d <- data.frame(dataset = tc, between = s$between_median, time = s$time_median,
                        technical = s$technical_median, residual = s$residual_median)
        d$repeatability <- s$repeatability_median; d$time_within <- s$time_within_median
        d$genes_time_gt10 <- s$share_time_gt10; d$subjects_repeated <- s$n_subjects_repeated
        d$n_genes <- s$n_genes; d$flags <- dataset_flags(s)
        return(d)
      }
      lg <- variance_long(con, tc, input$geneset)
      need_data(lg, "variance_categories")
      headline_table(lg, summ)
    })

    output$headline <- renderPlot({
      h <- headline()
      parts <- data.frame(dataset = rep(h$dataset, 4),
        part = factor(rep(c("between patients", "time", "technical", "residual"), each = nrow(h)),
                      levels = c("between patients", "time", "technical", "residual")),
        share = c(h$between, h$time, h$technical, h$residual))
      parts$dataset <- factor(parts$dataset, levels = h$dataset[order(h$between)])
      cols <- c(`between patients` = "#2C5F8A", time = "#E07B00", technical = "#7B5EA7", residual = "#BDBDBD")
      ggplot(parts, aes(share, dataset, fill = part)) + geom_col(width = 0.7, position = position_stack(reverse = TRUE)) +
        scale_fill_manual(values = cols, name = NULL) +
        scale_x_continuous(labels = scales::percent, expand = c(0, 0)) + coord_cartesian(xlim = c(0, 1)) +
        labs(x = "share of variance", y = NULL) + theme_minimal(base_size = 18) + theme(legend.position = "top", plot.margin = margin(5, 25, 5, 5))
    })
    output$headline_note <- renderUI({
      missing <- setdiff(time_course_parents(design), summ$dataset_id)
      if (length(missing) == 0) return(NULL)
      missing_notice(paste("No variance results for:", paste(missing, collapse = ", ")))
    })

    output$pc1_text <- renderUI({
      m <- composition_markers()
      p(class = "text-muted",
        "Spearman correlation of each dataset's PC1 (of the expression matrix) with its neutrophil, lymphocyte",
        "and neutrophil - lymphocyte marker scores. A large value means PC1 mostly tracks leukocyte composition.",
        " A score is the mean, over the markers present in the dataset, of each marker's expression z-scored across samples.",
        if (!is.null(m)) tagList(
          tags$br(), strong("Neutrophil markers: "), paste(m$neutrophil, collapse = ", "), ".",
          tags$br(), strong("Lymphocyte markers: "), paste(m$lymphocyte, collapse = ", "), ".",
          " (T cell: CD3D/E/G, CD2, CD5, CD6, CD247, LCK; B cell: CD19, CD79A/B, MS4A1, CD22.)"))
    })

    output$pc1_table <- DT::renderDT({
      d <- need_data(pc1_composition_table(con), "composition_correlations")
      d <- d[order(-abs(d$nl)), ]
      out <- data.frame(Dataset = d$dataset, Neutrophil = d$neutrophil, Lymphocyte = d$lymphocyte,
                        `Neutrophil - lymphocyte` = d$nl, `p (N - L)` = signif(d$p_nl, 2), check.names = FALSE)
      DT::datatable(out, rownames = FALSE, options = list(dom = "t", pageLength = 20)) |>
        DT::formatRound(c("Neutrophil", "Lymphocyte", "Neutrophil - lymphocyte"), 2) |>
        DT::formatStyle("Neutrophil - lymphocyte", fontWeight = DT::styleInterval(c(-0.5, 0.5), c("bold", "normal", "bold")))
    })

    output$recurrence <- DT::renderDT({
      g <- need_data(time_gene_recurrence(con), "time_gene_recurrence")
      names(g) <- c("Gene", "Datasets", "Which")
      DT::datatable(g, rownames = FALSE, options = list(dom = "ft", paging = FALSE, scrollY = "480px", scrollCollapse = TRUE))
    })

    output$headline_table <- DT::renderDT({
      h <- headline()
      out <- data.frame(Dataset = h$dataset, `Between patients` = h$between, Time = h$time, Technical = h$technical,
                        Residual = h$residual, Repeatability = h$repeatability, `Time / (time + residual)` = h$time_within,
                        `Genes with time > 10%` = h$genes_time_gt10, `Subjects sampled repeatedly` = h$subjects_repeated,
                        Flags = h$flags, check.names = FALSE)
      DT::datatable(out, rownames = FALSE, options = list(dom = "t", pageLength = 20)) |>
        DT::formatPercentage(c("Between patients", "Time", "Technical", "Residual", "Repeatability",
                            "Time / (time + residual)", "Genes with time > 10%"), 1)
    })

  })
}

mod_variance_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h4("Variance partition by dataset"),
    p(class = "text-muted",
      "The evidence behind the overview for the dataset chosen in the sidebar."),
    uiOutput(ns("notices")),
    uiOutput(ns("selected_header")),
    navset_card_tab(
      nav_panel("Design and model",
        uiOutput(ns("design_text")),
        vcard("Differential expression contrasts", DT::DTOutput(ns("de")))),
      nav_panel("Genes",
        layout_columns(col_widths = c(6, 6), fill = FALSE,
          vcard("Shares by expression and by variance rank", plotOutput(ns("bins"), height = "420px", fill = FALSE)),
          vcard("Model terms: share of variance per gene", DT::DTOutput(ns("terms")))),
        vcard("Genes with the largest time effect", DT::DTOutput(ns("time_genes")))),
      nav_panel("Within-person change",
        uiOutput(ns("within_text")),
        vcard("Trajectories on the first three expression PCs, by group",
              p(class = "text-muted",
                "Each thin line is one person across timepoints; thick lines are the group mean at each timepoint.",
                "Groups that diverge over time, or move in opposite directions, show as separate thick lines."),
              plotOutput(ns("traj_plot"), height = "520px", fill = FALSE)),
        layout_columns(col_widths = c(6, 6), fill = FALSE,
          vcard("Spread of each gene's time share", plotOutput(ns("time_hist"), height = "420px", fill = FALSE)),
          vcard("Change summary", DT::DTOutput(ns("within_table"))))),
      nav_panel("Sample checks",
        layout_columns(col_widths = c(6, 6), fill = FALSE,
          vcard("Sex: recorded vs inferred from Y and XIST genes", uiOutput(ns("sex_text")),
               plotOutput(ns("sex_plot"), height = "420px", fill = FALSE)),
          vcard("PC1 vs neutrophil - lymphocyte score", uiOutput(ns("nl_text")),
               plotOutput(ns("nl_plot"), height = "420px", fill = FALSE))),
        layout_columns(col_widths = c(6, 6), fill = FALSE,
          vcard("Is composition a person trait? Variance of sample-level traits",
               plotOutput(ns("trait_plot"), height = "420px", fill = FALSE)),
          vcard("Correlations between sample-level traits (Spearman)", plotOutput(ns("cor_plot"), height = "420px", fill = FALSE)))))
  )
}

mod_variance_server <- function(id, state, con) {
  moduleServer(id, function(input, output, session) {
    ready <- variance_ready(con)
    design <- if (ready) design_table(con) else NULL
    summ <- if (ready) variance_summary_table(con) else NULL
    ok <- !is.null(design) && !is.null(summ) && nrow(summ) > 0

    # the dataset whose results apply to the sidebar selection
    parent <- reactive({
      req(ok)
      d <- design[design$dataset_id == state$dataset(), , drop = FALSE]
      if (nrow(d) == 0) { note_missing("a repeated-measures design", paste("dataset", state$dataset())); return(NULL) }
      d$parent_dataset_id
    })
    has_summary <- reactive(!is.null(parent()) && parent() %in% summ$dataset_id)
    row <- reactive(summ[summ$dataset_id == parent(), , drop = FALSE])

    output$notices <- renderUI({
      if (!ready || !ok)
        return(missing_notice("Variance results.", "The database has no repeated-measures tables (schema 2); rerun the pipeline."))
      tagList(
        if (!is.null(parent()) && !has_summary()) missing_notice(
          paste0(state$dataset(), " has no variance results."),
          "It has no design (no sample metadata, or a cross-sectional dataset)."),
        if (!is.null(parent()) && parent() != state$dataset()) p(class = "text-muted small",
          state$dataset(), " is a subset; the detail tabs show its parent, ", parent(), "."))
    })

    output$selected_header <- renderUI({
      req(has_summary())
      h5(class = "mt-3", "Evidence for ", strong(parent()),
         if (parent() != state$dataset()) span(class = "text-muted", " (parent of ", state$dataset(), ")"))
    })

    output$bins <- renderPlot({
      req(has_summary())
      b <- need_data(variance_bins_for(con, parent()), "variance_bins")
      long <- rbind(
        data.frame(binning = b$binning, bin = b$bin, part = "between patients", share = b$between_subject),
        data.frame(binning = b$binning, bin = b$bin, part = "time", share = b$time),
        data.frame(binning = b$binning, bin = b$bin, part = "residual", share = b$residual))
      long$binning <- ifelse(long$binning == "expression_quartile", "by expression (Q1 lowest)", "by variance rank")
      ggplot(long, aes(bin, share, fill = part)) + geom_col(position = "dodge") + facet_wrap(~binning, scales = "free_x") +
        scale_fill_manual(values = c(`between patients` = "#2C5F8A", time = "#E07B00", residual = "#BDBDBD"), name = NULL) +
        scale_y_continuous(labels = scales::percent) + labs(x = NULL, y = "median share per gene") +
        theme_minimal(base_size = 16) + theme(legend.position = "top")
    })
    output$terms <- DT::renderDT({
      req(has_summary())
      t <- need_data(variance_terms_for(con, parent()), "variance_terms")
      t <- t[, c("term", "category", "median", "q90", "share_gt05", "share_gt10")]
      names(t) <- c("Term", "Category", "Median", "90th pct", "Genes > 5%", "Genes > 10%")
      DT::datatable(t, rownames = FALSE, options = list(dom = "t", pageLength = 20)) |>
        DT::formatPercentage(c("Median", "90th pct", "Genes > 5%", "Genes > 10%"), 1)
    })
    output$time_genes <- DT::renderDT({
      req(has_summary())
      g <- need_data(time_genes_for(con, parent()), "variance_time_genes", "no gene has time above 10%")
      names(g) <- c("Gene", "Time share", "Patient share")
      DT::datatable(g, rownames = FALSE, options = list(dom = "tp", pageLength = 8)) |> DT::formatPercentage(2:3, 1)
    })

    output$within_text <- renderUI({
      req(has_summary()); r <- row()
      if (is.na(r$n_subjects_repeated)) return(missing_notice("No subject was sampled more than once in this dataset."))
      p("Each later sample minus the same person's first sample: ", strong(r$n_changes), " changes in ",
        strong(r$n_subjects_repeated), " people. If time were a shared programme, the first principal component of these",
        " changes would dominate; here it carries ", strong(sprintf("%.0f%%", 100 * r$change_pc1_share)),
        " of their variance, and tracks the neutrophil-lymphocyte change with correlation ",
        strong(sprintf("%.2f", r$change_pc1_cor_nl)), ".")
    })
    output$traj_plot <- renderPlot({
      req(has_summary())
      d <- need_data(samples(), "samples")
      d <- d[!is.na(d$timepoint) & !is.na(d$subject), , drop = FALSE]
      need_data(d, "samples with a subject and timepoint")
      lv <- tryCatch(unlist(jsonlite::fromJSON(design$time_order[design$dataset_id == parent()])), error = function(e) NULL)
      lv <- if (is.null(lv)) sort(unique(d$timepoint)) else intersect(lv, unique(d$timepoint))
      d$time <- factor(d$timepoint, levels = lv)
      d$group <- ifelse(is.na(d$grp), "no group", d$grp)
      long <- do.call(rbind, lapply(c("pc1", "pc2", "pc3"), function(k)
        data.frame(subject = d$subject, time = d$time, group = d$group, pc = toupper(k), value = d[[k]])))
      long <- long[!is.na(long$value) & !is.na(long$time), ]
      need_data(long, "PC scores for this dataset")
      mean_df <- aggregate(value ~ time + group + pc, long, mean)
      ggplot(long, aes(time, value, color = group)) +
        geom_line(aes(group = subject), alpha = 0.15, linewidth = 0.4) +
        geom_line(data = mean_df, aes(group = group), linewidth = 1.6) +
        geom_point(data = mean_df, size = 3) +
        facet_wrap(~pc, scales = "free_y", ncol = 1) +
        labs(x = NULL, y = "score", color = NULL) + theme_minimal(base_size = 16) +
        theme(legend.position = "top", axis.text.x = element_text(angle = 30, hjust = 1))
    })
    output$time_hist <- renderPlot({
      req(has_summary())
      g <- need_data(gene_variance_for(con, parent()), "variance_partition artifact")
      g <- g[g$factor_gene %in% TRUE, , drop = FALSE]
      ggplot(g, aes(time_all)) + geom_histogram(bins = 50, fill = "#E07B00") + scale_x_continuous(labels = scales::percent) +
        scale_y_sqrt() + labs(x = "time + group x time share of a gene's variance", y = "genes (sqrt scale)") + theme_minimal(base_size = 16)
    })
    output$within_table <- DT::renderDT({
      req(has_summary()); r <- row()
      v <- c(`Share of change on its first PC` = r$change_pc1_share, `Correlation of that PC with the N/L change` = r$change_pc1_cor_nl,
             `Median R2 of gene change vs composition change (factor genes)` = r$median_r2_nl_factor_genes,
             `Median R2 of gene change vs technical change (factor genes)` = r$median_r2_technical_factor_genes)
      DT::datatable(data.frame(Measure = names(v), Value = round(unname(v), 3)), rownames = FALSE, options = list(dom = "t"))
    })

    samples <- reactive({ req(has_summary()); samples_for(con, parent()) })
    output$sex_text <- renderUI({
      req(has_summary()); r <- row()
      if (is.na(r$sex_agreement)) return(missing_notice("No recorded sex for this dataset.", "Inferred sex is shown without a comparison."))
      cls <- if (r$sex_agreement < 0.9) "text-danger" else ""
      p(class = cls, "Agreement ", strong(sprintf("%.1f%%", 100 * r$sex_agreement)), "; ",
        strong(r$n_subjects_sex_discordant), " discordant subjects; separation ", sprintf("%.1f", r$sex_separation), ".",
        if (r$sex_agreement < 0.9) " The annotations do not match the arrays; the dataset is kept and flagged.")
    })
    output$sex_plot <- renderPlot({
      s <- need_data(samples(), "samples")
      need_data(s$xist[!is.na(s$xist)], "sex-gene expression for this dataset")
      s$recorded <- ifelse(is.na(s$recorded_sex), "not recorded", s$recorded_sex)
      ggplot(s, aes(xist, chrY_score, color = recorded, shape = inferred_sex)) + geom_point(alpha = 0.7) +
        labs(x = "XIST expression", y = "chrY score", color = "recorded", shape = "inferred") + theme_minimal(base_size = 16)
    })
    output$nl_text <- renderUI({
      req(has_summary()); r <- row()
      p("Spearman rho of PC1 with the score: ", strong(sprintf("%.2f", r$rho_pc1_nl)),
        if (!is.na(r$nl_subject_fraction)) paste0("; the score is ", sprintf("%.0f%%", 100 * r$nl_subject_fraction),
          " patient, ", sprintf("%.0f%%", 100 * r$nl_time_fraction), " time."),
        if (abs(r$rho_pc1_nl) > 0.5) " PC1 is largely leukocyte composition.")
    })
    output$nl_plot <- renderPlot({
      s <- need_data(samples(), "samples"); need_data(s$nl_score[!is.na(s$nl_score)], "composition scores")
      ggplot(s, aes(nl_score, pc1, color = timepoint)) + geom_point(alpha = 0.7) +
        labs(x = "neutrophil - lymphocyte score", y = "PC1", color = NULL) + theme_minimal(base_size = 16)
    })
    output$trait_plot <- renderPlot({
      req(has_summary())
      t <- need_data(trait_variance_for(con, parent()), "trait_variance")
      t <- t[t$trait %in% c("nl_score", "neutrophil_score", "lymphocyte_score", "PC1", "PC2", "PC3", "expr_mean", "expr_sd") &
               t$category %in% names(VARIANCE_COLORS), ]
      need_data(t, "trait_variance")
      t$category <- factor(t$category, levels = names(VARIANCE_COLORS))
      ggplot(t, aes(fraction, trait, fill = category)) + geom_col(position = position_stack(reverse = TRUE)) +
        scale_fill_manual(values = VARIANCE_COLORS, labels = VARIANCE_LABELS[names(VARIANCE_COLORS)], drop = TRUE, name = NULL) +
        scale_x_continuous(labels = scales::percent) + labs(x = "share of the trait's variance", y = NULL) + theme_minimal(base_size = 16)
    })
    output$cor_plot <- renderPlot({
      req(has_summary())
      cc <- need_data(composition_correlations_for(con, parent()), "composition_correlations")
      lv <- unique(c(cc$x, cc$y))
      ggplot(cc, aes(x, y, fill = rho)) + geom_tile() + scale_fill_gradient2(limits = c(-1, 1), low = "#3B6FB6", high = "#C0392B", name = "rho") +
        scale_x_discrete(limits = lv) + scale_y_discrete(limits = rev(lv)) +
        labs(x = NULL, y = NULL) + theme_minimal(base_size = 16) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    })

    output$design_text <- renderUI({
      req(has_summary())
      d <- design[design$dataset_id == parent(), ]; r <- row()
      tagList(
        p(strong(d$n_samples), " samples, ", strong(d$n_subjects), " subjects (", d$n_dropped, " dropped for missing model variables); platform ", d$platform, "."),
        p("Variance model: ", tags$code(gsub("\\s+", " ", d$vp_formula))),
        p("Differential expression model: ", tags$code(gsub("\\s+", " ", d$dream_formula))),
        p(class = "text-muted small", sprintf("%.0f%% of genes had a singular fit; maximum collinearity %.2f; genome-wide subject correlation (dream) %.2f.",
          100 * r$share_singular, r$max_colinearity, r$dupcor_consensus)),
        if (!is.na(r$control_factor_fraction)) p(class = "text-muted small", sprintf(
          "Array control probes enter as a covariate (their first PC, %.0f%% of the controls' variance).", 100 * r$control_factor_fraction)))
    })
    output$de <- DT::renderDT({
      req(has_summary())
      d <- need_data(de_contrasts_for(con, parent()), "de_contrasts")
      names(d) <- c("Contrast", "Kind", "Comparison", "Genes tested", "FDR < 5%")
      DT::datatable(d, rownames = FALSE, options = list(dom = "tp", pageLength = 10))
    })
  })
}

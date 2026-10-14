# Remaining analyses: Biplots, Trajectory across time, Projection comparison,
# and Pattern drivers. `focus` and `state` as in mod_inspect.R.

component_choices <- function(con, focus) {
  fac <- fit_factors(con, focus$fit_id, focus$n_factors)
  setNames(seq_len(nrow(fac)), fac$label)
}

## ---- Biplots -----------------------------------------------------------------------------------

mod_biplot_ui <- function(id) {
  ns <- NS(id)
  tagList(
    layout_columns(col_widths = c(3, 3, 3, 3), fill = FALSE, selectInput(ns("x"), "x axis", choices = NULL), selectInput(ns("y"), "y axis", choices = NULL),
                   sliderInput(ns("top_n"), "Genes shown", 5, 50, 15, step = 5), uiOutput(ns("note"))),
    navset_card_tab(
      nav_panel("Genes",
        p(class = "text-muted", "Samples at their scores on the two factors, with arrows for the genes that load most on them: direction and relative length matter, not absolute scale."),
        plotOutput(ns("genes"), height = "560px")),
      nav_panel("Metadata",
        p(class = "text-muted", "The same sample cloud with every usable metadata field: numeric fields as a Spearman-correlation vector, categorical fields as a centroid per level."),
        plotOutput(ns("meta"), height = "560px"))))
}

mod_biplot_server <- function(id, state, focus, con) {
  moduleServer(id, function(input, output, session) {
    scores <- reactive(focus_scores(con, focus()))
    observe({
      ch <- component_choices(con, focus()); cur <- isolate(input$x)
      updateSelectInput(session, "x", choices = ch, selected = if (!is.null(cur) && cur %in% ch) cur else ch[1])
      updateSelectInput(session, "y", choices = ch, selected = if (length(ch) > 1) ch[2] else ch[1])
    })
    output$note <- renderUI(p(class = "text-muted mt-4", focus_label(focus())))
    axes <- reactive({ req(input$x, input$y); S <- scores(); validate(need(!is.null(S), "No scores for this fit.")); c(as.integer(input$x), as.integer(input$y)) })
    output$genes <- renderPlot({
      f <- focus(); validate(need(f$method != "wgcna", "WGCNA modules have no loadings, so there are no gene arrows; the Metadata tab works on module eigengenes."))
      ax <- axes(); S <- scores(); L <- load_loadings(con, f$fit_id, f$n_factors); validate(need(!is.null(L), "No loadings for this fit."))
      d <- biplot_data(L, S, ax[1], ax[2], input$top_n); d$arrows$gene <- label_genes(gene_labels(con, state$dataset()), d$arrows$gene)
      ggplot() + geom_point(data = d$scores, aes(x, y), alpha = 0.5, colour = "grey40") +
        geom_segment(data = d$arrows, aes(0, 0, xend = x, yend = y), colour = "#B2182B", arrow = arrow(length = unit(0.2, "cm"))) +
        ggrepel::geom_text_repel(data = d$arrows, aes(x, y, label = gene), colour = "#B2182B", size = 3.5, max.overlaps = Inf) +
        labs(x = colnames(S)[ax[1]], y = colnames(S)[ax[2]]) + theme_minimal(base_size = 12)
    })
    output$meta <- renderPlot({
      ax <- axes(); S <- scores(); m <- dataset_sample_metadata(con, state$dataset()); validate(need(!is.null(m), "No sample metadata."))
      d <- biplot_metadata_data(S, m, ax[1], ax[2]); validate(need(!is.null(d), "Too few samples shared with the metadata."))
      p <- ggplot() + geom_point(data = d$scores, aes(x, y), alpha = 0.4, colour = "grey60")
      if (nrow(d$vectors) > 0) p <- p + geom_segment(data = d$vectors, aes(0, 0, xend = x, yend = y), colour = "#2166AC", arrow = arrow(length = unit(0.2, "cm"))) +
        ggrepel::geom_text_repel(data = d$vectors, aes(x, y, label = field), colour = "#2166AC", size = 3.5, max.overlaps = Inf)
      if (nrow(d$centroids) > 0) p <- p + geom_point(data = d$centroids, aes(x, y, colour = field), size = 3) +
        ggrepel::geom_text_repel(data = d$centroids, aes(x, y, label = level, colour = field), size = 3.5, max.overlaps = Inf, show.legend = FALSE)
      p + labs(x = colnames(S)[ax[1]], y = colnames(S)[ax[2]], colour = "categorical field") + theme_minimal(base_size = 12)
    })
  })
}

## ---- Trajectory across time ---------------------------------------------------------------------

mod_trajectory_ui <- function(id) {
  ns <- NS(id)
  tagList(
    layout_columns(col_widths = c(3, 3, 6), fill = FALSE, selectInput(ns("component"), "Factor", choices = NULL),
                   selectInput(ns("colour"), "Colour subjects by", choices = "none"), uiOutput(ns("note"))),
    p(class = "text-muted", "The fit's basis projected (projectR) onto each sibling per-timepoint dataset of the family; each line is one subject's weight on the factor across the time course.",
      "Only available for datasets whose family has per-timepoint siblings."),
    plotOutput(ns("traj"), height = "520px"))
}

mod_trajectory_server <- function(id, state, focus, con) {
  moduleServer(id, function(input, output, session) {
    observe({
      ch <- component_choices(con, focus()); cur <- isolate(input$component)
      updateSelectInput(session, "component", choices = ch, selected = if (!is.null(cur) && cur %in% ch) cur else ch[1])
    })
    fields <- reactive(subject_level_fields(con, state$dataset()))
    observe({
      f <- fields(); updateSelectInput(session, "colour", choices = c("none", if (!is.null(f)) setdiff(names(f), "subject_id")))
    })
    output$note <- renderUI({
      n <- nrow(within_projections(con, focus()$fit_id))
      p(class = "text-muted mt-4", focus_label(focus()), ": ", n, " within-family projection", if (n != 1) "s")
    })
    output$traj <- renderPlot({
      f <- focus(); req(input$component)
      validate(need(f$method != "wgcna", "WGCNA modules are not projected onto other datasets, so there is no trajectory. Choose another method in the sidebar."))
      validate(need(nrow(within_projections(con, f$fit_id)) > 0, "This fit has no within-family projections (needs a dataset with per-timepoint siblings)."))
      label <- names(component_choices(con, f))[as.integer(input$component)]
      d <- trajectory_data(con, f$fit_id, label)
      validate(need(!is.null(d), "No usable projections: the timepoint datasets need a subject id column."))
      lev <- unique(d$dataset_id); d$dataset_id <- factor(d$dataset_id, levels = lev[order_timepoints(lev)])
      fl <- fields(); col <- input$colour
      if (!is.null(fl) && !identical(col, "none") && col %in% names(fl)) d$colour <- fl[[col]][match(d$subject_id, fl$subject_id)]
      p <- ggplot(d, aes(dataset_id, weight, group = subject_id))
      p <- if ("colour" %in% names(d)) p + geom_line(aes(colour = if (is.numeric(colour)) colour else factor(colour)), alpha = 0.5) + geom_point(aes(colour = if (is.numeric(colour)) colour else factor(colour)), alpha = 0.5, size = 1.5) + labs(colour = col)
           else p + geom_line(alpha = 0.3, colour = "#4C78A8") + geom_point(alpha = 0.4, colour = "#4C78A8", size = 1.5)
      p + labs(x = NULL, y = paste(label, "projected weight"),
               subtitle = paste(length(unique(d$subject_id)), "subjects across", length(lev), "timepoints")) + theme_minimal(base_size = 12)
    })
  })
}

## ---- Projection comparison -------------------------------------------------------------------------

mod_projection_ui <- function(id) {
  ns <- NS(id)
  tagList(
    layout_columns(col_widths = c(3, 3, 6), fill = FALSE, selectInput(ns("target"), "Project onto dataset", choices = NULL),
                   selectInput(ns("fit_b"), "Compare with (Fit B)", choices = NULL), uiOutput(ns("note"))),
    navset_card_tab(
      nav_panel("R² by method",
        p(class = "text-muted", "How much of the target's variance each method's basis explains when projected in (mean per-sample R²). Grey: every fit; coloured: representative fits. sPCA runs lower by construction: a sparse basis reconstructs new data worse."),
        plotOutput(ns("r2"), height = "420px")),
      nav_panel("Sample-space alignment",
        p(class = "text-muted", "Correlation between Fit A's and Fit B's projected per-sample weights on the target, samples matched by id."),
        plotOutput(ns("samples"), height = "480px")),
      nav_panel("Gene-space alignment",
        p(class = "text-muted", "Fit A's loadings projected (least squares, as projectR does) onto Fit B's basis on shared genes. Axes are reordered by best one-to-one match, since factor numbering is arbitrary for most methods."),
        plotOutput(ns("genes"), height = "520px"))))
}

mod_projection_server <- function(id, state, focus, con) {
  moduleServer(id, function(input, output, session) {
    observe({
      t <- projection_targets(con, state$dataset()); cur <- isolate(input$target)
      updateSelectInput(session, "target", choices = t, selected = if (!is.null(cur) && cur %in% t) cur else t[1])
    })
    pair <- reactive({ req(input$target); d <- projections_for_pair(con, state$dataset(), input$target); validate(need(nrow(d) > 0, "No projections of this dataset's fits onto that target.")); d })
    reps <- reactive({
      d <- pair(); d <- d[d$representative == 1, , drop = FALSE]; d <- d[order(match(d$method, METHOD_ORDER), d$rank), ]
      d$label <- paste(method_label(d$method), "rank", d$rank); d
    })
    observe({
      d <- reps(); d <- d[d$source_fit_id != focus()$fit_id, , drop = FALSE]; ch <- setNames(d$source_fit_id, d$label)
      cur <- isolate(input$fit_b)
      updateSelectInput(session, "fit_b", choices = ch, selected = if (!is.null(cur) && cur %in% ch) cur else ch[1])
    })
    output$note <- renderUI(p(class = "mt-4", "Fit A: ", tags$strong(state$dataset(), " / ", focus_label(focus())), " (sidebar). Both fits are from this dataset."))
    output$r2 <- renderPlot({
      d <- pair(); d$m <- factor(method_label(d$method), levels = method_label(METHOD_ORDER))
      ggplot(d, aes(m, mean_r_squared)) + geom_boxplot(outlier.shape = NA, fill = "grey90") +
        geom_jitter(data = d[d$representative != 1, ], width = 0.15, alpha = 0.3, colour = "grey40") +
        geom_jitter(data = d[d$representative == 1, ], aes(colour = method), width = 0.15, size = 2) +
        scale_colour_manual(values = METHOD_COLORS, guide = "none") +
        labs(x = NULL, y = "mean R²", subtitle = paste(state$dataset(), "→", input$target)) + theme_minimal(base_size = 12)
    })
    proj_of <- function(id, n = NULL) {
      d <- pair(); p <- d$path[d$source_fit_id == as.integer(id)][1]
      P <- tryCatch(readRDS(resolve_artifact(p))$projection, error = function(e) NULL)
      if (!is.null(P) && !is.null(n)) P <- P[seq_len(min(n, nrow(P))), , drop = FALSE]
      P
    }
    fit_a <- reactive({
      f <- focus(); validate(need(f$method != "wgcna", "WGCNA modules are not projected; choose another method in the sidebar."))
      validate(need(f$fit_id %in% pair()$source_fit_id, "The sidebar fit has no projection onto this target."))
      f
    })
    output$samples <- renderPlot({
      req(input$fit_b); fa <- fit_a(); validate(need(fa$fit_id != as.integer(input$fit_b), "Pick two different fits."))
      Pa <- proj_of(fa$fit_id, if (fa$method == "pca") fa$n_factors); Pb <- proj_of(input$fit_b); validate(need(!is.null(Pa) && !is.null(Pb), "Could not read a projection artifact."))
      Sa <- t(Pa); Sb <- t(Pb); shared <- intersect(rownames(Sa), rownames(Sb)); validate(need(length(shared) >= 3, "Too few shared samples."))
      cm <- cor(Sa[shared, , drop = FALSE], Sb[shared, , drop = FALSE])
      la <- focus_label(fa); lb <- reps()$label[reps()$source_fit_id == as.integer(input$fit_b)]
      heatmap_plot(cm, c(-1, 1), "correlation") + labs(x = paste(lb, "(B)"), y = paste(la, "(A)"), subtitle = paste(length(shared), "shared samples"))
    })
    output$genes <- renderPlot({
      req(input$fit_b); fa <- fit_a(); validate(need(fa$fit_id != as.integer(input$fit_b), "Pick two different fits."))
      La <- load_loadings(con, fa$fit_id, if (fa$method == "pca") fa$n_factors); Lb <- load_loadings(con, as.integer(input$fit_b))
      validate(need(!is.null(La) && !is.null(Lb), "One fit has no loadings."))
      pr <- loadings_projection(La, Lb); validate(need(!is.null(pr), "Too few shared genes."))
      P <- pr$projection; idx <- hungarian_match_abs(P); ord <- order(idx[, "a"])
      P <- P[c(idx[ord, "a"], setdiff(seq_len(nrow(P)), idx[, "a"])), c(idx[ord, "b"], setdiff(seq_len(ncol(P)), idx[, "b"])), drop = FALSE]
      la <- focus_label(fa); lb <- reps()$label[reps()$source_fit_id == as.integer(input$fit_b)]
      heatmap_plot(P, c(-1, 1), "coefficient") + labs(x = paste(la, "(A, projected)"), y = paste(lb, "(B, basis)"), subtitle = paste(pr$n_genes, "shared genes"))
    })
  })
}

## ---- Pattern drivers ------------------------------------------------------------------------------------

mod_drivers_ui <- function(id) {
  ns <- NS(id)
  tagList(
    layout_columns(col_widths = c(3, 5, 2, 2), fill = FALSE, selectInput(ns("factor"), "Factor", choices = NULL), selectInput(ns("cmp"), "Comparison", choices = NULL),
                   selectInput(ns("mode"), "Test", choices = c("Confidence intervals" = "CI", "Welch tests" = "PV")), uiOutput(ns("note"))),
    p(class = "text-muted", "projectR's projectionDriveR: which genes differ between two sample groups AND carry weight on the factor. Only the combinations computed by the pipeline are listed (factors of the representative fits, categorical fields with few levels, combinations with at least one shared significant gene)."),
    layout_columns(col_widths = c(7, 5), plotOutput(ns("plot"), height = "520px"), DTOutput(ns("table"))))
}

mod_drivers_server <- function(id, state, focus, con) {
  moduleServer(id, function(input, output, session) {
    drivers <- reactive({ f <- focus(); fit_drivers(con, f$fit_id, f$n_factors) })
    output$note <- renderUI(p(class = "text-muted mt-4", focus_label(focus())))
    observe({
      d <- drivers(); fac <- unique(d[c("factor_index", "label")])
      cur <- isolate(input$factor)
      updateSelectInput(session, "factor", choices = setNames(fac$factor_index, fac$label), selected = if (!is.null(cur) && cur %in% fac$factor_index) cur else fac$factor_index[1])
    })
    observe({
      d <- drivers(); req(input$factor); d <- d[d$factor_index == as.numeric(input$factor) & d$mode == input$mode, ]
      if (nrow(d) == 0) return(updateSelectInput(session, "cmp", choices = character(0)))
      key <- paste(d$grouping_col, d$group1_level, d$group2_level, sep = "|")
      updateSelectInput(session, "cmp", choices = setNames(key, paste0(d$grouping_col, ": ", d$group1_level, " vs ", d$group2_level)))
    })
    result <- reactive({
      d <- drivers()
      validate(need(nrow(d) > 0, if (focus()$method == "wgcna") "No pattern drivers are computed for WGCNA modules. Choose another method in the sidebar." else "No pattern drivers stored for this fit."))
      req(input$factor, input$cmp)
      p <- strsplit(input$cmp, "|", fixed = TRUE)[[1]]
      r <- d[d$factor_index == as.numeric(input$factor) & d$mode == input$mode & d$grouping_col == p[1] & d$group1_level == p[2] & d$group2_level == p[3], ]
      validate(need(nrow(r) >= 1, "No stored result for this combination."))
      res <- load_driver(r$path[1]); validate(need(!is.null(res), "The result file is missing."))
      list(res = res, row = r[1, ])
    })
    labels <- reactive(gene_labels(con, state$dataset()))
    shared_genes <- function(res) res$sig_genes$significant_shared_genes %||% res$sig_genes$PV_significant_shared_genes
    output$plot <- renderPlot({
      x <- result(); res <- x$res; lab <- labels()
      if (input$mode == "CI") {
        sg <- shared_genes(res); validate(need(length(sg) > 0, "No shared significant genes."))
        ci <- res$weighted_mean_ci; ci <- ci[ci$gene %in% sg, ]
        w <- res$normalized_weights[ci$gene]; ci$mid <- (ci$low + ci$high) / 2
        ci <- head(ci[order(-abs(ci$mid)), ], 25); ci$g <- factor(label_genes(lab, ci$gene), levels = rev(label_genes(lab, ci$gene)))
        ggplot(ci, aes(y = g)) + geom_linerange(aes(xmin = low, xmax = high), colour = "#4C78A8", linewidth = 1) + geom_point(aes(mid), size = 2) +
          geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
          labs(x = "weighted mean difference (confidence interval)", y = NULL, subtitle = paste(x$row$group1_level, "vs", x$row$group2_level, "-", length(sg), "shared significant genes")) +
          theme_minimal(base_size = 12)
      } else {
        ms <- res$weighted_mean_stats; sg <- shared_genes(res)
        ms$sig <- ms$gene %in% sg; ms$logp <- -log10(pmax(ms$welch_padj, 1e-300))
        ggplot(ms, aes(mean_diff, logp, colour = sig)) + geom_point(alpha = 0.5) +
          ggrepel::geom_text_repel(data = ms[ms$sig, ], aes(label = label_genes(lab, gene)), size = 3.5, max.overlaps = Inf, show.legend = FALSE) +
          scale_colour_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey60"), labels = c(`TRUE` = "shared significant", `FALSE` = "other"), name = NULL) +
          labs(x = "weighted mean difference", y = expression(-log[10]~adjusted~p)) + theme_minimal(base_size = 12) + theme(legend.position = "top")
      }
    })
    output$table <- renderDT({
      x <- result(); res <- x$res; lab <- labels(); sg <- shared_genes(res); validate(need(length(sg) > 0, "No shared significant genes."))
      d <- data.frame(Gene = label_genes(lab, sg), Weight = round(res$normalized_weights[sg], 3))
      datatable(d[order(-abs(d$Weight)), ], rownames = FALSE, selection = "none", options = list(dom = "ftp", pageLength = 12))
    })
  })
}

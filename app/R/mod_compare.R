# The Compare page: one tab per question, all reading the sidebar state.
#   1. optimal hyperparameter per method      (mod_hyperparam)
#   2. best method for the dataset            (mod_best_method)
#   3. similarity of optimal factors, methods (mod_cross_method)
#   4. similarity of factors across datasets  (mod_cross_dataset)
# `state` is a list of reactives built in app.R: dataset, metric, methods,
# curves (all present methods), optimal (one fit per selected method).

METHOD_COLORS <- c(pca = "#4C78A8", spca = "#72B7B2", ica = "#F58518",
                   nmf = "#54A24B", cogaps = "#B279A2", wgcna = "#E45756")
METHOD_LABELS <- c(pca = "PCA", spca = "sPCA", ica = "ICA", nmf = "NMF", cogaps = "CoGAPS", wgcna = "WGCNA")

method_label <- function(m) unname(METHOD_LABELS[m] %||% m)

heatmap_plot <- function(m, limits, name, current = NULL, digits = 2) {
  mark <- function(x) ifelse(x %in% current, paste0("\u25b6 ", x), x)
  dimnames(m) <- list(mark(rownames(m)), mark(colnames(m)))
  d <- as.data.frame(as.table(m), stringsAsFactors = FALSE)
  names(d) <- c("row", "col", "value")
  d$row <- factor(d$row, levels = rev(rownames(m))); d$col <- factor(d$col, levels = colnames(m))
  ggplot(d, aes(col, row, fill = value)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = ifelse(is.na(value), "", formatC(value, digits = digits, format = "f"))), size = 3) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = mean(limits),
                         limits = limits, oob = scales::squish, na.value = "grey92", name = name) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
}

## ---- 1. optimal hyperparameter -----------------------------------------------

mod_hyperparam_ui <- function(id) {
  ns <- NS(id)
  tagList(
    p(class = "text-muted",
      "Each panel is one method's representative fits across its hyperparameter. The ringed point is the",
      "optimal fit under the rule in the sidebar. Quality is comparable within a method, not across them."),
    plotOutput(ns("curves"), height = "420px"),
    DTOutput(ns("table")))
}

mod_hyperparam_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {
    curves <- reactive({
      cv <- state$curves()
      cv[cv$method %in% state$methods(), , drop = FALSE]
    })
    output$curves <- renderPlot({
      cv <- curves(); req(nrow(cv) > 0)
      opt <- state$optimal()
      long <- rbind(
        data.frame(cv[c("method", "rank_key", "fit_id")], panel = "in-sample quality", value = cv$in_sample_r2),
        data.frame(cv[c("method", "rank_key", "fit_id")], panel = paste("stability |", state$metric()), value = cv$stability))
      ringed <- long[paste(long$method, long$rank_key) %in% paste(opt$method, opt$rank_key), ]
      strip <- ifelse(long$method == "wgcna", "WGCNA (x = power)", method_label(long$method))
      lv <- ifelse(METHOD_ORDER == "wgcna", "WGCNA (x = power)", method_label(METHOD_ORDER))
      long$method <- factor(strip, levels = lv)
      ringed$method <- factor(ifelse(ringed$method == "wgcna", "WGCNA (x = power)", method_label(ringed$method)), levels = lv)
      ggplot(long, aes(rank_key, value, colour = method)) +
        geom_line(na.rm = TRUE) + geom_point(na.rm = TRUE, size = 1.8) +
        geom_point(data = ringed, shape = 21, size = 5, stroke = 1.2, colour = "black", fill = NA, na.rm = TRUE) +
        facet_grid(panel ~ method, scales = "free") +
        scale_colour_manual(values = setNames(METHOD_COLORS[METHOD_ORDER], lv), guide = "none") +
        labs(x = "rank", y = NULL) +
        theme_minimal(base_size = 12) + theme(strip.text = element_text(face = "bold"))
    })
    output$table <- renderDT({
      o <- state$optimal(); req(nrow(o) > 0)
      d <- data.frame(Method = method_label(o$method), `Optimal fit` = o$rank_label, Factors = o$n_factors,
                      `In-sample quality` = round(o$in_sample_r2, 3), Stability = round(o$stability, 3),
                      check.names = FALSE)
      datatable(d, rownames = FALSE, options = list(dom = "t", paging = FALSE), selection = "none")
    })
  })
}

## ---- 2. best method -----------------------------------------------------------

mod_best_method_ui <- function(id) {
  ns <- NS(id)
  tagList(
    radioButtons(ns("quality"), "Quality axis", inline = TRUE,
                 choices = c("In-sample" = "in_sample_r2", "Out-of-sample, within dataset family" = "oos_r2_within",
                             "Out-of-sample, other datasets" = "oos_r2_cross")),
    p(class = "text-muted", "Each point is a method's optimal fit. Top-right is best: high quality and reproducible across seeds.",
      "Out-of-sample quality is the median relative transfer R² of that method's fits on other samples (the share of what the target's own principal components explain, after removing gene means); WGCNA is never projected."),
    plotOutput(ns("scatter"), height = "380px"),
    DTOutput(ns("table")))
}

mod_best_method_server <- function(id, state, con) {
  moduleServer(id, function(input, output, session) {
    summary_df <- reactive({
      o <- state$optimal(); req(nrow(o) > 0)
      oos <- oos_r2_by_method(con, state$dataset())
      o$oos_r2_within <- oos$oos_r2_within[match(o$method, oos$method)]
      o$oos_r2_cross <- oos$oos_r2_cross[match(o$method, oos$method)]
      o
    })
    output$scatter <- renderPlot({
      d <- summary_df(); q <- input$quality
      d$q <- d[[q]]
      shown <- d[!is.na(d$q) & !is.na(d$stability), ]
      validate(need(nrow(shown) > 0, "No methods have both values on these axes for this dataset."))
      ggplot(shown, aes(q, stability, colour = method, label = method_label(method))) +
        geom_point(size = 4) + ggrepel::geom_text_repel(show.legend = FALSE, size = 4.2) +
        scale_colour_manual(values = METHOD_COLORS, guide = "none") +
        labs(x = names(which(c("In-sample" = "in_sample_r2", "Out-of-sample, within dataset family" = "oos_r2_within",
                               "Out-of-sample, other datasets" = "oos_r2_cross") == q)),
             y = paste("Stability |", state$metric())) +
        theme_minimal(base_size = 13)
    })
    output$table <- renderDT({
      d <- summary_df()
      out <- data.frame(Method = method_label(d$method), `Optimal fit` = d$rank_label,
                        `In-sample` = round(d$in_sample_r2, 3), Stability = round(d$stability, 3),
                        `Out-of-sample (family)` = round(d$oos_r2_within, 3),
                        `Out-of-sample (other)` = round(d$oos_r2_cross, 3), check.names = FALSE)
      datatable(out, rownames = FALSE, options = list(dom = "t", paging = FALSE), selection = "none")
    })
  })
}

## ---- 3. cross-method similarity ---------------------------------------------------

mod_cross_method_ui <- function(id) {
  ns <- NS(id)
  tagList(
    layout_columns(col_widths = c(6, 6),
      radioButtons(ns("view"), "Agreement measure", inline = TRUE,
                   choices = c("Gene sets (top-N Jaccard)" = "gene_jaccard", "Sample clusters (ARI)" = "sample_ari")),
      sliderInput(ns("top_n"), "Top genes per factor", min = 10, max = 200, value = 50, step = 10)),
    p(class = "text-muted",
      "Optimal fits compared pairwise. Gene-set agreement matches factors one-to-one (Hungarian, on |Jaccard|) and",
      "averages the matched overlaps; WGCNA modules are compared by full membership. Sample agreement clusters samples by each",
      "fit's scores (k = 4) and compares the partitions."),
    layout_columns(col_widths = c(5, 7),
      plotOutput(ns("matrix"), height = "380px", click = ns("pick")),
      tagList(
        uiOutput(ns("pair_title")),
        plotOutput(ns("pair_heat"), height = "300px"),
        DTOutput(ns("pair_table")))))
}

mod_cross_method_server <- function(id, state, con) {
  moduleServer(id, function(input, output, session) {
    agree <- reactive({
      o <- state$optimal()
      validate(need(nrow(o) >= 2, "Select at least two methods with fits for this dataset."))
      withProgress(message = "Comparing optimal fits", value = 0.5,
                   method_pairwise_agreement(con, o, input$top_n))
    }) |> bindCache(state$dataset(), state$optimal()$fit_id, input$top_n)
    output$matrix <- renderPlot({
      a <- agree(); m <- a[[input$view]]
      rownames(m) <- colnames(m) <- method_label(a$methods)
      heatmap_plot(m, c(0, 1), if (input$view == "gene_jaccard") "mean matched\nJaccard" else "ARI")
    })
    selected <- reactiveVal(NULL)
    observeEvent(input$pick, {
      a <- agree(); n <- length(a$methods)
      i <- round(input$pick$y); j <- round(input$pick$x)
      i <- n + 1 - i                      # rows are drawn top-down
      if (i >= 1 && i <= n && j >= 1 && j <= n && i != j) selected(sort(c(i, j)))
    })
    pair <- reactive({
      a <- agree(); s <- selected()
      if (is.null(s) && length(a$pairs) > 0) s <- match(strsplit(names(a$pairs)[1], "|", fixed = TRUE)[[1]], a$methods)
      req(s)
      key <- paste(a$methods[s], collapse = "|")
      list(methods = a$methods[s], ag = a$pairs[[key]])
    })
    output$pair_title <- renderUI({
      p <- pair()
      h5(paste(method_label(p$methods), collapse = " vs "), class = "mt-2")
    })
    output$pair_heat <- renderPlot({
      ag <- pair()$ag; validate(need(!is.null(ag), "No comparable factors for this pair."))
      heatmap_plot(ag$matrix, c(0, 1), "Jaccard")
    })
    output$pair_table <- renderDT({
      ag <- pair()$ag; req(ag)
      d <- ag$matches; d$jaccard <- round(d$jaccard, 3)
      names(d) <- c(method_label(pair()$methods[1]), method_label(pair()$methods[2]), "Jaccard")
      datatable(d, rownames = FALSE, options = list(dom = "tp", pageLength = 8), selection = "none")
    })
  })
}

## ---- 4. cross-dataset similarity --------------------------------------------------

mod_cross_dataset_ui <- function(id) {
  ns <- NS(id)
  tagList(
    radioButtons(ns("kind"), "Similarity", inline = TRUE,
                 choices = c("Sample space: relative transfer R²" = "r2", "Gene space: loadings agreement" = "gene_space")),
    p(class = "text-muted",
      "The sidebar's method, every pair of datasets, ordered by clustering; the selected dataset is marked \u25b6. relative R² is how much of the variation in the other",
      "dataset's samples (after removing each gene's mean) one dataset's factors explain, as a share of what that dataset's own top",
      "principal components explain; gene-space agreement is the mean |matched similarity| between the datasets'",
      "representative loadings on shared genes. Blank cells have no projection or comparison."),
    layout_columns(col_widths = c(7, 5),
      plotOutput(ns("heat"), height = "560px"),
      tagList(h5("Most similar to the selected dataset"), DTOutput(ns("nearest")))))
}

mod_cross_dataset_server <- function(id, state, focus, con) {
  moduleServer(id, function(input, output, session) {
    method <- reactive(focus()$method)
    long <- reactive({
      validate(need(method() != "wgcna", "WGCNA modules are not projected onto other datasets, so there is no cross-dataset similarity for them. Choose another method in the sidebar."))
      cross_dataset_similarity(con, method(), input$kind)
    }) |> bindCache(method(), input$kind)
    output$heat <- renderPlot({
      m <- similarity_matrix(long())
      validate(need(!is.null(m), "No cross-dataset results for this method."))
      lim <- c(0, 1)
      heatmap_plot(m, lim, if (input$kind == "r2") "relative R²" else "mean |diag|",
                   current = state$dataset(), digits = 1) +
        theme(axis.text = element_text(size = 8))
    })
    output$nearest <- renderDT({
      d <- long(); ds <- state$dataset()
      out <- d[d$source == ds, ]; inn <- d[d$target == ds, ]
      validate(need(nrow(out) + nrow(inn) > 0, "This dataset has no cross-dataset results for this method."))
      agg <- data.frame(dataset = union(out$target, inn$source), stringsAsFactors = FALSE)
      agg$to_it <- inn$value[match(agg$dataset, inn$source)]
      agg$from_it <- out$value[match(agg$dataset, out$target)]
      agg$mean <- rowMeans(agg[c("to_it", "from_it")], na.rm = TRUE)
      agg <- agg[order(-agg$mean), ]
      names(agg) <- c("Dataset", "Others → this", "This → others", "Mean")
      datatable(agg, rownames = FALSE, selection = "none", options = list(dom = "tp", pageLength = 12)) |>
        formatRound(2:4, 3)
    })
  })
}

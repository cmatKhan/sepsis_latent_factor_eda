# Pages that inspect one fit: Factor detail, Enrichment, Metadata associations,
# and the WGCNA module views. `focus` is a reactive list(method, fit_id,
# n_factors, rank_label) built in app.R from the sidebar's "Fit to inspect";
# `state` is the shared sidebar state (mod_compare.R).

focus_label <- function(f) paste(method_label(f$method), f$rank_label)

#' Sample scores of the focus fit with columns named by factor label (WGCNA's
#' eigengene columns are in a different order from factor_index).
focus_scores <- function(con, focus) {
  S <- load_scores(con, focus$fit_id)
  if (is.null(S)) return(NULL)
  fac <- fit_factors(con, focus$fit_id, focus$n_factors)
  if (focus$method == "wgcna") {
    keep <- intersect(fac$label, colnames(S))
    S[, keep, drop = FALSE]
  } else {
    S <- S[, seq_len(min(nrow(fac), ncol(S))), drop = FALSE]
    colnames(S) <- fac$label[seq_len(ncol(S))]
    S
  }
}

module_palette <- function(modules) {
  mods <- sort(unique(modules[modules != 0]))
  cols <- setNames(scales::hue_pal(l = 60, c = 80)(max(length(mods), 1))[seq_along(mods)], mods)
  c(`0` = "grey80", cols)
}

## ---- enrichment panel (shared by Factor detail and the Enrichment page) ----------------

enrichment_panel_ui <- function(id) {
  ns <- NS(id)
  tagList(
    layout_columns(col_widths = c(4, 4, 4),
      selectInput(ns("query"), "Analysis", choices = NULL),
      selectInput(ns("source"), "Gene-set collection", choices = "all"),
      sliderInput(ns("top"), "Terms shown", min = 5, max = 40, value = 15, step = 5)),
    plotOutput(ns("dots"), height = "380px"),
    DTOutput(ns("table")))
}

enrichment_panel_server <- function(id, con, focus, factor_id) {
  moduleServer(id, function(input, output, session) {
    queries <- reactive(enrichment_queries_for_fit(con, focus()$fit_id, focus()$n_factors))
    observe({
      q <- queries()
      if (nrow(q) == 0) { updateSelectInput(session, "query", choices = character(0)); return() }
      key <- paste(q$query_type, q$direction, sep = "|")
      lab <- paste0(QUERY_LABELS[q$query_type] %||% q$query_type, " (", q$direction, ")")
      cur <- isolate(input$query)
      updateSelectInput(session, "query", choices = setNames(key, lab), selected = if (!is.null(cur) && cur %in% key) cur else key[1])
    })
    terms <- reactive({
      req(input$query, factor_id())
      qd <- strsplit(input$query, "|", fixed = TRUE)[[1]]
      factor_enrichment(con, factor_id(), qd[1], qd[2])
    })
    observe({
      s <- sort(unique(terms()$source))
      updateSelectInput(session, "source", choices = c("all", s), selected = isolate(input$source) %||% "all")
    })
    shown <- reactive({
      d <- terms(); if (nzchar(input$source %||% "") && input$source != "all") d <- d[d$source == input$source, ]
      head(d, input$top)
    })
    output$dots <- renderPlot({
      d <- shown(); validate(need(nrow(d) > 0, "No significant terms stored for this factor and analysis."))
      d$term <- factor(substr(d$term_name, 1, 60), levels = rev(unique(substr(d$term_name, 1, 60))))
      d$logp <- -log10(pmax(d$p_value, 1e-300))
      ggplot(d, aes(logp, term, colour = source)) + geom_point(size = 3) +
        labs(x = expression(-log[10]~adjusted~p), y = NULL, colour = NULL) + theme_minimal(base_size = 12)
    })
    output$table <- renderDT({
      d <- shown(); req(nrow(d) > 0)
      d$p_value <- signif(d$p_value, 3); d$nes <- round(d$nes, 2)
      datatable(d[c("source", "term_name", "p_value", "nes", "intersection_size", "term_size")], rownames = FALSE,
                colnames = c("Collection", "Term", "adj. p", "NES", "Overlap", "Term size"),
                selection = "none", options = list(dom = "tp", pageLength = 10))
    })
  })
}

## ---- Factor detail ---------------------------------------------------------------------

mod_factor_ui <- function(id) {
  ns <- NS(id)
  tagList(
    layout_columns(col_widths = c(4, 8), fill = FALSE, selectInput(ns("factor"), "Factor", choices = NULL), uiOutput(ns("summary"))),
    uiOutput(ns("tabs")))
}

mod_factor_server <- function(id, state, focus, con) {
  moduleServer(id, function(input, output, session) {
    factors <- reactive({ f <- focus(); fit_factors(con, f$fit_id, f$n_factors) })
    output$tabs <- renderUI({
      ns <- session$ns
      panels <- list(
        nav_panel("Top genes",
          layout_columns(col_widths = c(3, 9), sliderInput(ns("top_n"), "Genes per side", 5, 50, 20, step = 5),
                         plotOutput(ns("genes"), height = "520px"))),
        nav_panel("Matches in other fits",
          p(class = "text-muted", "This factor's best match in every other fit of the method (|", textOutput(ns("metric_name"), inline = TRUE),
            "|). Same-rank points are seed replicates; other ranks show where the factor persists."),
          plotOutput(ns("matches"), height = "420px")),
        nav_panel("Enrichment", enrichment_panel_ui(ns("enrich"))),
        nav_panel("Sample scores",
          selectInput(ns("field"), "Colour / group by", choices = NULL),
          plotOutput(ns("scores"), height = "420px")))
      if (focus()$method == "wgcna") panels <- c(panels, list(
        nav_panel("Module quality",
          p(class = "text-muted", "kME is each gene's correlation with the module eigengene. A coherent module has members well above non-members and above the permutation null."),
          plotOutput(ns("kme_hist"), height = "460px")),
        nav_panel("Module genes",
          selectInput(ns("gs_field"), "Gene significance for", choices = NULL),
          layout_columns(col_widths = c(6, 6), plotOutput(ns("gs"), height = "400px"), DTOutput(ns("genes_table"))))))
      do.call(navset_card_tab, panels)
    })
    observe({
      req(focus()$method == "wgcna"); fl <- gene_significance_fields(con, state$dataset())
      updateSelectInput(session, "gs_field", choices = fl$field)
    })
    observe({
      fac <- factors(); cur <- isolate(input$factor)
      ch <- setNames(fac$factor_index, fac$label)
      updateSelectInput(session, "factor", choices = ch,
                        selected = if (!is.null(cur) && cur %in% fac$factor_index) cur else fac$factor_index[1])
    })
    row <- reactive({ req(input$factor); factors()[factors()$factor_index == as.numeric(input$factor), ] })
    output$metric_name <- renderText(state$metric())
    output$summary <- renderUI({
      r <- row(); st <- r[[paste0("stability_", state$metric())]]
      p(class = "mt-4", tags$strong(focus_label(focus()), " / ", r$label), " — ", r$n_genes, " genes",
        if (!is.na(st)) paste0("; stability ", round(st, 3)), if (!is.na(r$excess_kurtosis)) paste0("; excess kurtosis ", round(r$excess_kurtosis, 2)))
    })
    loadings <- reactive({ f <- focus(); if (f$method == "wgcna") NULL else load_loadings(con, f$fit_id, f$n_factors) })
    labels <- reactive(gene_labels(con, state$dataset()))
    output$genes <- renderPlot({
      f <- focus(); j <- as.numeric(input$factor); req(j)
      if (f$method == "wgcna") {
        gm <- wgcna_gene_modules(con, f$fit_id); gm <- gm[gm$module == j, ]
        validate(need(nrow(gm) > 0, "No genes in this module."))
        d <- head(gm[order(-gm$kme_own), ], input$top_n)
        d <- data.frame(gene = label_genes(labels(), d$gene), loading = d$kme_own, side = "positive")
        xl <- "kME (correlation with module eigengene)"
      } else {
        L <- loadings(); validate(need(!is.null(L), "No loadings artifact for this fit."))
        d <- top_loading_genes(L, j, input$top_n); d$gene <- label_genes(labels(), d$gene); xl <- "loading"
      }
      d$gene <- factor(make.unique(d$gene), levels = rev(make.unique(d$gene)))
      ggplot(d, aes(loading, gene, fill = side)) + geom_col() +
        scale_fill_manual(values = c(positive = "#B2182B", negative = "#2166AC"), name = NULL) +
        facet_wrap(~side, scales = "free", ncol = 2) + labs(x = xl, y = NULL) + theme_minimal(base_size = 12) +
        theme(legend.position = "none")
    })
    output$matches <- renderPlot({
      f <- focus(); r <- row()
      d <- factor_matches_everywhere(con, r$factor_id, if (f$method == "wgcna") "jaccard" else state$metric())
      validate(need(f$method != "pca", "PCA is one deterministic fit, so there are no other fits to match against."))
      validate(need(nrow(d) > 0, "No stored matches for this factor."))
      d$value <- abs(d$value); d$kind <- ifelse(d$same_rank == 1, "same rank", "other rank")
      d$x <- if (f$method == "wgcna") d$other_power else d$other_rank
      ggplot(d, aes(factor(x), value, colour = kind)) + geom_jitter(width = 0.2, height = 0, alpha = 0.7) +
        scale_colour_manual(values = c("same rank" = "#4C78A8", "other rank" = "#F58518"), name = NULL) +
        coord_cartesian(ylim = c(0, 1)) + labs(x = if (f$method == "wgcna") "power of the other fit" else "rank of the other fit",
           y = paste("matched |", if (f$method == "wgcna") "jaccard" else state$metric(), "|")) +
        theme_minimal(base_size = 12) + theme(legend.position = "top")
    })
    enrichment_panel_server("enrich", con, focus, reactive(row()$factor_id))

    # WGCNA only: the selected factor is a module (factor_index = module number)
    module <- reactive({ req(focus()$method == "wgcna", input$factor); as.numeric(input$factor) })
    module_pal <- reactive(module_palette(wgcna_gene_modules(con, focus()$fit_id)$module))
    output$kme_hist <- renderPlot({
      h <- wgcna_kme_histograms(con, focus()$fit_id, module()); validate(need(nrow(h) > 0, "No kME histogram stored for this module."))
      h$mid <- (h$bin_lo + h$bin_hi) / 2
      ggplot(h, aes(mid, count, fill = grp)) + geom_col(position = "identity", alpha = 0.6, width = diff(range(h$mid)) / length(unique(h$mid))) +
        facet_wrap(~grp, scales = "free_y", ncol = 1) +
        scale_fill_manual(values = c(member = "#B2182B", nonmember = "#4C78A8", permutation = "grey50"), guide = "none") +
        labs(x = "kME", y = "genes") + theme_minimal(base_size = 11)
    })
    module_genes <- reactive({
      g <- wgcna_gene_modules(con, focus()$fit_id); g <- g[g$module == module(), ]
      g$label <- label_genes(labels(), g$gene); g[order(-g$kme_own), ]
    })
    output$genes_table <- renderDT({
      g <- module_genes(); req(nrow(g) > 0)
      d <- data.frame(Gene = g$label, kME = round(g$kme_own, 3), `Next module` = g$next_module, `Next kME` = round(g$kme_next, 3),
                      `p (perm)` = signif(g$p_perm, 2), check.names = FALSE)
      datatable(d, rownames = FALSE, selection = "none", options = list(dom = "ftp", pageLength = 10))
    })
    output$gs <- renderPlot({
      req(input$gs_field); g <- module_genes(); gs <- gene_significance_for(con, state$dataset(), input$gs_field)
      d <- merge(g[c("gene", "kme_own")], gs, by = "gene"); validate(need(nrow(d) > 0, "No gene-significance values for this module."))
      ggplot(d, aes(kme_own, abs(statistic))) + geom_point(alpha = 0.5, colour = module_pal()[[as.character(module())]]) +
        geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "grey30") +
        labs(x = "kME", y = paste("|gene significance|:", input$gs_field),
             subtitle = "Genes central to the module that also track the field drive its association.") + theme_minimal(base_size = 12)
    })

    meta <- reactive(dataset_sample_metadata(con, state$dataset()))
    observe({
      m <- meta(); f <- if (is.null(m)) character(0) else setdiff(names(m), "sample_id")
      updateSelectInput(session, "field", choices = f)
    })
    output$scores <- renderPlot({
      S <- focus_scores(con, focus()); m <- meta()
      validate(need(!is.null(S) && !is.null(m), "No scores or sample metadata for this fit."))
      lab <- row()$label; validate(need(lab %in% colnames(S), "No scores for this factor."))
      req(input$field)
      d <- data.frame(sample_id = rownames(S), score = S[, lab])
      d$value <- m[[input$field]][match(d$sample_id, m$sample_id)]
      d <- d[!is.na(d$value), ]; validate(need(nrow(d) > 2, "No samples with this field."))
      p <- if (is.numeric(d$value)) ggplot(d, aes(value, score)) + geom_point(alpha = 0.7) + geom_smooth(method = "loess", formula = y ~ x, se = FALSE, colour = "grey30")
           else ggplot(d, aes(factor(value), score)) + geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.15, height = 0, alpha = 0.6)
      p + labs(x = input$field, y = paste(lab, "score")) + theme_minimal(base_size = 12) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    })
  })
}

## ---- Enrichment page ---------------------------------------------------------------------

mod_enrichment_ui <- function(id) {
  ns <- NS(id)
  navset_card_tab(
    nav_panel("All factors",
      selectInput(ns("query"), "Analysis", choices = NULL),
      p(class = "text-muted", "Number of significant terms per factor, and the strongest terms across the fit (top 4 per factor)."),
      layout_columns(col_widths = c(4, 8), plotOutput(ns("counts"), height = "460px"), plotOutput(ns("heat"), height = "460px"))),
    nav_panel("One factor",
      selectInput(ns("factor"), "Factor", choices = NULL),
      enrichment_panel_ui(ns("one"))))
}

mod_enrichment_server <- function(id, state, focus, con) {
  moduleServer(id, function(input, output, session) {
    factors <- reactive({ f <- focus(); fit_factors(con, f$fit_id, f$n_factors) })
    queries <- reactive(enrichment_queries_for_fit(con, focus()$fit_id, focus()$n_factors))
    observe({
      q <- queries(); key <- paste(q$query_type, q$direction, sep = "|")
      lab <- if (nrow(q) == 0) character(0) else paste0(QUERY_LABELS[q$query_type] %||% q$query_type, " (", q$direction, ")")
      cur <- isolate(input$query)
      updateSelectInput(session, "query", choices = setNames(key, lab), selected = if (!is.null(cur) && cur %in% key) cur else key[1])
      fac <- factors(); curf <- isolate(input$factor)
      updateSelectInput(session, "factor", choices = setNames(fac$factor_id, fac$label),
                        selected = if (!is.null(curf) && curf %in% fac$factor_id) curf else fac$factor_id[1])
    })
    terms <- reactive({
      req(input$query); qd <- strsplit(input$query, "|", fixed = TRUE)[[1]]
      f <- focus(); fit_enrichment(con, f$fit_id, qd[1], qd[2], f$n_factors)
    }) |> bindCache(focus()$fit_id, focus()$n_factors, input$query)
    output$counts <- renderPlot({
      d <- terms(); fac <- factors(); validate(need(nrow(d) > 0, "No stored enrichment for this fit and analysis."))
      n <- as.data.frame(table(label = d$label), stringsAsFactors = FALSE)
      n <- merge(data.frame(label = fac$label, stringsAsFactors = FALSE), n, all.x = TRUE); n$Freq[is.na(n$Freq)] <- 0
      n$label <- factor(n$label, levels = rev(fac$label))
      ggplot(n, aes(Freq, label)) + geom_col(fill = "#4C78A8") + labs(x = "significant terms", y = NULL) + theme_minimal(base_size = 12)
    })
    output$heat <- renderPlot({
      d <- terms(); validate(need(nrow(d) > 0, "No stored enrichment for this fit and analysis."))
      d <- d[order(d$p_value), ]
      top <- do.call(rbind, lapply(split(d, d$label), head, 4))
      keep <- d[d$term_name %in% top$term_name, ]
      keep$logp <- pmin(-log10(pmax(keep$p_value, 1e-300)), 30)
      keep$term <- substr(keep$term_name, 1, 45)
      tm <- tapply(keep$logp, list(keep$term, keep$label), max)
      tm[is.na(tm)] <- 0
      ord <- rownames(tm)[order(apply(tm, 1, which.max), -apply(tm, 1, max))]
      long <- as.data.frame(as.table(tm[ord, , drop = FALSE]), stringsAsFactors = FALSE); names(long) <- c("term", "factor", "logp")
      long$term <- factor(long$term, levels = rev(ord)); long$factor <- factor(long$factor, levels = factors()$label)
      ggplot(long, aes(factor, term, fill = logp)) + geom_tile(colour = "white") +
        scale_fill_gradient(low = "grey95", high = "#B2182B", name = "-log10 p") + labs(x = NULL, y = NULL) +
        theme_minimal(base_size = 11) + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
    })
    enrichment_panel_server("one", con, focus, reactive(as.numeric(req(input$factor))))
  })
}

## ---- Metadata associations -----------------------------------------------------------------

mod_metadata_ui <- function(id) {
  ns <- NS(id)
  navset_card_tab(
    nav_panel("Factor x field",
      p(class = "text-muted", "How strongly each factor's sample scores track each sample-metadata field: Spearman correlation for numeric fields,",
        "Kruskal-Wallis for the rest. Colour is -log10 of the BH-adjusted p-value over the whole grid; * marks adjusted p < 0.05;",
        "numeric fields show the correlation."),
      plotOutput(ns("grid"), height = "520px", click = ns("pick"))),
    nav_panel("One field",
      layout_columns(col_widths = c(4, 4, 4), selectInput(ns("field"), "Field", choices = NULL),
                     selectInput(ns("factor"), "Factor", choices = NULL),
                     uiOutput(ns("hint"))),
      plotOutput(ns("one"), height = "420px")))
}

mod_metadata_server <- function(id, state, focus, con) {
  moduleServer(id, function(input, output, session) {
    meta <- reactive(dataset_sample_metadata(con, state$dataset()))
    scores <- reactive(focus_scores(con, focus()))
    scan <- reactive({
      validate(need(!is.null(scores()) && !is.null(meta()), "This fit has no scores, or the dataset has no sample metadata."))
      association_scan(scores(), meta())
    }) |> bindCache(state$dataset(), focus()$fit_id, focus()$n_factors)
    observe({
      m <- meta(); S <- scores()
      f <- if (is.null(m)) character(0) else setdiff(names(m), "sample_id")
      updateSelectInput(session, "field", choices = f, selected = isolate(input$field) %||% f[1])
      updateSelectInput(session, "factor", choices = if (is.null(S)) character(0) else colnames(S))
    })
    output$grid <- renderPlot({
      d <- scan(); validate(need(nrow(d) > 0, "No usable fields."))
      d$logp <- pmin(-log10(pmax(d$padj, 1e-300)), 10)
      d$text <- ifelse(d$test == "spearman", sprintf("%.2f", d$statistic), "")
      d$text <- paste0(d$text, ifelse(d$padj < 0.05, "*", ""))
      lv <- colnames(scores()); d$component <- factor(d$component, levels = rev(lv))
      ggplot(d, aes(field, component, fill = logp)) + geom_tile(colour = "white") + geom_text(aes(label = text), size = 3) +
        scale_fill_gradient(low = "grey96", high = "#B2182B", limits = c(0, 10), name = "-log10\nadj. p") +
        labs(x = NULL, y = NULL) + theme_minimal(base_size = 12) +
        theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
    })
    observeEvent(input$pick, {
      d <- scan(); fields <- unique(d$field); comps <- colnames(scores())
      j <- round(input$pick$x); i <- length(comps) + 1 - round(input$pick$y)
      if (j >= 1 && j <= length(fields) && i >= 1 && i <= length(comps)) {
        updateSelectInput(session, "field", selected = fields[j]); updateSelectInput(session, "factor", selected = comps[i])
        showNotification(paste("Selected", comps[i], "x", fields[j], "-- see the One field tab"), duration = 3)
      }
    })
    output$hint <- renderUI(p(class = "text-muted mt-4", "Click a cell in the grid to jump here."))
    output$one <- renderPlot({
      req(input$field, input$factor); S <- scores(); m <- meta()
      validate(need(input$factor %in% colnames(S), "Pick a factor."))
      d <- data.frame(sample_id = rownames(S), score = S[, input$factor])
      d$value <- m[[input$field]][match(d$sample_id, m$sample_id)]
      d <- d[!is.na(d$value), ]; validate(need(nrow(d) > 2, "No samples with this field."))
      p <- if (is.numeric(d$value)) ggplot(d, aes(value, score)) + geom_point(alpha = 0.7) + geom_smooth(method = "loess", formula = y ~ x, se = FALSE, colour = "grey30")
           else ggplot(d, aes(factor(value), score)) + geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.15, height = 0, alpha = 0.6)
      p + labs(x = input$field, y = paste(input$factor, "score")) + theme_minimal(base_size = 12) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    })
  })
}

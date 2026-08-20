# Latent-factor stability explorer.
#
# Reads whatever the stability DB currently contains -- fully decoupled
# from ingest (R/ingest_results.R). DB path is parameterized via the
# STABILITY_DB environment variable, falling back to results/stability.sqlite
# relative to the project root.
#
# Drill-down levels (breadcrumb at top navigates back up):
#   0 dataset overview -> 1 method -> 2 rank/parameter -> 3 factor

library(shiny)
library(bslib)
library(DBI)
library(RSQLite)
library(ggplot2)
library(DT)

source(file.path("R", "db_helpers.R"), local = TRUE)

db_path <- Sys.getenv("STABILITY_DB", unset = "")
if (!nzchar(db_path)) {
  # fall back to the project-level default, whether launched from app/ or repo root
  cand <- c(file.path("..", "results", "stability.sqlite"),
            file.path("results", "stability.sqlite"))
  db_path <- cand[file.exists(cand)][1]
  if (is.na(db_path)) stop("No stability DB found -- set STABILITY_DB or run R/ingest_results.R first")
}
con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
options(stability.db_dir = normalizePath(dirname(db_path)))
onStop(function() DBI::dbDisconnect(con))

FACTORIZATION_METHODS <- c("pca", "nmf", "cogaps")

## ---- ui -----------------------------------------------------------------------

ui <- page_fillable(
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  padding = "1rem",
  layout_columns(
    col_widths = c(3, 3, 6),
    selectInput("dataset", NULL, choices = NULL, width = "100%"),
    radioButtons("metric", NULL, choices = METRICS, selected = "cosine", inline = TRUE),
    uiOutput("breadcrumb")
  ),
  uiOutput("level_ui")
)

## ---- server --------------------------------------------------------------------

server <- function(input, output, session) {

  nav <- reactiveValues(level = 0, method = NULL, rank = NULL,
                        fit = NULL, factor_index = NULL,
                        wgcna_fit = NULL, wto_a = NULL, wto_b = NULL)

  datasets <- list_datasets(con)
  updateSelectInput(session, "dataset", choices = datasets$dataset_id)

  metric <- reactive(check_metric(input$metric))
  ds <- reactive({ req(input$dataset); input$dataset })

  observeEvent(input$dataset, {
    nav$level <- 0; nav$method <- NULL; nav$rank <- NULL
    nav$fit <- NULL; nav$factor_index <- NULL
  })

  ## ---- breadcrumb ---------------------------------------------------------------

  output$breadcrumb <- renderUI({
    crumbs <- list(actionLink("bc_home", ds()))
    if (nav$level >= 1) crumbs <- c(crumbs, list(HTML("&nbsp;&raquo;&nbsp;"), actionLink("bc_method", toupper(nav$method))))
    if (nav$level >= 2 && nav$method %in% FACTORIZATION_METHODS) {
      crumbs <- c(crumbs, list(HTML("&nbsp;&raquo;&nbsp;"), actionLink("bc_rank", paste0("rank ", nav$rank))))
    }
    if (nav$level >= 2 && nav$method == "wgcna") {
      crumbs <- c(crumbs, list(HTML("&nbsp;&raquo;&nbsp;"), actionLink("bc_rank", paste0("power ", get_fit(con, nav$wgcna_fit)$power))))
    }
    if (nav$level >= 2 && nav$method == "wto") {
      crumbs <- c(crumbs, list(HTML("&nbsp;&raquo;&nbsp;"), actionLink("bc_rank", "run pair")))
    }
    if (nav$level >= 3) crumbs <- c(crumbs, list(HTML("&nbsp;&raquo;&nbsp;"), tags$b(paste0("factor ", nav$factor_index))))
    div(style = "font-size: 1.1rem; padding-top: 6px;", crumbs)
  })
  observeEvent(input$bc_home,   { nav$level <- 0 })
  observeEvent(input$bc_method, { nav$level <- 1 })
  observeEvent(input$bc_rank,   { nav$level <- 2 })

  ## ---- level dispatcher -----------------------------------------------------------

  output$level_ui <- renderUI({
    switch(as.character(nav$level),
      "0" = level0_ui(),
      "1" = level1_ui(),
      "2" = level2_ui(),
      "3" = level3_ui())
  })

  ## =====================================================================
  ## LEVEL 0 -- dataset overview
  ## =====================================================================

  level0_ui <- function() {
    ov <- method_overview(con, ds())
    methods <- sort(unique(ov$counts$method))
    cards <- lapply(methods, function(m) {
      cnt <- ov$counts[ov$counts$method == m, ]
      headline <- if (m == "wgcna") {
        if (nrow(ov$wgcna) > 0 && !is.na(ov$wgcna$ari)) sprintf("mean ARI %.2f", ov$wgcna$ari) else "--"
      } else if (m == "wto") {
        if (nrow(ov$wto) > 0 && !is.na(ov$wto$pearson)) sprintf("mean edge r %.2f", ov$wto$pearson) else "--"
      } else {
        s <- ov$stability[ov$stability$method == m, metric()]
        if (length(s) == 1 && !is.na(s)) sprintf("mean matched %s %.2f", metric(), s) else "--"
      }
      card(
        card_header(tags$b(toupper(m))),
        card_body(
          p(tags$span(style = "font-size:1.6rem;", headline)),
          p(sprintf("%d fits ok / %d failed / %d missing",
                    cnt$n_ok, cnt$n_failed, cnt$n_missing)),
          actionButton(paste0("open_", m), "Explore", class = "btn-primary btn-sm")
        )
      )
    })
    tagList(
      h4("Method overview -- cross-seed stability at a glance"),
      layout_column_wrap(width = 1 / max(length(cards), 1), !!!cards),
      h5("Ingested job families"),
      tableOutput("l0_ingests")
    )
  }
  output$l0_ingests <- renderTable(list_ingests(con, ds()))

  # one observer per possible method, registered exactly once (registering
  # inside observe() would stack duplicate handlers on every re-run)
  for (m in c(FACTORIZATION_METHODS, "wgcna", "wto")) {
    local({
      m_local <- m
      observeEvent(input[[paste0("open_", m_local)]], {
        nav$method <- m_local
        nav$level <- 1
      }, ignoreInit = TRUE)
    })
  }

  ## =====================================================================
  ## LEVEL 1 -- method view
  ## =====================================================================

  level1_ui <- function() {
    m <- nav$method
    if (m %in% FACTORIZATION_METHODS) {
      navset_card_tab(
        nav_panel("Seed stability vs rank",
          p("Distribution of matched factor similarity across all seed pairs, per rank. Click a rank to drill in."),
          plotOutput("l1_stability", click = "l1_stability_click", height = "420px")),
        nav_panel("Rank selection (masking-CV)",
          plotOutput("l1_maskcv", height = "420px")),
        nav_panel("Cross-rank persistence",
          p("Mean matched similarity between factors found at different ranks (all seed pairs)."),
          plotOutput("l1_crossrank", height = "380px"),
          layout_columns(col_widths = c(3, 9),
            selectInput("l1_ref_rank", "Track factors from rank:", choices = NULL),
            plotOutput("l1_trajectory", height = "320px")))
      )
    } else if (m == "wgcna") {
      navset_card_tab(
        nav_panel("Module stability (ARI)",
          p("Adjusted Rand Index between module assignments at each pair of soft-threshold powers. Click a diagonal-adjacent cell or use the selector to drill into a power."),
          plotOutput("l1_wgcna_ari", height = "420px"),
          selectInput("l1_wgcna_power", "Drill into power:", choices = NULL),
          actionButton("l1_wgcna_go", "Explore power", class = "btn-primary btn-sm")),
        nav_panel("Modules vs power",
          plotOutput("l1_wgcna_counts", height = "420px"))
      )
    } else if (m == "wto") {
      navset_card_tab(
        nav_panel("Run-pair stability",
          p("Edge-value correlation and significant-edge overlap between wTO runs (seed x bootstrap-n grid)."),
          plotOutput("l1_wto_heat", height = "420px"),
          layout_columns(col_widths = c(4, 4, 4),
            selectInput("l1_wto_a", "Run A:", choices = NULL),
            selectInput("l1_wto_b", "Run B:", choices = NULL),
            actionButton("l1_wto_go", "Compare edges", class = "btn-primary btn-sm",
                         style = "margin-top: 32px;")))
      )
    }
  }

  l1_stab_data <- reactive({
    req(nav$method %in% FACTORIZATION_METHODS)
    seed_stability_by_rank(con, ds(), nav$method)
  })

  output$l1_stability <- renderPlot({
    d <- l1_stab_data()
    validate(need(nrow(d) > 0,
      paste0("No same-rank seed pairs for ", toupper(nav$method),
             " -- PCA is deterministic (one fit per rank), so cross-seed stability ",
             "doesn't apply; see the masking-CV and cross-rank tabs instead.")))
    d$sim <- d[[metric()]]
    ggplot(d, aes(x = factor(rank), y = sim)) +
      geom_boxplot(outlier.size = 0.6, fill = "grey85") +
      labs(x = "rank (number of factors)", y = paste("matched", metric(), "similarity"),
           title = paste(toupper(nav$method), "-- cross-seed factor stability by rank")) +
      theme_minimal(base_size = 14)
  })
  observeEvent(input$l1_stability_click, {
    d <- l1_stab_data(); req(nrow(d) > 0)
    ranks <- sort(unique(d$rank))
    i <- round(input$l1_stability_click$x)
    if (i >= 1 && i <= length(ranks)) {
      nav$rank <- ranks[i]
      nav$level <- 2
    }
  })

  output$l1_maskcv <- renderPlot({
    d <- maskcv_curve(con, ds(), nav$method); req(nrow(d) > 0)
    if (all(is.na(d$alpha))) {
      ggplot(d, aes(rank, mse)) + geom_line() + geom_point(size = 2) +
        labs(title = paste(toupper(nav$method), "masking-CV: held-out reconstruction MSE"),
             x = "rank", y = "held-out MSE") +
        theme_minimal(base_size = 14)
    } else {
      ggplot(d, aes(factor(rank), factor(alpha), fill = mse)) +
        geom_tile() +
        geom_text(aes(label = ifelse(is.na(mse), "failed", sprintf("%.2f", mse))), size = 3) +
        scale_fill_viridis_c(na.value = "grey70", direction = -1) +
        labs(title = paste(toupper(nav$method), "masking-CV: held-out MSE by rank x alpha"),
             x = "rank", y = "alpha") +
        theme_minimal(base_size = 14)
    }
  })

  l1_crossrank_data <- reactive({
    req(nav$method %in% FACTORIZATION_METHODS)
    crossrank_matrix(con, ds(), nav$method)
  })
  output$l1_crossrank <- renderPlot({
    d <- l1_crossrank_data(); req(nrow(d) > 0)
    d$sim <- d[[metric()]]
    # symmetrize for display
    d2 <- rbind(d[, c("rank_a", "rank_b", "sim")],
                setNames(d[, c("rank_b", "rank_a", "sim")], c("rank_a", "rank_b", "sim")))
    ggplot(d2, aes(factor(rank_a), factor(rank_b), fill = sim)) +
      geom_tile() +
      scale_fill_viridis_c(limits = c(0, 1)) +
      labs(x = "rank A", y = "rank B", fill = metric(),
           title = "Mean matched similarity between ranks (factor persistence)") +
      theme_minimal(base_size = 14)
  })

  observe({
    req(nav$level == 1, nav$method %in% FACTORIZATION_METHODS)
    d <- l1_crossrank_data()
    ranks <- sort(unique(c(d$rank_a, d$rank_b)))
    updateSelectInput(session, "l1_ref_rank", choices = ranks)
  })

  output$l1_trajectory <- renderPlot({
    req(input$l1_ref_rank)
    ref_rank <- as.integer(input$l1_ref_rank)
    # per-factor best-match similarity from reference rank to every other rank,
    # averaged over seed pairs: pull matched cross-rank pairs touching ref rank
    d <- DBI::dbGetQuery(con,
      "SELECT fa.rank AS rank_a, fb.rank AS rank_b, fp.factor_a, fp.factor_b,
              fp.cosine, fp.pearson, fp.spearman
       FROM factor_pairs fp
       JOIN fits fa ON fa.fit_id = fp.fit_a
       JOIN fits fb ON fb.fit_id = fp.fit_b
       WHERE fp.matched = 1 AND fp.same_rank = 0
         AND fa.dataset_id = ?1 AND fa.method = ?2
         AND (fa.rank = ?3 OR fb.rank = ?3)",
      params = list(ds(), nav$method, ref_rank))
    req(nrow(d) > 0)
    # orient so the reference rank's factor is 'factor_ref'
    ref_side_a <- d$rank_a == ref_rank
    d$factor_ref <- ifelse(ref_side_a, d$factor_a, d$factor_b)
    d$other_rank <- ifelse(ref_side_a, d$rank_b, d$rank_a)
    d$sim <- d[[metric()]]
    agg <- aggregate(sim ~ factor_ref + other_rank, data = d, FUN = mean)
    ggplot(agg, aes(other_rank, sim, group = factor_ref, color = factor(factor_ref))) +
      geom_line() + geom_point(size = 1.4) +
      labs(x = "other rank", y = paste("mean best-match", metric()),
           color = paste0("factor @ rank ", ref_rank),
           title = paste0("Persistence of rank-", ref_rank, " factors across other ranks")) +
      coord_cartesian(ylim = c(min(0, min(agg$sim)), 1)) +
      theme_minimal(base_size = 14)
  })

  # WGCNA level 1
  output$l1_wgcna_ari <- renderPlot({
    d <- wgcna_ari_matrix(con, ds()); req(nrow(d) > 0)
    d2 <- rbind(d, setNames(d[, c("power_b", "power_a", "ari")], c("power_a", "power_b", "ari")))
    ggplot(d2, aes(factor(power_a), factor(power_b), fill = ari)) +
      geom_tile() + geom_text(aes(label = sprintf("%.2f", ari)), size = 3) +
      scale_fill_viridis_c(limits = c(0, 1)) +
      labs(x = "power", y = "power", title = "Module assignment ARI between soft-threshold powers") +
      theme_minimal(base_size = 14)
  })
  output$l1_wgcna_counts <- renderPlot({
    f <- wgcna_fits(con, ds()); req(nrow(f) > 0)
    ggplot(f, aes(factor(power), n_factors)) +
      geom_col(fill = "grey70") +
      labs(x = "soft-threshold power", y = "number of modules (excl. unassigned)",
           title = "Module count vs power") +
      theme_minimal(base_size = 14)
  })
  observe({
    req(nav$level == 1, nav$method == "wgcna")
    f <- wgcna_fits(con, ds())
    updateSelectInput(session, "l1_wgcna_power", choices = setNames(f$fit_id, paste0("power ", f$power)))
  })
  observeEvent(input$l1_wgcna_go, {
    nav$wgcna_fit <- as.integer(input$l1_wgcna_power)
    nav$level <- 2
  })

  # wTO level 1
  output$l1_wto_heat <- renderPlot({
    d <- wto_pair_matrix(con, ds()); req(nrow(d) > 0)
    lab <- function(n, s) paste0("n=", n, ",seed=", s)
    d$a <- lab(d$n_a, d$seed_a); d$b <- lab(d$n_b, d$seed_b)
    d2 <- rbind(d[, c("a", "b", "pearson", "jaccard_sig")],
                setNames(d[, c("b", "a", "pearson", "jaccard_sig")], c("a", "b", "pearson", "jaccard_sig")))
    long <- rbind(data.frame(a = d2$a, b = d2$b, value = d2$pearson, what = "edge wTO Pearson r"),
                  data.frame(a = d2$a, b = d2$b, value = d2$jaccard_sig, what = "significant-edge Jaccard"))
    ggplot(long, aes(a, b, fill = value)) +
      geom_tile() + facet_wrap(~what) +
      scale_fill_viridis_c(limits = c(0, 1)) +
      labs(x = NULL, y = NULL, title = "wTO run-pair stability") +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  observe({
    req(nav$level == 1, nav$method == "wto")
    f <- wto_fits(con, ds())
    ch <- setNames(f$fit_id, paste0("n=", f$n_boot, ", seed=", f$seed))
    updateSelectInput(session, "l1_wto_a", choices = ch)
    updateSelectInput(session, "l1_wto_b", choices = ch)
  })
  observeEvent(input$l1_wto_go, {
    nav$wto_a <- as.integer(input$l1_wto_a)
    nav$wto_b <- as.integer(input$l1_wto_b)
    nav$level <- 2
  })

  ## =====================================================================
  ## LEVEL 2 -- rank / parameter view
  ## =====================================================================

  level2_ui <- function() {
    m <- nav$method
    if (m %in% FACTORIZATION_METHODS) {
      f <- fits_at_rank(con, ds(), m, nav$rank)
      seed_choices <- setNames(f$fit_id, paste0("seed ", f$seed))
      navset_card_tab(
        nav_panel("Seed-pair matrix",
          p("Mean matched factor similarity for every seed pair at this rank."),
          plotOutput("l2_seedpair", height = "420px")),
        nav_panel("Per-factor stability",
          layout_columns(col_widths = c(3, 9),
            selectInput("l2_ref_fit", "Reference seed:", choices = seed_choices),
            p("Each factor's matched similarity to every other seed. Click a factor to drill in.")),
          plotOutput("l2_factor_stability", click = "l2_factor_click", height = "400px")),
        nav_panel("Factor x factor heatmap",
          layout_columns(col_widths = c(3, 3, 6),
            selectInput("l2_fit_a", "Fit A (seed):", choices = seed_choices),
            selectInput("l2_fit_b", "Fit B (seed):", choices = seed_choices,
                        selected = if (length(seed_choices) > 1) seed_choices[[2]] else seed_choices[[1]]),
            p("All factor-pair similarities; Hungarian matches outlined.")),
          plotOutput("l2_ff_heatmap", height = "440px"))
      )
    } else if (m == "wgcna") {
      f <- wgcna_fits(con, ds())
      navset_card_tab(
        nav_panel("Modules at this power",
          tableOutput("l2_wgcna_sizes"),
          h6("Best Hungarian-matched module at every other power (Jaccard):"),
          DTOutput("l2_wgcna_matches")),
        nav_panel("Module overlap vs another power",
          selectInput("l2_wgcna_other", "Compare with power:",
                      choices = setNames(f$fit_id, paste0("power ", f$power))),
          plotOutput("l2_wgcna_jaccard", height = "440px"))
      )
    } else if (m == "wto") {
      navset_card_tab(
        nav_panel("Edge scatter",
          p("Per-edge wTO values between the two selected runs (sampled if very large)."),
          plotOutput("l2_wto_scatter", height = "460px"))
      )
    }
  }

  output$l2_seedpair <- renderPlot({
    d <- seedpair_matrix(con, ds(), nav$method, nav$rank); req(nrow(d) > 0)
    d$sim <- d[[metric()]]
    d2 <- rbind(d[, c("seed_a", "seed_b", "sim")],
                setNames(d[, c("seed_b", "seed_a", "sim")], c("seed_a", "seed_b", "sim")))
    ggplot(d2, aes(factor(seed_a), factor(seed_b), fill = sim)) +
      geom_tile() + geom_text(aes(label = sprintf("%.2f", sim)), size = 3) +
      scale_fill_viridis_c(limits = c(0, 1)) +
      labs(x = "seed", y = "seed", fill = metric(),
           title = paste0(toupper(nav$method), " rank ", nav$rank, ": mean matched similarity per seed pair")) +
      theme_minimal(base_size = 14)
  })

  l2_factor_data <- reactive({
    req(input$l2_ref_fit)
    factor_stability_for_fit(con, as.integer(input$l2_ref_fit))
  })
  output$l2_factor_stability <- renderPlot({
    d <- l2_factor_data(); req(nrow(d) > 0)
    d$sim <- d[[metric()]]
    ggplot(d, aes(factor(factor_index), sim)) +
      geom_boxplot(fill = "grey85", outlier.shape = NA) +
      geom_jitter(width = 0.15, size = 1, alpha = 0.6) +
      ylim(0, 1) +
      labs(x = "factor", y = paste("matched", metric(), "to other seeds"),
           title = "Per-factor cross-seed stability (reference seed)") +
      theme_minimal(base_size = 14)
  })
  observeEvent(input$l2_factor_click, {
    d <- l2_factor_data(); req(nrow(d) > 0)
    idxs <- sort(unique(d$factor_index))
    i <- round(input$l2_factor_click$x)
    if (i >= 1 && i <= length(idxs)) {
      nav$fit <- as.integer(input$l2_ref_fit)
      nav$factor_index <- idxs[i]
      nav$level <- 3
    }
  })

  output$l2_ff_heatmap <- renderPlot({
    req(input$l2_fit_a, input$l2_fit_b)
    a <- as.integer(input$l2_fit_a); b <- as.integer(input$l2_fit_b)
    req(a != b)
    d <- factor_pair_heatmap_data(con, a, b); req(nrow(d) > 0)
    d$sim <- d[[metric()]]
    ggplot(d, aes(factor(factor_a), factor(factor_b), fill = sim)) +
      geom_tile() +
      geom_tile(data = d[d$matched == 1, ], color = "red", linewidth = 1, fill = NA) +
      geom_text(aes(label = sprintf("%.2f", sim)), size = 3) +
      scale_fill_viridis_c(limits = c(-1, 1)) +
      labs(x = "factor (fit A)", y = "factor (fit B)", fill = metric(),
           title = "Factor x factor similarity (red outline = Hungarian match)") +
      theme_minimal(base_size = 14)
  })

  # WGCNA level 2
  output$l2_wgcna_sizes <- renderTable({
    req(nav$wgcna_fit)
    wgcna_module_sizes(con, nav$wgcna_fit)
  })
  output$l2_wgcna_matches <- renderDT({
    req(nav$wgcna_fit)
    d <- wgcna_module_best_matches(con, ds(), nav$wgcna_fit)
    datatable(d[order(d$module, d$other_power), ], rownames = FALSE,
              options = list(pageLength = 15))
  })
  output$l2_wgcna_jaccard <- renderPlot({
    req(nav$wgcna_fit, input$l2_wgcna_other)
    other <- as.integer(input$l2_wgcna_other)
    req(other != nav$wgcna_fit)
    d <- wgcna_module_jaccard(con, nav$wgcna_fit, other); req(nrow(d) > 0)
    ggplot(d, aes(factor(module_a), factor(module_b), fill = jaccard)) +
      geom_tile() +
      geom_tile(data = d[d$matched == 1, ], color = "red", linewidth = 1, fill = NA) +
      geom_text(aes(label = sprintf("%.2f", jaccard)), size = 3) +
      scale_fill_viridis_c(limits = c(0, 1)) +
      labs(x = "module (this power)", y = "module (other power)",
           title = "Module membership Jaccard (red outline = Hungarian match)") +
      theme_minimal(base_size = 14)
  })

  # wTO level 2
  output$l2_wto_scatter <- renderPlot({
    req(nav$wto_a, nav$wto_b, nav$wto_a != nav$wto_b)
    fa <- get_fit(con, nav$wto_a); fb <- get_fit(con, nav$wto_b)
    ea <- arrow::read_parquet(resolve_artifact(fa$loadings_file))
    eb <- arrow::read_parquet(resolve_artifact(fb$loadings_file))
    ea$key <- paste(pmin(ea$node1, ea$node2), pmax(ea$node1, ea$node2), sep = "|")
    eb$key <- paste(pmin(eb$node1, eb$node2), pmax(eb$node1, eb$node2), sep = "|")
    m <- match(ea$key, eb$key); ok <- !is.na(m)
    d <- data.frame(wto_a = ea$wto[ok], wto_b = eb$wto[m[ok]],
                    sig_both = !is.na(ea$padj[ok]) & ea$padj[ok] < 0.05 &
                               !is.na(eb$padj[m[ok]]) & eb$padj[m[ok]] < 0.05)
    if (nrow(d) > 100000) d <- d[sample(nrow(d), 100000), ]
    ggplot(d, aes(wto_a, wto_b, color = sig_both)) +
      geom_point(size = 0.4, alpha = 0.4) +
      geom_abline(linetype = "dashed") +
      scale_color_manual(values = c(`FALSE` = "grey60", `TRUE` = "firebrick")) +
      labs(x = paste0("wTO (n=", fa$n_boot, ", seed=", fa$seed, ")"),
           y = paste0("wTO (n=", fb$n_boot, ", seed=", fb$seed, ")"),
           color = "sig. both runs",
           title = "Per-edge wTO stability between runs") +
      theme_minimal(base_size = 14)
  })

  ## =====================================================================
  ## LEVEL 3 -- factor view
  ## =====================================================================

  level3_ui <- function() {
    navset_card_tab(
      nav_panel("Loadings",
        sliderInput("l3_topn", "Top genes:", min = 10, max = 100, value = 25, step = 5),
        plotOutput("l3_loadings", height = "420px")),
      nav_panel("Matches everywhere",
        p("This factor's Hungarian-matched partner in every other fit (all seeds and ranks)."),
        DTOutput("l3_matches"),
        plotOutput("l3_match_scatter", height = "380px")),
      nav_panel("Enrichment",
        layout_columns(col_widths = c(3, 3, 3, 3),
          radioButtons("l3_enrich_type", "Query:", c("ORA" = "ora", "GSEA" = "gsea"), inline = TRUE),
          radioButtons("l3_enrich_dir", "Loadings:", c("positive" = "pos", "negative" = "neg"), inline = TRUE),
          sliderInput("l3_enrich_topn", "Top genes (ORA):", min = 50, max = 300, value = 100, step = 50),
          actionButton("l3_enrich_go", "Run / load enrichment", class = "btn-primary",
                       style = "margin-top: 24px;")),
        DTOutput("l3_enrichment"))
    )
  }

  l3_loadings <- reactive({
    req(nav$fit, nav$factor_index)
    L <- load_loadings(con, nav$fit); req(!is.null(L))
    L[, nav$factor_index]
  })

  output$l3_loadings <- renderPlot({
    v <- l3_loadings()
    topn <- input$l3_topn %||% 25
    ord <- order(abs(v), decreasing = TRUE)[seq_len(min(topn, length(v)))]
    d <- data.frame(gene = names(v)[ord], loading = v[ord])
    d$gene <- factor(d$gene, levels = rev(d$gene))
    ggplot(d, aes(loading, gene)) +
      geom_col(fill = "grey40") +
      labs(title = sprintf("%s rank %s seed %s -- factor %d: top %d loadings",
                           toupper(nav$method), nav$rank,
                           get_fit(con, nav$fit)$seed, nav$factor_index, nrow(d)),
           y = NULL) +
      theme_minimal(base_size = 13)
  })

  l3_matches <- reactive({
    req(nav$fit, nav$factor_index)
    d <- factor_matches_everywhere(con, nav$fit, nav$factor_index)
    d[order(-d$cosine), ]
  })
  output$l3_matches <- renderDT({
    d <- l3_matches()
    datatable(d, rownames = FALSE, selection = "single",
              options = list(pageLength = 10)) |>
      formatRound(c("cosine", "pearson", "spearman"), 3)
  })
  output$l3_match_scatter <- renderPlot({
    d <- l3_matches()
    sel <- input$l3_matches_rows_selected
    req(length(sel) == 1)
    row <- d[sel, ]
    v1 <- l3_loadings()
    L2 <- load_loadings(con, row$other_fit); req(!is.null(L2))
    v2 <- L2[, row$other_factor]
    shared <- intersect(names(v1), names(v2))
    dd <- data.frame(a = v1[shared], b = v2[shared])
    ggplot(dd, aes(a, b)) +
      geom_point(size = 0.6, alpha = 0.5) +
      labs(x = "this factor's loadings",
           y = sprintf("matched factor %d (rank %d, seed %s)",
                       row$other_factor, row$other_rank, row$other_seed),
           title = sprintf("Loading agreement (%s = %.3f)", metric(), row[[metric()]])) +
      theme_minimal(base_size = 13)
  })

  enrich_result <- eventReactive(input$l3_enrich_go, {
    req(nav$fit, nav$factor_index)
    factor_id <- get_factor_id(con, nav$fit, nav$factor_index)
    qtype <- input$l3_enrich_type
    direction <- if (nav$method == "pca") input$l3_enrich_dir else "pos"
    cached <- enrichment_cached(con, factor_id, qtype, direction)
    if (!is.null(cached)) return(cached)

    v <- l3_loadings()
    genes <- if (qtype == "ora") {
      n <- input$l3_enrich_topn
      if (direction == "neg") names(sort(v))[seq_len(min(n, length(v)))]
      else names(sort(v, decreasing = TRUE))[seq_len(min(n, length(v)))]
    } else {
      if (direction == "neg") names(sort(v)) else names(sort(v, decreasing = TRUE))
    }
    res <- withProgress(message = "Querying g:Profiler ...", {
      tryCatch(
        gprofiler2::gost(
          query = genes, organism = "hsapiens", significant = TRUE,
          ordered_query = (qtype == "gsea"),
          correction_method = "fdr",
          sources = c("GO:BP", "GO:MF", "REAC", "KEGG", "WP")),
        error = function(e) {
          showNotification(paste("g:Profiler query failed:", conditionMessage(e)), type = "error")
          NULL
        })
    })
    enrichment_store(con, factor_id, qtype, direction, res)
    enrichment_cached(con, factor_id, qtype, direction)
  })
  output$l3_enrichment <- renderDT({
    d <- enrich_result()
    if (is.null(d) || nrow(d) == 0) {
      return(datatable(data.frame(note = "No significant terms (or query returned nothing)"), rownames = FALSE))
    }
    datatable(d, rownames = FALSE, options = list(pageLength = 15)) |>
      formatSignif("p_value", 3)
  })
}

`%||%` <- function(a, b) if (is.null(a)) b else a

shinyApp(ui, server)

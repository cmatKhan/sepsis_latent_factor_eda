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
source(file.path("R", "metadata_helpers.R"), local = TRUE)
source(file.path("R", "comparison_helpers.R"), local = TRUE)
# prepare_loadings()/pair_similarities() -- reused unmodified for the
# cross-method "Compare to PCA" views (see app/R/comparison_helpers.R)
source(file.path("..", "R", "lib", "ingest", "similarity.R"), local = TRUE)

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

# In-app additive schema upgrade for DBs built by an older ingest version --
# mirrors the columns ensure_schema() in R/lib/ingest/db.R adds, so the app
# works against a DB that hasn't been re-ingested yet.
local({
  ensure_column <- function(con, table, col, decl) {
    info <- DBI::dbGetQuery(con, sprintf("PRAGMA table_info(%s)", table))
    if (!(col %in% info$name)) DBI::dbExecute(con, sprintf("ALTER TABLE %s ADD COLUMN %s %s", table, col, decl))
  }
  DBI::dbExecute(con,
    "CREATE TABLE IF NOT EXISTS dataset_metadata_sources (
       dataset_id TEXT NOT NULL, kind TEXT NOT NULL CHECK (kind IN ('sample','feature')),
       path TEXT NOT NULL, id_col TEXT NOT NULL, registered_at TEXT,
       UNIQUE(dataset_id, kind))")
  ensure_column(con, "fits", "scores_file", "TEXT")
  ensure_column(con, "enrichment_cache", "query_size", "INTEGER")
})

#' Larger, darker text for the per-factor (Level 3) plots -- gene/term
#' labels in particular tend to be numerous and easy to lose at the
#' default ggplot text size/color.
readable_factor_theme <- function(base_size = 15) {
  theme_minimal(base_size = base_size) +
    theme(
      axis.text   = element_text(size = base_size, color = "black"),
      axis.title  = element_text(size = base_size + 1, color = "black", face = "bold"),
      plot.title  = element_text(size = base_size + 2, color = "black", face = "bold"),
      strip.text  = element_text(size = base_size - 1, color = "black", face = "bold"),
      legend.text = element_text(size = base_size - 1, color = "black")
    )
}

FACTORIZATION_METHODS <- c("pca", "nmf", "cogaps")   # loadings/scores-capable (single-fit views)
STABILITY_METHODS     <- c("nmf", "cogaps")          # multi-seed -- stability views apply
COMPARE_TO_PCA_METHODS <- c("nmf", "cogaps")         # gene-loadings space comparable to PCA

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
      } else if (m == "pca") {
        # PCA is deterministic (one fit per rank) -- no stability metric
        # applies; masking-CV rank selection is what's still meaningful
        mc <- maskcv_curve(con, ds(), "pca")
        if (nrow(mc) > 0 && any(!is.na(mc$mse))) {
          sprintf("best rank (masking-CV): %d", mc$rank[which.min(mc$mse)])
        } else "--"
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
      card(
        class = "mb-3",
        card_body(
          style = "background-color: #f4f8fb;",
          tags$b("How to explore this app: "),
          "click ", tags$b("Explore"), " on a method card below to see its stability ",
          "overview (Level 1). From there, click a point on a plot (e.g. a rank in the ",
          "seed-stability plot, a power in the WGCNA ARI heatmap) -- or use a rank/power ",
          "selector where there's no plot to click -- to drill into that specific fit ",
          "(Level 2). Within a fit, click a factor in the per-factor stability plot to ",
          "see its top genes and run functional enrichment (Level 3). Use the ",
          tags$b("breadcrumb"), " at the top of the page to navigate back up at any level."
        )
      ),
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
    if (m == "pca") {
      # PCA is deterministic (one fit per rank) -- no seed-stability or
      # cross-rank views apply; masking-CV rank selection still does, plus
      # a plain rank-select (no click-through plot left to drill via)
      navset_card_tab(
        nav_panel("Rank selection (masking-CV)",
          plotOutput("l1_maskcv", height = "420px"),
          layout_columns(col_widths = c(6, 6),
            selectInput("l1_pca_rank", "Explore rank:", choices = NULL),
            actionButton("l1_pca_go", "Explore this rank", class = "btn-primary btn-sm",
                         style = "margin-top: 24px;")))
      )
    } else if (m %in% STABILITY_METHODS) {
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
    req(nav$method %in% STABILITY_METHODS)
    seed_stability_by_rank(con, ds(), nav$method)
  })

  output$l1_stability <- renderPlot({
    d <- l1_stab_data()
    validate(need(nrow(d) > 0,
      paste0("No same-rank seed pairs for ", toupper(nav$method),
             " at any rank -- see the masking-CV and cross-rank tabs instead.")))
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
    req(nav$method %in% STABILITY_METHODS)
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
    req(nav$level == 1, nav$method %in% STABILITY_METHODS)
    d <- l1_crossrank_data()
    ranks <- sort(unique(c(d$rank_a, d$rank_b)))
    updateSelectInput(session, "l1_ref_rank", choices = ranks)
  })

  # PCA's reduced Level 1: plain rank-select + button (no click-through
  # plot left to drill via, since the stability tabs are gone)
  observe({
    req(nav$level == 1, nav$method == "pca")
    mc <- maskcv_curve(con, ds(), "pca")
    ranks <- sort(unique(mc$rank))
    if (length(ranks) == 0) ranks <- distinct_ranks(con, ds(), "pca")
    updateSelectInput(session, "l1_pca_rank", choices = ranks)
  })
  observeEvent(input$l1_pca_go, {
    req(input$l1_pca_rank)
    nav$rank <- as.integer(input$l1_pca_rank)
    nav$level <- 2
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
      panels <- list()
      if (m %in% STABILITY_METHODS) {
        panels <- c(panels, list(
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
        ))
      }
      panels <- c(panels, list(
        nav_panel("Sample scores",
          layout_columns(col_widths = c(3, 9),
            selectInput("l2_scores_fit", "Fit (seed):", choices = seed_choices),
            selectInput("l2_scores_sort", "Sort/group by metadata field:", choices = NULL)),
          DTOutput("l2_scores_table")),
        nav_panel("Factor correlation (this fit)",
          p("How redundant are this fit's own factors with each other, at the sample-score level?"),
          selectInput("l2_corr_fit", "Fit (seed):", choices = seed_choices),
          plotOutput("l2_factor_corr", height = "420px")),
        nav_panel("Metadata associations",
          p("Spearman correlation (numeric fields) / Kruskal-Wallis (categorical fields) between each factor's sample scores and every sample-metadata column. p-values BH-adjusted across the whole grid shown."),
          selectInput("l2_assoc_fit", "Fit (seed):", choices = seed_choices),
          plotOutput("l2_assoc_heatmap", height = "440px"))
      ))
      if (m %in% COMPARE_TO_PCA_METHODS) {
        panels <- c(panels, list(
          nav_panel("Compare to PCA",
            p("Is ", toupper(m), " estimating latent factors better than PCA at this rank, and how do its factors relate to PCA's?"),
            selectInput("l2_cmp_fit", "Fit (seed):", choices = seed_choices),
            h6("Reconstruction quality (masking-CV, held-out MSE)"),
            p(em("One fixed mask draw per rank -- these are point estimates, not a tested difference.")),
            tableOutput("l2_cmp_mse"),
            h6("Factor similarity (signed cosine; negative = matches PCA's opposite-signed side)"),
            plotOutput("l2_cmp_heatmap", height = "380px"),
            DTOutput("l2_cmp_table"),
            layout_columns(col_widths = c(6, 6),
              actionButton("l2_cmp_view_this", "View selected factor (this method)", class = "btn-outline-primary btn-sm"),
              actionButton("l2_cmp_view_pca", "View matched PCA factor", class = "btn-outline-primary btn-sm")),
            h6("Biological interpretability -- already-queried enrichment (Level 3), side by side"),
            p("Reflects only what's already been run via Level 3's Enrichment tab; nothing is queried automatically here."),
            layout_columns(col_widths = c(6, 6),
              tagList(strong(toupper(m)), DTOutput("l2_cmp_enrich_this")),
              tagList(strong("PCA"), DTOutput("l2_cmp_enrich_pca"))))
        ))
      }
      do.call(navset_card_tab, panels)
    } else if (m == "wgcna") {
      f <- wgcna_fits(con, ds())
      power_choices <- setNames(f$fit_id, paste0("power ", f$power))
      pca_fits_all <- DBI::dbGetQuery(con,
        "SELECT fit_id, rank FROM fits WHERE dataset_id = ? AND method = 'pca' AND status = 'ok' ORDER BY rank",
        params = list(ds()))
      pca_choices <- setNames(pca_fits_all$fit_id, paste0("PCA rank ", pca_fits_all$rank))
      navset_card_tab(
        nav_panel("Modules at this power",
          tableOutput("l2_wgcna_sizes"),
          h6("Best Hungarian-matched module at every other power (Jaccard):"),
          DTOutput("l2_wgcna_matches")),
        nav_panel("Module overlap vs another power",
          selectInput("l2_wgcna_other", "Compare with power:", choices = power_choices),
          plotOutput("l2_wgcna_jaccard", height = "440px")),
        nav_panel("Module eigengenes",
          p("Eigengene scores require re-running wgcna_grid with the current R/methods/wgcna.R (records sample ids) and re-ingesting with --overwrite; older fits show no data here."),
          layout_columns(col_widths = c(3, 9),
            selectInput("l2_wgcna_scores_fit", "Power:", choices = power_choices),
            selectInput("l2_wgcna_scores_sort", "Sort/group by metadata field:", choices = NULL)),
          DTOutput("l2_wgcna_scores_table")),
        nav_panel("Metadata associations",
          selectInput("l2_wgcna_assoc_fit", "Power:", choices = power_choices),
          plotOutput("l2_wgcna_assoc_heatmap", height = "440px")),
        nav_panel("Compare to PCA",
          p("Similarity only -- WGCNA has no masking-CV-style reconstruction MSE (ARI-based module stability isn't the same kind of quantity), and there's no existing per-module enrichment feature yet to compare biological interpretability against."),
          layout_columns(col_widths = c(6, 6),
            selectInput("l2_wgcna_cmp_fit", "Power:", choices = power_choices),
            selectInput("l2_wgcna_cmp_pca", "Compare with:", choices = pca_choices)),
          plotOutput("l2_wgcna_cmp_heatmap", height = "420px"))
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

  # ---- generic sample-scores / metadata-association helpers, shared by
  # factorization methods (factor scores) and WGCNA (module eigengenes) ----

  sample_meta_reactive <- reactive({
    req(input$dataset)
    dataset_metadata(con, ds(), "sample")
  })

  #' Populate a metadata-field selectInput generically from whatever
  #' columns the dataset's sample metadata table has (none if unregistered).
  update_meta_field_choices <- function(session, input_id) {
    m <- sample_meta_reactive()
    choices <- if (is.null(m)) character(0) else setdiff(names(m), "sample_id")
    updateSelectInput(session, input_id, choices = choices)
  }

  scores_table_output <- function(fit_id, sort_field) {
    L <- load_scores(con, fit_id)
    validate(need(!is.null(L), "No sample scores available for this fit (older WGCNA fits predate eigengene capture -- re-run and re-ingest with --overwrite)."))
    d <- as.data.frame(L)
    d <- cbind(sample_id = rownames(L), d)
    m <- sample_meta_reactive()
    if (!is.null(m) && nzchar(sort_field %||% "") && sort_field %in% names(m)) {
      d <- merge(d, m[, c("sample_id", sort_field)], by = "sample_id", all.x = TRUE)
      d <- d[order(d[[sort_field]]), ]
    }
    datatable(d, rownames = FALSE, options = list(pageLength = 15)) |>
      formatRound(setdiff(names(d), c("sample_id", sort_field)), 3)
  }

  assoc_heatmap_plot <- function(fit_id, title_prefix) {
    L <- load_scores(con, fit_id)
    validate(need(!is.null(L), "No sample scores available for this fit."))
    m <- sample_meta_reactive()
    validate(need(!is.null(m), "No sample metadata registered for this dataset (see dataset config sample_metadata_path/sample_id_col)."))
    d <- generic_association_scan(L, m, "sample_id")
    validate(need(nrow(d) > 0, "Not enough overlap between sample scores and sample metadata to test."))
    d$neg_log10_padj <- -log10(pmax(d$padj, 1e-300))
    ggplot(d, aes(field, factor(component), fill = neg_log10_padj)) +
      geom_tile() +
      geom_text(aes(label = ifelse(padj < 0.05, "*", "")), size = 5) +
      scale_fill_gradient(low = "white", high = "firebrick", name = "-log10(BH p)") +
      labs(x = NULL, y = NULL,
           title = paste(title_prefix, "-- factor x metadata association (* padj < 0.05)")) +
      theme_minimal(base_size = 13) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }

  # factorization methods
  observe({
    req(nav$level == 2, nav$method %in% FACTORIZATION_METHODS)
    update_meta_field_choices(session, "l2_scores_sort")
  })
  output$l2_scores_table <- renderDT({
    req(input$l2_scores_fit)
    scores_table_output(as.integer(input$l2_scores_fit), input$l2_scores_sort)
  })
  output$l2_factor_corr <- renderPlot({
    req(input$l2_corr_fit)
    L <- load_scores(con, as.integer(input$l2_corr_fit))
    validate(need(!is.null(L), "No sample scores available for this fit."))
    cm <- cor(L)
    d <- as.data.frame(as.table(cm))
    names(d) <- c("factor_a", "factor_b", "cor")
    ggplot(d, aes(factor_a, factor_b, fill = cor)) +
      geom_tile() + geom_text(aes(label = sprintf("%.2f", cor)), size = 3) +
      scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick", limits = c(-1, 1)) +
      labs(x = NULL, y = NULL, title = "Within-fit factor correlation (sample scores)") +
      theme_minimal(base_size = 14)
  })
  output$l2_assoc_heatmap <- renderPlot({
    req(input$l2_assoc_fit)
    assoc_heatmap_plot(as.integer(input$l2_assoc_fit), toupper(nav$method))
  })

  # ---- "Compare to PCA" (nmf/cogaps) ----

  l2_cmp_data <- reactive({
    req(input$l2_cmp_fit, nav$rank)
    compare_loadings_to_pca(con, ds(), as.integer(input$l2_cmp_fit), nav$rank)
  })
  output$l2_cmp_mse <- renderTable({
    req(input$l2_cmp_fit, nav$rank)
    this_mse <- maskcv_best_mse(con, ds(), nav$method, nav$rank)
    pca_mse  <- maskcv_best_mse(con, ds(), "pca", nav$rank)
    data.frame(
      method = c(toupper(nav$method), "PCA"),
      `held-out MSE` = c(this_mse, pca_mse),
      check.names = FALSE
    )
  })
  output$l2_cmp_heatmap <- renderPlot({
    d <- l2_cmp_data()
    validate(need(!is.null(d), sprintf("No PCA fit at rank %s to compare against.", nav$rank)))
    pca_fit_id <- attr(d, "pca_fit_id")
    L <- load_loadings(con, as.integer(input$l2_cmp_fit))
    Lp <- load_loadings(con, pca_fit_id)
    prep_a <- prepare_loadings(L); prep_b <- prepare_loadings(Lp)
    sims <- pair_similarities(prep_a, prep_b)
    grid <- expand.grid(component = seq_len(nrow(sims$cosine)), pca_component = seq_len(ncol(sims$cosine)))
    grid$cosine <- as.vector(sims$cosine)
    grid$matched <- as.integer(mapply(function(a, b) any(d$component == a & d$pca_component == b),
                                       grid$component, grid$pca_component))
    ggplot(grid, aes(factor(component), factor(pca_component), fill = cosine)) +
      geom_tile() +
      geom_tile(data = grid[grid$matched == 1, ], color = "red", linewidth = 1, fill = NA) +
      geom_text(aes(label = sprintf("%.2f", cosine)), size = 3) +
      scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick", limits = c(-1, 1)) +
      labs(x = paste(toupper(nav$method), "factor"), y = "PCA component", fill = "cosine",
           title = "Signed factor similarity (red outline = best |match|, negative = opposite-signed side)") +
      theme_minimal(base_size = 13)
  })
  output$l2_cmp_table <- renderDT({
    d <- l2_cmp_data()
    validate(need(!is.null(d), "No PCA fit at this rank."))
    datatable(d, rownames = FALSE, selection = "single", options = list(pageLength = 10)) |>
      formatRound(c("cosine", "pearson", "spearman"), 3)
  })
  observeEvent(input$l2_cmp_view_this, {
    d <- l2_cmp_data(); sel <- input$l2_cmp_table_rows_selected
    req(!is.null(d), length(sel) == 1)
    nav$fit <- as.integer(input$l2_cmp_fit)
    nav$factor_index <- d$component[sel]
    nav$level <- 3
  })
  observeEvent(input$l2_cmp_view_pca, {
    d <- l2_cmp_data(); sel <- input$l2_cmp_table_rows_selected
    req(!is.null(d), length(sel) == 1)
    pca_fit_id <- attr(d, "pca_fit_id")
    nav$method <- "pca"
    nav$fit <- pca_fit_id
    nav$factor_index <- d$pca_component[sel]
    nav$level <- 3
  })
  output$l2_cmp_enrich_this <- renderDT({
    req(input$l2_cmp_fit)
    d <- cached_enrichment_summary(con, as.integer(input$l2_cmp_fit))
    if (nrow(d) == 0) return(datatable(data.frame(note = "Nothing queried yet"), rownames = FALSE))
    datatable(d, rownames = FALSE, options = list(pageLength = 5, dom = "tp")) |> formatSignif("min_p_value", 3)
  })
  output$l2_cmp_enrich_pca <- renderDT({
    d <- l2_cmp_data()
    validate(need(!is.null(d), ""))
    pca_fit_id <- attr(d, "pca_fit_id")
    dd <- cached_enrichment_summary(con, pca_fit_id)
    if (nrow(dd) == 0) return(datatable(data.frame(note = "Nothing queried yet"), rownames = FALSE))
    datatable(dd, rownames = FALSE, options = list(pageLength = 5, dom = "tp")) |> formatSignif("min_p_value", 3)
  })

  # ---- "Compare to PCA" (wgcna) ----

  output$l2_wgcna_cmp_heatmap <- renderPlot({
    req(input$l2_wgcna_cmp_fit, input$l2_wgcna_cmp_pca)
    d <- compare_scores_to_pca(con, as.integer(input$l2_wgcna_cmp_fit), as.integer(input$l2_wgcna_cmp_pca))
    validate(need(!is.null(d), "Eigengenes unavailable for this power (see note above) or no PCA fit selected."))
    Lo <- load_scores(con, as.integer(input$l2_wgcna_cmp_fit))
    Lp <- load_scores(con, as.integer(input$l2_wgcna_cmp_pca))
    shared <- intersect(rownames(Lo), rownames(Lp))
    cm <- suppressWarnings(cor(Lo[shared, , drop = FALSE], Lp[shared, , drop = FALSE], method = "pearson"))
    grid <- expand.grid(module = seq_len(nrow(cm)), pca_component = seq_len(ncol(cm)))
    grid$cor <- as.vector(cm)
    grid$matched <- as.integer(mapply(function(a, b) any(d$component == a & d$pca_component == b),
                                       grid$module, grid$pca_component))
    ggplot(grid, aes(factor(module), factor(pca_component), fill = cor)) +
      geom_tile() +
      geom_tile(data = grid[grid$matched == 1, ], color = "red", linewidth = 1, fill = NA) +
      geom_text(aes(label = sprintf("%.2f", cor)), size = 3) +
      scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick", limits = c(-1, 1)) +
      labs(x = "WGCNA module (eigengene)", y = "PCA component (score)", fill = "pearson r",
           title = "Module eigengene vs PCA score correlation (red outline = best |match|)") +
      theme_minimal(base_size = 13)
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
  observe({
    req(nav$level == 2, nav$method == "wgcna")
    update_meta_field_choices(session, "l2_wgcna_scores_sort")
  })
  output$l2_wgcna_scores_table <- renderDT({
    req(input$l2_wgcna_scores_fit)
    scores_table_output(as.integer(input$l2_wgcna_scores_fit), input$l2_wgcna_scores_sort)
  })
  output$l2_wgcna_assoc_heatmap <- renderPlot({
    req(input$l2_wgcna_assoc_fit)
    assoc_heatmap_plot(as.integer(input$l2_wgcna_assoc_fit), "WGCNA")
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
        layout_columns(col_widths = c(6, 6),
          sliderInput("l3_topn", "Top genes:", min = 10, max = 100, value = 25, step = 5),
          selectInput("l3_label_col", "Label genes by:", choices = "feature_id")),
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
        navset_card_tab(
          nav_panel("Table", DTOutput("l3_enrichment")),
          nav_panel("Mirror bar",
            p("Both the positive- and negative-loading direction must be run/loaded (via the button above) for a full mirror; PCA only -- NMF/CoGAPS weights are non-negative, so this shows the single positive side."),
            plotOutput("l3_enrich_mirror", height = "560px")),
          nav_panel("Dot plot", plotOutput("l3_enrich_dot", height = "500px"))
        ))
    )
  }

  l3_loadings <- reactive({
    req(nav$fit, nav$factor_index)
    L <- load_loadings(con, nav$fit); req(!is.null(L))
    L[, nav$factor_index]
  })

  feature_meta_reactive <- reactive({
    req(input$dataset)
    dataset_metadata(con, ds(), "feature")
  })
  observe({
    req(nav$level == 3)
    m <- feature_meta_reactive()
    choices <- if (is.null(m)) "feature_id" else c("feature_id", setdiff(names(m), "feature_id"))
    updateSelectInput(session, "l3_label_col", choices = choices)
  })

  output$l3_loadings <- renderPlot({
    v <- l3_loadings()
    topn <- input$l3_topn %||% 25
    ord <- order(abs(v), decreasing = TRUE)[seq_len(min(topn, length(v)))]
    d <- data.frame(feature_id = names(v)[ord], loading = v[ord], stringsAsFactors = FALSE)

    label_col <- input$l3_label_col %||% "feature_id"
    m <- feature_meta_reactive()
    if (!is.null(m) && label_col %in% names(m) && label_col != "feature_id") {
      d <- merge(d, m[, c("feature_id", label_col)], by = "feature_id", all.x = TRUE)
      d$label <- ifelse(is.na(d[[label_col]]) | !nzchar(as.character(d[[label_col]])),
                        d$feature_id, as.character(d[[label_col]]))
    } else {
      d$label <- d$feature_id
    }
    d <- d[order(abs(d$loading), decreasing = TRUE), ]
    d$label <- factor(d$label, levels = rev(d$label))
    ggplot(d, aes(loading, label)) +
      geom_col(fill = "grey40") +
      labs(title = sprintf("%s rank %s seed %s -- factor %d: top %d loadings",
                           toupper(nav$method), nav$rank,
                           get_fit(con, nav$fit)$seed, nav$factor_index, nrow(d)),
           y = NULL) +
      readable_factor_theme(15)
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
      readable_factor_theme(15)
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

  # combines whatever's cached for BOTH loading directions (PCA only has a
  # meaningful negative side; NMF/CoGAPS weights are non-negative so this
  # degenerates to just "pos") -- independent of which single direction the
  # radio button currently has selected, so a mirror plot can show both
  # sides once each has been run at least once
  enrich_both <- reactive({
    input$l3_enrich_go  # re-check cache whenever a run/load happens
    req(nav$fit, nav$factor_index)
    factor_id <- get_factor_id(con, nav$fit, nav$factor_index)
    qtype <- input$l3_enrich_type
    dirs <- if (nav$method == "pca") c("pos", "neg") else "pos"
    res <- lapply(dirs, function(dir) {
      cc <- enrichment_cached(con, factor_id, qtype, dir)
      if (is.null(cc) || nrow(cc) == 0) return(NULL)
      cc$direction <- dir
      cc
    })
    res <- Filter(Negate(is.null), res)
    if (length(res) == 0) return(NULL)
    do.call(rbind, res)
  })

  output$l3_enrich_mirror <- renderPlot({
    d <- enrich_both()
    validate(need(!is.null(d) && nrow(d) > 0,
      "No cached enrichment yet -- click 'Run / load enrichment' for each direction you want to see."))
    d$log10p <- -log10(pmax(d$p_value, 1e-300))
    d$log10p <- ifelse(d$direction == "neg", -d$log10p, d$log10p)
    d <- do.call(rbind, lapply(split(d, d$source), function(s) {
      s[order(-abs(s$log10p)), ][seq_len(min(20, nrow(s))), , drop = FALSE]
    }))
    # a term can appear in both the positive- and negative-loading rows
    # (mirror bar) -- dedupe before building factor levels, or factor()
    # errors on the duplicate ("factor level [...] is duplicated")
    d$term_name <- factor(d$term_name, levels = unique(d$term_name[order(d$log10p)]))
    ggplot(d, aes(log10p, term_name, fill = direction)) +
      geom_col() +
      geom_vline(xintercept = 0, linewidth = 0.3) +
      scale_fill_manual(values = c(pos = "#d73027", neg = "#4575b4")) +
      facet_wrap(~source, scales = "free_y") +
      labs(x = "-log10(p)  [negative direction flipped]", y = NULL, fill = "loadings",
           title = sprintf("%s factor %d -- %s (mirror bar)",
                           toupper(nav$method), nav$factor_index, toupper(input$l3_enrich_type))) +
      readable_factor_theme(14)
  })

  output$l3_enrich_dot <- renderPlot({
    d <- enrich_both()
    validate(need(!is.null(d) && nrow(d) > 0,
      "No cached enrichment yet -- click 'Run / load enrichment' first."))
    d <- d[!is.na(d$query_size) & d$query_size > 0, , drop = FALSE]
    validate(need(nrow(d) > 0,
      "Gene-ratio dot plot needs query_size, which pre-existing cached rows (from before this feature) don't have -- click 'Run / load enrichment' again to refresh."))
    d$gene_ratio <- d$intersection_size / d$query_size
    d$log10p <- -log10(pmax(d$p_value, 1e-300))
    d <- do.call(rbind, lapply(split(d, interaction(d$source, d$direction, drop = TRUE)), function(s) {
      s[order(-s$log10p), ][seq_len(min(20, nrow(s))), , drop = FALSE]
    }))
    ggplot(d, aes(gene_ratio, reorder(term_name, gene_ratio), color = log10p, size = intersection_size)) +
      geom_point() +
      scale_color_viridis_c(name = "-log10(p)") +
      facet_wrap(~source + direction, scales = "free_y") +
      labs(x = "Gene ratio", y = NULL,
           title = sprintf("%s factor %d -- %s (dot plot)",
                           toupper(nav$method), nav$factor_index, toupper(input$l3_enrich_type))) +
      readable_factor_theme(14)
  })
}

`%||%` <- function(a, b) if (is.null(a)) b else a

shinyApp(ui, server)

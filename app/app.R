# Latent-factor stability explorer.
#
# Reads whatever the stability DB currently contains -- fully decoupled
# from ingest (R/ingest_results.R). DB path is parameterized via the
# STABILITY_DB environment variable, falling back to results/stability.sqlite
# relative to the project root.
#
# Drill-down levels (breadcrumb at top navigates back up):
#   0 dataset overview -> 1 method -> 2 rank/parameter -> 3 factor
# Each breadcrumb segment past "Home" is a dropdown of sibling options at
# that level, so switching e.g. rank or factor jumps straight there; an
# "up one level" link next to it collapses back up without changing values.

library(shiny)
library(bslib)
library(DBI)
library(RSQLite)
library(ggplot2)
library(DT)

source(file.path("R", "db_helpers.R"), local = TRUE)
source(file.path("R", "metadata_helpers.R"), local = TRUE)
source(file.path("R", "comparison_helpers.R"), local = TRUE)
# prepare_loadings()/pair_similarities() -- reused unmodified for
# the standalone "Compare methods" screen (see app/R/comparison_helpers.R)
source(file.path("..", "R", "lib", "ingest", "similarity.R"), local = TRUE)
# build_ensembl_map()/remap_to_ensembl() -- pure functions, safe to reuse
# unmodified (see app/R/metadata_helpers.R::ensembl_map_for_dataset(),
# which supplies build_ensembl_map()'s input from the DB rather than
# config/*.yml, keeping this app decoupled from ingest).
source(file.path("..", "R", "lib", "ingest", "symbol_mapping.R"), local = TRUE)

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
  # dataset.ensembl_col -- see R/lib/ingest/db.R's matching migration and
  # app/R/metadata_helpers.R::ensembl_map_for_dataset(), which reads this
  # to build the app's on-demand enrichment Ensembl remap.
  ensure_column(con, "dataset_metadata_sources", "ensembl_col", "TEXT")
  ensure_column(con, "dataset_metadata_sources", "symbol_col", "TEXT")
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

FACTORIZATION_METHODS <- c("pca", "nmf", "cogaps", "spca", "ica", "cp", "tucker")   # loadings/scores-capable
STABILITY_METHODS     <- c("nmf", "cogaps", "ica")   # multi-seed -- stability views apply
# PCA/sPCA are both deterministic given (rank [+ sPCA's para]) -- no
# cross-seed stability question -- so Level 1 gets a scree-style
# reconstruction-error-by-rank plot + "Explore rank" dropdown (see
# level1_ui()) rather than either the seed-stability tabs
# (STABILITY_METHODS) or the flat fit list (DIRECT_FIT_METHODS below).
# (A masking-CV rank-selection family was planned here at one point but
# never actually wired into the config-driven method-registry system --
# no `_maskcv` job has ever been generated for any dataset -- so it was
# removed entirely rather than left as a permanently-empty panel.)
SCREE_RANK_METHODS    <- c("pca", "spca")
DIRECT_FIT_METHODS    <- c("cp", "tucker")           # no rank-selection UI at all -- Level 1 lists fits directly (like WGCNA's power list)
TENSOR_METHODS        <- c("cp", "tucker")           # third (time-mode) loading matrix; "scores" is SUBJECT-mode, not sample-mode
# Methods eligible as a side in the standalone "Compare methods" screen's
# gene-loadings (Jaccard/cosine) comparison -- excludes "wgcna" from this
# list specifically (it has no continuous loadings) but wgcna still
# participates via its own module-gene-set path (see gene_sets_from_wgcna()).
LOADINGS_METHODS      <- c("pca", "nmf", "cogaps", "spca", "ica", "cp", "tucker")
# Methods eligible on the sample-level (scores-based) side of that screen --
# excludes CP/Tucker (subject-mode scores, not sample-mode; see the
# TENSOR_METHODS note throughout this file).
SAMPLE_SCORE_METHODS  <- c("pca", "nmf", "cogaps", "spca", "ica", "wgcna")
# Methods with signed loadings -- both loading directions are meaningful
# to query separately (NMF/CoGAPS weights are non-negative, so only "pos"
# applies there). Must match R/ingest_jobs/fgsea_job.R's own `dirs <-
# if (method %in% c(...)) c("pos","neg") else "pos"` exactly, or ingested
# neg-direction ORA/GSEA rows for ica/spca are computed but unreachable
# from the app.
NEG_DIRECTION_METHODS <- c("pca", "ica", "spca")
# Methods R/create_ingest_slurm_bundle.R's --stage enrichment stages
# fgsea_grid for -- every loadings-bearing method, including cp/tucker
# (representative_fit_ids() falls through to "every ok fit" for those,
# same as sPCA used to before it got its own per-K collapse). WGCNA is
# NOT in this list -- it has no loadings to rank, so it gets its own
# separate batch job family instead (R/ingest_jobs/wgcna_ora_job.R,
# whole-module ORA only, no GSEA).
FGSEA_GRID_METHODS    <- c("pca", "nmf", "cogaps", "spca", "ica", "cp", "tucker")

#' Enrichment query-type radio choices for one method -- ALL of "ora"/
#' "gsea"/"fgsea"/"cogaps_fora" are computed exclusively by the cluster
#' ingest pipeline now (R/ingest_jobs/fgsea_job.R + R/ingest_enrichment_
#' results.R) -- this app never computes enrichment live, only reads
#' whatever's already cached (see enrichment_cached()). "fgsea" (batch
#' Hallmark GSEA) and "gsea"/"ora" apply to every method fgsea_grid
#' covers; "cogaps_fora" (batch CoGAPS marker-gene ORA) only exists for
#' CoGAPS fits (see R/ingest_jobs/fgsea_job.R).
enrich_query_choices <- function(method) {
  if (!(method %in% FGSEA_GRID_METHODS)) return(c("ORA" = "ora", "GSEA" = "gsea"))
  choices <- c("ORA" = "ora", "GSEA" = "gsea", "Hallmark GSEA (batch)" = "fgsea")
  if (identical(method, "cogaps")) choices <- c(choices, "CoGAPS marker ORA (batch)" = "cogaps_fora")
  choices
}

## ---- ui ------------------------------------------------------------------

ui <- page_fillable(
  theme = bs_theme(bootswatch = "flatly"),
  padding = c(8, 12, 8, 12),
  gap = "0.4rem",
  layout_columns(
    col_widths = c(3, 2, 2, 5),
    # fill = FALSE keeps this control strip sized to its own (short)
    # content instead of splitting the page's height 50/50 with
    # uiOutput("level_ui") below it -- the latter is what was pushing
    # all the actual tab content (plots, rank selectors, etc.) down
    # into the bottom half of the screen with a wall of empty space
    # above it.
    fill = FALSE,
    selectInput("dataset", NULL, choices = NULL, width = "100%"),
    radioButtons("metric", NULL, choices = METRICS, selected = "cosine", inline = TRUE),
    uiOutput("mode_toggle"),
    uiOutput("breadcrumb")
  ),
  # flex: 1 1 auto + min-height: 0 makes this div claim all remaining
  # vertical space in the page_fillable layout (now that the control
  # strip above is fill = FALSE) while still allowing its own content to
  # scroll internally when taller than the viewport -- min-height: 0 is
  # required for that: a flex child otherwise refuses to shrink below its
  # content's natural height, which silently defeats overflow-y: auto.
  div(style = "flex: 1 1 auto; min-height: 0; overflow-y: auto;",
      uiOutput("level_ui"))
)

## ---- server ---------------------------------------------------------------

server <- function(input, output, session) {

  nav <- reactiveValues(level = 0, method = NULL, rank = NULL,
                         fit = NULL, factor_index = NULL,
                         wgcna_fit = NULL, mode = "explore")

  datasets <- list_datasets(con)
  updateSelectInput(session, "dataset", choices = datasets$dataset_id)

  metric <- reactive(check_metric(input$metric))
  ds <- reactive({ req(input$dataset); input$dataset })

  observeEvent(input$dataset, {
    nav$level <- 0; nav$method <- NULL; nav$rank <- NULL
    nav$fit <- NULL; nav$factor_index <- NULL; nav$mode <- "explore"
  })

  ## ---- explore/compare mode toggle -----------------------------------------

  output$mode_toggle <- renderUI({
    if (nav$mode == "explore") {
      actionButton("goto_compare", "Compare methods", class = "btn-outline-primary btn-sm")
    } else {
      actionButton("goto_explore", "Back to explore", class = "btn-outline-secondary btn-sm")
    }
  })
  observeEvent(input$goto_compare, { nav$mode <- "compare" })
  observeEvent(input$goto_explore, { nav$mode <- "explore" })

  ## ---- breadcrumb (dropdowns for lateral switching, "up" to collapse) ------

  output$breadcrumb <- renderUI({
    if (nav$mode == "compare") {
      return(div(style = "font-size: 1.1rem; padding-top: 6px;", tags$b(paste(ds(), "\u00bb Compare methods"))))
    }
    sep <- HTML("&nbsp;&raquo;&nbsp;")
    crumbs <- list(actionLink("bc_home", ds()))

    if (nav$level >= 1) {
      present <- sort(unique(method_overview(con, ds())$counts$method))
      method_choices <- intersect(c(FACTORIZATION_METHODS, "wgcna"), present)
      crumbs <- c(crumbs, list(sep, selectInput("bc_method", NULL,
        choices = setNames(method_choices, toupper(method_choices)),
        selected = nav$method, width = "130px")))
    }
    if (nav$level >= 2 && nav$method %in% DIRECT_FIT_METHODS) {
      f <- direct_fits(con, ds(), nav$method)
      fit_labels <- vapply(f$fit_id, function(id) fit_descriptor(con, nav$method, id), character(1))
      crumbs <- c(crumbs, list(sep, selectInput("bc_rank", NULL,
        choices = setNames(f$fit_id, fit_labels), selected = nav$fit, width = "220px")))
    } else if (nav$level >= 2 && nav$method %in% FACTORIZATION_METHODS) {
      ranks <- distinct_ranks(con, ds(), nav$method)
      crumbs <- c(crumbs, list(sep, selectInput("bc_rank", NULL,
        choices = setNames(ranks, paste0("rank ", ranks)), selected = nav$rank, width = "130px")))
    }
    if (nav$level >= 2 && nav$method == "wgcna") {
      f <- wgcna_fits(con, ds())
      crumbs <- c(crumbs, list(sep, selectInput("bc_rank", NULL,
        choices = setNames(f$fit_id, paste0("power ", f$power)), selected = nav$wgcna_fit, width = "150px")))
    }
    if (nav$level >= 3) {
      if (nav$method == "wgcna") {
        idxs <- wgcna_module_sizes(con, nav$wgcna_fit)$module
        unit <- "module"
      } else {
        L <- load_loadings(con, nav$fit)
        idxs <- if (!is.null(L)) seq_len(ncol(L)) else nav$factor_index
        unit <- "factor"
      }
      crumbs <- c(crumbs, list(sep, selectInput("bc_factor", NULL,
        choices = setNames(idxs, paste(unit, idxs)), selected = nav$factor_index, width = "120px")))
    }
    if (nav$level > 0) {
      crumbs <- c(crumbs, list(HTML("&nbsp;&nbsp;"), actionLink("bc_up", HTML("&uarr; up one level"))))
    }
    div(style = "font-size: 1.1rem; display: flex; align-items: center; gap: 2px; flex-wrap: wrap;", crumbs)
  })

  observeEvent(input$bc_home, { nav$level <- 0 })
  observeEvent(input$bc_up, { nav$level <- max(0, nav$level - 1) })

  observeEvent(input$bc_method, {
    req(input$bc_method)
    if (!identical(input$bc_method, nav$method)) {
      nav$method <- input$bc_method
      nav$rank <- NULL; nav$fit <- NULL; nav$factor_index <- NULL; nav$wgcna_fit <- NULL
      nav$level <- 1
    }
  })

  observeEvent(input$bc_rank, {
    req(input$bc_rank, nav$method)
    if (nav$method %in% DIRECT_FIT_METHODS) {
      fid <- as.integer(input$bc_rank)
      if (!identical(fid, nav$fit)) {
        nav$fit <- fid
        nav$rank <- get_fit(con, fid)$rank
        nav$factor_index <- NULL
        nav$level <- 2
      }
    } else if (nav$method == "wgcna") {
      fid <- as.integer(input$bc_rank)
      if (!identical(fid, nav$wgcna_fit)) {
        nav$wgcna_fit <- fid
        nav$factor_index <- NULL
        nav$level <- 2
      }
    } else {
      rk <- as.integer(input$bc_rank)
      if (!identical(rk, nav$rank)) {
        nav$rank <- rk
        nav$fit <- NULL; nav$factor_index <- NULL
        nav$level <- 2
      }
    }
  })

  observeEvent(input$bc_factor, {
    req(input$bc_factor)
    fi <- as.integer(input$bc_factor)
    if (!identical(fi, nav$factor_index)) {
      nav$factor_index <- fi
      nav$level <- 3
    }
  })

  ## ---- level dispatcher ------------------------------------------------------

  output$level_ui <- renderUI({
    if (nav$mode == "compare") return(compare_ui())
    switch(as.character(nav$level),
      "0" = level0_ui(),
      "1" = level1_ui(),
      "2" = level2_ui(),
      "3" = level3_ui())
  })

  ## =============================================================================
  ## COMPARE METHODS -- standalone cross-method comparison (not part of
  ## the dataset > method > rank > factor drill-down; see nav$mode)
  ## =============================================================================

  compare_ui <- function() {
    navset_card_tab(
      nav_panel("Setup",
        p("Pick any two fits -- same method or different -- to compare. ",
          "WGCNA participates as a gene-set side (its modules) in the ",
          "Jaccard comparison, but not in the sample-level comparison ",
          "(module eigengenes work fine there too, actually -- see below)."),
        layout_columns(col_widths = c(6, 6),
          card(card_header(tags$b("Side A")), card_body(
            selectInput("cmp_method_a", "Method:", choices = NULL),
            selectInput("cmp_fit_a", "Fit:", choices = NULL))),
          card(card_header(tags$b("Side B")), card_body(
            selectInput("cmp_method_b", "Method:", choices = NULL),
            selectInput("cmp_fit_b", "Fit:", choices = NULL)))
        )),
      nav_panel("Gene-set (Jaccard) similarity",
        p("Jaccard overlap of the top-N most strongly weighted genes per factor -- a different, more scale-robust metric than the continuous cosine similarity used elsewhere in this app. For WGCNA, \"top genes\" is a module's full membership (N/sign controls below don't apply to that side)."),
        layout_columns(col_widths = c(4, 4, 4),
          sliderInput("cmp_topn", "Top N genes (ignored for WGCNA side):", min = 10, max = 300, value = 50, step = 10),
          selectInput("cmp_sign_a", "Side A sign:", c("both", "positive", "negative")),
          selectInput("cmp_sign_b", "Side B sign:", c("both", "positive", "negative"))),
        plotOutput("cmp_jaccard_heatmap", height = "480px"),
        h6("Best (Hungarian, by |Jaccard|) match per Side-A factor:"),
        DTOutput("cmp_jaccard_table"),
        layout_columns(col_widths = c(6, 6),
          actionButton("cmp_view_a", "View selected factor (Side A)", class = "btn-outline-primary btn-sm"),
          actionButton("cmp_view_b", "View selected factor (Side B)", class = "btn-outline-primary btn-sm"))),
      nav_panel("Sample-level comparison",
        p("How do samples project onto each side's factors? Correlation between the two score matrices, plus independent clustering of each side's samples compared via Adjusted Rand Index. CP/Tucker excluded (subject-mode scores, not sample-mode)."),
        plotOutput("cmp_scores_heatmap", height = "420px"),
        sliderInput("cmp_k", "Number of clusters (k), applied to both sides:", min = 2, max = 10, value = 4),
        verbatimTextOutput("cmp_ari"),
        h6("Cluster A x cluster B contingency table:"),
        tableOutput("cmp_crosstab")),
      nav_panel("Enrichment cross-reference",
        p("Reflects whatever enrichment has already been computed for these fits by the cluster ingest pipeline (R/ingest_jobs/fgsea_job.R + R/ingest_enrichment_results.R) -- all enrichment computation happens there now, never live in this app."),
        layout_columns(col_widths = c(6, 6),
          tagList(strong("Side A"), DTOutput("cmp_enrich_a")),
          tagList(strong("Side B"), DTOutput("cmp_enrich_b"))))
    )
  }

  observe({
    req(nav$mode == "compare", input$dataset)
    present <- sort(unique(method_overview(con, ds())$counts$method))
    choices <- intersect(c(LOADINGS_METHODS, "wgcna"), present)
    updateSelectInput(session, "cmp_method_a", choices = choices)
    updateSelectInput(session, "cmp_method_b", choices = choices,
                       selected = if (length(choices) > 1) choices[2] else choices[1])
  })
  observe({
    req(nav$mode == "compare", input$cmp_method_a)
    ids <- all_fits_for_compare(con, ds(), input$cmp_method_a)
    labels <- vapply(ids, function(id) fit_descriptor(con, input$cmp_method_a, id), character(1))
    updateSelectInput(session, "cmp_fit_a", choices = setNames(ids, labels))
  })
  observe({
    req(nav$mode == "compare", input$cmp_method_b)
    ids <- all_fits_for_compare(con, ds(), input$cmp_method_b)
    labels <- vapply(ids, function(id) fit_descriptor(con, input$cmp_method_b, id), character(1))
    updateSelectInput(session, "cmp_fit_b", choices = setNames(ids, labels))
  })

  cmp_method_a <- reactive({ req(input$cmp_method_a); input$cmp_method_a })
  cmp_method_b <- reactive({ req(input$cmp_method_b); input$cmp_method_b })
  cmp_fit_a <- reactive({ req(input$cmp_fit_a); as.integer(input$cmp_fit_a) })
  cmp_fit_b <- reactive({ req(input$cmp_fit_b); as.integer(input$cmp_fit_b) })

  # ---- gene-set (Jaccard) comparison ----

  cmp_sets_a <- reactive({
    comparator_gene_sets(con, cmp_method_a(), cmp_fit_a(), input$cmp_topn %||% 50, input$cmp_sign_a %||% "both")
  })
  cmp_sets_b <- reactive({
    comparator_gene_sets(con, cmp_method_b(), cmp_fit_b(), input$cmp_topn %||% 50, input$cmp_sign_b %||% "both")
  })
  cmp_jaccard <- reactive({
    a <- cmp_sets_a(); b <- cmp_sets_b()
    validate(need(!is.null(a) && !is.null(b) && length(a$sets) > 0 && length(b$sets) > 0,
                  "No gene sets available for one or both sides (missing loadings/module artifact?)."))
    jm <- jaccard_matrix(setNames(a$sets, a$label), setNames(b$sets, b$label))
    list(matrix = jm, a = a, b = b)
  })
  cmp_best_matches <- reactive({
    cj <- cmp_jaccard()
    match_idx <- hungarian_match_abs(cj$matrix)
    data.frame(
      side_a = cj$a$label[match_idx[, "a"]],
      side_b = cj$b$label[match_idx[, "b"]],
      jaccard = cj$matrix[match_idx],
      a_factor_index = cj$a$factor_index[match_idx[, "a"]],
      b_factor_index = cj$b$factor_index[match_idx[, "b"]]
    )
  })
  output$cmp_jaccard_heatmap <- renderPlot({
    cj <- cmp_jaccard()
    jm <- cj$matrix
    bm <- cmp_best_matches()
    grid <- expand.grid(a = seq_len(nrow(jm)), b = seq_len(ncol(jm)))
    grid$jaccard <- as.vector(jm)
    grid$a_lab <- factor(rownames(jm)[grid$a], levels = rownames(jm))
    grid$b_lab <- factor(colnames(jm)[grid$b], levels = colnames(jm))
    grid$matched <- as.integer(mapply(function(al, bl) any(bm$side_a == al & bm$side_b == bl),
                                       as.character(grid$a_lab), as.character(grid$b_lab)))
    ggplot(grid, aes(a_lab, b_lab, fill = jaccard)) +
      geom_tile() +
      geom_tile(data = grid[grid$matched == 1, ], color = "red", linewidth = 1, fill = NA) +
      geom_text(aes(label = sprintf("%.2f", jaccard)), size = 3) +
      scale_fill_viridis_c(limits = c(0, 1)) +
      labs(x = toupper(cmp_method_a()), y = toupper(cmp_method_b()), fill = "Jaccard",
           title = "Top-N gene-set Jaccard overlap (red outline = best |match| per Side-A factor)") +
      theme_minimal(base_size = 13) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  output$cmp_jaccard_table <- renderDT({
    d <- cmp_best_matches()
    datatable(d[, c("side_a", "side_b", "jaccard")], rownames = FALSE, selection = "single",
              options = list(pageLength = 10)) |> formatRound("jaccard", 3)
  })
  observeEvent(input$cmp_view_a, {
    d <- cmp_best_matches(); sel <- input$cmp_jaccard_table_rows_selected
    req(length(sel) == 1)
    nav$method <- cmp_method_a(); nav$fit <- cmp_fit_a()
    nav$factor_index <- d$a_factor_index[sel]
    nav$mode <- "explore"; nav$level <- 3
  })
  observeEvent(input$cmp_view_b, {
    d <- cmp_best_matches(); sel <- input$cmp_jaccard_table_rows_selected
    req(length(sel) == 1)
    nav$method <- cmp_method_b(); nav$fit <- cmp_fit_b()
    nav$factor_index <- d$b_factor_index[sel]
    nav$mode <- "explore"; nav$level <- 3
  })

  # ---- sample-level comparison (correlation + clustering/ARI) ----

  cmp_scores_pair <- reactive({
    req(cmp_method_a() %in% SAMPLE_SCORE_METHODS, cmp_method_b() %in% SAMPLE_SCORE_METHODS)
    compare_scores(con, cmp_fit_a(), cmp_fit_b())
  })
  output$cmp_scores_heatmap <- renderPlot({
    cs <- cmp_scores_pair()
    validate(need(!is.null(cs),
      "No sample-score overlap (or CP/Tucker involved -- excluded here, subject-mode scores)."))
    cm <- cs$matrix
    grid <- expand.grid(a = seq_len(nrow(cm)), b = seq_len(ncol(cm)))
    grid$cor <- as.vector(cm)
    ggplot(grid, aes(factor(a), factor(b), fill = cor)) +
      geom_tile() + geom_text(aes(label = sprintf("%.2f", cor)), size = 3) +
      scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick", limits = c(-1, 1)) +
      labs(x = paste(toupper(cmp_method_a()), "factor"), y = paste(toupper(cmp_method_b()), "factor"),
           fill = "pearson r", title = "Sample-score correlation between the two fits") +
      theme_minimal(base_size = 13)
  })
  cmp_cluster <- reactive({
    req(cmp_method_a() %in% SAMPLE_SCORE_METHODS, cmp_method_b() %in% SAMPLE_SCORE_METHODS, input$cmp_k)
    cluster_and_ari(con, cmp_fit_a(), cmp_fit_b(), k = input$cmp_k)
  })
  output$cmp_ari <- renderText({
    ca <- cmp_cluster()
    validate(need(!is.null(ca), "Not enough shared samples to cluster."))
    sprintf("Adjusted Rand Index between the two methods' %d-cluster sample groupings: %.3f (n = %d shared samples)",
            input$cmp_k, ca$ari, ca$n_shared)
  })
  output$cmp_crosstab <- renderTable({
    ca <- cmp_cluster(); req(!is.null(ca))
    cbind(`cluster (A)` = rownames(ca$crosstab), ca$crosstab)
  }, rownames = FALSE)

  # ---- enrichment cross-reference ----
  # Purely a read of whatever the cluster ingest pipeline has already
  # computed (R/ingest_jobs/fgsea_job.R + R/ingest_enrichment_results.R,
  # or R/ingest_jobs/wgcna_ora_job.R for WGCNA) -- no live compute, no
  # refresh trigger needed; these just re-render whenever the selected
  # fits change.

  output$cmp_enrich_a <- renderDT({
    d <- cached_enrichment_summary(con, cmp_fit_a())
    if (nrow(d) == 0) return(datatable(data.frame(note = "Nothing computed yet"), rownames = FALSE))
    datatable(d, rownames = FALSE, options = list(pageLength = 8, dom = "tp")) |> formatSignif("min_p_value", 3)
  })
  output$cmp_enrich_b <- renderDT({
    d <- cached_enrichment_summary(con, cmp_fit_b())
    if (nrow(d) == 0) return(datatable(data.frame(note = "Nothing computed yet"), rownames = FALSE))
    datatable(d, rownames = FALSE, options = list(pageLength = 8, dom = "tp")) |> formatSignif("min_p_value", 3)
  })

  ## =============================================================================
  ## LEVEL 0 -- dataset overview
  ## =============================================================================

  level0_ui <- function() {
    ov <- method_overview(con, ds())
    methods <- sort(unique(ov$counts$method))
    cards <- lapply(methods, function(m) {
      cnt <- ov$counts[ov$counts$method == m, ]
      headline <- if (m == "wgcna") {
        if (nrow(ov$wgcna) > 0 && !is.na(ov$wgcna$ari)) sprintf("mean ARI %.2f", ov$wgcna$ari) else "--"
      } else if (m %in% SCREE_RANK_METHODS) {
        # PCA/sPCA are deterministic given their rank (+ sPCA's para) --
        # no cross-seed stability metric applies; lowest in-sample
        # reconstruction error (the scree-plot minimum) is what's still
        # meaningful
        mc <- scree_mse_by_rank(con, ds(), m)
        if (nrow(mc) > 0 && any(!is.na(mc$mse))) {
          sprintf("best rank (lowest MSE): %d", mc$rank[which.min(mc$mse)])
        } else "--"
      } else if (m %in% DIRECT_FIT_METHODS) {
        # CP/Tucker: deterministic, no masking-CV family at all --
        # headline off the ordinary (in-sample) reconstruction MSE instead
        bf <- best_direct_fit(con, ds(), m)
        if (!is.null(bf)) sprintf("best fit (lowest MSE): %s", fit_descriptor(con, m, bf$fit_id)) else "--"
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
          "see its top genes and run functional enrichment (Level 3), or use a fit's ",
          "\"Enrichment overview\" tab at Level 2 to run every factor at once. Use the ",
          tags$b("breadcrumb"), " at the top of the page to jump straight to a different ",
          "method/rank/factor via its dropdown, or the adjacent \"up one level\" link to ",
          "step back without changing the current selection."
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
  for (m in c(FACTORIZATION_METHODS, "wgcna")) {
    local({
      m_local <- m
      observeEvent(input[[paste0("open_", m_local)]], {
        nav$method <- m_local
        nav$level <- 1
      }, ignoreInit = TRUE)
    })
  }

  ## =============================================================================
  ## LEVEL 1 -- method view
  ## =============================================================================

  level1_ui <- function() {
    m <- nav$method
    if (m %in% SCREE_RANK_METHODS) {
      # PCA/sPCA are deterministic given their rank (+ sPCA's para) -- no
      # seed-stability or cross-rank views apply; a scree-style
      # reconstruction-error-by-rank plot (in-sample MSE, already
      # populated per rank in `fits` for every pca/spca fit) stands in for
      # rank selection, plus a plain rank-select (no click-through plot
      # left to drill via)
      navset_card_tab(
        nav_panel("Reconstruction error by rank",
          plotOutput("l1_scree", height = "420px"),
          layout_columns(col_widths = c(6, 6),
            selectInput("l1_scree_rank", "Explore rank:", choices = NULL),
            actionButton("l1_scree_go", "Explore this rank", class = "btn-primary btn-sm",
                         style = "margin-top: 24px;")))
      )
    } else if (m %in% STABILITY_METHODS) {
      navset_card_tab(
        nav_panel("Seed stability vs rank",
          p("Distribution of matched factor similarity across all seed pairs, per rank. Click a rank to drill in."),
          plotOutput("l1_stability", click = "l1_stability_click", height = "420px")),
        nav_panel("Reconstruction error across seeds",
          p("Mean reconstruction error (in-sample, seed-sweep family) +/- 1 SD per rank -- a stability question distinct from masking-CV's single fixed mask draw."),
          plotOutput("l1_seed_mse", height = "420px")),
        nav_panel("Cross-rank persistence",
          p("Mean matched similarity between factors found at different ranks (all seed pairs)."),
          plotOutput("l1_crossrank", height = "380px"),
          layout_columns(col_widths = c(3, 9),
            selectInput("l1_ref_rank", "Track factors from rank:", choices = NULL),
            plotOutput("l1_trajectory", height = "320px")))
      )
    } else if (m %in% DIRECT_FIT_METHODS) {
      # sPCA/CP/Tucker: no seed dimension to sweep -- every parameter
      # combination is a single deterministic fit, so Level 1 just lists
      # them directly (same idiom as WGCNA's power list) rather than the
      # rank-then-seed drill-down the other factorization methods use.
      navset_card_tab(
        nav_panel("All fits",
          p("No seed sweep for this method -- each parameter combination below is a single deterministic fit."),
          DTOutput("l1_direct_table"),
          layout_columns(col_widths = c(6, 6),
            selectInput("l1_direct_fit", "Explore fit:", choices = NULL),
            actionButton("l1_direct_go", "Explore this fit", class = "btn-primary btn-sm",
                         style = "margin-top: 24px;")))
      )
    } else if (m == "wgcna") {
      navset_card_tab(
        nav_panel("Module stability (ARI)",
          p("Adjusted Rand Index between module assignments at each pair of soft-threshold powers. Click a diagonal-adjacent cell or use the selector to drill into a power."),
          plotOutput("l1_wgcna_ari", height = "420px"),
          selectInput("l1_wgcna_power", "Drill into power:", choices = NULL),
          actionButton("l1_wgcna_go", "Explore power", class = "btn-primary btn-sm")),
        nav_panel("Modules vs power",
          plotOutput("l1_wgcna_counts", height = "420px")),
        nav_panel("Scale-free topology fit",
          p("WGCNA::pickSoftThreshold()'s diagnostic: does each candidate power actually produce a scale-free network? A signed R² near/above the usual 0.9 reference line (dashed) is the conventional justification for a power choice -- this pipeline otherwise only judges power by the module-stability views in the other two tabs."),
          plotOutput("l1_wgcna_sft", height = "420px"))
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

  output$l1_scree <- renderPlot({
    d <- scree_mse_by_rank(con, ds(), nav$method); req(nrow(d) > 0)
    if (all(is.na(d$alpha))) {
      ggplot(d, aes(rank, mse)) + geom_line() + geom_point(size = 2) +
        labs(title = paste(toupper(nav$method), "-- in-sample reconstruction error by rank (scree plot)"),
             x = "rank", y = "reconstruction MSE") +
        theme_minimal(base_size = 14)
    } else {
      # sPCA: `alpha` stores the `para` sparsity penalty crossed with rank
      # (see fits_at_rank()'s doc) -- facet/color by it since multiple
      # values exist per rank.
      ggplot(d, aes(factor(rank), factor(alpha), fill = mse)) +
        geom_tile() +
        geom_text(aes(label = ifelse(is.na(mse), "failed", sprintf("%.2f", mse))), size = 3) +
        scale_fill_viridis_c(na.value = "grey70", direction = -1) +
        labs(title = paste(toupper(nav$method), "-- in-sample reconstruction MSE by rank x alpha (para)"),
             x = "rank", y = "alpha (para)") +
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

  # PCA/sPCA's reduced Level 1: plain rank-select + button (no
  # click-through plot left to drill via, since the stability tabs are
  # gone). Defaults to the lowest-reconstruction-MSE rank -- always
  # overridable via the same control.
  observe({
    req(nav$level == 1, nav$method %in% SCREE_RANK_METHODS)
    mc <- scree_mse_by_rank(con, ds(), nav$method)
    ranks <- sort(unique(mc$rank))
    if (length(ranks) == 0) ranks <- distinct_ranks(con, ds(), nav$method)
    best <- scree_best_rank(con, ds(), nav$method)
    updateSelectInput(session, "l1_scree_rank", choices = ranks,
                       selected = if (!is.na(best) && best %in% ranks) best else ranks[1])
  })
  observeEvent(input$l1_scree_go, {
    req(input$l1_scree_rank)
    nav$rank <- as.integer(input$l1_scree_rank)
    nav$level <- 2
  })

  output$l1_seed_mse <- renderPlot({
    d <- seed_sweep_mse_by_rank(con, ds(), nav$method); req(nrow(d) > 0)
    ggplot(d, aes(rank, mean_mse)) +
      geom_errorbar(aes(ymin = mean_mse - sd_mse, ymax = mean_mse + sd_mse), width = 0.3, color = "grey50") +
      geom_line() + geom_point(size = 2) +
      labs(x = "rank", y = "in-sample reconstruction MSE (mean +/- 1 SD across seeds)",
           title = paste(toupper(nav$method), "-- reconstruction error variability across seeds, by rank")) +
      theme_minimal(base_size = 14)
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
  output$l1_wgcna_sft <- renderPlot({
    d <- wgcna_sft_fit(con, ds())
    validate(need(nrow(d) > 0,
      "No scale-free-topology fit computed yet for this dataset -- re-ingest (R/ingest_results.R or R/cache_dataset_matrices.R) to populate it."))
    # signed R^2 (WGCNA's own diagnostic-plot convention): flips sign so
    # the curve rises toward 1 for a well-behaved (positive-connectivity-
    # correlated) fit instead of just reporting the unsigned R^2.
    d$signed_r_sq <- -sign(d$slope) * d$sft_r_sq
    long <- rbind(
      data.frame(power = d$power, panel = "Scale independence (signed R²)", value = d$signed_r_sq),
      data.frame(power = d$power, panel = "Mean connectivity", value = d$mean_k)
    )
    ggplot(long, aes(power, value)) +
      geom_line() + geom_point(size = 2) +
      geom_hline(data = data.frame(panel = "Scale independence (signed R²)", yint = 0.9),
                 aes(yintercept = yint), color = "red", linetype = "dashed") +
      facet_wrap(~panel, scales = "free_y") +
      labs(x = "soft-threshold power", y = NULL,
           title = "WGCNA::pickSoftThreshold() diagnostic (dashed line: conventional R² = 0.9 reference)") +
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

  # sPCA/CP/Tucker level 1 -- direct fit list (see level1_ui())
  output$l1_direct_table <- renderDT({
    req(nav$method %in% DIRECT_FIT_METHODS)
    f <- direct_fits(con, ds(), nav$method); req(nrow(f) > 0)
    f$fit <- vapply(f$fit_id, function(id) fit_descriptor(con, nav$method, id), character(1))
    datatable(f[, c("fit", "mse", "n_factors")], rownames = FALSE, options = list(pageLength = 15)) |>
      formatSignif("mse", 4)
  })
  observe({
    req(nav$level == 1, nav$method %in% DIRECT_FIT_METHODS)
    f <- direct_fits(con, ds(), nav$method)
    labels <- vapply(f$fit_id, function(id) fit_descriptor(con, nav$method, id), character(1))
    updateSelectInput(session, "l1_direct_fit", choices = setNames(f$fit_id, labels))
  })
  observeEvent(input$l1_direct_go, {
    req(input$l1_direct_fit)
    nav$fit <- as.integer(input$l1_direct_fit)
    nav$rank <- get_fit(con, nav$fit)$rank   # NA for Tucker (no single rank) -- harmless
    nav$level <- 2
  })

  ## =============================================================================
  ## LEVEL 2 -- rank / parameter view
  ## =============================================================================

  level2_ui <- function() {
    m <- nav$method
    if (m %in% DIRECT_FIT_METHODS) {
      # CP/Tucker: single fit already chosen at Level 1 (nav$fit) -- no
      # seed selector needed, unlike the pca/spca/nmf/cogaps/ica branch
      # below. "scores" are SUBJECT-mode, not sample-mode -- sample-
      # metadata association is deferred (would need subject-id
      # reconciliation, see R/README.md); factor correlation still works
      # generically regardless of what the rows represent. (sPCA used to
      # live in this branch too -- it moved to the rank-selector
      # FACTORIZATION_METHODS branch below once it got its own scree-style
      # rank view; see SCREE_RANK_METHODS.)
      panels <- list(
        nav_panel("Factor correlation (this fit)",
          p("How redundant are this fit's own components with each other -- computed on the SUBJECT-mode scores (not sample-mode; sample-metadata association isn't available yet for tensor methods)."),
          plotOutput("l2_direct_factor_corr", height = "420px")),
        nav_panel("Time profile",
          p("Each component's value across timepoint levels -- the third, time-mode loading matrix CP/Tucker produce that no other method has."),
          plotOutput("l2_time_profile", height = "420px"))
      )
      panels <- c(panels, list(
        nav_panel("Enrichment overview",
          p("Whatever functional enrichment the cluster ingest pipeline has already computed for every factor in this fit -- ORA and GSEA, positive direction (this method has no meaningful negative side). Jump straight to any one factor's full detail via Level 3's Enrichment tab."),
          DTOutput("l2_direct_enrich_overview_table"),
          layout_columns(col_widths = c(6, 6),
            radioButtons("l2_direct_enrich_overview_type", "Query shown:", c("ORA" = "ora", "GSEA" = "gsea"), inline = TRUE),
            actionButton("l2_direct_enrich_overview_view", "View selected factor at Level 3", class = "btn-outline-primary btn-sm",
                         style = "margin-top: 24px;")))
      ))
      do.call(navset_card_tab, panels)
    } else if (m %in% FACTORIZATION_METHODS) {
      f <- fits_at_rank(con, ds(), m, nav$rank)
      # sPCA has no real seed -- multiple fits at the same rank are
      # disambiguated by `para` (stored in fits.alpha, see
      # extract_result()'s spca branch) instead.
      seed_choices <- if (m == "spca") {
        setNames(f$fit_id, paste0("para ", f$alpha))
      } else {
        setNames(f$fit_id, paste0("seed ", f$seed))
      }
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
            selectInput("l2_scores_fit", "Fit:", choices = seed_choices),
            selectInput("l2_scores_sort", "Sort/group by metadata field:", choices = NULL)),
          DTOutput("l2_scores_table")),
        nav_panel("Factor correlation (this fit)",
          p("How redundant are this fit's own factors with each other, at the sample-score level?"),
          selectInput("l2_corr_fit", "Fit:", choices = seed_choices),
          plotOutput("l2_factor_corr", height = "420px")),
        nav_panel("Metadata associations",
          p("Spearman correlation (numeric fields) / Kruskal-Wallis (categorical fields) between each factor's sample scores and every sample-metadata column. p-values BH-adjusted across the whole grid shown."),
          selectInput("l2_assoc_fit", "Fit:", choices = seed_choices),
          plotOutput("l2_assoc_heatmap", height = "440px")),
        nav_panel("Enrichment overview",
          p("Whatever functional enrichment the cluster ingest pipeline has already computed for every factor in this fit -- ORA and GSEA, and (for PCA/ICA/sPCA) both loading directions. Jump straight to any one factor's full detail via Level 3's Enrichment tab."),
          selectInput("l2_enrich_all_fit", "Fit:", choices = seed_choices),
          DTOutput("l2_enrich_overview_table"),
          layout_columns(col_widths = c(4, 4, 4),
            radioButtons("l2_enrich_overview_type", "Query shown:", enrich_query_choices(m), inline = TRUE),
            radioButtons("l2_enrich_overview_dir", "Direction shown:", c("positive" = "pos", "negative" = "neg"), inline = TRUE),
            actionButton("l2_enrich_overview_view", "View selected factor at Level 3", class = "btn-outline-primary btn-sm",
                         style = "margin-top: 24px;")))
      ))
      do.call(navset_card_tab, panels)
    } else if (m == "wgcna") {
      f <- wgcna_fits(con, ds())
      power_choices <- setNames(f$fit_id, paste0("power ", f$power))
      navset_card_tab(
        nav_panel("Modules at this power",
          tableOutput("l2_wgcna_sizes"),
          layout_columns(col_widths = c(6, 6),
            selectInput("l2_wgcna_module_view", "View module:", choices = NULL),
            actionButton("l2_wgcna_module_go", "View module genes / enrichment", class = "btn-primary btn-sm",
                         style = "margin-top: 24px;")),
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
        nav_panel("Enrichment overview",
          p("ORA already computed by the cluster ingest pipeline (R/ingest_jobs/wgcna_ora_job.R) for every module in this power (WGCNA has no GSEA/direction side -- see the Level 3 module view). Jump straight to any one module's full detail."),
          DTOutput("l2_wgcna_enrich_overview_table"),
          actionButton("l2_wgcna_enrich_overview_view", "View selected module at Level 3", class = "btn-outline-primary btn-sm"))
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

  # ---- CP/Tucker Level 2 (direct fit, no seed selector) ----

  output$l2_direct_factor_corr <- renderPlot({
    req(nav$fit)
    L <- load_scores(con, nav$fit)
    validate(need(!is.null(L), "No sample/subject scores available for this fit."))
    cm <- cor(L)
    d <- as.data.frame(as.table(cm))
    names(d) <- c("factor_a", "factor_b", "cor")
    ggplot(d, aes(factor_a, factor_b, fill = cor)) +
      geom_tile() + geom_text(aes(label = sprintf("%.2f", cor)), size = 3) +
      scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick", limits = c(-1, 1)) +
      labs(x = NULL, y = NULL, title = "Within-fit factor correlation") +
      theme_minimal(base_size = 14)
  })
  output$l2_time_profile <- renderPlot({
    req(nav$fit)
    tl <- load_time_loadings(con, nav$fit)
    validate(need(!is.null(tl), "No time-mode loadings available for this fit."))
    d <- data.frame(
      time = factor(rep(rownames(tl), ncol(tl)), levels = rownames(tl)),
      component = factor(rep(colnames(tl), each = nrow(tl)), levels = colnames(tl))
    )
    d$value <- as.vector(tl)
    ggplot(d, aes(time, value, group = component, color = component)) +
      geom_line() + geom_point(size = 2) +
      labs(x = "timepoint", y = "loading", color = "component",
           title = paste(toupper(nav$method), "-- component trajectories across timepoints")) +
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
  observe({
    req(nav$level == 2, nav$method == "wgcna")
    sizes <- wgcna_module_sizes(con, nav$wgcna_fit)
    updateSelectInput(session, "l2_wgcna_module_view",
                       choices = setNames(sizes$module, paste0("module ", sizes$module, " (", sizes$n_genes, " genes)")))
  })
  observeEvent(input$l2_wgcna_module_go, {
    req(input$l2_wgcna_module_view, nav$wgcna_fit)
    nav$fit <- nav$wgcna_fit
    nav$factor_index <- as.integer(input$l2_wgcna_module_view)
    nav$level <- 3
  })

  # ---- Level 2: "Enrichment overview" -- browse whatever the cluster
  # ingest pipeline has already computed for every factor of a fit, and
  # jump to any one factor's full detail. No live compute anywhere in this
  # app -- see R/ingest_jobs/fgsea_job.R / R/ingest_jobs/wgcna_ora_job.R +
  # R/ingest_enrichment_results.R for where enrichment actually gets
  # computed now. ----

  #' Cache lookup for one factor's enrichment that also recognizes a
  #' numerically-identical loading vector already cached under a
  #' DIFFERENT fit_id. PCA gets a distinct fit_id per rank, but
  #' prcomp(rank. = k) just truncates the same underlying SVD (see
  #' R/methods/pca.R) -- so PC 3 at rank 7 and PC 3 at rank 10 are
  #' byte-identical. fgsea_grid's representative_fit_ids() only ever
  #' computes enrichment for PCA's single max-rank fit (R/lib/ingest/
  #' redundancy.R), so every OTHER PCA rank's factors would otherwise show
  #' nothing at all -- this is a pure read-time lookup (no DB write) that
  #' finds the equivalent factor's already-computed result instead.
  #' Restricted to PCA since that's the only method with this guarantee --
  #' NMF/CoGAPS/sPCA/CP/Tucker fits are never numerically identical to one
  #' another, so the search below would just waste time looking.
  find_or_reuse_enrichment <- function(fit_id, factor_index, qtype, direction) {
    factor_id <- get_factor_id(con, fit_id, factor_index)
    cached <- enrichment_cached(con, factor_id, qtype, direction)
    if (!is.null(cached)) return(cached)

    fit <- get_fit(con, fit_id)
    if (nrow(fit) == 0 || fit$method != "pca") return(NULL)

    v <- load_loadings(con, fit_id); if (is.null(v)) return(NULL)
    v <- v[, factor_index]
    siblings <- DBI::dbGetQuery(con,
      "SELECT fit_id FROM fits
       WHERE dataset_id = ? AND method = 'pca' AND family = 'seed_sweep'
         AND fit_id != ? AND status = 'ok'",
      params = list(fit$dataset_id, fit_id))$fit_id
    for (other_fit in siblings) {
      L2 <- load_loadings(con, other_fit)
      if (is.null(L2)) next
      for (fi2 in seq_len(ncol(L2))) {
        v2 <- L2[, fi2]
        if (length(v) != length(v2) || !identical(names(v), names(v2))) next
        if (!isTRUE(all.equal(as.numeric(v), as.numeric(v2), tolerance = 1e-8))) next
        other_factor_id <- get_factor_id(con, other_fit, fi2)
        hit <- enrichment_cached(con, other_factor_id, qtype, direction)
        if (!is.null(hit)) return(hit)
      }
    }
    NULL
  }

  #' Core work for "run every (factor, query type, direction) combination
  #' for a fit that isn't already cached (or reusable from an equivalent
  #' PCA factor at a different rank)" -- i.e. everything a user could
  #' otherwise trigger one factor/direction/query-type at a time from
  #' Level 3. `method` is taken as an explicit argument rather than read
  #' off `nav$method` so this is safely callable from contexts where
  #' nav$method isn't (yet) the target method -- e.g. the Level 0
  #' "best seed per rank" button below, which never navigates at all.
  #' Cheap (cache-only, no computation) per-factor summary for one query
  #' type/direction -- shared by both the seed-selector and direct-fit
  #' "Enrichment overview" panels. Purely a read of whatever
  #' R/ingest_jobs/fgsea_job.R + R/ingest_enrichment_results.R have
  #' already computed on the cluster -- all enrichment computation happens
  #' there now, never live in this app.
  enrich_overview_table <- function(fit_id, qtype, direction) {
    L <- load_loadings(con, fit_id); req(!is.null(L))
    rows <- lapply(seq_len(ncol(L)), function(fi) {
      factor_id <- get_factor_id(con, fit_id, fi)
      cc <- enrichment_cached(con, factor_id, qtype, direction)
      if (is.null(cc)) {
        data.frame(factor = fi, status = "not computed yet", n_terms = NA_integer_,
                   top_term = NA_character_, min_p_value = NA_real_)
      } else if (nrow(cc) == 0) {
        data.frame(factor = fi, status = "computed -- no significant terms", n_terms = 0L,
                   top_term = NA_character_, min_p_value = NA_real_)
      } else {
        data.frame(factor = fi, status = "ok", n_terms = nrow(cc),
                   top_term = cc$term_name[1], min_p_value = cc$p_value[1])
      }
    })
    do.call(rbind, rows)
  }

  # factorization methods (pca/nmf/cogaps/spca/ica) -- fit chosen via seed selector
  l2_enrich_overview_data <- reactive({
    req(input$l2_enrich_all_fit)
    enrich_overview_table(as.integer(input$l2_enrich_all_fit),
                           input$l2_enrich_overview_type %||% "ora",
                           input$l2_enrich_overview_dir %||% "pos")
  })
  output$l2_enrich_overview_table <- renderDT({
    d <- l2_enrich_overview_data()
    datatable(d, rownames = FALSE, selection = "single", options = list(pageLength = 15)) |>
      formatSignif("min_p_value", 3)
  })
  observeEvent(input$l2_enrich_overview_view, {
    d <- l2_enrich_overview_data()
    sel <- input$l2_enrich_overview_table_rows_selected
    req(length(sel) == 1)
    nav$fit <- as.integer(input$l2_enrich_all_fit)
    nav$factor_index <- d$factor[sel]
    nav$level <- 3
  })

  # sPCA/CP/Tucker -- single fit already chosen at Level 1 (nav$fit)
  l2_direct_enrich_overview_data <- reactive({
    req(nav$fit)
    enrich_overview_table(nav$fit, input$l2_direct_enrich_overview_type %||% "ora", "pos")
  })
  output$l2_direct_enrich_overview_table <- renderDT({
    d <- l2_direct_enrich_overview_data()
    datatable(d, rownames = FALSE, selection = "single", options = list(pageLength = 15)) |>
      formatSignif("min_p_value", 3)
  })
  observeEvent(input$l2_direct_enrich_overview_view, {
    d <- l2_direct_enrich_overview_data()
    sel <- input$l2_direct_enrich_overview_table_rows_selected
    req(length(sel) == 1)
    nav$factor_index <- d$factor[sel]
    nav$level <- 3
  })

  # WGCNA -- ORA only, whole module membership (no topN/direction, matching
  # the single-module enrichment at Level 3). Purely a read of whatever
  # R/ingest_jobs/wgcna_ora_job.R + R/ingest_enrichment_results.R have
  # already computed on the cluster -- no live compute in this app.
  l2_wgcna_enrich_overview_data <- reactive({
    req(nav$wgcna_fit)
    sizes <- wgcna_module_sizes(con, nav$wgcna_fit)
    rows <- lapply(sizes$module, function(mod) {
      factor_id <- get_factor_id(con, nav$wgcna_fit, mod)
      cc <- if (length(factor_id) == 1) enrichment_cached(con, factor_id, "ora", "pos") else NULL
      if (is.null(cc)) {
        data.frame(module = mod, status = "not computed yet", n_terms = NA_integer_,
                   top_term = NA_character_, min_p_value = NA_real_)
      } else if (nrow(cc) == 0) {
        data.frame(module = mod, status = "computed -- no significant terms", n_terms = 0L,
                   top_term = NA_character_, min_p_value = NA_real_)
      } else {
        data.frame(module = mod, status = "ok", n_terms = nrow(cc),
                   top_term = cc$term_name[1], min_p_value = cc$p_value[1])
      }
    })
    do.call(rbind, rows)
  })
  output$l2_wgcna_enrich_overview_table <- renderDT({
    d <- l2_wgcna_enrich_overview_data()
    datatable(d, rownames = FALSE, selection = "single", options = list(pageLength = 15)) |>
      formatSignif("min_p_value", 3)
  })
  observeEvent(input$l2_wgcna_enrich_overview_view, {
    d <- l2_wgcna_enrich_overview_data()
    sel <- input$l2_wgcna_enrich_overview_table_rows_selected
    req(length(sel) == 1)
    nav$fit <- nav$wgcna_fit
    nav$factor_index <- d$module[sel]
    nav$level <- 3
  })

  ## =============================================================================
  ## LEVEL 3 -- factor view
  ## =============================================================================

  level3_ui <- function() {
    if (nav$method == "wgcna") {
      return(navset_card_tab(
        nav_panel(
          "Module genes",
          p("Member genes of this module -- no continuous ranking available (WGCNA doesn't store intramodular connectivity/kME today), so this is the full membership list."),
          selectInput("l3_label_col", "Label genes by:", choices = "feature_id"),
          DTOutput("l3_wgcna_genes")
        ),
        nav_panel(
          "Enrichment (ORA)",
          p("GSEA doesn't apply here -- a WGCNA module is an unranked gene set, not a continuous loading to order by. Shows whatever the cluster ingest pipeline (R/ingest_jobs/wgcna_ora_job.R) has already computed."),
          DTOutput("l3_wgcna_enrichment")
        )
      ))
    }
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
        p("Shows whatever the cluster ingest pipeline (R/ingest_jobs/fgsea_job.R + R/ingest_enrichment_results.R) has already computed for this factor -- no live compute in this app."),
        layout_columns(col_widths = c(4, 4, 4),
          radioButtons("l3_enrich_type", "Query:", enrich_query_choices(nav$method), inline = TRUE),
          radioButtons("l3_enrich_dir", "Loadings:", c("positive" = "pos", "negative" = "neg"), inline = TRUE)),
        navset_card_tab(
          nav_panel("Table", DTOutput("l3_enrichment")),
          nav_panel("Mirror bar",
            p("Both the positive- and negative-loading direction need to already be computed for a full mirror; PCA/ICA/sPCA only -- NMF/CoGAPS weights are non-negative, so this shows the single positive side."),
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
    # Default to gene symbol (DISPLAY ONLY -- see build_display_map()'s
    # header) rather than feature_id, when a symbol column is
    # available/resolvable for this dataset.
    default_sel <- if (!is.null(m)) resolve_symbol_col(con, ds(), names(m)) %||% "feature_id" else "feature_id"
    updateSelectInput(session, "l3_label_col", choices = choices, selected = default_sel)
  })

  output$l3_loadings <- renderPlot({
    v <- l3_loadings()
    topn <- input$l3_topn %||% 25
    ord <- order(abs(v), decreasing = TRUE)[seq_len(min(topn, length(v)))]
    d <- data.frame(feature_id = names(v)[ord], loading = v[ord], stringsAsFactors = FALSE)

    # DISPLAY ONLY -- see build_display_map()'s header. Fallback chain:
    # chosen label_col (defaults to gene symbol) -> Ensembl gene id ->
    # feature_id itself, never left blank.
    label_col <- input$l3_label_col %||% "feature_id"
    if (label_col == "feature_id") {
      d$label <- d$feature_id
    } else {
      disp <- build_display_map(con, ds(), label_col)
      d$label <- if (!is.null(disp)) unname(disp[d$feature_id]) else d$feature_id
      d$label[is.na(d$label)] <- d$feature_id[is.na(d$label)]
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

  #' Pure read of whatever the cluster ingest pipeline has already
  #' computed for this factor -- no live compute anywhere in this app (see
  #' R/ingest_jobs/fgsea_job.R / R/ingest_enrichment_results.R). A plain
  #' reactive (not eventReactive/a button) so it just updates whenever the
  #' selected factor/query type/direction changes.
  enrich_result <- reactive({
    req(nav$fit, nav$factor_index)
    qtype <- input$l3_enrich_type
    direction <- if (nav$method %in% NEG_DIRECTION_METHODS) input$l3_enrich_dir else "pos"
    factor_id <- get_factor_id(con, nav$fit, nav$factor_index)

    if (qtype %in% c("fgsea", "cogaps_fora")) {
      # Hallmark GSEA / CoGAPS marker-gene ORA are only ever computed for
      # representative fits at ingest time -- no PCA-cross-rank reuse
      # concept applies to these two.
      return(enrichment_cached(con, factor_id, qtype, direction))
    }
    find_or_reuse_enrichment(nav$fit, nav$factor_index, qtype, direction)
  })
  output$l3_enrichment <- renderDT({
    d <- enrich_result()
    if (is.null(d) || nrow(d) == 0) {
      return(datatable(data.frame(
        note = "Not yet computed for this factor (all enrichment is computed by the cluster ingest pipeline -- see R/ingest_jobs/fgsea_job.R + R/ingest_enrichment_results.R), or no significant terms were found."
      ), rownames = FALSE))
    }
    datatable(d, rownames = FALSE, options = list(pageLength = 15)) |>
      formatSignif("p_value", 3)
  })

  # combines whatever's cached for BOTH loading directions (PCA only has a
  # meaningful negative side; NMF/CoGAPS weights are non-negative so this
  # degenerates to just "pos") -- independent of which single direction the
  # radio button currently has selected, so a mirror plot can show both
  # sides once each has been computed
  enrich_both <- reactive({
    req(nav$fit, nav$factor_index)
    factor_id <- get_factor_id(con, nav$fit, nav$factor_index)
    qtype <- input$l3_enrich_type
    dirs <- if (nav$method %in% NEG_DIRECTION_METHODS) c("pos", "neg") else "pos"
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

  # ---- WGCNA module view (Level 3 analog) ----

  l3_wgcna_genes_reactive <- reactive({
    req(nav$method == "wgcna", nav$fit, nav$factor_index)
    DBI::dbGetQuery(con, "SELECT gene FROM wgcna_modules WHERE fit_id = ? AND module = ?",
                    params = list(nav$fit, nav$factor_index))$gene
  })
  output$l3_wgcna_genes <- renderDT({
    genes <- l3_wgcna_genes_reactive()
    # DISPLAY ONLY -- see build_display_map()'s header. `gene` (the raw,
    # native feature_id) is kept as its own column alongside the
    # human-readable label, not replaced by it.
    label_col <- input$l3_label_col %||% "feature_id"
    label <- if (label_col == "feature_id") {
      genes
    } else {
      disp <- build_display_map(con, ds(), label_col)
      out <- if (!is.null(disp)) unname(disp[genes]) else genes
      out[is.na(out)] <- genes[is.na(out)]
      out
    }
    datatable(data.frame(gene = genes, label = label), rownames = FALSE, options = list(pageLength = 20))
  })
  #' Pure read of whatever R/ingest_jobs/wgcna_ora_job.R +
  #' R/ingest_enrichment_results.R have already computed for this module --
  #' no live compute in this app.
  wgcna_enrich_result <- reactive({
    req(nav$method == "wgcna", nav$fit, nav$factor_index)
    factor_id <- get_factor_id(con, nav$fit, nav$factor_index)
    validate(need(length(factor_id) == 1,
      "This WGCNA fit predates module-level enrichment support -- re-ingest this family with --overwrite (see R/ingest_results.R) to enable it."))
    enrichment_cached(con, factor_id, "ora", "pos")
  })
  output$l3_wgcna_enrichment <- renderDT({
    d <- wgcna_enrich_result()
    if (is.null(d) || nrow(d) == 0) {
      return(datatable(data.frame(note = "Not yet computed for this module (see R/ingest_jobs/wgcna_ora_job.R), or no significant terms were found."), rownames = FALSE))
    }
    datatable(d, rownames = FALSE, options = list(pageLength = 15)) |>
      formatSignif("p_value", 3)
  })
}

`%||%` <- function(a, b) if (is.null(a)) b else a

shinyApp(ui, server)

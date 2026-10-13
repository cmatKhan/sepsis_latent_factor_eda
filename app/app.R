# Latent-factor stability explorer.
#
# Reads the stability DB the targets pipeline writes (R/targets/db_write.R,
# docs/database.qmd), read-only. DB path: the STABILITY_DB environment
# variable, else `db_path` in config/pipeline.yml, relative to the project root.
#
# Layout: a sidebar holds the choices every page shares (dataset, methods,
# similarity metric, how the "optimal" fit is picked); the main panel is the
# Compare page, whose four tabs answer
#   1. what is the optimal hyperparameter for each method on this dataset
#   2. which method is best / most stable on this dataset
#   3. how similar are the optimal factors across methods
#   4. how similar are the factors across datasets (same method)
# Further pages: Explore fits (one method across its hyperparameter), Inspect fit
# (factor detail, enrichment, metadata, biplots, pattern drivers), and Across datasets (trajectory across time, projection comparison).

library(shiny)
library(bslib)
library(DBI)
library(RSQLite)
library(ggplot2)
library(DT)

source(file.path("R", "db_helpers.R"), local = TRUE)
source(file.path("R", "comparison_helpers.R"), local = TRUE)
source(file.path("R", "mod_compare.R"), local = TRUE)
source(file.path("R", "mod_explore.R"), local = TRUE)
source(file.path("R", "data_helpers.R"), local = TRUE)
source(file.path("R", "mod_inspect.R"), local = TRUE)
source(file.path("R", "mod_extra.R"), local = TRUE)
source(file.path("R", "mod_time.R"), local = TRUE)
# prepare_loadings() / hungarian_match_abs(): pure functions shared with ingestion
source(file.path("..", "R", "lib", "ingest", "similarity.R"), local = TRUE)

db_path <- Sys.getenv("STABILITY_DB", unset = "")
if (!nzchar(db_path)) {
  root <- if (file.exists("config/pipeline.yml")) "." else ".."
  rel <- yaml::read_yaml(file.path(root, "config", "pipeline.yml"))$db_path
  db_path <- file.path(root, rel)
}
if (!file.exists(db_path)) stop("No stability DB at ", db_path, " -- set STABILITY_DB or run the targets pipeline first")
con <- DBI::dbConnect(RSQLite::SQLite(), db_path, flags = RSQLite::SQLITE_RO)
options(stability.db_dir = normalizePath(dirname(db_path)))
onStop(function() DBI::dbDisconnect(con))

DATASETS <- list_datasets(con)


ui <- page_navbar(
  title = "Latent-factor stability",
  theme = bs_theme(bootswatch = "flatly"),
  sidebar = sidebar(
    width = 300,
    selectInput("dataset", "Dataset", choices = DATASETS),
    selectInput("fit_method", "Method", choices = NULL),
    selectInput("fit_rank", "Rank", choices = NULL),
    radioButtons("metric", "Similarity metric", choices = METRICS, inline = TRUE),
    accordion(
      open = FALSE,
      accordion_panel(
        "Compare settings",
        helpText("Used only by the Compare page, which sets several methods side by side."),
        checkboxGroupInput("methods", "Methods to compare", choices = NULL),
        radioButtons("rule", "Optimal fit rule", choices = c(
          "Highest quality that is still stable" = "stable",
          "Highest in-sample quality" = "max_r2")),
        conditionalPanel("input.rule == 'stable'",
          sliderInput("threshold", "Minimum stability", min = 0.5, max = 1, value = 0.9, step = 0.01)),
        helpText("Quality rises with rank, so the max-quality rule picks the largest rank.",
                 "The stable rule picks the largest rank whose factors still reproduce across seeds.")))),
  nav_panel("Time structure", mod_time_ui("time")),
  nav_panel("Compare",
    uiOutput("coverage"),
    navset_card_tab(
      id = "compare_tab",
      nav_panel("1. Optimal hyperparameter", mod_hyperparam_ui("hyper")),
      nav_panel("2. Best method", mod_best_method_ui("best")),
      nav_panel("3. Across methods", mod_cross_method_ui("xmethod")),
      nav_panel("4. Across datasets", mod_cross_dataset_ui("xdata")))),
  nav_panel("Explore fits", mod_explore_ui("explore")),
  nav_menu("Inspect fit",
    nav_panel("Factor detail", mod_factor_ui("factor")),
    nav_panel("Enrichment", mod_enrichment_ui("enrich")),
    nav_panel("Metadata", mod_metadata_ui("meta")),
    nav_panel("Biplots", mod_biplot_ui("biplot")),
    nav_panel("Pattern drivers", mod_drivers_ui("drivers"))),
  nav_menu("Across datasets",
    nav_panel("Trajectory across time", mod_trajectory_ui("traj")),
    nav_panel("Projection comparison", mod_projection_ui("proj"))))

server <- function(input, output, session) {
  present <- reactive(methods_for_dataset(con, input$dataset))

  observeEvent(present(), {
    p <- present()
    updateCheckboxGroupInput(session, "methods", choices = setNames(p, method_label(p)), selected = p, inline = TRUE)
  })

  curves <- reactive({
    req(input$dataset)
    dataset_curves(con, input$dataset, present(), input$metric)
  }) |> bindCache(input$dataset, input$metric)

  optimal <- reactive({
    cv <- curves()
    cv <- cv[cv$method %in% input$methods, , drop = FALSE]
    optimal_fits(cv, input$rule, input$threshold)
  })

  state <- list(dataset = reactive(input$dataset), metric = reactive(input$metric),
                methods = reactive(input$methods), curves = curves, optimal = optimal)

  # The fit every page looks at: the sidebar's method at a rank (default: that
  # method's optimal fit). Pages that need a second fit or dataset have their
  # own selector for just that.
  observe({
    p <- present(); req(length(p) > 0)
    cur <- isolate(input$fit_method)
    updateSelectInput(session, "fit_method", choices = setNames(p, method_label(p)),
                      selected = if (!is.null(cur) && cur %in% p) cur else p[1])
  })
  observeEvent(list(input$fit_method, curves(), optimal()), {
    req(input$fit_method)
    cv <- curves(); cv <- cv[cv$method == input$fit_method, , drop = FALSE]
    opt <- optimal(); o <- opt$rank_key[opt$method == input$fit_method]
    updateSelectInput(session, "fit_rank", label = if (input$fit_method == "wgcna") "Power" else "Rank",
                      choices = setNames(cv$rank_key, cv$rank_label),
                      selected = if (length(o) == 1) o else cv$rank_key[1])
  })
  focus <- reactive({
    req(input$fit_method, input$fit_rank)
    cv <- curves()
    r <- cv[cv$method == input$fit_method & cv$rank_key == as.numeric(input$fit_rank), , drop = FALSE]
    req(nrow(r) == 1)
    list(method = r$method, fit_id = r$fit_id, n_factors = r$n_factors, rank_label = r$rank_label, rank_key = r$rank_key)
  })

  output$coverage <- renderUI({
    cov <- fit_coverage(con, input$dataset)
    miss <- setdiff(METHOD_ORDER, cov$method)
    if (length(miss) == 0) return(NULL)
    p(class = "text-muted small", "No fits for this dataset: ", paste(method_label(miss), collapse = ", "))
  })

  mod_hyperparam_server("hyper", state)
  mod_best_method_server("best", state, con)
  mod_cross_method_server("xmethod", state, con)
  mod_cross_dataset_server("xdata", state, focus, con)
  mod_time_server("time", state, focus, con)
  mod_explore_server("explore", state, focus, con)
  mod_factor_server("factor", state, focus, con)
  mod_enrichment_server("enrich", state, focus, con)
  mod_metadata_server("meta", state, focus, con)
  mod_biplot_server("biplot", state, focus, con)
  mod_drivers_server("drivers", state, focus, con)
  mod_trajectory_server("traj", state, focus, con)
  mod_projection_server("proj", state, focus, con)
}

shinyApp(ui, server)

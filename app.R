library(shiny)
library(shinythemes)
library(fontawesome)
library(DT)
library(plotly)
library(ggplot2)
library(limma)
library(ggrepel)
library(shinycssloaders)
library(fgsea)
library(msigdbr)
library(dplyr)
library(tidyr)

source("helper_functions.R")

# Load CHARM object once when app starts
Charm.object <- readRDS("data/Charm.object.RDS")
sh_effect_vector <- readRDS("data/shRNA_Efficiency.Rds")
Charm.object_K562 <- readRDS("data/Charm.object_K562.RDS")
sh_effect_vector_K562 <- readRDS("data/shRNA_Efficiency_K562.Rds")
Charm.object_HEPG2 <- readRDS("data/Charm.object_HEPG2.RDS")
sh_effect_vector_HEPG2 <- readRDS("data/shRNA_Efficiency_HEPG2.Rds")

#Similarity stuff
similar_expression_all <- readRDS("data/RBPs.t_All.RDS")
similar_expression_K562 <- readRDS("data/RBPs.t_K562.RDS")
similar_expression_HEPG2 <- readRDS("data/RBPs.t_HEPG2.RDS")
similar_gsea_all <- readRDS("data/RBPs.gsea_All.RDS")
similar_gsea_K562 <- readRDS("data/RBPs.gsea_K562.RDS")
similar_gsea_HEPG2 <- readRDS("data/RBPs.gsea_HEPG2.RDS")

#SearchBarPopulation
GenesBoth <- readRDS("data/AvailableGenes_both.RDS")
GenesK562 <- readRDS("data/AvailableGenes_K562.RDS")
GenesHEPG2 <- readRDS("data/AvailableGenes_HEPG2.RDS")
EventsBoth <- readRDS("data/AvailableEvents_Both.RDS")
EventsK562 <- readRDS("data/AvailableEvents_K562.RDS")
EventsHEPG2 <- readRDS("data/AvailableEvents_HEPG2.RDS")

ui <- fluidPage(
  theme = shinytheme("flatly"),

  # Custom CSS for bold tabs
  tags$head(
    tags$style(HTML("
      .search-box {
        border-radius: 20px;
        padding: 6px 12px;
        border: 1px solid #ccc;
        box-shadow: inset 0 1px 3px rgba(0,0,0,0.1);
        font-size: 14px;
      }
      .search-box:focus {
        border-color: #2c3e50;
        outline: none;
        box-shadow: 0 0 5px rgba(44,62,80,0.4);
      }
    ")),
    tags$script(HTML("
      Shiny.addCustomMessageHandler('toggleCursor', function(state) {
        document.body.style.cursor = state ? 'wait' : 'default';
      });
    "))
  ),

  # Tabset with 5 tabs
  tabsetPanel(
    type = "tabs",
    id = "main_tabs",

    # Home tab
    tabPanel(
      tagList(fa("home", fill = "black", height = "1em"), " Home"),
      fluidPage(
        fluidRow(
          column(
            4,
            div(
              style = "display: flex; align-items: center; justify-content: center; margin: 30px 0;",
              div(style = "margin-right: 15px;", fa("gem", fill = "black", height = "4em")),
              div(
                style = "text-align: left;",
                tags$h2("CHARM", style = "margin: 0; font-weight: bold;"),
                tags$h4("Comprehensive Hub of Alternative Regulatory Mapping", style = "margin: 0; font-weight: bold;")
              )
            )
          )
        ),
        fluidRow(
          column(
            6,
            div(
              style = "background-color: #EAEBEB; color: black; border: 2px solid #2c3e50; border-radius: 20px; padding: 40px; text-align: center; margin: 15px;",
              div(style = "margin-bottom: 10px;", fa("globe", fill = "black", height = "3em")),
              tags$h3("Explore", style = "font-weight: bold;"),
              tags$p("Investigate how a known RBP affects expression, splicing, and binding, based on ENCODE's gene silencing series and eCLIP data. Find out which genes/gene sets are more influenced by RBPs")
            )
          ),
          column(
            6,
            div(
              style = "background-color: #EAEBEB; color: black; border: 2px solid #2c3e50; border-radius: 20px; padding: 40px; text-align: center; margin: 15px;",
              div(style = "margin-bottom: 10px;", fa("map", fill = "black", height = "3em")),
              tags$h3("Discover", style = "font-weight: bold;"),
              tags$p("Users can input their own expression, splicing, or binding data to discover which RBPs are more likely to be altered on your biological system.")
            )
          )
        )
      )
    ),

    # Expression tab
    tabPanel(
      tagList(fa("dna", fill = "black", height = "1em"), " Expression"),
      fluidPage(
        fluidRow(
          # Left sidebar
          # Left sidebar (width = 3)
          column(
            width = 3,
            wellPanel(

              # --- Mode selection ---
              radioButtons(
                inputId = "expr_mode",
                label = "Select mode:",
                choices = c("Explore Mode", "Discovery Mode"),
                selected = "Explore Mode"
              ),

              # --- Explore Mode ---
              conditionalPanel(
                condition = "input.expr_mode == 'Explore Mode'",

                # Dataset selection
                selectInput(
                  inputId = "expr_dataset",
                  label = "Select option:",
                  choices = c("Both Cells", "K562", "HEPG2", "Similar RBPs")
                ),

                # ---- Search by RBP ----
                div(
                  style = "margin-top: 15px;",
                  tags$h5("Search by RBP"),
                  tags$p("Find out how RBPs influence expression and gene set enrichment",
                         style = "font-size: 12px; color: #666; margin-top: -5px;"),
                  div(
                    style = "display: flex; align-items: center;",
                    selectizeInput(
                      inputId = "expr_search_rbp",
                      label = NULL,
                      choices = NULL,   # dynamically populated in server
                      multiple = FALSE,
                      options = list(placeholder = "RBP"),
                      width = "400px"
                    ),
                    conditionalPanel(
                      condition = "input.expr_dataset != 'Similar RBPs'",
                      div(
                        style = "display: flex; align-items: center; margin-left: 10px;",
                        actionButton("search_btn_rbp", tagList(fa("search"), " Search"),
                                     class = "btn btn-primary", style = "margin-right: 10px; border-radius: 20px;"),
                        actionButton("reset_btn_rbp", tagList(fa("redo"), " Reset"),
                                     class = "btn btn-secondary", style = "border-radius: 20px;")
                      )
                    )
                  )
                ),

                # ---- Search by Gene (hidden for Similar RBPs) ----
                conditionalPanel(
                  condition = "input.expr_dataset != 'Similar RBPs'",
                  div(
                    style = "margin-top: 25px;",
                    tags$h5("Search by Gene"),
                    selectizeInput(
                      inputId = "expr_search_gene",
                      label = NULL,
                      choices = NULL,   # dynamically populated in server
                      multiple = FALSE,
                      options = list(placeholder = "Gene"),
                      width = "400px"
                    ),
                    div(
                      style = "margin-top: 5px;",
                      actionButton("search_btn_gene", tagList(fa("search"), " Search"),
                                   class = "btn btn-primary", style = "margin-right: 10px; border-radius: 20px;"),
                      actionButton("reset_btn_gene", tagList(fa("redo"), " Reset"),
                                   class = "btn btn-secondary", style = "border-radius: 20px;")
                    )
                  )
                ),

                # ---- Search by Hallmark Gene Set (hidden for Similar RBPs) ----
                conditionalPanel(
                  condition = "input.expr_dataset != 'Similar RBPs'",
                  div(
                    style = "margin-top: 25px;",
                    tags$h5("Search by Hallmark Gene Set"),
                    selectizeInput(
                      inputId = "expr_search_hallmark",
                      label = NULL,
                      choices = NULL,   # dynamically populated in server
                      multiple = FALSE,
                      options = list(placeholder = "Hallmark Gene Set"),
                      width = "400px"
                    ),
                    div(
                      style = "margin-top: 5px;",
                      actionButton("search_btn_hallmark", tagList(fa("search"), " Search"),
                                   class = "btn btn-primary", style = "margin-right: 10px; border-radius: 20px;"),
                      actionButton("reset_btn_hallmark", tagList(fa("redo"), " Reset"),
                                   class = "btn btn-secondary", style = "border-radius: 20px;")
                    )
                  )
                ),

                # ---- Similar RBPs options (single, no duplication) ----

                conditionalPanel(
                  condition = "input.expr_dataset == 'Similar RBPs'",

                  # Correlation type
                  radioButtons(
                    inputId = "similar_mode",
                    label = "Select correlation type:",
                    choices = c("By Gene Expression" = "expr", "By Gene Set Enrichment" = "gsea"),
                    selected = "expr",
                    inline = TRUE
                  ),

                  # Expression correlation inputs - ONLY show when expr mode is selected
                  conditionalPanel(
                    condition = "input.similar_mode == 'expr'",
                    selectizeInput(
                      inputId = "similar_rbps_select_expr",
                      label = "Compare with specific RBPs (optional)",
                      choices = NULL,
                      multiple = TRUE,
                      options = list(placeholder = "Select one or more RBPs")
                    ),
                    numericInput("correl_num_expr", "Show top N correlated RBPs (optional):", value = NA, min = 1),
                    numericInput("n_pos_expr", "Show top N positive correlations (optional):", value = NA, min = 1),
                    numericInput("n_neg_expr", "Show top N negative correlations (optional):", value = NA, min = 1)
                  ),

                  # GSEA correlation inputs - ONLY show when gsea mode is selected
                  conditionalPanel(
                    condition = "input.similar_mode == 'gsea'",
                    selectizeInput(
                      inputId = "similar_rbps_select_gsea",
                      label = "Compare with specific RBPs (optional)",
                      choices = NULL,
                      multiple = TRUE,
                      options = list(placeholder = "Select one or more RBPs")
                    ),
                    numericInput("correl_num_gsea", "Show top N correlated RBPs (optional):", value = NA, min = 1),
                    numericInput("n_pos_gsea", "Show top N positive correlations (optional):", value = NA, min = 1),
                    numericInput("n_neg_gsea", "Show top N negative correlations (optional):", value = NA, min = 1)
                  ),

                  # Plot / Reset buttons
                  div(
                    style = "display: flex; align-items: center; margin-top: 15px;",
                    actionButton("similar_plot_btn", tagList(fa("chart-line"), " Plot"),
                                 class = "btn btn-primary", style = "margin-right: 10px; border-radius: 20px;"),
                    actionButton("similar_reset_btn", tagList(fa("redo"), " Reset"),
                                 class = "btn btn-secondary", style = "border-radius: 20px;")
                  )
                )
              ),

              # --- Discovery Mode: User File Upload ---
              conditionalPanel(
                condition = "input.expr_mode == 'Discovery Mode'",
                hr(),
                tags$p("Alternatively, upload your own table with differential expression values. Data must have no header, and the first column must be the HGNC gene symbol, and the second column the t-stats."),
                fileInput("user_file_expr", "Upload your file:", accept = c(".txt")),
                uiOutput("file_warning_expr"),
                uiOutput("user_file_options")
              )
            )
          )
          ,
          # Right content area
          column(
            width = 9,
            tags$h3("Expression Data Results"),

            # --- Explore Mode (default mode) ---
            conditionalPanel(
              condition = "input.expr_mode == 'Explore Mode' && input.expr_dataset != 'Similar RBPs'",
              # --- Warning message ---
              div(
                "⚠ Please press reset after every plot!",
                style = "border: 2px solid #f0ad4e;
               background-color: #fff3cd;
               padding: 8px;
               border-radius: 6px;
               font-weight: bold;
               color: #856404;"
              ),

              # --- Gene search plot (appear at top when searched) ---
              conditionalPanel(
                condition = "input.search_btn_gene > 0",
                fluidRow(
                  column(
                    width = 12,
                    plotOutput("expr_gene_plot", height = "600px")
                  )
                )
              ),

              # --- Hallmark gene set search plot (appear at top when searched) ---
              conditionalPanel(
                condition = "input.search_btn_hallmark > 0",
                fluidRow(
                  column(
                    width = 12,
                    plotOutput("expr_hallmark_plot", height = "600px")
                  )
                )
              ),

              # --- Default RBP plots ---
              fluidRow(
                column(width = 6, plotOutput("expr_violin", height = "400px")),
                column(
                  width = 6,
                  plotOutput("shrna_plot"),
                  uiOutput("shrna_warning")
                )
              ),
              fluidRow(
                column(width = 6, plotlyOutput("volcano_plot", height = "450px")),
                column(width = 6, DTOutput("volcano_table"))
              ),
              fluidRow(
                column(width = 6, plotOutput("gsea_plot", height = "600px")),
                column(width = 6, DTOutput("geneset_table"))
              )
            ),

            # --- Similar RBPs (Explore Mode, special layout) ---
            conditionalPanel(
              condition = "input.expr_mode == 'Explore Mode' && input.expr_dataset == 'Similar RBPs'",
              div(
                "⚠ Please press reset after every plot!",
                style = "border: 2px solid #f0ad4e;
               background-color: #fff3cd;
               padding: 8px;
               border-radius: 6px;
               font-weight: bold;
               color: #856404;
               margin-bottom: 10px;"
              ),
              uiOutput("similar_expr_plots")
            ),

            # --- Discovery Mode ---
            conditionalPanel(
              condition = "input.expr_mode == 'Discovery Mode'",
              div(
                "⚠ Please press reset after every plot!",
                style = "border: 2px solid #f0ad4e;
               background-color: #fff3cd;
               padding: 8px;
               border-radius: 6px;
               font-weight: bold;
               color: #856404;
               margin-bottom: 10px;"
              ),
              uiOutput("userfilesimilar_expr"),
              uiOutput("userfilesimilar_gsea")
            )
          )
        )
      )
    ),

    # Splicing tab
    tabPanel(
      tagList(fa("scissors", fill = "black", height = "1em"), " Splicing"),
      fluidPage(
        fluidRow(
          column(
            width = 3,
            wellPanel(
              radioButtons("splice_mode", "Select mode:", choices = c("Explore Mode", "Discovery Mode"), selected = "Explore Mode"),

              conditionalPanel(
                condition = "input.splice_mode == 'Explore Mode'",
                selectInput(
                  "splice_dataset", "Select option:",
                  choices = c("Both Cells", "K562", "HEPG2", "Similar RBPs")
                ),

                # ---- Search by RBP (always visible) ----
                div(
                  style = "margin-top: 15px;",
                  tags$h5("Search by RBP"),
                  tags$p(
                    "Visualize how splicing of an RBP impacts alternative splicing across datasets.",
                    style = "font-size: 12px; color: #666; margin-top: -5px;"
                  ),
                  div(
                    style = "display: flex; align-items: center;",
                    selectizeInput(
                      inputId = "splice_search",
                      label = NULL,
                      choices = NULL,
                      multiple = FALSE,
                      options = list(placeholder = "RBP"),
                      width = "400px"
                    ),
                    conditionalPanel(
                      condition = "input.splice_dataset != 'Similar RBPs'",
                      div(
                        style = "display: flex; align-items: center; margin-left: 10px;",
                        actionButton(
                          "splice_search_btn",
                          tagList(fa("search"), " Search"),
                          class = "btn btn-primary",
                          style = "margin-right: 10px; border-radius: 20px;"
                        ),
                        actionButton(
                          "splice_reset_btn",
                          tagList(fa("redo"), " Reset"),
                          class = "btn btn-secondary",
                          style = "border-radius: 20px;"
                        )
                      )
                    )
                  )
                ),

                # ---- Search by Event ID (hidden when Similar RBPs is selected) ----
                conditionalPanel(
                  condition = "input.splice_dataset != 'Similar RBPs'",
                  div(
                    style = "margin-top: 25px;",
                    tags$h5("Search by Event ID"),
                    tags$p(
                      list(
                        "Find out which RBPs most strongly regulate a specific splicing event. Event IDs are obtained from ",
                        tags$a(
                          href = "https://vastdb.crg.eu/",
                          "VastDB",
                          target = "_blank",
                          style = "color: #007bff; text-decoration: none;"
                        ),
                        "."
                      ),
                      style = "font-size: 12px; color: #666; margin-top: -5px;"
                    ),
                    div(
                      style = "display: flex; align-items: center;",
                      selectizeInput(
                        inputId = "splice_search_event",   # <-- consistent ID
                        label = NULL,
                        choices = NULL,
                        multiple = FALSE,
                        options = list(placeholder = "Event ID"),
                        width = "400px"
                      ),
                      div(
                        style = "display: flex; align-items: center; margin-left: 10px;",
                        actionButton(
                          "search_btn_event",
                          tagList(fa("search"), " Search"),
                          class = "btn btn-primary",
                          style = "margin-right: 10px; border-radius: 20px;"
                        )
                      )
                    )
                  )
                ),
                conditionalPanel(
                  condition = "input.splice_dataset == 'Similar RBPs'",
                  selectizeInput("similar_rbps_select_splice", "Compare with specific RBPs (optional)", choices = NULL, multiple = TRUE, options = list(placeholder = "Select one or more RBPs")),
                  numericInput("correl_num_splice", "Show top N correlated RBPs (optional):", value = NA, min = 1),
                  numericInput("n_pos_splice", "Show top N positive correlations (optional):", value = NA, min = 1),
                  numericInput("n_neg_splice", "Show top N negative correlations (optional):", value = NA, min = 1),
                  div(
                    style = "display: flex; align-items: center; margin-top: 15px;",
                    actionButton("similar_plot_btn_splice", tagList(fa("chart-line"), " Plot"), class = "btn btn-primary", style = "margin-right: 10px; border-radius: 20px;"),
                    actionButton("similar_reset_btn_splice", tagList(fa("redo"), " Reset"), class = "btn btn-secondary", style = "border-radius: 20px;")
                  )
                )
              ),
              conditionalPanel(
                condition = "input.splice_mode == 'Discovery Mode'",
                hr(),
                tags$p(
                  list(
                    "Alternatively, upload your own table with differential splicing values. Table must be directly obtained from ",
                    tags$a(
                      href = "https://compbio.imm.medicina.ulisboa.pt/app/betAS",
                      "betAS",
                      target = "_blank",                # opens in a new tab
                      style = "color: #007bff; text-decoration: none;" # optional styling
                    ),
                    "."
                  )
                ),
                fileInput("user_file_splice", "Upload your file:", accept = c(".txt")),
                uiOutput("file_warning_splice"),
                uiOutput("user_file_options_splice")
              )
            )
          ),
          #Right Side, where plots will appear
          # Right Side, where plots will appear
          column(
            width = 9,
            tags$h3("Splicing Data Results"),

            # --- Explore Mode (default mode, NOT Similar RBPs) ---
            conditionalPanel(
              condition = "input.splice_mode == 'Explore Mode' && input.splice_dataset != 'Similar RBPs'",
              # --- Warning message ---
              div(
                "⚠ Please press reset after every plot!",
                style = "border: 2px solid #f0ad4e;
               background-color: #fff3cd;
               padding: 8px;
               border-radius: 6px;
               font-weight: bold;
               color: #856404;
               margin-bottom: 10px;"
              ),

              # --- Default Explore-Mode Plots ---
              fluidRow(
                column(width = 6, plotOutput("violin_splice_plot", height = "400px")),
                column(width = 6, plotOutput("plot_shrna_effect", height = "400px"))
              ),
              fluidRow(
                column(width = 6, plotlyOutput("plot_splice_volcano", height = "450px")),
                column(width = 6, DTOutput("splice_volcano_table"))
              ),
              fluidRow(
                column(width = 12, plotOutput("heatmap_splicing_dpsi", height = "600px"))
              )
            ),

            # --- Explore Mode: Similar RBPs ---
            conditionalPanel(
              condition = "input.splice_mode == 'Explore Mode' && input.splice_dataset == 'Similar RBPs'",
              div(
                "⚠ Please press reset after every plot!",
                style = "border: 2px solid #f0ad4e;
               background-color: #fff3cd;
               padding: 8px;
               border-radius: 6px;
               font-weight: bold;
               color: #856404;
               margin-bottom: 10px;"
              ),
              uiOutput("similar_splice_plots")
            ),

            # --- Discovery Mode ---
            conditionalPanel(
              condition = "input.splice_mode == 'Discovery Mode'",
              div(
                "⚠ Please press reset after every plot!",
                style = "border: 2px solid #f0ad4e;
               background-color: #fff3cd;
               padding: 8px;
               border-radius: 6px;
               font-weight: bold;
               color: #856404;
               margin-bottom: 10px;"
              ),
              shinycssloaders::withSpinner(
                plotOutput("user_file_initial_plot_splice", height = "420px"),
                type = 6
              ),
              br(),
              uiOutput("similar_splice_plots_file")
            )
          )
        )
      )
    ),

    # Binding tab
    tabPanel(
      tagList(fa("link", fill = "black", height = "1em"), " Binding"),
      fluidPage(
        fluidRow(
          column(
            width = 3,
            wellPanel(
              radioButtons("binding_mode", "Select mode:", choices = c("Explore Mode", "Discovery Mode"), selected = "Explore Mode"),
              conditionalPanel(
                condition = "input.binding_mode == 'Explore Mode'",
                selectInput("binding_dataset", "Select option:", choices = c("Both Cells", "K562", "HEPG2", "Similar RBPs")),
                div(
                  style = "display: flex; align-items: center; margin-top: 15px;",
                  selectizeInput("binding_search", NULL, choices = NULL, multiple = FALSE, options = list(placeholder = "RBP"), width = "400px"),
                  actionButton("search_btn", tagList(fa("search"), " Search"), class = "btn btn-primary", style = "margin-left: 10px; border-radius: 20px;"),
                  actionButton("reset_btn", tagList(fa("redo"), " Reset"), class = "btn btn-secondary", style = "margin-left: 10px; border-radius: 20px;")
                )
              ),
              conditionalPanel(
                condition = "input.binding_mode == 'Discovery Mode'",
                hr(),
                tags$p("Alternatively, upload your own table with differential splicing values."),
                fileInput("user_file_binding", "Upload your file:", accept = c(".txt")),
                uiOutput("file_warning_binding")
              )
            )
          ),
          column(width = 9, tags$h3("Binding Data Results"), tags$p("Filtered results will appear here."))
        )
      )
    ),

    # Network tab
    tabPanel(
      tagList(fa("project-diagram", fill = "black", height = "1em"), " Network"),
      fluidPage(
        fluidRow(
          column(
            width = 3,
            wellPanel(
              radioButtons("network_mode", "Select mode:", choices = c("Explore Mode", "Discovery Mode"), selected = "Explore Mode"),
              conditionalPanel(
                condition = "input.network_mode == 'Explore Mode'",
                selectInput("network_dataset", "Select option:", choices = c("Both Cells", "K562", "HEPG2", "Similar RBPs")),
                div(
                  style = "display: flex; align-items: center; margin-top: 15px;",
                  selectizeInput("network_search", NULL, choices = NULL, multiple = FALSE, options = list(placeholder = "RBP"), width = "400px"),
                  actionButton("search_btn", tagList(fa("search"), " Search"), class = "btn btn-primary", style = "margin-left: 10px; border-radius: 20px;"),
                  actionButton("reset_btn", tagList(fa("redo"), " Reset"), class = "btn btn-secondary", style = "margin-left: 10px; border-radius: 20px;")
                )
              )
            )
          ),
          column(width = 9, tags$h3("Network Results"), tags$p("Network-related results will appear here."))
        )
      )
    )
  ),

  # GitHub link floating across ALL tabs
  tags$div(
    style = "
      position: fixed;
      bottom: 20px;
      right: 20px;
      z-index: 1000;
      background-color: #EAEBEB;
      padding: 8px 12px;
      border-radius: 10px;
      border: 1px solid #2c3e50;
      text-align: center;
    ",
    tags$a(
      href = "https://github.com/DiseaseTranscriptomicsLab/CHARM",
      target = "_blank",
      style = "color: black; text-decoration: none;",
      fa("github", fill = "black", height = "1.5em"),
      " GitHub"
    )
  )
)

# ---- Server ----

server <- function(input, output, session) {
  result_splice <- reactiveVal(NULL)

  # ---- Helper functions to pick correct dataset ----
  # Expression tab
  current_charm_expr <- reactive({
    req(input$expr_dataset)
    switch(input$expr_dataset,
           "K562"  = Charm.object_K562,
           "HEPG2" = Charm.object_HEPG2,
           Charm.object)
  })

  # Splicing tab
  current_charm_splice <- reactive({
    req(input$splice_dataset)
    switch(input$splice_dataset,
           "K562"  = Charm.object_K562,
           "HEPG2" = Charm.object_HEPG2,
           Charm.object)
  })

  current_shrna_expr <- reactive({
    req(input$splice_dataset)
    switch(input$splice_dataset,
           "K562"  = sh_effect_vector_K562,
           "HEPG2" = sh_effect_vector_HEPG2,
           sh_effect_vector)
  })

  current_shrna_splice <- reactive({
    req(input$splice_dataset)
    switch(input$splice_dataset,
           "K562"  = sh_effect_vector_K562,
           "HEPG2" = sh_effect_vector_HEPG2,
           sh_effect_vector)
  })


  # Reactive storage
  display_table <- reactiveVal(NULL)
  rbp_current   <- reactiveVal(NULL)

  # ---- Search bar population (Expression tab) ----
  observe({
    req(current_charm_expr(), input$expr_dataset)

    # 1️⃣ RBP search bar (depends on selected dataset)
    rbp_choices <- names(current_charm_expr())
    updateSelectizeInput(
      session,
      "expr_search_rbp",
      choices = rbp_choices,
      server = TRUE
    )

    # 2️⃣ Gene search bar (depends on selected dataset)
    gene_choices <- switch(
      input$expr_dataset,
      "K562"  = GenesK562,
      "HEPG2" = GenesHEPG2,
      GenesBoth # default for "Both Cells" or anything else
    )
    updateSelectizeInput(
      session,
      "expr_search_gene",
      choices = sort(unique(gene_choices)),
      server = TRUE
    )

    # 3️⃣ Hallmark Gene Set search bar (always from RBFOX2)
    hallmark_choices <- Charm.object[["RBFOX2"]]$GSEA$pathway
    updateSelectizeInput(
      session,
      "expr_search_hallmark",
      choices = sort(unique(hallmark_choices)),
      server = TRUE
    )
  })


  # ---- Search bar population (Splicing tab) ----
  observe({
    req(current_charm_splice())

    rbp_choices_splice <- names(current_charm_splice())
    updateSelectizeInput(
      session,
      "splice_search",
      choices = rbp_choices_splice,
      server = TRUE
    )

    rbp_choices_similar <- names(current_charm_splice())
    updateSelectizeInput(session, "similar_rbps_select_splice", choices = rbp_choices_similar)
  })

  # ---- Populate Event ID choices based on dataset ----
  observe({
    req(current_charm_splice())
    req(input$splice_dataset)

    event_choices <- switch(
      input$splice_dataset,
      "Both Cells" = EventsBoth,
      "K562"       = EventsK562,
      "HEPG2"      = EventsHEPG2,
      NULL
    )

    updateSelectizeInput(
      session,
      "splice_search_event",  # must match UI
      choices = event_choices,
      server = TRUE
    )
  })

  # ---- Event ID Search button ----
  observeEvent(input$search_btn_event, {
    req(input$splice_search_event)
    event_id <- input$splice_search_event

    charm_obj <- current_charm_splice()

    output$heatmap_splicing_dpsi <- renderPlot({
      plot_event_dpsi_barplot(charm_obj, event_id)
    })
  })

  ###EXPRESSION (user File)
  upload_ok <- reactiveVal(FALSE)
  user_expr_df <- reactiveVal(NULL)

  observeEvent(input$user_file_expr, {
    req(input$user_file_expr)
    file_path <- input$user_file_expr$datapath

    df <- tryCatch(
      read.table(file_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE),
      error = function(e) NULL
    )

    if (is.null(df)) {
      upload_ok(FALSE)
      user_expr_df(NULL)
      error_msg <- "Could not read the file. Make sure it is tab-delimited."
    } else if (ncol(df) != 2) {
      upload_ok(FALSE)
      user_expr_df(NULL)
      error_msg <- "File must have exactly 2 columns."
    } else if (!is.character(df[[1]])) {
      upload_ok(FALSE)
      user_expr_df(NULL)
      error_msg <- "First column must be character (gene names)."
    } else if (!is.numeric(df[[2]])) {
      upload_ok(FALSE)
      user_expr_df(NULL)
      error_msg <- "Second column must be numeric (t-statistics)."
    } else {
      colnames(df) <- c("Gene", "t")
      upload_ok(TRUE)
      user_expr_df(df)   # save dataframe globally
      error_msg <- NULL
    }

    output$file_warning_expr <- renderUI({
      if (upload_ok()) {
        div(style="color:green;font-weight:bold;margin-top:10px;", "Upload complete!")
      } else {
        div(style="color:red;font-weight:bold;margin-top:10px;", paste("Upload failed:", error_msg))
      }
    })
  })

  # Only show extra options if upload is successful
  output$user_file_options <- renderUI({
    if (!upload_ok()) return(NULL)

    tagList(
      radioButtons(
        inputId = "user_file_mode_expr",
        label = "Select correlation type:",
        choices = c("By Gene Expression" = "expr",
                    "By Gene Set Enrichment" = "gsea"),
        selected = "expr",
        inline = TRUE
      ),

      # Gene Expression options
      conditionalPanel(
        condition = "input.user_file_mode_expr == 'expr'",
        selectizeInput(
          inputId = "user_file_compare_expr",
          label = "Compare with specific RBPs (optional)",
          choices = names(Charm.object),
          multiple = TRUE,
          options = list(placeholder = "Select one or more RBPs")
        ),
        numericInput(
          "user_file_topN_expr",
          "Show top N correlated RBPs (optional):",
          value = NA, min = 1
        ),
        helpText("Tip: By selecting a number here, this will take precedence over the two inputs below."),
        numericInput(
          "user_file_n_pos_expr",
          "Show top N positive correlations (optional):",
          value = NA, min = 1
        ),
        numericInput(
          "user_file_n_neg_expr",
          "Show top N negative correlations (optional):",
          value = NA, min = 1
        ),
        helpText("Tip: You may use one or both of the numeric inputs above."),
        # Plot / Reset buttons
        div(
          style = "display: flex; align-items: center; margin-top: 15px;",
          actionButton(
            "user_file_plot_btn",
            tagList(fa("chart-line"), " Plot"),
            class = "btn btn-primary",
            style = "margin-right: 10px; border-radius: 20px;"
          ),
          actionButton(
            "user_file_reset_btn",
            tagList(fa("redo"), " Reset"),
            class = "btn btn-secondary",
            style = "border-radius: 20px;"
          )
        )
      ),

      # GSEA options
      conditionalPanel(
        condition = "input.user_file_mode_expr == 'gsea'",
        selectizeInput(
          inputId = "user_file_compare_gsea",
          label = "Compare with specific RBPs (optional)",
          choices = names(Charm.object),
          multiple = TRUE,
          options = list(placeholder = "Select one or more RBPs")
        ),
        numericInput(
          "user_file_topN_gsea",
          "Show top N correlated RBPs (optional):",
          value = NA, min = 1
        ),
        helpText("Tip: By selecting a number here, this will take precedence over the two inputs below."),
        numericInput(
          "user_file_n_pos_gsea",
          "Show top N positive correlations (optional):",
          value = NA, min = 1
        ),
        numericInput(
          "user_file_n_neg_gsea",
          "Show top N negative correlations (optional):",
          value = NA, min = 1
        ),
        helpText("Tip: You may use one or both of the numeric inputs above."),
        # Plot / Reset buttons
        div(
          style = "display: flex; align-items: center; margin-top: 15px;",
          actionButton(
            "user_file_plot_btn",
            tagList(fa("chart-line"), " Plot"),
            class = "btn btn-primary",
            style = "margin-right: 10px; border-radius: 20px;"
          ),
          actionButton(
            "user_file_reset_btn",
            tagList(fa("redo"), " Reset"),
            class = "btn btn-secondary",
            style = "border-radius: 20px;"
          )
        )
      )
    )
  })

  # ---- User File Similarity Plots ----
  observeEvent(input$user_file_plot_btn, {
    req(user_expr_df())
    mode <- input$user_file_mode_expr

    # Choose datasets and plotting functions
    if (mode == "expr") {
      datasets <- list(
        "Both Cells" = similar_expression_all,
        "K562"       = similar_expression_K562,
        "HEPG2"      = similar_expression_HEPG2
      )
      selected_rbps <- input$user_file_compare_expr
      scatter_fun <- correl_exp_rbp_plotly
      heatmap_fun <- exp_correl
      target_ui <- "userfilesimilar_expr"
    } else {  # mode == "gsea"
      datasets <- list(
        "Both Cells" = similar_gsea_all,
        "K562"       = similar_gsea_K562,
        "HEPG2"      = similar_gsea_HEPG2
      )
      selected_rbps <- input$user_file_compare_gsea
      scatter_fun <- correl_scatter_gsea_plotly
      heatmap_fun <- gsea_correl
      target_ui <- "userfilesimilar_gsea"
    }

    if (length(selected_rbps) == 0) selected_rbps <- NULL

    output[[target_ui]] <- renderUI({
      tagList(
        tags$div(
          "Generating plots, please wait...",
          style = "font-weight:bold;color:#A10702;margin-bottom:15px;"
        ),

        # Generate all dataset plots dynamically
        lapply(names(datasets), function(ds_name) {
          plotname <- paste0("userfile_plot_", ds_name)

          # Dynamic output creation
          if (!is.null(selected_rbps) && length(selected_rbps) == 1) {
            output[[plotname]] <- renderPlotly({
              scatter_fun(
                datasets[[ds_name]],
                user_expr_df(),
                selected_rbps[1],
                plot_title = ds_name
              )
            })
            plot_ui <- plotlyOutput(plotname, height = "500px")
          } else {
            output[[plotname]] <- renderPlot({
              heat_res <- heatmap_fun(
                datasets[[ds_name]],
                user_expr_df(),
                correl_num = if (mode == "expr") input$user_file_topN_expr else input$user_file_topN_gsea,
                n_pos      = if (mode == "expr") input$user_file_n_pos_expr else input$user_file_n_pos_gsea,
                n_neg      = if (mode == "expr") input$user_file_n_neg_expr else input$user_file_n_neg_gsea,
                other_rbps = if (!is.null(selected_rbps) && length(selected_rbps) > 1) selected_rbps else NULL
              )
              heat_res$heatmap   # ✅ FIXED
            })
            plot_ui <- plotOutput(plotname, height = "500px")
          }

          # Return UI element for each dataset
          column(
            width = 12,
            tags$h4(ds_name, style = "text-align:center;"),
            shinycssloaders::withSpinner(plot_ui)
          )
        }) %>% tagList()  # flatten
      )
    })
  })
  # ---- Reset button ----
  observeEvent(input$user_file_reset_btn, {
    output$userfilesimilar_expr <- renderUI(NULL)
    output$userfilesimilar_gsea <- renderUI(NULL)
    updateSelectizeInput(session, "user_file_compare_expr", selected = "")
    updateSelectizeInput(session, "user_file_compare_gsea", selected = "")
  })

  ### SPICING (user File)
  upload_ok_splice <- reactiveVal(FALSE)
  user_splice_df <- reactiveVal(NULL)

  # ---- File upload ----
  observeEvent(input$user_file_splice, {
    req(input$user_file_splice)
    file_path <- input$user_file_splice$datapath

    df <- tryCatch(
      read.table(file_path, header = TRUE, sep = " ", stringsAsFactors = FALSE),
      error = function(e) NULL
    )

    if (is.null(df)) {
      upload_ok_splice(FALSE)
      user_splice_df(NULL)
      error_msg <- "Could not read the file. Make sure it is space-delimited."
    } else if (ncol(df) != 4) {
      upload_ok_splice(FALSE)
      user_splice_df(NULL)
      error_msg <- "File must have exactly 4 columns."
    } else if (!is.character(df[[1]])) {
      upload_ok_splice(FALSE)
      user_splice_df(NULL)
      error_msg <- "First column must be character (Event IDs)."
    } else if (!is.numeric(df[[3]])) {
      upload_ok_splice(FALSE)
      user_splice_df(NULL)
      error_msg <- "Third column must be numeric (dPSI values)."
    } else {
      colnames(df) <- c("Event.ID", "Gene", "dPSI", "PDiff")
      upload_ok_splice(TRUE)
      user_splice_df(df)
      error_msg <- NULL

    }

    output$file_warning_splice <- renderUI({
      if (upload_ok_splice()) {
        div(style = "color:green;font-weight:bold;margin-top:10px;",
            "Upload complete!")
      } else {
        div(style = "color:red;font-weight:bold;margin-top:10px;",
            paste("Upload failed:", error_msg))
      }
    })
  })
  # ---- Only show extra options if upload is successful ----
  output$user_file_options_splice <- renderUI({
    if (!upload_ok_splice()) return(NULL)

    tagList(
      selectizeInput(
        inputId = "user_file_compare_splice",
        label = "Compare with specific RBPs (optional)",
        choices = names(Charm.object),
        multiple = TRUE,
        options = list(placeholder = "Select one or more RBPs")
      ),
      numericInput(
        "user_file_topN_splice",
        "Show top N correlated RBPs (optional):",
        value = NA, min = 1
      ),
      helpText("Tip: You may use one or both of the numeric inputs below."),
      numericInput(
        "user_file_n_pos_splice",
        "Show top N positive correlations (optional):",
        value = NA, min = 1
      ),
      numericInput(
        "user_file_n_neg_splice",
        "Show top N negative correlations (optional):",
        value = NA, min = 1
      ),
      div(
        style = "display: flex; align-items: center; margin-top: 15px;",
        actionButton(
          "user_file_plot_btn_splice",
          tagList(fa("chart-line"), " Plot"),
          class = "btn btn-primary",
          style = "margin-right: 10px; border-radius: 20px;"
        ),
        actionButton(
          "user_file_reset_btn_splice",
          tagList(fa("redo"), " Reset"),
          class = "btn btn-secondary",
          style = "border-radius: 20px;"
        )
      )
    )
  })

  # Auto-preview violin plot for uploaded file
  output$user_file_initial_plot_splice <- renderPlot({
    req(upload_ok_splice())          # only plot after successful upload
    df <- user_splice_df()
    req(df)

    # derive Type from Event.ID (same logic as your function)
    dpsi_table <- df %>%
      mutate(Type = case_when(
        startsWith(Event.ID, "HsaEX")  ~ "ES",
        startsWith(Event.ID, "HsaINT") ~ "IR",
        TRUE                           ~ NA_character_
      )) %>%
      filter(!is.na(Type))

    # If no ES/IR events present, show a message plot
    if (nrow(dpsi_table) == 0) {
      plot.new()
      title(main = "No ES or IR events detected in uploaded file")
      return(invisible(NULL))
    }

    # Calculate average ΔPSI by type (safe)
    avg_vals <- dpsi_table %>%
      group_by(Type) %>%
      summarise(mean_dPSI = mean(dPSI, na.rm = TRUE)) %>%
      pivot_wider(names_from = Type, values_from = mean_dPSI)

    # Build subtitle text safely (handle missing types)
    mean_ES <- ifelse("ES" %in% names(avg_vals), sprintf("%.3f", avg_vals$ES), "NA")
    mean_IR <- ifelse("IR" %in% names(avg_vals), sprintf("%.3f", avg_vals$IR), "NA")
    subtitle_text <- paste0("Mean ΔPSI — ES: ", mean_ES, " | IR: ", mean_IR)

    # Violin + jitter plot similar to your function
    ggplot(dpsi_table, aes(x = Type, y = dPSI)) +
      geom_jitter(width = 0.2, alpha = 0.6, size = 1.5) +
      geom_violin(alpha = .7) +
      theme_minimal() +
      ylab("ΔPSI (shRNA - CTRL)") +
      xlab("Event Type") +
      ggtitle("Uploaded data", subtitle = subtitle_text) +
      coord_flip() +
      geom_hline(yintercept = 0, linetype = "dashed") +
      theme(
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20, family = "Arial MS"),
        plot.subtitle = element_text(hjust = 0.5)
      )
  })

  # ---- Render similarity plots ----
  observeEvent(input$user_file_plot_btn_splice, {
    req(user_splice_df())
    selected_rbps <- input$user_file_compare_splice

    datasets <- list(
      "Both Cells" = Charm.object,
      "K562" = Charm.object_K562,
      "HEPG2" = Charm.object_HEPG2
    )

    output$similar_splice_plots_file <- renderUI({
      tagList(
        tags$div("Generating plots, please wait...", style = "font-weight:bold;color:#A10702;margin-bottom:15px;"),
        lapply(names(datasets), function(ds_name) {
          plotname <- paste0("userfile_plot_splice_", ds_name)

          if (!is.null(selected_rbps) && length(selected_rbps) == 1) {
            output[[plotname]] <- renderPlotly({
              correl_splicing_rbp_plotly(datasets[[ds_name]], user_splice_df(), selected_rbps)
            })
            plot_ui <- plotlyOutput(plotname, height = "500px")
          } else {
            output[[plotname]] <- renderPlot({
              heat_res <- splicing_correl(
                datasets[[ds_name]],
                user_splice_df(),
                correl_num = input$user_file_topN_splice,
                n_pos = input$user_file_n_pos_splice,
                n_neg = input$user_file_n_neg_splice,
                other_rbps = if (!is.null(selected_rbps) && length(selected_rbps) > 1) selected_rbps else NULL
              )
              heat_res$heatmap  # <-- make sure to return $heatmap, not $plot
            })
            plot_ui <- plotOutput(plotname, height = "500px")
          }

          column(
            width = 12,
            tags$h4(ds_name, style = "text-align:center;"),
            shinycssloaders::withSpinner(plot_ui)
          )
        }) %>% tagList()
      )
    })
  })
  # ---- Reset button ----
  observeEvent(input$user_file_reset_btn_splice, {
    # Clear the UI for user_file plots
    output$similar_splice_plots_file <- renderUI(NULL)

    # Reset the selectize input
    updateSelectizeInput(session, "user_file_compare_splice", selected = "")

    # Reset numeric inputs if any
    updateNumericInput(session, "user_file_topN_splice", value = NA)
    updateNumericInput(session, "user_file_n_pos_splice", value = NA)
    updateNumericInput(session, "user_file_n_neg_splice", value = NA)
  })

  #BINDING
  observe({
    if (is.null(input$user_file_binding)) return(NULL)
    file_path <- input$user_file_binding$datapath
    if (!file.exists(file_path)) return(NULL)

    df <- tryCatch(
      read.table(file_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE),
      error = function(e) NULL
    )

    upload_ok <- FALSE
    error_msg <- NULL

    if (is.null(df)) {
      error_msg <- "Could not read the file. Make sure it is tab-delimited."
    } else if (ncol(df) != 2) {
      error_msg <- "File must have exactly 2 columns."
    } else if (!is.character(df[[1]])) {
      error_msg <- "First column must be character (gene names)."
    } else if (!is.numeric(df[[2]])) {
      error_msg <- "Second column must be numeric (t-statistics)."
    } else {
      upload_ok <- TRUE
    }

    output$file_warning_binding <- renderUI({
      if (upload_ok) {
        div(style="color:green;font-weight:bold;margin-top:10px;", "Upload complete!")
      } else {
        div(style="color:red;font-weight:bold;margin-top:10px;", paste("Upload failed:", error_msg))
      }
    })
  })

  # ---- Volcano plot helper ----
  make_volcano_plot <- function(tbl, rbp) {
    if (is.null(tbl) || nrow(tbl) == 0) {
      ggplot() + annotate("text", x = 0, y = 0, label = paste("No results for", rbp)) + theme_void()
    } else {
      ggplot(tbl, aes(x = logFC, y = B,
                      text = paste0("Gene: ", gene,
                                    "<br>logFC: ", round(logFC,2),
                                    "<br>B: ", round(B,2),
                                    "<br>P: ", signif(P.Value,3)))) +
        geom_point(aes(color = highlight), alpha = 0.7) +
        scale_color_manual(values = c("None"="#CCCCCC",
                                      "RBP"="#A10702",
                                      "Selected"="#008057")) +
        theme_bw() +
        labs(title=paste(rbp,"KD"), x="Log2 Fold-Change", y="B-statistic") +
        theme(legend.position="none", plot.title=element_text(hjust=0.5)) +
        theme(
          axis.line = element_line(colour = "black"),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          text = element_text(size = 15, face = "bold"),
          legend.position = "none"
        )
    }
  }

  # ---- React to clicking a point in the volcano ----
  observeEvent(event_data("plotly_click", source = "volcano"), {
    ed <- event_data("plotly_click", source = "volcano")
    tbl <- display_table()
    if (!is.null(ed) && nrow(tbl) > 0) {
      # Find closest point
      clicked_gene <- tbl$gene[which.min(abs(tbl$logFC - ed$x) + abs(tbl$B - ed$y))]

      # Update highlights
      tbl$highlight <- ifelse(tbl$gene == clicked_gene, "Selected",
                              ifelse(tbl$gene == rbp_current(), "RBP", "None"))

      # Reorder table to bring selected gene on top
      tbl <- rbind(tbl[tbl$gene == clicked_gene, ], tbl[tbl$gene != clicked_gene, ])
      display_table(tbl)

      # Re-render volcano plot
      output$volcano_plot <- renderPlotly({
        ggplotly(make_volcano_plot(display_table(), rbp_current()), tooltip = "text", source = "volcano") %>%
          event_register("plotly_click")
      })
    }
  })

  # ---- React to selecting a row in the table ----
  observeEvent(input$volcano_table_rows_selected, {
    sel_row <- input$volcano_table_rows_selected
    if (is.null(sel_row)) return()
    tbl <- display_table()
    sel_gene <- tbl$gene[sel_row]

    # Update highlights
    tbl$highlight <- ifelse(tbl$gene == sel_gene, "Selected",
                            ifelse(tbl$gene == rbp_current(), "RBP", "None"))
    display_table(tbl)

    # Re-render volcano plot
    output$volcano_plot <- renderPlotly({
      ggplotly(make_volcano_plot(display_table(), rbp_current()), tooltip = "text", source = "volcano") %>%
        event_register("plotly_click")
    })
  })

  # ---- Main: search button triggers all plots ----
  observeEvent(input$search_btn_rbp, {
    req(input$expr_search_rbp)
    rbp_sel <- input$expr_search_rbp
    rbp_current(rbp_sel)

    charm_obj <- current_charm_expr()
    shrna_obj <- current_shrna_expr()

    exp_list <- names(charm_obj)
    if (is.null(exp_list) || !(rbp_sel %in% exp_list)) {
      showModal(modalDialog(
        title = "Warning",
        paste("This RBP has no information available for", input$cell_line),
        easyClose = TRUE,
        footer = NULL
      ))
      return(NULL)
    }

    session$sendCustomMessage("toggleCursor", TRUE)

    output$expr_violin <- renderPlot({ violinplotter(charm_obj, rbp_sel) })

    output$shrna_plot <- renderPlot({ plot_shRNA_effect(shrna_obj, rbp_sel) })
    output$shrna_warning <- renderUI({
      stats <- shrna_obj[rbp_sel,,drop=FALSE]
      if (is.null(stats) || nrow(stats)==0) return(NULL)
      logFC <- stats$logFC
      pval <- stats$P.Value
      if (pval > 0.05 || logFC > -0.5) {
        div(style="margin-top:10px;font-size:16px;font-weight:bold;color:red;",
            "WARNING: The efficiency of this knockdown is uncertain. Proceed with caution.")
      } else NULL
    })

    result <- plot_rbp_volcano(charm_obj, rbp_sel)
    tbl <- as.data.frame(result$top_table)
    if (!"gene" %in% colnames(tbl)) tbl$gene <- rownames(tbl)
    tbl <- tbl[, c("gene", setdiff(colnames(tbl),"gene"))]
    tbl$highlight <- ifelse(tbl$gene==rbp_sel,"RBP","None")
    display_table(tbl)

    output$volcano_table <- renderDT({
      datatable(display_table(),
                filter = "none",
                selection = "single",
                rownames = FALSE,
                options = list(pageLength = 10, scrollX = TRUE, searching = TRUE))
    })

    output$volcano_plot <- renderPlotly({
      ggplotly(make_volcano_plot(display_table(), rbp_sel), tooltip="text", source="volcano") %>%
        event_register("plotly_click")
    })

    gsea_result <- plot_gsea(charm_obj, rbp_sel, thresh = 0.05)
    output$gsea_plot <- renderPlot({ gsea_result$gsea_plot })
    output$geneset_table <- renderDT({
      datatable(gsea_result$geneset_table,
                filter = "none",
                rownames = FALSE,
                options = list(pageLength = 10, scrollX = TRUE))
    })

    session$sendCustomMessage("toggleCursor", FALSE)
  })
  # ---- Gene Search ----
  observeEvent(input$search_btn_gene, {
    req(input$expr_search_gene)
    gene <- input$expr_search_gene
    charm_obj <- current_charm_expr()

    output$expr_gene_plot <- renderPlot({
      validate(
        need(gene %in% unlist(lapply(charm_obj, function(x) rownames(x$DEGenes))),
             paste("Gene", gene, "not found in any RBP DEGenes tables."))
      )
      plot_gene_logFC_barplot(charm_obj, gene)
    })

    session$sendCustomMessage("toggleCursor", FALSE) # turn cursor back
  })

  # ---- Hallmark Gene Set Search ----
  observeEvent(input$search_btn_hallmark, {
    req(input$expr_search_hallmark)
    geneset <- input$expr_search_hallmark
    charm_obj <- current_charm_expr()

    output$expr_hallmark_plot <- renderPlot({
      validate(
        need(geneset %in% unlist(lapply(charm_obj, function(x) x$GSEA$pathway)),
             paste0("Geneset '", geneset, "' not found in any RBP GSEA tables."))
      )
      plot_hallmark_nes_barplot(charm_obj, geneset)
    })

    session$sendCustomMessage("toggleCursor", FALSE)
  })

  # Reactive values to store what the user selected
  similar_plot_inputs <- reactiveVal(list(
    rbp1 = NULL,
    selected_rbps = NULL,
    mode = NULL
  ))

  # ---- Populate Similar RBPs selectize inputs ----
  observe({
    req(Charm.object)  # make sure Charm.object exists
    rbp_names <- names(Charm.object)

    # Expression correlation selectize
    updateSelectizeInput(
      session,
      "similar_rbps_select_expr",
      choices = rbp_names,
      server = TRUE
    )

    # GSEA correlation selectize
    updateSelectizeInput(
      session,
      "similar_rbps_select_gsea",
      choices = rbp_names,
      server = TRUE
    )
  })
  # ---- Similar RBPs: Expression/GSEA correlation ----
  observeEvent(input$similar_plot_btn, {
    req(input$expr_dataset)
    if (input$expr_dataset != "Similar RBPs") return(NULL)
    req(input$expr_search_rbp)
    rbp1 <- input$expr_search_rbp
    mode <- input$similar_mode

    # Datasets
    datasets <- list(
      "Both Cells" = if(mode == "expr") similar_expression_all else similar_gsea_all,
      "K562"       = if(mode == "expr") similar_expression_K562 else similar_gsea_K562,
      "HEPG2"      = if(mode == "expr") similar_expression_HEPG2 else similar_gsea_HEPG2
    )

    # Functions
    scatter_fun <- if(mode == "expr") correl_exp_rbp_plotly else correl_scatter_gsea_plotly
    heatmap_fun <- if(mode == "expr") exp_correl else gsea_correl

    # Inputs
    selected_rbps <- if(mode == "expr") input$similar_rbps_select_expr else input$similar_rbps_select_gsea
    correl_num <- if(mode == "expr") input$correl_num_expr else input$correl_num_gsea
    n_pos      <- if(mode == "expr") input$n_pos_expr else input$n_pos_gsea
    n_neg      <- if(mode == "expr") input$n_neg_expr else input$n_neg_gsea

    if(length(selected_rbps) == 0) selected_rbps <- NULL

    # Render plots UI
    output$similar_expr_plots <- renderUI({
      # Use lapply and wrap columns in tagList() properly
      tagList(
        lapply(names(datasets), function(ds_name) {
          plotname <- paste0("similar_plot_", ds_name)

          if (!is.null(selected_rbps) && length(selected_rbps) == 1) {
            output[[plotname]] <- renderPlotly({
              scatter_fun(datasets[[ds_name]], rbp1, selected_rbps[1], plot_title = ds_name)
            })
            plot_ui <- plotlyOutput(plotname, height = "500px")
          } else {
            output[[plotname]] <- renderPlot({
              heat_res <- heatmap_fun(
                datasets[[ds_name]], rbp1,
                correl_num = correl_num,
                n_pos = n_pos,
                n_neg = n_neg,
                other_rbps = if(!is.null(selected_rbps) && length(selected_rbps) > 1) selected_rbps else NULL
              )
              # always return ggplot object
              heat_res$heatmap
            })
            plot_ui <- plotOutput(plotname, height = "500px")
          }

          # Return the column for renderUI
          column(
            width = 12,
            tags$h4(ds_name, style = "text-align:center;"),
            shinycssloaders::withSpinner(plot_ui)
          )
        }) %>% tagList()  # <-- flatten the list of columns
      )
    })
  })

  # ---- Main: Splicing Explore ----
  observeEvent(input$splice_search_btn, {
    req(input$splice_search)
    rbp_sel <- input$splice_search
    rbp_current(rbp_sel)

    charm_obj <- current_charm_splice()
    shrna_obj <- current_shrna_splice()

    splice_list <- names(charm_obj)
    if (is.null(splice_list) || !(rbp_sel %in% splice_list)) {
      showModal(modalDialog(
        title = "Warning",
        paste("This RBP has no splicing information available."),
        easyClose = TRUE,
        footer = NULL
      ))
      return(NULL)
    }

    session$sendCustomMessage("toggleCursor", TRUE)



    # ---- Top row: Violin + shRNA knockdown plots ----
    output$violin_splice_plot <- renderPlot({
      violin_splice_plot(charm_obj, rbp_sel)
    })

    output$plot_shrna_effect <- renderPlot({
      plot_shRNA_effect(shrna_obj, rbp_sel)
    })

    output$shrna_warning_splice <- renderUI({
      stats <- shrna_obj[rbp_sel,,drop=FALSE]
      if (is.null(stats) || nrow(stats) == 0) return(NULL)
      logFC <- stats$logFC
      pval <- stats$P.Value
      if (pval > 0.05 || logFC > -0.5) {
        div(style="margin-top:10px;font-size:16px;font-weight:bold;color:red;",
            "WARNING: The efficiency of this knockdown is uncertain. Proceed with caution.")
      } else NULL
    })

    # ---- Second row: Volcano plot ----
    splice_res <- plot_splice_volcano(charm_obj, rbp_sel)
    if (is.null(splice_res)) {
      output$plot_splice_volcano <- renderPlot({
        ggplot() + annotate("text", x = 0, y = 0, label = paste("No results for", rbp_sel)) + theme_void()
      })
      output$splice_volcano_table <- DT::renderDataTable(NULL)
      return()
    }

    # Store results in reactiveVal
    result_splice(splice_res)

    # ---- Helper to render volcano ----
    render_volcano <- function(tbl) {
      p <- ggplot(tbl, aes(x = dPSI, y = Pdiff, key = Event.ID, color = highlight)) +
        geom_point(alpha = 0.7) +
        scale_color_manual(values = c("None"="#CCCCCC","RBP"="#A10702","Selected"="#008057")) +
        theme_bw() +
        theme(legend.position = "none") +
        labs(title = paste("Volcano Plot:", rbp_sel), x = "ΔPSI (shRNA - CTRL)", y = "PDiff")
      ggplotly(p, tooltip = "key", source = "splice_volcano") %>%
        event_register("plotly_click")
    }

    # ---- Render volcano plot ----
    output$plot_splice_volcano <- renderPlotly({
      req(result_splice())
      render_volcano(result_splice()$top_table)
    })

    # ---- Render table ----
    output$splice_volcano_table <- DT::renderDataTable({
      req(result_splice())
      result_splice()$top_table %>% select(-text)
    }, rownames = FALSE, selection = "single", options = list(pageLength = 10))

    # ---- Shared function to update highlights ----
    update_highlight <- function(selected_event) {
      if (is.null(selected_event) || length(selected_event) == 0) return(NULL)
      tbl <- result_splice()$top_table
      if (!selected_event %in% tbl$Event.ID) return(NULL)

      tbl$highlight <- ifelse(tbl$Event.ID == selected_event, "Selected",
                              ifelse(tbl$Event.ID %in% input$splice_search, "RBP", "None"))

      result_splice(list(
        top_table = tbl,
        volcano_plot = result_splice()$volcano_plot
      ))

      # Sync table selection
      proxy <- DT::dataTableProxy("splice_volcano_table")
      sel_idx <- which(tbl$highlight == "Selected")
      DT::selectRows(proxy, sel_idx)
    }

    # ---- React to selecting a row in the table ----
    observeEvent(input$splice_volcano_table_rows_selected, {
      sel_row <- input$splice_volcano_table_rows_selected
      if (is.null(sel_row)) return()
      selected_event <- result_splice()$top_table$Event.ID[sel_row]
      update_highlight(selected_event)
    })

    # ---- React to clicking a point in the volcano plot ----
    observeEvent(event_data("plotly_click", source = "splice_volcano"), {
      click <- event_data("plotly_click", source = "splice_volcano")
      if (is.null(click) || is.null(click$key)) return()
      clicked_event <- click$key
      update_highlight(clicked_event)
    })

    session$sendCustomMessage("toggleCursor", FALSE)
  })

  # ---- Event ID Search ----
  observeEvent(input$splice_event_btn, {
    req(input$splice_event_search)
    event_id <- input$splice_event_search
    charm_obj <- current_charm_splice()

    output$heatmap_splicing_dpsi <- renderPlot({
      plot_event_dpsi_barplot(charm_obj, event_id)
    })
  })

  # ---- Similar RBPs: Splicing correlation ----
  observeEvent(input$similar_plot_btn_splice, {
    req(input$splice_dataset)
    if (input$splice_dataset != "Similar RBPs") return(NULL)
    req(input$splice_search)
    rbp1 <- input$splice_search

    # Only one mode currently
    mode <- "splice"

    # Choose datasets
    datasets <- list(
      "Both Cells" = Charm.object,
      "K562"       = Charm.object_K562,
      "HEPG2"      = Charm.object_HEPG2
    )

    scatter_fun <- correl_splicing_rbp_plotly
    heatmap_fun <- splicing_correl

    selected_rbps <- input$similar_rbps_select_splice
    correl_num <- input$correl_num_splice
    n_pos      <- input$n_pos_splice
    n_neg      <- input$n_neg_splice

    if (length(selected_rbps) == 0) selected_rbps <- NULL

    output$similar_splice_plots <- renderUI({
      tagList(
        lapply(names(datasets), function(ds_name) {
          plotname <- paste0("similar_splice_plot_", ds_name)

          if (!is.null(selected_rbps) && length(selected_rbps) == 1) {
            output[[plotname]] <- renderPlotly({
              scatter_fun(datasets[[ds_name]], rbp1, selected_rbps)
            })
            plot_ui <- plotlyOutput(plotname, height = "500px")
          } else {
            output[[plotname]] <- renderPlot({
              heat_res <- heatmap_fun(
                datasets[[ds_name]], rbp1,
                correl_num = correl_num,
                n_pos = n_pos,
                n_neg = n_neg,
                other_rbps = if (!is.null(selected_rbps) && length(selected_rbps) > 1) selected_rbps else NULL
              )
              heat_res$heatmap   # ✅ FIXED LINE
            })
            plot_ui <- plotOutput(plotname, height = "500px")
          }

          column(
            width = 12,
            tags$h4(ds_name, style = "text-align:center;"),
            shinycssloaders::withSpinner(plot_ui)
          )
        }) %>% tagList()  # flatten
      )
    })
  })


  # ---- Reset Similar RBPs Splicing plots ----
  observeEvent(input$similar_reset_btn_splice, {
    output$similar_splice_plots <- renderUI(NULL)
    updateSelectizeInput(session, "similar_rbps_select_splice", selected = "")
    updateNumericInput(session, "correl_num_splice", value = NA)
    updateNumericInput(session, "n_pos_splice", value = NA)
    updateNumericInput(session, "n_neg_splice", value = NA)
  })

  # ---- Reset buttons Expression ----
  observeEvent(input$similar_reset_btn, {
    output$similar_expr_plots <- renderUI(NULL)
    updateSelectizeInput(session, "similar_rbps_select_expr", selected = "")
    updateSelectizeInput(session, "similar_rbps_select_gsea", selected = "")
  })

  observeEvent(input$reset_btn_rbp, {
    output$expr_violin    <- renderPlot(NULL)
    output$shrna_plot     <- renderPlot(NULL)
    output$shrna_warning  <- renderUI(NULL)
    output$volcano_plot   <- renderPlotly(NULL)
    output$volcano_table  <- renderDT(NULL)
    output$gsea_plot      <- renderPlot(NULL)
    output$geneset_table  <- renderDT(NULL)
    display_table(NULL)
    rbp_current(NULL)
    updateSelectizeInput(session, "expr_search", selected = "")
  })
  observeEvent(input$reset_btn_gene, {
    output$expr_gene_plot <- renderPlot(NULL)
  })

  observeEvent(input$reset_btn_hallmark, {
    output$expr_hallmark_plot <- renderPlot(NULL)
  })
  # ---- Reset buttons for Splicing tab ----

  # 1️⃣ Similar RBPs reset
  observeEvent(input$similar_reset_btn_splice, {
    output$similar_splice_plots <- renderUI(NULL)
    updateSelectizeInput(session, "similar_rbps_select_splice", selected = "")
  })

  # 2️⃣ RBP search reset
  observeEvent(input$splice_reset_btn, {
    output$violin_splice_plot    <- renderPlot(NULL)
    output$plot_shrna_effect     <- renderPlot(NULL)
    output$shrna_warning_splice  <- renderUI(NULL)
    output$plot_splice_volcano   <- renderPlotly(NULL)
    output$splice_volcano_table  <- renderDT(NULL)
    result_splice(NULL)  # if using reactiveVal to store results
    rbp_current(NULL)    # if using reactiveVal for selected RBP
    updateSelectizeInput(session, "splice_search", selected = "")
  })

  # 3️⃣ Event ID search reset
  # Reset Event heatmap
  observeEvent(input$splice_event_reset_btn, {
    output$heatmap_splicing_dpsi <- renderPlot(NULL)  # clear plot
    updateSelectizeInput(session, "splice_event_search", selected = "")  # clear search
  })
}


# end server


shinyApp(ui, server)

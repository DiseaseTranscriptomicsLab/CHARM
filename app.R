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


###### UI

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
          column(
            width = 3,
            wellPanel(
              radioButtons(
                inputId = "expr_mode",
                label = "Select mode:",
                choices = c("Explore Mode", "Discovery Mode"),
                selected = "Explore Mode"
              ),
              conditionalPanel(
                condition = "input.expr_mode == 'Explore Mode'",
                selectInput(
                  inputId = "expr_dataset",
                  label = "Select option:",
                  choices = c("Both Cells", "K562", "HEPG2", "Similar RBPs")
                ),
                div(
                  style = "display: flex; align-items: center; margin-top: 15px;",
                  selectizeInput(
                    inputId = "expr_search",
                    label = NULL,
                    choices = NULL,
                    multiple = FALSE,
                    options = list(placeholder = "RBP"),
                    width = "400px"
                  ),
                  conditionalPanel(
                    condition = "input.expr_dataset != 'Similar RBPs'",
                    div(
                      style = "display: flex; align-items: center; margin-left: 10px;",
                      actionButton("search_btn", tagList(fa("search"), " Search"), class = "btn btn-primary", style = "margin-right: 10px; border-radius: 20px;"),
                      actionButton("reset_btn", tagList(fa("redo"), " Reset"), class = "btn btn-secondary", style = "border-radius: 20px;")
                    )
                  )
                ),
                conditionalPanel(
                  condition = "input.expr_dataset == 'Similar RBPs'",
                  radioButtons(
                    inputId = "similar_mode",
                    label = "Select correlation type:",
                    choices = c("By Gene Expression" = "expr", "By Gene Set Enrichment" = "gsea"),
                    selected = "expr",
                    inline = TRUE
                  ),
                  conditionalPanel(
                    condition = "input.similar_mode == 'expr'",
                    selectizeInput("similar_rbps_select_expr", "Compare with specific RBPs (optional)", choices = NULL, multiple = TRUE, options = list(placeholder = "Select one or more RBPs")),
                    helpText("• If you select one RBP, a scatter plot will be shown.\n• If you select multiple RBPs, correlation values will be summarised."),
                    numericInput("correl_num_expr", "Show top N correlated RBPs (optional):", value = NA, min = 1),
                    numericInput("n_pos_expr", "Show top N positive correlations (optional):", value = NA, min = 1),
                    numericInput("n_neg_expr", "Show top N negative correlations (optional):", value = NA, min = 1)
                  ),
                  conditionalPanel(
                    condition = "input.similar_mode == 'gsea'",
                    selectizeInput("similar_rbps_select_gsea", "Compare with specific RBPs (optional)", choices = NULL, multiple = TRUE, options = list(placeholder = "Select one or more RBPs")),
                    numericInput("correl_num_gsea", "Show top N correlated RBPs (optional):", value = NA, min = 1),
                    numericInput("n_pos_gsea", "Show top N positive correlations (optional):", value = NA, min = 1),
                    numericInput("n_neg_gsea", "Show top N negative correlations (optional):", value = NA, min = 1)
                  ),
                  div(
                    style = "display: flex; align-items: center; margin-top: 15px;",
                    actionButton("similar_plot_btn", tagList(fa("chart-line"), " Plot"), class = "btn btn-primary", style = "margin-right: 10px; border-radius: 20px;"),
                    actionButton("similar_reset_btn", tagList(fa("redo"), " Reset"), class = "btn btn-secondary", style = "border-radius: 20px;")
                  )
                )
              ),
              conditionalPanel(
                condition = "input.expr_mode == 'Discovery Mode'",
                hr(),
                tags$p("Alternatively, upload your own table with differential expression values:"),
                fileInput("user_file_expr", "Upload your file:", accept = c(".txt")),
                uiOutput("file_warning_expr"),
                uiOutput("user_file_options")
              )
            )
          ),
          # Right content area
          column(
            width = 9,
            tags$h3("Expression Data Results"),
            div(
              "⚠ Please press reset after every plot!",
              style = "border: 2px solid #f0ad4e; background-color: #fff3cd; padding: 8px; border-radius: 6px; font-weight: bold; color: #856404;"
            ),
            conditionalPanel(condition = "input.expr_dataset == 'Similar RBPs'", uiOutput("similar_expr_plots")),
            conditionalPanel(
              condition = "input.expr_dataset != 'Similar RBPs'",
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
                selectInput("splice_dataset", "Select option:", choices = c("Both Cells", "K562", "HEPG2", "Similar RBPs")),
                div(
                  style = "display: flex; align-items: center; margin-top: 15px;",
                  selectizeInput("splice_search", NULL, choices = NULL, multiple = FALSE, options = list(placeholder = "RBP"), width = "400px"),
                  conditionalPanel(
                    condition = "input.splice_dataset != 'Similar RBPs'",
                    div(
                      style = "display: flex; align-items: center; margin-left: 10px;",
                      actionButton("splice_search_btn", tagList(fa("search"), " Search"), class = "btn btn-primary", style = "margin-right: 10px; border-radius: 20px;"),
                      actionButton("splice_reset_btn", tagList(fa("redo"), " Reset"), class = "btn btn-secondary", style = "border-radius: 20px;")
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
                tags$p("Alternatively, upload your own table with differential splicing values:"),
                fileInput("user_file_splice", "Upload your file:", accept = c(".txt")),
                uiOutput("file_warning_splice"),
                uiOutput("user_file_options_splice")
              )
            )
          ),
          column(
            width = 9,
            tags$h3("Splicing Data Results"),
            div(
              "⚠ Please press reset after every plot!",
              style = "border: 2px solid #f0ad4e; background-color: #fff3cd; padding: 8px; border-radius: 6px; font-weight: bold; color: #856404;"
            ),
            conditionalPanel(condition = "input.splice_dataset == 'Similar RBPs'", uiOutput("similar_splice_plots")),
            conditionalPanel(
              condition = "input.splice_dataset != 'Similar RBPs'",
              fluidRow(
                column(width = 6, plotOutput("violin_splice_plot", height = "400px")),
                column(width = 6, plotOutput("plot_shrna_effect", height = "400px"))
              ),
              fluidRow(
                column(width = 6, plotlyOutput("plot_splice_volcano", height = "450px")),
                column(width = 6, DTOutput("splice_volcano_table"))
              )
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

  # ---- Search bar population ----
  observe({
    req(Charm.object)
    updateSelectizeInput(session, "expr_search",
                         choices = names(Charm.object),
                         server = TRUE)
  })

  observe({
    req(Charm.object)
    updateSelectizeInput(session, "splice_search",
                         choices = names(Charm.object),
                         server = TRUE)
  })

  observe({
    rbp_choices_expr <- names(current_charm_expr())

    updateSelectizeInput(session, "similar_rbps_select_expr", choices = rbp_choices_expr)
    updateSelectizeInput(session, "similar_rbps_select_gsea", choices = rbp_choices_expr)
  })

  observe({
    rbp_choices_splice <- names(current_charm_splice())

    updateSelectizeInput(session, "similar_rbps_select_splice", choices = rbp_choices_splice)
  })

  ###EXPRESSION
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
        tags$div("Generating plots, please wait...",
                 style = "font-weight:bold;color:#A10702;margin-bottom:15px;"),

        lapply(names(datasets), function(ds_name) {
          plotname <- paste0("userfile_plot_", ds_name)

          output[[plotname]] <- if (!is.null(selected_rbps) && length(selected_rbps) == 1) {
            renderPlotly({
              scatter_fun(datasets[[ds_name]], user_expr_df(), selected_rbps, plot_title = ds_name)
            })
          } else {
            renderPlot({
              heat_res <- heatmap_fun(
                datasets[[ds_name]], user_expr_df(),
                correl_num = if (mode == "expr") input$user_file_topN_expr else input$user_file_topN_gsea,
                n_pos      = if (mode == "expr") input$user_file_n_pos_expr else input$user_file_n_pos_gsea,
                n_neg      = if (mode == "expr") input$user_file_n_neg_expr else input$user_file_n_neg_gsea,
                other_rbps = if (!is.null(selected_rbps) && length(selected_rbps) > 1) selected_rbps else NULL
              )
              heat_res$heatmap
            })
          }

          column(
            width = 12,
            tags$h4(ds_name, style="text-align:center;"),
            if (!is.null(selected_rbps) && length(selected_rbps) == 1) {
              shinycssloaders::withSpinner(plotlyOutput(plotname, height="500px"))
            } else {
              shinycssloaders::withSpinner(plotOutput(plotname, height="500px"))
            }
          )
        })
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

  #SPLICING
  observe({
    if (is.null(input$user_file_splicing)) return(NULL)
    file_path <- input$user_file_splicing$datapath
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

    output$file_warning_splicing <- renderUI({
      if (upload_ok) {
        div(style="color:green;font-weight:bold;margin-top:10px;", "Upload complete!")
      } else {
        div(style="color:red;font-weight:bold;margin-top:10px;", paste("Upload failed:", error_msg))
      }
    })
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
  observeEvent(input$search_btn, {
    req(input$expr_search)
    rbp_sel <- input$expr_search
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

  # ---- Similar RBPs: Expression/GSEA correlation ----
  observeEvent(input$similar_plot_btn, {
    req(input$expr_search)
    rbp1 <- input$expr_search
    mode <- input$similar_mode

    # Choose datasets and plotting functions
    if (mode == "expr") {
      datasets <- list(
        "Both Cells" = similar_expression_all,
        "K562"       = similar_expression_K562,
        "HEPG2"      = similar_expression_HEPG2
      )
      selected_rbps <- input$similar_rbps_select_expr
      scatter_fun <- correl_exp_rbp_plotly
      heatmap_fun <- exp_correl
    } else {  # mode == "gsea"
      datasets <- list(
        "Both Cells" = similar_gsea_all,
        "K562"       = similar_gsea_K562,
        "HEPG2"      = similar_gsea_HEPG2
      )
      selected_rbps <- input$similar_rbps_select_gsea
      scatter_fun <- correl_scatter_gsea_plotly
      heatmap_fun <- gsea_correl
    }

    if (length(selected_rbps) == 0) selected_rbps <- NULL

    output$similar_expr_plots <- renderUI({
      tagList(
        tags$div("Generating plots, please wait...",
                 style = "font-weight:bold;color:#A10702;margin-bottom:15px;"),

        lapply(names(datasets), function(ds_name) {
          plotname <- paste0("similar_plot_", ds_name)

          # Render scatter if 1 RBP, otherwise heatmap
          output[[plotname]] <- if (!is.null(selected_rbps) && length(selected_rbps) == 1) {
            renderPlotly({
              scatter_fun(datasets[[ds_name]], rbp1, selected_rbps, plot_title = ds_name)
            })
          } else {
            renderPlot({
              heat_res <- heatmap_fun(
                datasets[[ds_name]], rbp1,
                correl_num = if (mode == "expr") input$correl_num_expr else input$correl_num_gsea,
                n_pos      = if (mode == "expr") input$n_pos_expr      else input$n_pos_gsea,
                n_neg      = if (mode == "expr") input$n_neg_expr      else input$n_neg_gsea,
                other_rbps = if (!is.null(selected_rbps) && length(selected_rbps) > 1) selected_rbps else NULL
              )
              heat_res$heatmap
            })
          }

          column(
            width = 12,   # full-width vertical layout
            tags$h4(ds_name, style="text-align:center;"),
            if (!is.null(selected_rbps) && length(selected_rbps) == 1) {
              shinycssloaders::withSpinner(plotlyOutput(plotname, height="500px"))
            } else {
              shinycssloaders::withSpinner(plotOutput(plotname, height="500px"))
            }
          )
        })
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

  # ---- Reset buttons ----
  observeEvent(input$similar_reset_btn, {
    output$similar_expr_plots <- renderUI(NULL)
    updateSelectizeInput(session, "similar_rbps_select_expr", selected = "")
    updateSelectizeInput(session, "similar_rbps_select_gsea", selected = "")
  })

  observeEvent(input$reset_btn, {
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

}


# end server


shinyApp(ui, server)

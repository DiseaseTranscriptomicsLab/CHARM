library(shiny)
library(shinythemes)
library(fontawesome)
library(DT)
library(plotly)
library(ggplot2)
library(limma)
library(ggrepel)
library(shinycssloaders)

source("helper_functions.R")

# Load CHARM object once when app starts
Charm.object <- readRDS("data/Charm.object.RDS")
sh_effect_vector <- readRDS("data/shRNA_Efficiency.Rds")


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
      tagList(
        fa("home", fill = "black", height = "1em"),
        " Home"
      ),
      fluidPage(
        # Title / Logo row
        fluidRow(
          column(
            4,
            div(
              style = "display: flex; align-items: center; justify-content: center; margin: 30px 0;",
              div(style = "margin-right: 15px;", fa("gem", fill = "black", height = "4em")),
              div(
                style = "text-align: left;",
                tags$h2("CHARM", style = "margin: 0; font-weight: bold;"),
                tags$h4("Comprehensive Hub of Alternative Regulatory Mapping",
                        style = "margin: 0; font-weight: bold;")
              )
            )
          )
        ),

        # Explore & Discover boxes
        fluidRow(
          column(
            6,
            div(
              style = "background-color: #EAEBEB; color: black;
                       border: 2px solid #2c3e50; border-radius: 20px;
                       padding: 40px; text-align: center; margin: 15px;",
              div(style = "margin-bottom: 10px;", fa("globe", fill = "black", height = "3em")),
              tags$h3("Explore", style = "font-weight: bold;"),
              tags$p("Investigate how a known RBP affects expression, splicing, and binding, based on ENCODE's gene silencing series and eCLIP data.")
            )
          ),
          column(
            6,
            div(
              style = "background-color: #EAEBEB; color: black;
                       border: 2px solid #2c3e50; border-radius: 20px;
                       padding: 40px; text-align: center; margin: 15px;",
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

        # Main layout: sidebar (left) + results (right)
        fluidRow(

          # Left sidebar for mode selection + search
          column(
            width = 3,
            wellPanel(
              # Choose mode: Explore or Discovery
              radioButtons(
                inputId = "expr_mode",
                label = "Select mode:",
                choices = c("Explore Mode", "Discovery Mode"),
                selected = "Explore Mode"
              ),

              # Conditional: show dropdown only in Explore Mode
              conditionalPanel(
                condition = "input.expr_mode == 'Explore Mode'",
                selectInput(
                  inputId = "expr_dataset",
                  label = "Select option:",
                  choices = c("Both Cells", "K562", "HEPG2", "Similar RBPs")
                ),

                # Search bar + buttons below the dropdown
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

                  actionButton(
                    "search_btn",
                    tagList(fa("search"), " Search"),
                    class = "btn btn-primary",
                    style = "margin-left: 10px; border-radius: 20px;"
                  ),

                  actionButton(
                    "reset_btn",
                    tagList(fa("redo"), " Reset"),
                    class = "btn btn-secondary",
                    style = "margin-left: 10px; border-radius: 20px;"
                  )
                )
              ),

              # Conditional: show file upload only in Discovery Mode
              conditionalPanel(
                condition = "input.expr_mode == 'Discovery Mode'",
                hr(),
                tags$p("Alternatively, upload your own table with differential expression values:"),
                tags$p(
                  "CHARM currently accepts tables where the first column are gene names, ",
                  "and the second column are differential expression values (t-statistics). ",
                  "This table must be in .txt format and must not contain a header."
                ),

                fileInput(
                  inputId = "user_file",
                  label = "Upload your file:",
                  accept = c(".txt")
                ),

                # Placeholder for inline warning
                uiOutput("file_warning")
              )
            )
          ),

          # Right content area
          column(
            width = 9,
            tags$h3("Expression Data Results"),

            # Top plots: Violin + shRNA
            fluidRow(
              column(
                width = 6,
                plotOutput("expr_violin"),
                uiOutput("expr_note")
              ),
              column(
                width = 6,
                plotOutput("shrna_plot"),
                uiOutput("shrna_warning")
              )
            ),

            # Add vertical space before volcano + table
            div(style = "margin-top: 50px;",
                fluidRow(
                  column(
                    width = 6,
                    plotlyOutput("volcano_plot", height = "600px")
                  ),
                  column(
                    width = 6,
                    DTOutput("volcano_table")
                  )
                )
            ),

            uiOutput("expr_placeholder")
          )

        )
      )
    )
    , # end tabPanel
    # Splicing tab
    tabPanel(
      tagList(fa("scissors", fill = "black", height = "1em"), " Splicing"),
      fluidPage(

        # Main layout: sidebar (left) + results (right)
        fluidRow(
          # Left sidebar for mode selection + search
          column(
            width = 3,
            wellPanel(
              # Choose mode: Explore or Discovery
              radioButtons(
                inputId = "expr_mode",
                label = "Select mode:",
                choices = c("Explore Mode", "Discovery Mode"),
                selected = "Explore Mode"
              ),

              # Conditional: show dropdown only in Explore Mode
              conditionalPanel(
                condition = "input.expr_mode == 'Explore Mode'",
                selectInput(
                  inputId = "expr_dataset",
                  label = "Select option:",
                  choices = c("Both Cells", "K562", "HEPG2", "Similar RBPs")
                ),

                # Search bar + buttons below the dropdown
                div(
                  style = "display: flex; align-items: center; margin-top: 15px;",

                  selectizeInput(
                    inputId = "expr_search",
                    label = NULL,
                    choices = NULL,
                    multiple = FALSE,
                    options = list(placeholder = "RBP"),
                    width = "400px"   # ⬅️ try 400px or "100%" if you want full row
                  ),

                  actionButton(
                    "search_btn",
                    tagList(fa("search"), " Search"),
                    class = "btn btn-primary",
                    style = "margin-left: 10px; border-radius: 20px;"
                  ),

                  actionButton(
                    "reset_btn",
                    tagList(fa("redo"), " Reset"),
                    class = "btn btn-secondary",
                    style = "margin-left: 10px; border-radius: 20px;"
                  )
                )
              ),

              # Conditional: show file upload only in Discovery Mode
              conditionalPanel(
                condition = "input.expr_mode == 'Discovery Mode'",
                hr(),
                tags$p("Alternatively, upload your own table with differential splicing values:"),
                tags$p(
                  "CHARM currently accepts tables where the first column are gene names, ",
                  "and the second column are differential splicing values (dPSI), preferably obtained after processing with betAS.",
                  "This table must be in .txt format and must not contain a header."
                ),

                fileInput(
                  inputId = "user_file",
                  label = "Upload your file:",
                  accept = c(".txt")
                ),

                # Placeholder for inline warning
                uiOutput("file_warning")
              )
            )
          ),

          # Right content area
          column(
            width = 9,
            tags$h3("Splicing Data Results"),
            tags$p("This is where filtered results will appear based on the selected option and search input or uploaded file.")
          )
        )
      )
    ),
    # Binding tab
    tabPanel(
      tagList(fa("link", fill = "black", height = "1em"), " Binding"),
      fluidPage(

        # Main layout: sidebar (left) + results (right)
        fluidRow(
          # Left sidebar for mode selection + search
          column(
            width = 3,
            wellPanel(
              # Choose mode: Explore or Discovery
              radioButtons(
                inputId = "expr_mode",
                label = "Select mode:",
                choices = c("Explore Mode", "Discovery Mode"),
                selected = "Explore Mode"
              ),

              # Conditional: show dropdown only in Explore Mode
              conditionalPanel(
                condition = "input.expr_mode == 'Explore Mode'",
                selectInput(
                  inputId = "expr_dataset",
                  label = "Select option:",
                  choices = c("Both Cells", "K562", "HEPG2", "Similar RBPs")
                ),

                # Search bar + buttons below the dropdown
                div(
                  style = "display: flex; align-items: center; margin-top: 15px;",

                  selectizeInput(
                    inputId = "expr_search",
                    label = NULL,
                    choices = NULL,
                    multiple = FALSE,
                    options = list(placeholder = "RBP"),
                    width = "400px"   # ⬅️ try 400px or "100%" if you want full row
                  ),

                  actionButton(
                    "search_btn",
                    tagList(fa("search"), " Search"),
                    class = "btn btn-primary",
                    style = "margin-left: 10px; border-radius: 20px;"
                  ),

                  actionButton(
                    "reset_btn",
                    tagList(fa("redo"), " Reset"),
                    class = "btn btn-secondary",
                    style = "margin-left: 10px; border-radius: 20px;"
                  )
                )
              ),

              # Conditional: show file upload only in Discovery Mode
              conditionalPanel(
                condition = "input.expr_mode == 'Discovery Mode'",
                hr(),
                tags$p("Alternatively, upload your own table with differential splicing values. Charm's integrated tool, eCLIPSE, will process the data:"),
                tags$p(
                  "CHARM currently accepts tables where the first column are gene names, ",
                  "and the second column are differential splicing values (dPSI), preferably obtained after processing with betAS.",
                  "This table must be in .txt format and must not contain a header. WARNING: THIS TAKES A LONG TIME. Go have yourself a coffee."
                ),

                fileInput(
                  inputId = "user_file",
                  label = "Upload your file:",
                  accept = c(".txt")
                ),

                # Placeholder for inline warning
                uiOutput("file_warning")
              )
            )
          ),

          # Right content area
          column(
            width = 9,
            tags$h3("Binding Data Results"),
            tags$p("This is where filtered results will appear based on the selected option and search input or uploaded file.")
          )
        )
      )
    ),

    # Network tab
    tabPanel(
      tagList(fa("project-diagram", fill = "black", height = "1em"), " Network"),
      fluidPage(

        # Main layout: sidebar (left) + results (right)
        fluidRow(
          # Left sidebar for mode selection + search
          column(
            width = 3,
            wellPanel(
              # Choose mode: Explore or Discovery
              radioButtons(
                inputId = "expr_mode",
                label = "Select mode:",
                choices = c("Explore Mode", "Discovery Mode"),
                selected = "Explore Mode"
              ),

              # Conditional: show dropdown only in Explore Mode
              conditionalPanel(
                condition = "input.expr_mode == 'Explore Mode'",
                selectInput(
                  inputId = "expr_dataset",
                  label = "Select option:",
                  choices = c("Both Cells", "K562", "HEPG2", "Similar RBPs")
                ),

                # Search bar + buttons below the dropdown
                div(
                  style = "display: flex; align-items: center; margin-top: 15px;",

                  selectizeInput(
                    inputId = "expr_search",
                    label = NULL,
                    choices = NULL,
                    multiple = FALSE,
                    options = list(placeholder = "RBP"),
                    width = "400px"   # ⬅️ try 400px or "100%" if you want full row
                  ),

                  actionButton(
                    "search_btn",
                    tagList(fa("search"), " Search"),
                    class = "btn btn-primary",
                    style = "margin-left: 10px; border-radius: 20px;"
                  ),

                  actionButton(
                    "reset_btn",
                    tagList(fa("redo"), " Reset"),
                    class = "btn btn-secondary",
                    style = "margin-left: 10px; border-radius: 20px;"
                  )
                )
              ),

              # Conditional: show file upload only in Discovery Mode
              conditionalPanel(
                condition = "input.expr_mode == 'Discovery Mode'",
                hr(),
                tags$p("Alternatively, upload your own table with expression values:"),
                tags$p(
                  "CHARM currently accepts tables where the first column are gene names, ",
                  "and the second column are differential expression values (t-statistics). ",
                  "This table must be in .txt format and must not contain a header."
                ),

                fileInput(
                  inputId = "user_file",
                  label = "Upload your file:",
                  accept = c(".txt")
                ),

                # Placeholder for inline warning
                uiOutput("file_warning")
              )
            )
          ),

          # Right content area
          column(
            width = 9,
            tags$h3("Network Results"),
            tags$p("This is where network related results will appear based on the selected option and search input.")
          )
        )
      )
    )),

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

  # reactive storage
  display_table <- reactiveVal(NULL)
  rbp_current <- reactiveVal(NULL)

  # ---- Search bar population ----
  observe({
    req(Charm.object)
    updateSelectizeInput(session, "expr_search",
                         choices = unique(Charm.object$Experiment),
                         server = TRUE)
  })

  # ---- File upload validation (unchanged) ----
  observeEvent(input$user_file, {
    req(input$user_file)
    file_path <- input$user_file$datapath
    df <- tryCatch(
      read.table(file_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE),
      error = function(e) NULL
    )

    upload_ok <- TRUE
    error_msg <- NULL
    if (is.null(df)) {
      upload_ok <- FALSE
      error_msg <- "Could not read the file. Make sure it is tab-delimited."
    } else if (ncol(df) != 2) {
      upload_ok <- FALSE
      error_msg <- "File must have exactly 2 columns."
    } else if (!is.character(df[[1]])) {
      upload_ok <- FALSE
      error_msg <- "First column must be character."
    } else if (!is.numeric(df[[2]])) {
      upload_ok <- FALSE
      error_msg <- "Second column must be numeric."
    }

    output$file_warning <- renderUI({
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
        scale_color_manual(values = c("None"="#CCCCCC","RBP"="#A10702","Selected"="#008057")) +
        theme_bw() +
        labs(title=paste(rbp,"KD"), x="Log Fold-Change", y="B-statistic") +
        theme(legend.position="none", plot.title=element_text(hjust=0.5))
    }
  }

  # ---- Main: search button triggers all plots ----
  observeEvent(input$search_btn, {
    req(input$expr_search)
    rbp_sel <- input$expr_search
    rbp_current(rbp_sel)

    session$sendCustomMessage("toggleCursor", TRUE)  # start wait cursor

    # Violin plot
    output$expr_violin <- renderPlot({ violinplotter(Charm.object, rbp_sel) })
    output$expr_note <- renderUI({
      div(style="margin-top:10px;font-size:16px;font-weight:bold;color:red;",
          "Red dots correspond to the controls from paired gene silencing experiment.")
    })

    # shRNA effect
    output$shrna_plot <- renderPlot({ plot_shRNA_effect(sh_effect_vector, rbp_sel) })
    output$shrna_warning <- renderUI({
      stats <- sh_effect_vector[rbp_sel,,drop=FALSE]
      if (is.null(stats) || nrow(stats)==0) return(NULL)
      logFC <- stats$logFC
      pval <- stats$P.Value
      if (pval > 0.05 || logFC > -0.5) {
        div(style="margin-top:10px;font-size:16px;font-weight:bold;color:red;",
            "WARNING: The efficiency of this knockdown is uncertain. Proceed with caution.")
      } else NULL
    })

    # Volcano computation
    result <- plot_rbp_volcano(Charm.object, rbp_sel)
    tbl <- as.data.frame(result$top_table)
    if (!"gene" %in% colnames(tbl)) tbl$gene <- rownames(tbl)
    tbl <- tbl[, c("gene", setdiff(colnames(tbl),"gene"))]
    if (!"highlight" %in% colnames(tbl)) tbl$highlight <- ifelse(tbl$gene==rbp_sel,"RBP","None")
    else tbl$highlight <- ifelse(tbl$gene==rbp_sel,"RBP","None")
    display_table(tbl)

    # render table
    output$volcano_table <- renderDT({
      datatable(
        display_table(),
        filter = "none",       # no per-column filters
        selection = "single",
        rownames = FALSE,
        options = list(
          pageLength = 10,
          scrollX = TRUE,
          searching = TRUE     # global search enabled
        )
      )
    })
    # render volcano
    output$volcano_plot <- renderPlotly({
      ggplotly(make_volcano_plot(display_table(), rbp_sel), tooltip="text", source="volcano") %>%
        event_register("plotly_click")
    })

    session$sendCustomMessage("toggleCursor", FALSE)  # end wait cursor
  })


  # ---- React to clicking a point in the volcano ----
  observeEvent(event_data("plotly_click", source = "volcano"), {
    ed <- event_data("plotly_click", source = "volcano")
    tbl <- display_table()   # get current table
    if (!is.null(ed) && nrow(tbl) > 0) {
      # Find closest gene to click
      clicked_gene <- tbl$gene[which.min(abs(tbl$logFC - ed$x) + abs(tbl$B - ed$y))]

      # Reorder table: clicked gene on top
      tbl <- rbind(tbl[tbl$gene == clicked_gene, ], tbl[tbl$gene != clicked_gene, ])
      display_table(tbl)  # update reactiveVal

      # Highlight selected gene in the table for the plot
      tmp_tbl <- tbl
      tmp_tbl$highlight <- ifelse(tmp_tbl$gene == clicked_gene, "Selected", tmp_tbl$highlight)

      # Re-render volcano plot
      plt <- ggplot(tmp_tbl, aes(x = logFC, y = B,
                                 text = paste0("Gene: ", gene,
                                               "<br>logFC: ", round(logFC, 2),
                                               "<br>B: ", round(B, 2),
                                               "<br>P: ", signif(P.Value, 3)))) +
        geom_point(aes(color = highlight), alpha = 0.7) +
        scale_color_manual(values = c(
          "None" = "#CCCCCC",
          "Other" = "#23586C",
          "RBP" = "#A10702",
          "Selected" = "#008057"
        )) +
        theme_bw() +
        labs(title = paste("Volcano plot:", rbp_current()),
             x = "Log Fold-Change", y = "B-statistic") +
        theme(legend.position = "none",
              plot.title = element_text(hjust = 0.5))

      output$volcano_plot <- renderPlotly({
        ggplotly(plt, tooltip = "text", source = "volcano")
      })
    }
  })

  # ---- React to selecting a row in the table ----
  observeEvent(input$volcano_table_rows_selected, {
    sel_row <- input$volcano_table_rows_selected
    if (is.null(sel_row)) return()
    tbl <- display_table()
    sel_gene <- tbl$gene[sel_row]
    tbl$highlight <- ifelse(tbl$gene==sel_gene,"Selected", ifelse(tbl$gene==rbp_current(),"RBP","None"))
    display_table(tbl)

    # re-render volcano
    output$volcano_plot <- renderPlotly({
      ggplotly(make_volcano_plot(display_table(), rbp_current()), tooltip="text", source="volcano") %>%
        event_register("plotly_click")
    })
  })

} # end server


shinyApp(ui, server)

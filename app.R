library(shiny)
library(shinythemes)
library(fontawesome)

ui <- fluidPage(
  theme = shinytheme("darkly"),

  # Tabset with Home as default
  tabsetPanel(
    type = "tabs",
    id = "main_tabs",

    # Home tab
    tabPanel(
      "Home",
      fluidPage(
        # Title / Logo row
        fluidRow(
          column(
            4,
            div(
              style = "display: flex; align-items: center; justify-content: center; margin: 30px 0;",
              div(style = "margin-right: 15px;", fa("gem", fill = "white", height = "4em")),
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
              style = "background-color: #f8f9fa; color: black;
                       border: 1px solid #2c3e50; border-radius: 20px;
                       padding: 40px; text-align: center; margin: 15px;",
              div(style = "margin-bottom: 10px;", fa("globe", fill = "black", height = "3em")),
              tags$h3("Explore", style = "font-weight: bold;"),
              tags$p("Investigate how a known RBP affects expression, splicing, and binding, based on ENCODE's gene silencing series and eCLIP data.")
            )
          ),
          column(
            6,
            div(
              style = "background-color: #f8f9fa; color: black;
                       border: 1px solid #2c3e50; border-radius: 20px;
                       padding: 40px; text-align: center; margin: 15px;",
              div(style = "margin-bottom: 10px;", fa("map", fill = "black", height = "3em")),
              tags$h3("Discover", style = "font-weight: bold;"),
              tags$p("Users can input their own expression, splicing, or binding data to discover which RBPs are more likely to be altered on your biological system.")
            )
          )
        ),

        # GitHub link at bottom-right
        tags$div(
          style = "
            position: fixed;
            bottom: 20px;
            right: 20px;
            z-index: 1000;
            background-color: #f8f9fa;
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
    ),

    # Expression tab
    tabPanel(
      "Expression",
      fluidPage(
        tags$h2("Expression Data"),
        tags$p("This is where users can explore expression data.")
      )
    ),

    # Splicing tab
    tabPanel(
      "Splicing",
      fluidPage(
        tags$h2("Splicing Data"),
        tags$p("This is where users can explore alternative splicing data.")
      )
    ),

    # Binding tab
    tabPanel(
      "Binding",
      fluidPage(
        tags$h2("Binding Data"),
        tags$p("This is where users can explore RBP binding data.")
      )
    ),

    # Network tab
    tabPanel(
      "Network",
      fluidPage(
        tags$h2("Network"),
        tags$p("This is where users can explore RBP interaction networks.")
      )
    )
  )
)

server <- function(input, output, session) {}

shinyApp(ui, server)

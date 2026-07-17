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
library(qs2)
library(ggpubr)
library(png)
library(grid)
library(base64enc)
logo_b64 <- base64enc::dataURI(file = "www/Charm_logo.png", mime = "image/png")

source("helper_functions.R")

# --- Data directory -----------------------------------------------------
# Points at the qs2-serialized, dtype-optimized data store produced by
# optimize_data.R. Falls back to the legacy QS_Files/ store (old `qs`
# format) if the optimized store hasn't been generated yet, so the app
# never fails to start -- NOTE: that fallback path requires the archived
# `qs` package to still be installed, since qs2 cannot read old .qs files.
# Once data/QS2_Files/ has been validated in production, delete
# data/QS_Files/ and this fallback can be removed.
DATA_DIR <- if (dir.exists("data/QS2_Files")) {
  "data/QS2_Files"
} else {
  "data/QS_Files"
}
DATA_EXT <- if (DATA_DIR == "data/QS2_Files") ".qs2" else ".qs"

qload <- function(name) {
  path <- file.path(DATA_DIR, paste0(name, DATA_EXT))
  if (DATA_EXT == ".qs2") qs2::qs_read(path) else qs::qread(path)
}

# --- Lazy, single-instance dataset-family caches ---------------------------
# Charm.object (3 cell-line variants: Both/K562/HEPG2) and the Binding tables
# (6 variants: event type x cell line) are the two largest pieces of the data
# store. Eagerly loading every variant at startup (the old behavior) pulled
# everything into memory regardless of which single variant a given session
# was actually looking at -- ~6.7GB compressed on disk was ballooning to
# ~35GB in R's memory. Since every session runs in its own dedicated
# container (confirmed: no pooled/shared R process across users), it's safe
# to keep this cache at app/global scope -- there is exactly one session per
# process, so evicting the previous variant on a dataset switch can never
# pull data out from under a different concurrent user.
#
# Trade-off: switching cell lines/event types now costs a real (if brief)
# qs2 decompression instead of being instant. And a handful of features that
# genuinely compare *across* cell lines in one view (the "Similar RBPs"
# comparison-with-uploaded-file panels, and the RBP similarity network
# builder) will still materialize more than one variant at once for the
# duration of that specific action -- that's inherent to what those features
# do, not a gap in this cache.
charm_cache <- new.env(parent = emptyenv())

get_charm_object <- function(dataset) {
  key <- switch(dataset %||% "Both", "K562" = "K562", "HEPG2" = "HEPG2", "Both")
  if (is.null(charm_cache$key) || !identical(charm_cache$key, key)) {
    file <- switch(key, "K562" = "Charm.object_K562", "HEPG2" = "Charm.object_HEPG2", "Charm.object")
    charm_cache$obj <- qload(file)
    charm_cache$key <- key
    gc(FALSE)  # release the evicted variant's memory right away, not at the next scheduled GC
  }
  charm_cache$obj
}

binding_cache <- new.env(parent = emptyenv())

get_binding_object <- function(eventtype, dataset) {
  is_ir <- grepl("Intron Retention", eventtype %||% "", ignore.case = TRUE)
  cell  <- switch(dataset %||% "Both Cells", "K562" = "K562", "HEPG2" = "HEPG2", "Both Cells")
  key   <- paste(if (is_ir) "IR" else "ES", cell)
  if (is.null(binding_cache$key) || !identical(binding_cache$key, key)) {
    file <- if (is_ir) {
      switch(cell, "K562" = "Binding_IR_K562", "HEPG2" = "Binding_IR_HEPG2", "Binding_IR_both")
    } else {
      switch(cell, "K562" = "Binding_K562", "HEPG2" = "Binding_HEPG2", "Binding_both")
    }
    binding_cache$obj <- qload(file)
    binding_cache$key <- key
    gc(FALSE)
  }
  binding_cache$obj
}

`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a

# Load small, always-needed objects once when the app starts (all under
# ~90MB combined -- not worth swapping).
sh_effect_vector       <- qload("shRNA_Efficiency")
sh_effect_vector_K562  <- qload("shRNA_Efficiency_K562")
sh_effect_vector_HEPG2 <- qload("shRNA_Efficiency_HEPG2")

# Warm the default ("Both Cells", Exon Skipping) variants so the app's
# default view loads with no extra delay, and capture the RBP/pathway name
# lists that several static UI dropdowns need regardless of which dataset
# variant ends up active later. Name sets are identical across all
# Charm.object / Binding_* variants (confirmed via direct inspection during
# the qs2 migration), so these lists stay valid even after the cache above
# evicts the object they were read from.
.default_charm <- get_charm_object("Both")
ALL_RBP_NAMES <- names(.default_charm)
HALLMARK_CHOICES <- sort(unique(.default_charm[["RBFOX2"]]$GSEA$pathway))

.default_binding_es <- get_binding_object("Exon Skipping", "Both Cells")
ALL_BINDING_RBP_NAMES_ES <- names(.default_binding_es)
# Intron Retention isn't the default event type, so this read gets evicted
# again the moment a session actually requests ES data -- that's fine, it's
# a one-time startup cost just to capture the name list.
ALL_BINDING_RBP_NAMES_IR <- names(qload("Binding_IR_both"))
rm(.default_charm, .default_binding_es)

#Similarity stuff
similar_expression_all    <- qload("RBPs.t_All")
similar_expression_K562   <- qload("RBPs.t_K562")
similar_expression_HEPG2  <- qload("RBPs.t_HEPG2")
similar_gsea_all          <- qload("RBPs.gsea_All")
similar_gsea_K562         <- qload("RBPs.gsea_K562")
similar_gsea_HEPG2        <- qload("RBPs.gsea_HEPG2")

#For eCLIPSE
eCLIPSE_Intron_full  <- qload("eCLIPSE_INTRONRET")
eCLIPSE_BigExon_full <- qload("eCLIPSE_BIGEXON")

#SearchBarPopulation
GenesBoth   <- qload("AvailableGenes_both")
GenesK562   <- qload("AvailableGenes_K562")
GenesHEPG2  <- qload("AvailableGenes_HEPG2")
EventsBoth  <- qload("AvailableEvents_Both")
EventsK562  <- qload("AvailableEvents_K562")
EventsHEPG2 <- qload("AvailableEvents_HEPG2")

#Binding similarity stuff — loaded lazily on first Network tab use (see server)
# similar_binding_* objects are loaded via the binding_sim_cache reactive below


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
      Shiny.addCustomMessageHandler('clickButton', function(msg) {
        setTimeout(function() {
          var btn = document.getElementById(msg.id);
          if (btn) btn.click();
        }, msg.delay || 300);
      });
      Shiny.addCustomMessageHandler('showStatus', function(msg) {
        var bar = document.getElementById('charm-status-bar');
        if (!bar) return;
        if (msg.show) {
          bar.innerText = msg.text || 'Loading\u2026';
          bar.style.display = 'block';
        } else {
          bar.style.display = 'none';
        }
      });
    ")),
    
    # Persistent status bar — shown during long computations
    tags$div(
      id = "charm-status-bar",
      style = "display:none; position:fixed; top:0; left:0; right:0; z-index:9999;
               background:#2c3e50; color:white; font-weight:bold; font-size:14px;
               padding:10px 20px; text-align:center;
               box-shadow:0 2px 8px rgba(0,0,0,0.3);"
    )
  ),
  
  # Tabset with 5 tabs
  tabsetPanel(
    type = "tabs",
    id = "main_tabs",
    
    # Home tab
    tabPanel(
      tagList(fa("home", fill = "black", height = "1em"), " Home"),
      value = "Home",
      fluidPage(
        fluidRow(
          column(
            4,
            div(
              style = "display: flex; align-items: center; justify-content: center; margin: 30px 0;",
              div(style = "margin-right: 15px;", 
                  tags$img(src = logo_b64, height = "200px")),  # <-- logo_b64, not a file path
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
            4,
            tags$div(
              style = "cursor:pointer;",
              onclick = "Shiny.setInputValue('home_go_explore', Math.random());",
              tags$div(
                style = "background-color:#EAEBEB; color:black; border:2px solid #2c3e50; border-radius:20px; padding:40px; text-align:center; margin:15px;",
                div(style = "margin-bottom:10px;", fa("globe", fill = "black", height = "3em")),
                tags$h3("Explore", style = "font-weight:bold;"),
                tags$p("Investigate how a known RBP affects expression, splicing, and binding, based on ENCODE's gene silencing series and eCLIP data.")
              )
            )
          ),
          column(
            4,
            tags$div(
              style = "cursor:pointer;",
              onclick = "Shiny.setInputValue('home_go_discover', Math.random());",
              tags$div(
                style = "background-color:#EAEBEB; color:black; border:2px solid #2c3e50; border-radius:20px; padding:40px; text-align:center; margin:15px;",
                div(style = "margin-bottom:10px;", fa("map", fill = "black", height = "3em")),
                tags$h3("Discover", style = "font-weight:bold;"),
                tags$p("Upload your own expression, splicing, or binding data to discover which RBPs are most likely altered in your biological system.")
              )
            )
          ),
          column(
            4,
            tags$div(
              style = "cursor:pointer;",
              onclick = "Shiny.setInputValue('home_go_network', Math.random());",
              tags$div(
                style = "background-color:#EAEBEB; color:black; border:2px solid #2c3e50; border-radius:20px; padding:40px; text-align:center; margin:15px;",
                div(style = "margin-bottom:10px;", fa("project-diagram", fill = "black", height = "3em")),
                tags$h3("Network", style = "font-weight:bold;"),
                tags$p("Explore multi-omic RBP similarity via an integrated network combining expression, splicing, and binding layers.")
              )
            )
          )
        )
      )
    ),
    
    # Expression tab
    tabPanel(
      tagList(fa("dna", fill = "black", height = "1em"), " Expression"),
      value = "Expression",
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
                choices = c("Explore", "Discovery"),
                selected = "Explore"
              ),
              
              # --- Explore Mode ---
              conditionalPanel(
                condition = "input.expr_mode == 'Explore'",
                
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
                condition = "input.expr_mode == 'Discovery'",
                hr(),
                tags$p("Alternatively, upload your own table with differential expression values. Data must have a header, and the first column must be the HGNC gene symbol, and the second column the t-stats."),
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
            tags$h3("Expression Results"),
            uiOutput("expr_nav_bar"),
            
            # --- Explore Mode (default mode) ---
            conditionalPanel(
              condition = "input.expr_mode == 'Explore' && input.expr_dataset != 'Similar RBPs'",
              
              # --- Top panel: shows gene OR hallmark plot, whichever was searched most recently ---
              # (server-controlled visibility, see expr_top_panel_ui in server)
              uiOutput("expr_top_panel_ui"),
              
              # --- Default RBP plots ---
              fluidRow(
                column(width = 6,
                  plotOutput("expr_violin", height = "400px"),
                  uiOutput("dl_ui_expr_violin")
                ),
                column(
                  width = 6,
                  plotOutput("shrna_plot"),
                  uiOutput("dl_ui_shrna"),
                  uiOutput("shrna_warning")
                )
              ),
              fluidRow(
                column(width = 6,
                  plotlyOutput("volcano_plot", height = "450px"),
                  uiOutput("dl_ui_volcano")
                ),
                column(width = 6,
                  DTOutput("volcano_table"),
                  uiOutput("dl_ui_volcano_table")
                )
              ),
              fluidRow(
                column(width = 6,
                  plotOutput("gsea_plot", height = "600px"),
                  uiOutput("dl_ui_gsea")
                ),
                column(width = 6,
                  DTOutput("geneset_table"),
                  uiOutput("dl_ui_geneset_table")
                )
              )
            ),
            
            # --- Similar RBPs (Explore Mode, special layout) ---
            conditionalPanel(
              condition = "input.expr_mode == 'Explore' && input.expr_dataset == 'Similar RBPs'",
              uiOutput("similar_expr_plots")
            ),
            
            # --- Discovery Mode ---
            conditionalPanel(
              condition = "input.expr_mode == 'Discovery'",
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
      value = "Splicing",
      fluidPage(
        fluidRow(
          column(
            width = 3,
            wellPanel(
              radioButtons("splice_mode", "Select mode:", choices = c("Explore", "Discovery"), selected = "Explore"),
              
              conditionalPanel(
                condition = "input.splice_mode == 'Explore'",
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
                        ),
                        actionButton(
                          "splice_event_reset_btn",
                          tagList(fa("redo"), " Reset"),
                          class = "btn btn-secondary",
                          style = "border-radius: 20px;"
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
                condition = "input.splice_mode == 'Discovery'",
                hr(),
                tags$p(
                  list(
                    "Alternatively, upload your own table with differential splicing values. Table must be directly obtained from ",
                    tags$a(
                      href   = "https://compbio.imm.medicina.ulisboa.pt/app/betAS",
                      "betAS",
                      target = "_blank",
                      style  = "color: #007bff; text-decoration: none;"
                    ),
                    "."
                  )
                ),
                fileInput("user_file_splice", "Upload your file:", accept = c(".txt")),
                uiOutput("file_warning_splice"),
                uiOutput("user_file_options_splice"),
                
                # ---- eCLIPSE binding map section (appears after successful upload) ----
                conditionalPanel(
                  condition = "input.splice_mode == 'Discovery'",
                  uiOutput("eclipse_raw_options_ui")
                )
              )
            )
          ),
          #Right Side, where plots will appear
          # Right Side, where plots will appear
          column(
            width = 9,
            tags$h3("Splicing Data Results"),
            uiOutput("splice_nav_bar"),
            
            # --- Explore Mode (default mode, NOT Similar RBPs) ---
            conditionalPanel(
              condition = "input.splice_mode == 'Explore' && input.splice_dataset != 'Similar RBPs'",
              
              # --- Top panel: Event ID plot, shown when searched (server-controlled) ---
              uiOutput("splice_top_panel_ui"),
              
              # --- Default Explore-Mode Plots ---
              fluidRow(
                column(width = 6,
                  plotOutput("violin_splice_plot", height = "400px"),
                  uiOutput("dl_ui_splice_violin")
                ),
                column(width = 6,
                  plotOutput("plot_shrna_effect", height = "400px"),
                  uiOutput("dl_ui_shrna_splice")
                )
              ),
              fluidRow(
                column(width = 6,
                  plotlyOutput("plot_splice_volcano", height = "450px"),
                  uiOutput("dl_ui_splice_volcano")
                ),
                column(width = 6,
                  DTOutput("splice_volcano_table"),
                  uiOutput("dl_ui_splice_volcano_table")
                )
              )
            ),
            
            # --- Explore Mode: Similar RBPs ---
            conditionalPanel(
              condition = "input.splice_mode == 'Explore' && input.splice_dataset == 'Similar RBPs'",
              uiOutput("similar_splice_plots")
            ),
            
            # --- Discovery Mode ---
            conditionalPanel(
              condition = "input.splice_mode == 'Discovery'",
              shinycssloaders::withSpinner(
                plotOutput("user_file_initial_plot_splice", height = "420px"),
                type = 6
              ),
              br(),
              uiOutput("similar_splice_plots_file"),
              br(),
              uiOutput("eclipse_raw_plot_ui")   # <-- new
            )
          )
        )
      )
    ),
    
    # Binding tab
    tabPanel(
      tagList(fa("link", fill = "black", height = "1em"), " Binding"),
      value = "Binding",
      fluidPage(
        fluidRow(
          column(
            width = 3,
            wellPanel(
              radioButtons("binding_eventtype", "Event Type:",
                           choices = c("Exon Skipping", "Intron Retention"),
                           selected = "Exon Skipping"),
              
              radioButtons("binding_mode", "Select mode:",
                           choices = c("Explore", "Discovery"),
                           selected = "Explore"),
              
              # -------------------------
              # EXPLORE MODE
              # -------------------------
              conditionalPanel(
                condition = "input.binding_mode == 'Explore'",
                
                selectInput("binding_dataset", "Select option:",
                            choices = c("Both Cells", "K562", "HEPG2", "Similar Profiles")),
                
                div(
                  style = "display: flex; align-items: center; margin-top: 15px;",
                  selectizeInput("binding_search", NULL,
                                 choices = NULL, multiple = FALSE,
                                 options = list(placeholder = "RBP"),
                                 width = "400px"),
                  actionButton("search_btn",
                               tagList(fa("search"), " Search"),
                               class = "btn btn-primary",
                               style = "margin-left: 10px; border-radius: 20px;"),
                  actionButton("reset_btn",
                               tagList(fa("redo"), " Reset"),
                               class = "btn btn-secondary",
                               style = "margin-left: 10px; border-radius: 20px;")
                ),
                
                # Show only when not "Similar RBPs"
                conditionalPanel(
                  condition = "input.binding_dataset != 'Similar Profiles'",
                  
                  
                  hr(),
                  tags$h5("eCLIPSE Binding Map"),
                  tags$p("Investigate where a chosen RBP binds relative to splicing events regulated by the shRNA knockdown above. Select one target for a full binding map, or multiple targets for a heatmap overview.",
                         style = "font-size: 12px; color: #666; margin-top: -5px;"),
                  
                  uiOutput("binding_target_ui"),
                  uiOutput("binding_target_warning"),
                  uiOutput("binding_dpsi_ui"),
                  
                  selectInput("binding_metric", "Metric:",
                              choices = c("FDR", "EffectSize"),
                              selected = "FDR")
                ),
                
                # Show only when "Similar Profiles" is selected
                conditionalPanel(
                  condition = "input.binding_dataset == 'Similar Profiles'",
                  
                  hr(),
                  tags$h5("Similar Profiles"),
                  tags$p(
                    "Find RBP\u2013Target binding profiles that are most similar to a chosen query pair, in terms of increased or decreased events.",
                    style = "font-size: 12px; color: #666; margin-top: -5px;"
                  ),
                  
                  uiOutput("binding_explore_sim_target_ui"),
                  
                  radioButtons(
                    "binding_explore_sim_direction",
                    "Direction:",
                    choices  = c("Increased" = "inc", "Decreased" = "dec"),
                    selected = "inc",
                    inline   = TRUE
                  ),
                  
                  numericInput(
                    "binding_explore_sim_topN",
                    "Show top N most similar profiles (optional):",
                    value = NA, min = 1
                  ),
                  helpText("Tip: selecting a number here takes precedence over the two below."),
                  numericInput(
                    "binding_explore_sim_n_close",
                    "Show top N closest profiles (optional):",
                    value = NA, min = 1
                  ),
                  numericInput(
                    "binding_explore_sim_n_far",
                    "Show top N most dissimilar profiles (optional):",
                    value = NA, min = 1
                  ),
                  
                  div(
                    style = "display: flex; align-items: center; margin-top: 15px;",
                    actionButton(
                      "binding_explore_sim_plot_btn",
                      tagList(fa("chart-line"), " Plot"),
                      class = "btn btn-primary",
                      style = "margin-right: 10px; border-radius: 20px;"
                    ),
                    actionButton(
                      "binding_explore_sim_reset_btn",
                      tagList(fa("redo"), " Reset"),
                      class = "btn btn-secondary",
                      style = "border-radius: 20px;"
                    )
                  )
                )
              ),
              
              # -------------------------
              # DISCOVERY MODE
              # -------------------------
              conditionalPanel(
                condition = "input.binding_mode == 'Discovery'",
                hr(),
                tags$p(
                  list(
                    "Alternatively, upload your own table with differential splicing values. Table must be directly obtained from ",
                    tags$a(
                      href   = "https://compbio.imm.medicina.ulisboa.pt/app/betAS",
                      "betAS",
                      target = "_blank",
                      style  = "color: #007bff; text-decoration: none;"
                    ),
                    "."
                  )
                ),
                fileInput("user_file_binding", "Upload your file:", accept = c(".txt")),
                uiOutput("file_warning_binding"),
                uiOutput("binding_discovery_options"),
                uiOutput("binding_similar_profiles_options")
              )
            )
          ),
          
          # =========================================================
          # RIGHT SIDE: RESULTS + PLOT AREA
          # =========================================================
          column(
            width = 9,
            
            tags$h3("Binding Data Results"),
            uiOutput("binding_nav_bar"),
            
            div(
              tagList(
                fa("clock", fill = "#555", height = "0.9em"),
                " Selecting many targets (or \"All\") may take a moment — larger selections require computing binding profiles for each target."
              ),
              style = "border: 1px solid #c8cbcf;
        background-color: #f8f9fa;
        padding: 7px 10px;
        border-radius: 6px;
        font-size: 12px;
        color: #555;
        margin-bottom: 12px;"
            ),
            
            # Explore Mode plot
            conditionalPanel(
              condition = "input.binding_mode == 'Explore' && input.binding_dataset != 'Similar Profiles'",
              uiOutput("binding_plot_ui")
            ),
            
            # Explore Mode — Similar Profiles plot
            conditionalPanel(
              condition = "input.binding_mode == 'Explore' && input.binding_dataset == 'Similar Profiles'",
              uiOutput("binding_explore_similar_plot_ui")
            ),
            
            # Discovery Mode plot
            conditionalPanel(
              condition = "input.binding_mode == 'Discovery'",
              uiOutput("binding_discovery_plot_ui"),
              uiOutput("binding_similar_profiles_plot_ui")
            )
            
          )
        )
      )
    ),
    
    # Network tab
    tabPanel(
      tagList(fa("project-diagram", fill = "black", height = "1em"), " Network"),
      value = "Network",
      fluidPage(
        fluidRow(
          # ── Left sidebar ──────────────────────────────────────────────────
          column(
            width = 3,
            wellPanel(
              # Cell line
              radioButtons(
                "network_cellline",
                "Cell line:",
                choices  = c("Both", "K562", "HEPG2"),
                selected = "Both"
              ),

              hr(),

              # ── Plot type ──
              tags$h5("Plot type:", style = "font-weight: bold;"),
              radioButtons(
                "network_plot_type",
                label    = NULL,
                choices  = c("MDS", "PCA", "t-SNE", "Dendrogram"),
                selected = "MDS"
              ),

              hr(),

              # Data layers
              tags$h5("Data layers:", style = "font-weight: bold;"),
              checkboxGroupInput(
                "network_layers",
                label    = NULL,
                choices  = c("Expression", "Splicing", "Binding"),
                selected = "Binding"
              ),

              # Slow-layer warning
              conditionalPanel(
                condition = "input.network_layers.indexOf('Expression') !== -1 ||
                             input.network_layers.indexOf('Splicing') !== -1",
                tags$div(
                  style = "background:#fff3cd; border:1px solid #ffc107;
                           border-radius:6px; padding:8px 10px; margin-top:6px;
                           font-size:11px; color:#856404;",
                  tags$b("\u26a0 Note:"),
                  " Expression and Splicing layers build large feature matrices
                    and may take 10\u201330 seconds to compute depending on the
                    number of RBPs. Please be patient after clicking Plot."
                )
              ),

              # ── Binding sub-options ──
              conditionalPanel(
                condition = "input.network_layers.indexOf('Binding') !== -1",
                tags$h6("Event type:", style = "margin-top: 10px; font-weight: bold;"),
                checkboxGroupInput(
                  "network_binding_type",
                  label    = NULL,
                  choices  = c("Intron Retention", "Exon Skipping"),
                  selected = c("Intron Retention", "Exon Skipping")
                ),
                tags$h6("Direction:", style = "margin-top: 8px; font-weight: bold;"),
                checkboxGroupInput(
                  "network_binding_dir",
                  label    = NULL,
                  choices  = c("Increased", "Decreased"),
                  selected = c("Increased", "Decreased")
                )
              ),

              hr(),

              # ── Layer weights ──
              tags$h5("Layer weights:", style = "font-weight: bold;"),
              tags$p(
                "Weights are automatically normalised to sum to 1.",
                style = "font-size: 11px; color: #888; margin-top: -5px;"
              ),

              conditionalPanel(
                condition = "input.network_layers.indexOf('Expression') !== -1",
                sliderInput("network_weight_expr",   "Expression weight:",
                            min = 0, max = 1, value = 0.33, step = 0.05)
              ),
              conditionalPanel(
                condition = "input.network_layers.indexOf('Splicing') !== -1",
                sliderInput("network_weight_splice", "Splicing weight:",
                            min = 0, max = 1, value = 0.33, step = 0.05)
              ),
              conditionalPanel(
                condition = "input.network_layers.indexOf('Binding') !== -1",
                sliderInput("network_weight_binding","Binding weight:",
                            min = 0, max = 1, value = 0.34, step = 0.05)
              ),

              hr(),

              # ── RBP highlight search ──
              tags$h5("Highlight RBPs:", style = "font-weight: bold;"),
              tags$p("Selected RBPs will be highlighted in red on the plot.",
                     style = "font-size: 11px; color: #888; margin-top: -5px;"),
              selectizeInput(
                "network_highlight_rbps",
                label   = NULL,
                choices = sort(ALL_RBP_NAMES),
                multiple = TRUE,
                options  = list(placeholder = "Type to search RBPs...")
              ),

              hr(),

              actionButton(
                "network_plot_btn",
                tagList(fa("chart-line"), " Plot"),
                class = "btn btn-primary",
                style = "border-radius: 20px; width: 100%;"
              )
            )
          ),

          # ── Right plot area ───────────────────────────────────────────────
          column(
            width = 9,
            tags$h3("RBP Similarity Network"),
            tags$p(
              "Multi-omic ordination of RBP similarity. Each point is an RBP;
               proximity reflects similarity across the selected data layers
               (Expression: t-statistics across shared genes; Splicing: dPSI
               across shared events; Binding: profile shape across regulated
               splice sites). Layer weights are normalised to sum to 1.",
              style = "color: #555; font-size: 13px; margin-bottom: 15px;"
            ),
            shinycssloaders::withSpinner(
              plotOutput("network_mds_plot", height = "600px"),
              type = 6
            ),
            uiOutput("dl_ui_network")
          )
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
  
  # ── Downloadable plot stores ──────────────────────────────────────────────
  # Each reactiveVal holds the last ggplot rendered for that output,
  # so downloadHandler can call ggsave() on it at any time.
  dl_store <- list(
    expr_violin    = reactiveVal(NULL),
    shrna          = reactiveVal(NULL),
    volcano        = reactiveVal(NULL),
    gsea           = reactiveVal(NULL),
    expr_gene      = reactiveVal(NULL),
    expr_hallmark  = reactiveVal(NULL),
    splice_violin  = reactiveVal(NULL),
    shrna_splice   = reactiveVal(NULL),
    splice_volcano = reactiveVal(NULL),
    binding        = reactiveVal(NULL),
    network        = reactiveVal(NULL)
  )
  
  # ── Status bar helpers ────────────────────────────────────────────────────
  show_status <- function(msg = "Loading\u2026") {
    session$sendCustomMessage("showStatus", list(show = TRUE, text = msg))
    session$sendCustomMessage("toggleCursor", TRUE)
  }
  hide_status <- function() {
    session$sendCustomMessage("showStatus", list(show = FALSE))
    session$sendCustomMessage("toggleCursor", FALSE)
  }
  
  # ── Download handlers ─────────────────────────────────────────────────────
  # Helper: build a plot downloadHandler from a dl_store key
  make_dl_handler <- function(store_key, filename_prefix,
                              width = 10, height = 7, dpi = 300) {
    downloadHandler(
      filename = function() {
        rbp <- if (!is.null(input$expr_search_rbp) && nchar(input$expr_search_rbp) > 0)
                 input$expr_search_rbp
               else if (!is.null(input$splice_search) && nchar(input$splice_search) > 0)
                 input$splice_search
               else if (!is.null(input$binding_search) && nchar(input$binding_search) > 0)
                 input$binding_search
               else "plot"
        paste0(filename_prefix, "_", rbp, "_", Sys.Date(), ".png")
      },
      content = function(file) {
        p <- dl_store[[store_key]]()
        req(p)
        ggplot2::ggsave(file, plot = p, device = "png",
                        width = width, height = height, dpi = dpi, bg = "white")
      }
    )
  }
  
  # Helper: build a table downloadHandler from a reactiveVal holding a data.frame
  make_table_dl_handler <- function(data_fn, filename_prefix) {
    downloadHandler(
      filename = function() paste0(filename_prefix, "_", Sys.Date(), ".csv"),
      content  = function(file) {
        df <- data_fn()
        req(!is.null(df) && nrow(df) > 0)
        write.csv(df, file, row.names = FALSE)
      }
    )
  }
  
  # Helper: render a download button only once the corresponding store is non-NULL
  make_dl_ui <- function(store_key, dl_id, label = "\u2193 Download plot") {
    output[[paste0("dl_ui_", store_key)]] <- renderUI({
      req(!is.null(dl_store[[store_key]]()))
      downloadButton(dl_id, label,
                     class = "btn btn-sm btn-default",
                     style = "margin-top:5px;")
    })
  }
  
  # ── Plot download handlers ─────────────────────────────────────────────────
  output$dl_expr_violin   <- make_dl_handler("expr_violin",   "expression_violin",   width = 8,  height = 6)
  output$dl_shrna         <- make_dl_handler("shrna",         "shrna_efficiency",    width = 8,  height = 6)
  output$dl_volcano       <- make_dl_handler("volcano",       "expression_volcano",  width = 8,  height = 7)
  output$dl_gsea          <- make_dl_handler("gsea",          "expression_gsea",     width = 10, height = 8)
  output$dl_expr_gene     <- make_dl_handler("expr_gene",     "expression_gene",     width = 8,  height = 7)
  output$dl_expr_hallmark <- make_dl_handler("expr_hallmark", "expression_hallmark", width = 8,  height = 7)
  output$dl_splice_violin <- make_dl_handler("splice_violin", "splicing_violin",     width = 8,  height = 6)
  output$dl_shrna_splice  <- make_dl_handler("shrna_splice",  "shrna_efficiency_splice", width = 8, height = 6)
  output$dl_splice_volcano<- make_dl_handler("splice_volcano","splicing_volcano",    width = 8,  height = 7)
  output$dl_binding       <- make_dl_handler("binding",       "binding_map",         width = 12, height = 9)
  output$dl_network       <- make_dl_handler("network",       "network",             width = 12, height = 9)
  
  # ── Conditional download button UIs (appear only after plot is rendered) ──
  make_dl_ui("expr_violin",   "dl_expr_violin")
  make_dl_ui("shrna",         "dl_shrna")
  make_dl_ui("volcano",       "dl_volcano")
  make_dl_ui("gsea",          "dl_gsea")
  make_dl_ui("expr_gene",     "dl_expr_gene")
  make_dl_ui("expr_hallmark", "dl_expr_hallmark")
  make_dl_ui("splice_violin", "dl_splice_violin")
  make_dl_ui("shrna_splice",  "dl_shrna_splice")
  make_dl_ui("splice_volcano","dl_splice_volcano")
  make_dl_ui("binding",       "dl_binding")
  make_dl_ui("network",       "dl_network")
  
  # ── Table download handlers ────────────────────────────────────────────────
  output$dl_volcano_table        <- make_table_dl_handler(display_table,  "expression_volcano_table")
  output$dl_geneset_table        <- make_table_dl_handler(gsea_table_rv,  "expression_gsea_table")
  output$dl_splice_volcano_table <- make_table_dl_handler(
    reactive({ res <- result_splice(); if (!is.null(res)) res$top_table else NULL }),
    "splicing_volcano_table"
  )
  
  # Table download button UIs
  output$dl_ui_volcano_table <- renderUI({
    req(!is.null(display_table()) && nrow(display_table()) > 0)
    downloadButton("dl_volcano_table", "\u2193 Download table",
                   class = "btn btn-sm btn-default", style = "margin-top:5px;")
  })
  output$dl_ui_geneset_table <- renderUI({
    req(!is.null(gsea_table_rv()) && nrow(gsea_table_rv()) > 0)
    downloadButton("dl_geneset_table", "\u2193 Download table",
                   class = "btn btn-sm btn-default", style = "margin-top:5px;")
  })
  output$dl_ui_splice_volcano_table <- renderUI({
    res <- result_splice()
    req(!is.null(res) && !is.null(res$top_table) && nrow(res$top_table) > 0)
    downloadButton("dl_splice_volcano_table", "\u2193 Download table",
                   class = "btn btn-sm btn-default", style = "margin-top:5px;")
  })
  
  # ── Lazy-load all 12 similar_binding_* objects on first Network/Binding use ──
  # Using a reactive means the 12 qs::qread() calls only happen the first time
  # this reactive is evaluated, and the result is cached for all subsequent calls.
  # This shaves the 12 file reads off app startup time.
  binding_sim_cache <- reactive({
    keys <- c(
      "ES_both_inc",  "ES_both_dec",
      "ES_K562_inc",  "ES_K562_dec",
      "ES_HEPG2_inc", "ES_HEPG2_dec",
      "IR_both_inc",  "IR_both_dec",
      "IR_K562_inc",  "IR_K562_dec",
      "IR_HEPG2_inc", "IR_HEPG2_dec"
    )
    objs <- lapply(keys, function(k) {
      qload(paste0("similar_binding_", k))
    })
    names(objs) <- paste0("similar_binding_", keys)
    objs
  })
  observeEvent(input$home_go_explore, {
    updateTabsetPanel(session, "main_tabs", selected = "Expression")
  })
  observeEvent(input$home_go_discover, {
    updateTabsetPanel(session, "main_tabs", selected = "Expression")
    updateRadioButtons(session, "expr_mode", selected = "Discovery")
  })
  observeEvent(input$home_go_network, {
    updateTabsetPanel(session, "main_tabs", selected = "Network")
  })
  
  
  image_to_ggplot <- function(path) {
    img <- png::readPNG(path)
    g <- grid::rasterGrob(img, interpolate = TRUE)
    ggplot() + annotation_custom(g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
      theme_void()
  }
  
  
  
  result_splice <- reactiveVal(NULL)
  
  # ---- Helper functions to pick correct dataset ----
  # Expression tab
  current_charm_expr <- reactive({
    req(input$expr_dataset)
    get_charm_object(input$expr_dataset)
  })

  # Splicing tab
  current_charm_splice <- reactive({
    req(input$splice_dataset)
    get_charm_object(input$splice_dataset)
  })
  
  current_shrna_expr <- reactive({
    req(input$expr_dataset)
    switch(input$expr_dataset,
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
  gsea_table_rv <- reactiveVal(NULL)
  
  # ── Cross-tab RBP navigation ──────────────────────────────────────────────
  # nav_rbp holds the RBP to jump to, plus the source tab so we don't
  # fire the observer back on the originating tab.
  nav_rbp <- reactiveVal(list(rbp = NULL, from = NULL))
  binding_nav_trigger <- reactiveVal(NULL)
  binding_nav_target  <- reactiveVal(NULL)   # stores the pre-chosen target for nav
  
  # When binding_nav_trigger is set, watch for binding_search to update
  # (which happens synchronously when updateSelectizeInput fires), then
  # select the target and click Search. We watch binding_search rather than
  # binding_target because binding_target is NULL until the user interacts
  # with it, so its change never fires when navigating programmatically.
  observeEvent(input$binding_search, {
    trigger_rbp <- binding_nav_trigger()
    req(!is.null(trigger_rbp))
    req(input$binding_search == trigger_rbp)
    
    data <- binding_data()
    req(!is.null(data))
    available_targets <- names(data[[trigger_rbp]])
    
    target_to_use <- if (trigger_rbp %in% available_targets) {
      trigger_rbp
    } else if (length(available_targets) > 0) {
      available_targets[1]
    } else {
      NULL
    }
    
    req(!is.null(target_to_use))
    binding_nav_trigger(NULL)          # clear so this doesn't re-fire
    binding_nav_target(target_to_use)  # binding_target_ui reads this to pre-select
                                        # itself at creation time (see renderUI below);
                                        # results() also reads it as a fallback.

    # Click Search after a delay so results() fires with binding_nav_target set
    session$sendCustomMessage("clickButton", list(id = "search_btn", delay = 500))
  }, ignoreNULL = TRUE, ignoreInit = TRUE)
  
  # Helper: render the navigation bar for a given tab
  make_nav_bar <- function(rbp, current_tab) {
    other_tabs <- setdiff(c("Expression", "Splicing", "Binding"), current_tab)
    links <- lapply(other_tabs, function(tab) {
      btn_id <- paste0("nav_to_", tab, "_from_", current_tab)
      icon_name <- switch(tab,
        "Expression" = "dna",
        "Splicing"   = "scissors",
        "Binding"    = "link"
      )
      actionButton(
        btn_id,
        tagList(fa(icon_name, height = "0.9em"), paste0(" View in ", tab, " \u2192")),
        class = "btn btn-sm btn-outline-secondary",
        style = "margin-right:8px; border-radius:20px;"
      )
    })
    div(
      style = "background:#f0f4f8; border:1px solid #d0d8e4; border-radius:8px;
               padding:8px 14px; margin-bottom:12px; display:flex; align-items:center;",
      tags$span(
        style = "font-weight:bold; margin-right:12px; color:#2c3e50;",
        paste0("Viewing: ", rbp)
      ),
      tagList(links)
    )
  }
  
  # Expression nav bar — shown when an RBP is active
  output$expr_nav_bar <- renderUI({
    rbp <- rbp_current()
    req(!is.null(rbp))
    make_nav_bar(rbp, "Expression")
  })
  
  # Splicing nav bar
  output$splice_nav_bar <- renderUI({
    rbp <- rbp_current()
    req(!is.null(rbp))
    make_nav_bar(rbp, "Splicing")
  })
  
  # Binding nav bar — shown when binding_search is set and a plot has been made
  output$binding_nav_bar <- renderUI({
    rbp <- input$binding_search
    req(!is.null(rbp) && nchar(rbp) > 0 && show_plot())
    make_nav_bar(rbp, "Binding")
  })
  
  # ── Navigation button handlers ─────────────────────────────────────────────
  # Expression → Splicing
  observeEvent(input$nav_to_Splicing_from_Expression, {
    rbp <- rbp_current(); req(!is.null(rbp))
    updateTabsetPanel(session, "main_tabs", selected = "Splicing")
    updateRadioButtons(session, "splice_mode", selected = "Explore")
    updateSelectizeInput(session, "splice_search", selected = rbp)
    nav_rbp(list(rbp = rbp, from = "Expression", to = "Splicing"))
  })
  
  # Expression → Binding
  observeEvent(input$nav_to_Binding_from_Expression, {
    rbp <- rbp_current(); req(!is.null(rbp))
    updateTabsetPanel(session, "main_tabs", selected = "Binding")
    updateRadioButtons(session, "binding_mode", selected = "Explore")
    updateSelectizeInput(session, "binding_search", selected = rbp)
    nav_rbp(list(rbp = rbp, from = "Expression", to = "Binding"))
  })
  
  # Splicing → Expression
  observeEvent(input$nav_to_Expression_from_Splicing, {
    rbp <- rbp_current(); req(!is.null(rbp))
    updateTabsetPanel(session, "main_tabs", selected = "Expression")
    updateRadioButtons(session, "expr_mode", selected = "Explore")
    updateSelectizeInput(session, "expr_search_rbp", selected = rbp)
    nav_rbp(list(rbp = rbp, from = "Splicing", to = "Expression"))
  })
  
  # Splicing → Binding
  observeEvent(input$nav_to_Binding_from_Splicing, {
    rbp <- rbp_current(); req(!is.null(rbp))
    updateTabsetPanel(session, "main_tabs", selected = "Binding")
    updateRadioButtons(session, "binding_mode", selected = "Explore")
    updateSelectizeInput(session, "binding_search", selected = rbp)
    nav_rbp(list(rbp = rbp, from = "Splicing", to = "Binding"))
  })
  
  # Binding → Expression
  observeEvent(input$nav_to_Expression_from_Binding, {
    rbp <- input$binding_search; req(!is.null(rbp) && nchar(rbp) > 0)
    updateTabsetPanel(session, "main_tabs", selected = "Expression")
    updateRadioButtons(session, "expr_mode", selected = "Explore")
    updateSelectizeInput(session, "expr_search_rbp", selected = rbp)
    nav_rbp(list(rbp = rbp, from = "Binding", to = "Expression"))
  })
  
  # Binding → Splicing
  observeEvent(input$nav_to_Splicing_from_Binding, {
    rbp <- input$binding_search; req(!is.null(rbp) && nchar(rbp) > 0)
    updateTabsetPanel(session, "main_tabs", selected = "Splicing")
    updateRadioButtons(session, "splice_mode", selected = "Explore")
    updateSelectizeInput(session, "splice_search", selected = rbp)
    nav_rbp(list(rbp = rbp, from = "Binding", to = "Splicing"))
  })
  
  # ── Auto-fire search when landing via navigation ───────────────────────────
  # Each observer fires only when nav_rbp targets *its* tab, so the originating
  # tab's search is never re-triggered.
  observeEvent(nav_rbp(), {
    nav <- nav_rbp()
    req(!is.null(nav$rbp), !is.null(nav$to))
    
    if (nav$to == "Expression") {
      rbp_sel   <- nav$rbp
      charm_obj <- current_charm_expr()
      shrna_obj <- current_shrna_expr()
      req(rbp_sel %in% names(charm_obj))
      rbp_current(rbp_sel)
      expr_top_panel(NULL)
      output$expr_violin    <- renderPlot(NULL)
      output$shrna_plot     <- renderPlot(NULL)
      output$volcano_plot   <- renderPlotly(NULL)
      output$volcano_table  <- renderDT(NULL)
      output$gsea_plot      <- renderPlot(NULL)
      output$geneset_table  <- renderDT(NULL)
      display_table(NULL)
      output$expr_violin <- renderPlot({
        p <- violinplotter(charm_obj, rbp_sel)
        dl_store$expr_violin(p); p
      })
      output$shrna_plot <- renderPlot({
        p <- plot_shRNA_effect(shrna_obj, rbp_sel)
        dl_store$shrna(p); p
      })
      result <- plot_rbp_volcano(charm_obj, rbp_sel)
      tbl <- as.data.frame(result$top_table)
      if (!"gene" %in% colnames(tbl)) tbl$gene <- rownames(tbl)
      tbl <- tbl[, c("gene", setdiff(colnames(tbl), "gene"))]
      tbl$highlight <- ifelse(tbl$gene == rbp_sel, "RBP", "None")
      display_table(tbl)
      output$volcano_plot <- renderPlotly({
        p <- make_volcano_plot(display_table(), rbp_sel)
        dl_store$volcano(p)
        ggplotly(p, tooltip = "text", source = "volcano") %>% event_register("plotly_click")
      })
      output$volcano_table <- renderDT({
        datatable(display_table(), filter = "none", selection = "single",
                  rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE, searching = TRUE))
      })
      gsea_result <- plot_gsea(charm_obj, rbp_sel, thresh = 0.05)
      gsea_table_rv(gsea_result$geneset_table)
      output$gsea_plot <- renderPlot({
        p <- gsea_result$gsea_plot; dl_store$gsea(p); p
      })
      output$geneset_table <- renderDT({
        datatable(gsea_result$geneset_table, filter = "none", rownames = FALSE,
                  options = list(pageLength = 10, scrollX = TRUE))
      })
    }
    
    if (nav$to == "Splicing") {
      rbp_sel   <- nav$rbp
      charm_obj <- current_charm_splice()
      shrna_obj <- current_shrna_splice()
      req(rbp_sel %in% names(charm_obj))
      rbp_current(rbp_sel)
      splice_top_panel(NULL)
      output$violin_splice_plot  <- renderPlot(NULL)
      output$plot_shrna_effect   <- renderPlot(NULL)
      output$plot_splice_volcano <- renderPlotly(NULL)
      output$splice_volcano_table <- renderDT(NULL)
      output$violin_splice_plot <- renderPlot({
        p <- violin_splice_plot(charm_obj, rbp_sel)
        dl_store$splice_violin(p); p
      })
      output$plot_shrna_effect <- renderPlot({
        p <- plot_shRNA_effect(shrna_obj, rbp_sel)
        dl_store$shrna_splice(p); p
      })
      splice_res <- plot_splice_volcano(charm_obj, rbp_sel)
      req(!is.null(splice_res))
      result_splice(splice_res)
      output$plot_splice_volcano <- renderPlotly({
        req(result_splice())
        tbl <- result_splice()$top_table
        p <- ggplot(tbl, aes(x = dPSI, y = Pdiff, key = Event.ID, color = highlight)) +
          geom_point(alpha = 0.7) +
          scale_color_manual(values = c("None" = "#CCCCCC", "RBP" = "#A10702", "Selected" = "#008057")) +
          theme_bw() + theme(legend.position = "none") +
          labs(title = paste("Volcano Plot:", rbp_sel), x = "dPSI (shRNA - CTRL)", y = "PDiff")
        dl_store$splice_volcano(p)
        ggplotly(p, tooltip = "key", source = "splice_volcano") %>% event_register("plotly_click")
      })
      output$splice_volcano_table <- DT::renderDataTable({
        req(result_splice())
        result_splice()$top_table %>% select(-text)
      }, rownames = FALSE, selection = "single", options = list(pageLength = 10))
    }
    
    if (nav$to == "Binding") {
      rbp <- nav$rbp
      show_status(paste0("Loading binding data for ", rbp, " \u2014 please wait\u2026"))
      binding_nav_trigger(rbp)
    }
    
    nav_rbp(list(rbp = NULL, from = NULL, to = NULL))  # reset to avoid re-firing
  }, ignoreNULL = TRUE, ignoreInit = TRUE)
  
  # Tracks which panel (if any) should appear pinned at the top of the
  # Expression Explore page: "gene", "hallmark", or NULL (none).
  expr_top_panel <- reactiveVal(NULL)
  
  output$expr_top_panel_ui <- renderUI({
    panel <- expr_top_panel()
    if (is.null(panel)) return(NULL)
    if (panel == "gene") {
      tagList(
        fluidRow(column(width = 12, plotOutput("expr_gene_plot", height = "600px"))),
        uiOutput("dl_ui_expr_gene")
      )
    } else if (panel == "hallmark") {
      tagList(
        fluidRow(column(width = 12, plotOutput("expr_hallmark_plot", height = "600px"))),
        uiOutput("dl_ui_expr_hallmark")
      )
    } else {
      NULL
    }
  })
  
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
    
    # 3️⃣ Hallmark Gene Set search bar (always from RBFOX2; precomputed at
    # startup as HALLMARK_CHOICES so this doesn't force-reload the "Both"
    # Charm.object variant every time the dataset selector changes)
    updateSelectizeInput(
      session,
      "expr_search_hallmark",
      choices = HALLMARK_CHOICES,
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
  
  # Tracks whether the Event ID panel should appear pinned at the top of the
  # Splicing Explore page: "event" or NULL (none).
  splice_top_panel <- reactiveVal(NULL)
  
  output$splice_top_panel_ui <- renderUI({
    panel <- splice_top_panel()
    if (is.null(panel)) return(NULL)
    if (panel == "event") {
      fluidRow(column(width = 12, plotOutput("heatmap_splicing_dpsi", height = "600px")))
    } else {
      NULL
    }
  })
  
  # ---- Event ID Search button ----
  observeEvent(input$search_btn_event, {
    req(input$splice_search_event)
    event_id <- input$splice_search_event
    
    charm_obj <- current_charm_splice()
    
    splice_top_panel("event")
    
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
      read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE),
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
        df_prev <- user_expr_df()
        tagList(
          div(style = "color:green; font-weight:bold; margin-top:10px;",
              sprintf("\u2713 Upload successful — %d genes detected.", nrow(df_prev))),
          tags$small(style = "color:#666;",
                     paste("Columns:", paste(colnames(df_prev), collapse = ", "))),
          br(),
          DT::renderDataTable(
            head(df_prev, 6),
            options = list(dom = "t", ordering = FALSE, pageLength = 6),
            rownames = FALSE,
            class    = "compact"
          )
        )
      } else {
        div(style = "color:red; font-weight:bold; margin-top:10px;",
            paste("Upload failed:", error_msg))
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
          choices = ALL_RBP_NAMES,
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
          choices = ALL_RBP_NAMES,
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
      read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE),
      error = function(e) NULL
    )
    
    if (is.null(df)) {
      upload_ok_splice(FALSE)
      user_splice_df(NULL)
      error_msg <- "Could not read the file. Make sure it is tab-delimited."
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
        df_prev <- user_splice_df()
        tagList(
          div(style = "color:green; font-weight:bold; margin-top:10px;",
              sprintf("\u2713 Upload successful — %d events detected.", nrow(df_prev))),
          tags$small(style = "color:#666;",
                     paste("Columns:", paste(colnames(df_prev), collapse = ", "))),
          br(),
          DT::renderDataTable(
            head(df_prev, 6),
            options = list(dom = "t", ordering = FALSE, pageLength = 6),
            rownames = FALSE,
            class    = "compact"
          )
        )
      } else {
        div(style = "color:red; font-weight:bold; margin-top:10px;",
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
        choices = ALL_RBP_NAMES,
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
  
  # ---- eCLIPSE binding map controls (Discovery Mode, shown after upload) ----
  output$eclipse_raw_options_ui <- renderUI({
    if (!upload_ok_splice()) return(NULL)
    
    tagList(
      hr(),
      tags$h5("eCLIPSE Binding Map"),
      tags$p("Visualise where the RBP binds relative to splicing events in your uploaded data.",
             style = "font-size: 12px; color: #666; margin-top: -5px;"),
      
      selectInput(
        "eclipse_raw_event_type",
        "Event type:",
        choices  = c("Exon Skipping", "Intron Retention"),
        selected = "Exon Skipping"
      ),
      
      selectizeInput(
        "eclipse_raw_rbp",
        "RBP to visualise binding for:",
        choices = sort(unique(eCLIPSE_BigExon_full$RBP)),  # or a shared RBP list
        options = list(placeholder = "Select an RBP")
      ),
      
      numericInput(
        "eclipse_raw_psi_thresh",
        "dPSI threshold:",
        value = 0.05, min = 0, max = 1, step = 0.01
      ),
      
      selectInput(
        "eclipse_raw_metric",
        "Metric:",
        choices  = c("FDR", "EffectSize"),
        selected = "FDR"
      ),
      
      div(
        style = "display: flex; align-items: center; margin-top: 15px;",
        actionButton(
          "eclipse_raw_plot_btn",
          tagList(fa("chart-line"), " Plot"),
          class = "btn btn-primary",
          style = "margin-right: 10px; border-radius: 20px;"
        ),
        actionButton(
          "eclipse_raw_reset_btn",
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
    
    # Calculate average dPSI by type (safe)
    avg_vals <- dpsi_table %>%
      group_by(Type) %>%
      summarise(mean_dPSI = mean(dPSI, na.rm = TRUE)) %>%
      pivot_wider(names_from = Type, values_from = mean_dPSI)
    
    # Build subtitle text safely (handle missing types)
    mean_ES <- ifelse("ES" %in% names(avg_vals), sprintf("%.3f", avg_vals$ES), "NA")
    mean_IR <- ifelse("IR" %in% names(avg_vals), sprintf("%.3f", avg_vals$IR), "NA")
    subtitle_text <- paste0("Mean dPSI — ES: ", mean_ES, " | IR: ", mean_IR)
    
    # Violin + jitter plot similar to your function
    ggplot(dpsi_table, aes(x = Type, y = dPSI)) +
      geom_jitter(width = 0.2, alpha = 0.6, size = 1.5) +
      geom_violin(alpha = .7) +
      theme_minimal() +
      ylab("dPSI (shRNA - CTRL)") +
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
    
    # NOTE: this feature compares all 3 cell lines side by side, so it
    # necessarily materializes all 3 Charm.object variants at once for the
    # duration of this action -- that's inherent to what it does, not
    # something the single-instance cache above can avoid.
    datasets <- list(
      "Both Cells" = get_charm_object("Both"),
      "K562" = get_charm_object("K562"),
      "HEPG2" = get_charm_object("HEPG2")
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
  
  # ---- eCLIPSE raw: plot button ----
  observeEvent(input$eclipse_raw_plot_btn, {
    req(upload_ok_splice(), user_splice_df())
    
    event_type   <- input$eclipse_raw_event_type
    psi_thresh   <- input$eclipse_raw_psi_thresh
    metric       <- input$eclipse_raw_metric
    
    # Pick the correct pre-loaded eCLIPSE map
    rnamapfile <- if (grepl("Intron Retention", event_type, ignore.case = TRUE)) {
      eCLIPSE_Intron_full
    } else {
      eCLIPSE_BigExon_full
    }
    
    output$eclipse_raw_plot_ui <- renderUI({
      shinycssloaders::withSpinner(
        plotOutput("eclipse_raw_plot", height = "600px"),
        type = 6
      )
    })
    
    output$eclipse_raw_plot <- renderPlot({
      withProgress(message = "Generating eCLIPSE binding map...", value = 0.1, {
        eCLIPSE_raw_user(
          rnamapfile   = rnamapfile,
          ASfile       = user_splice_df(),
          rnaBP        = "user",          # label for the map title
          event_type   = event_type,
          PSIthreshold = psi_thresh,
          metric       = metric,
          plot         = TRUE,
          title        = "Uploaded data"
        )
      })
    })
  })
  
  # ---- eCLIPSE raw: reset button ----
  observeEvent(input$eclipse_raw_reset_btn, {
    output$eclipse_raw_plot_ui <- renderUI(NULL)
  })
  
  #BINDING (Discovery Mode)
  upload_ok_binding <- reactiveVal(FALSE)
  user_binding_df   <- reactiveVal(NULL)

  observeEvent(input$user_file_binding, {
    req(input$user_file_binding)
    file_path <- input$user_file_binding$datapath

    df <- tryCatch(
      read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE),
      error = function(e) NULL
    )

    error_msg <- NULL

    if (is.null(df)) {
      error_msg <- "Could not read the file. Make sure it is tab-delimited."
    } else if (ncol(df) != 4) {
      error_msg <- "File must have exactly 4 columns (Event.ID, Gene, dPSI, PDiff)."
    } else if (!is.character(df[[1]])) {
      error_msg <- "First column must be character (Event IDs)."
    } else if (!is.numeric(df[[3]])) {
      error_msg <- "Third column must be numeric (dPSI values)."
    } else {
      colnames(df) <- c("Event.ID", "Gene", "dPSI", "PDiff")
    }

    if (is.null(error_msg)) {
      upload_ok_binding(TRUE)
      user_binding_df(df)
    } else {
      upload_ok_binding(FALSE)
      user_binding_df(NULL)
    }

    output$file_warning_binding <- renderUI({
      if (upload_ok_binding()) {
        df_prev <- user_binding_df()
        tagList(
          div(style = "color:green; font-weight:bold; margin-top:10px;",
              sprintf("\u2713 Upload successful — %d events detected.", nrow(df_prev))),
          tags$small(style = "color:#666;",
                     paste("Columns:", paste(colnames(df_prev), collapse = ", "))),
          br(),
          DT::renderDataTable(
            head(df_prev, 6),
            options = list(dom = "t", ordering = FALSE, pageLength = 6),
            rownames = FALSE,
            class    = "compact"
          )
        )
      } else {
        div(style = "color:red; font-weight:bold; margin-top:10px;",
            paste("Upload failed:", error_msg))
      }
    })
  })

  # ---- Binding Discovery Mode: dynamic options shown after successful upload ----
  output$binding_discovery_options <- renderUI({
    if (!upload_ok_binding()) return(NULL)
    tagList(
      hr(),
      tags$h5("eCLIPSE Binding Map"),
      tags$p("Visualise where an RBP (or multiple RBPs) binds relative to splicing events in your uploaded data. Select one target for a full binding map, or multiple targets for a heatmap overview.",
             style = "font-size: 12px; color: #666; margin-top: -5px;"),
      uiOutput("binding_disc_target_ui"),
      selectInput(
        "binding_disc_psi_thresh",
        "\u0394PSI threshold:",
        choices  = c("0.01", "0.05", "0.10", "0.15", "0.20", "0.25", "0.30"),
        selected = "0.05"
      ),
      selectInput(
        "binding_disc_metric",
        "Metric:",
        choices  = c("FDR", "EffectSize"),
        selected = "FDR"
      ),
      div(
        style = "display: flex; align-items: center; margin-top: 15px;",
        actionButton(
          "binding_disc_plot_btn",
          tagList(fa("chart-line"), " Plot"),
          class = "btn btn-primary",
          style = "margin-right: 10px; border-radius: 20px;"
        ),
        actionButton(
          "binding_disc_reset_btn",
          tagList(fa("redo"), " Reset"),
          class = "btn btn-secondary",
          style = "border-radius: 20px;"
        )
      )
    )
  })

  # Target (RBP whose binding profile to visualise) — Discovery Mode only needs this
  output$binding_disc_target_ui <- renderUI({
    req(upload_ok_binding(), input$binding_eventtype)
    target_choices <- if (grepl("Intron Retention", input$binding_eventtype, ignore.case = TRUE)) {
      sort(ALL_BINDING_RBP_NAMES_IR)
    } else {
      sort(ALL_BINDING_RBP_NAMES_ES)
    }
    selectizeInput(
      "binding_disc_target",
      "Target(s) — RBP binding profile to visualise:",
      choices  = c("All", target_choices),
      multiple = TRUE,
      options  = list(placeholder = "Select one or more targets, or \'All\'")
    )
  })

  
  # ---- Binding Discovery Mode: plot button ----
  observeEvent(input$binding_disc_plot_btn, {
    req(upload_ok_binding(), user_binding_df(),
        input$binding_disc_target,
        input$binding_eventtype)

    event_type   <- input$binding_eventtype
    psi_thresh   <- as.numeric(input$binding_disc_psi_thresh)
    metric       <- input$binding_disc_metric
    targets_sel  <- input$binding_disc_target

    rnamapfile <- if (grepl("Intron Retention", event_type, ignore.case = TRUE)) {
      eCLIPSE_Intron_full
    } else {
      eCLIPSE_BigExon_full
    }

    # Resolve "All" to every available target
    all_rbps <- if (grepl("Intron Retention", event_type, ignore.case = TRUE)) {
      ALL_BINDING_RBP_NAMES_IR
    } else {
      ALL_BINDING_RBP_NAMES_ES
    }
    if ("All" %in% targets_sel) targets_sel <- all_rbps

    n_targets <- length(targets_sel)

    if (n_targets == 1) {
      # Single target: full eCLIPSE line-plot
      output$binding_discovery_plot_ui <- renderUI({
        shinycssloaders::withSpinner(
          plotOutput("binding_disc_plot", height = "600px"),
          type = 6
        )
      })

      output$binding_disc_plot <- renderPlot({
        withProgress(message = "Generating eCLIPSE binding map...", value = 0.1, {
          eCLIPSE_raw_user(
            rnamapfile   = rnamapfile,
            ASfile       = user_binding_df(),
            rnaBP        = targets_sel,
            event_type   = event_type,
            PSIthreshold = psi_thresh,
            metric       = metric,
            plot         = TRUE,
            title        = targets_sel
          )
        })
      })

    } else {
      # Multiple targets: heatmap
      height_px <- min(900 + (n_targets - 1) * 10, 2000)

      output$binding_discovery_plot_ui <- renderUI({
        shinycssloaders::withSpinner(
          plotOutput("binding_disc_plot", height = paste0(height_px, "px")),
          type = 6
        )
      })

      output$binding_disc_plot <- renderPlot({
        withProgress(message = "Generating heatmap...", value = 0.1, {
          build_binding_heatmap_user(
            rnamapfile   = rnamapfile,
            ASfile       = user_binding_df(),
            targets      = targets_sel,
            event_type   = event_type,
            PSIthreshold = psi_thresh,
            metric       = metric
          )
        })
      })
    }
  })
  # ---- Binding Discovery Mode: reset button ----
  observeEvent(input$binding_disc_reset_btn, {
    output$binding_discovery_plot_ui <- renderUI(NULL)
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
    
    # Auto-clear previous outputs before rendering new ones —
    # eliminates the need for the "please reset" warning banner
    output$expr_violin    <- renderPlot(NULL)
    output$shrna_plot     <- renderPlot(NULL)
    output$shrna_warning  <- renderUI(NULL)
    output$volcano_plot   <- renderPlotly(NULL)
    output$volcano_table  <- renderDT(NULL)
    output$gsea_plot      <- renderPlot(NULL)
    output$geneset_table  <- renderDT(NULL)
    display_table(NULL)
    
    # Searching the main RBP plot brings it back to the top: hide gene/hallmark panel
    expr_top_panel(NULL)
    
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
    
    show_status("Generating expression plots \u2014 please wait\u2026")
    
    output$expr_violin <- renderPlot({
      p <- violinplotter(charm_obj, rbp_sel)
      dl_store$expr_violin(p)
      p
    })
    
    output$shrna_plot <- renderPlot({
      p <- plot_shRNA_effect(shrna_obj, rbp_sel)
      dl_store$shrna(p)
      p
    })
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
      p <- make_volcano_plot(display_table(), rbp_sel)
      dl_store$volcano(p)
      ggplotly(p, tooltip = "text", source = "volcano") %>%
        event_register("plotly_click")
    })
    
    gsea_result <- plot_gsea(charm_obj, rbp_sel, thresh = 0.05)
    gsea_table_rv(gsea_result$geneset_table)
    output$gsea_plot <- renderPlot({
      p <- gsea_result$gsea_plot
      dl_store$gsea(p)
      hide_status()
      p
    })
    output$geneset_table <- renderDT({
      datatable(gsea_result$geneset_table,
                filter = "none",
                rownames = FALSE,
                options = list(pageLength = 10, scrollX = TRUE))
    })
  })
  # ---- Gene Search ----
  observeEvent(input$search_btn_gene, {
    req(input$expr_search_gene)
    gene <- input$expr_search_gene
    charm_obj <- current_charm_expr()
    show_status("Generating gene plot \u2014 please wait\u2026")
    expr_top_panel("gene")
    
    output$expr_gene_plot <- renderPlot({
      validate(
        need(gene %in% unlist(lapply(charm_obj, function(x) rownames(x$DEGenes))),
             paste("Gene", gene, "not found in any RBP DEGenes tables."))
      )
      p <- plot_gene_logFC_barplot(charm_obj, gene)
      dl_store$expr_gene(p)
      hide_status()
      p
    })
    
  })
  
  # ---- Hallmark Gene Set Search ----
  observeEvent(input$search_btn_hallmark, {
    req(input$expr_search_hallmark)
    geneset <- input$expr_search_hallmark
    charm_obj <- current_charm_expr()
    show_status("Generating hallmark plot \u2014 please wait\u2026")
    expr_top_panel("hallmark")
    
    output$expr_hallmark_plot <- renderPlot({
      validate(
        need(geneset %in% unlist(lapply(charm_obj, function(x) x$GSEA$pathway)),
             paste0("Geneset '", geneset, "' not found in any RBP GSEA tables."))
      )
      p <- plot_hallmark_nes_barplot(charm_obj, geneset)
      dl_store$expr_hallmark(p)
      hide_status()
      p
    })
    
  })
  
  # Reactive values to store what the user selected
  similar_plot_inputs <- reactiveVal(list(
    rbp1 = NULL,
    selected_rbps = NULL,
    mode = NULL
  ))
  
  # ---- Populate Similar RBPs selectize inputs ----
  observe({
    req(ALL_RBP_NAMES)  # make sure the RBP name list exists
    rbp_names <- ALL_RBP_NAMES
    
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
    
    # Auto-clear previous output
    output$similar_expr_plots <- renderUI(NULL)
    
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
    
    # Auto-clear previous outputs before rendering new ones
    output$violin_splice_plot   <- renderPlot(NULL)
    output$plot_shrna_effect    <- renderPlot(NULL)
    output$shrna_warning_splice <- renderUI(NULL)
    output$plot_splice_volcano  <- renderPlotly(NULL)
    output$splice_volcano_table <- renderDT(NULL)
    
    # Searching the main RBP plot brings it back to the top: hide Event ID panel
    splice_top_panel(NULL)
    
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
    
    show_status("Generating splicing plots \u2014 please wait\u2026")
    
    
    
    # ---- Top row: Violin + shRNA knockdown plots ----
    output$violin_splice_plot <- renderPlot({
      p <- violin_splice_plot(charm_obj, rbp_sel)
      dl_store$splice_violin(p)
      hide_status()
      p
    })
    
    output$plot_shrna_effect <- renderPlot({
      p <- plot_shRNA_effect(shrna_obj, rbp_sel)
      dl_store$shrna_splice(p)
      p
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
        labs(title = paste("Volcano Plot:", rbp_sel), x = "dPSI (shRNA - CTRL)", y = "PDiff")
      dl_store$splice_volcano(p)
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
    
    # Auto-clear previous output
    output$similar_splice_plots <- renderUI(NULL)
    
    # Only one mode currently
    mode <- "splice"
    
    # Choose datasets. NOTE: compares all 3 cell lines side by side, so this
    # necessarily materializes all 3 Charm.object variants at once for the
    # duration of this action -- inherent to the feature, not a cache gap.
    datasets <- list(
      "Both Cells" = get_charm_object("Both"),
      "K562"       = get_charm_object("K562"),
      "HEPG2"      = get_charm_object("HEPG2")
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
    gsea_table_rv(NULL)
    dl_store$expr_violin(NULL); dl_store$shrna(NULL)
    dl_store$volcano(NULL);     dl_store$gsea(NULL)
    updateSelectizeInput(session, "expr_search", selected = "")
  })
  observeEvent(input$reset_btn_gene, {
    output$expr_gene_plot <- renderPlot(NULL)
    if (identical(expr_top_panel(), "gene")) expr_top_panel(NULL)
  })
  
  observeEvent(input$reset_btn_hallmark, {
    output$expr_hallmark_plot <- renderPlot(NULL)
    if (identical(expr_top_panel(), "hallmark")) expr_top_panel(NULL)
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
    result_splice(NULL)
    rbp_current(NULL)
    dl_store$splice_violin(NULL); dl_store$shrna_splice(NULL)
    dl_store$splice_volcano(NULL)
    updateSelectizeInput(session, "splice_search", selected = "")
  })
  
  # 3️⃣ Event ID search reset
  # Reset Event heatmap
  observeEvent(input$splice_event_reset_btn, {
    output$heatmap_splicing_dpsi <- renderPlot(NULL)  # clear plot
    updateSelectizeInput(session, "splice_event_search", selected = "")  # clear search
  })
  
  # =========================
  # 1. Binding dataset selector
  # =========================
  binding_data <- reactive({
    req(input$binding_eventtype, input$binding_dataset)

    if (input$binding_eventtype %in% c("Exon Skipping", "Intron Retention")) {
      return(get_binding_object(input$binding_eventtype, input$binding_dataset))
    }

    NULL
  })
  
  
  
  # =========================
  # 2. RBP selector
  # =========================
  observeEvent(binding_data(), {
    data <- binding_data()
    if (!is.null(data)) {
      updateSelectizeInput(
        session, "binding_search",
        choices = names(data),
        server = TRUE
      )
    }
  })
  
  # =========================
  # 3. Target selector
  # =========================
  output$binding_target_ui <- renderUI({
    req(input$binding_search)
    data <- binding_data()
    rbp <- input$binding_search
    all_targets <- sort(names(data[[rbp]]))

    # If we're arriving here via cross-tab navigation, binding_nav_target()
    # already holds the target to pre-select. Bake it into `selected` at
    # creation time instead of updating it afterward via
    # updateSelectizeInput() -- that update message races the client-side
    # creation of this very widget and gets silently dropped if it arrives
    # first, which left binding_target (and the binding_dpsi UI that depends
    # on it) unset and results() never firing.
    nav_tgt <- binding_nav_target()
    preselected <- if (!is.null(nav_tgt) && nav_tgt %in% all_targets) nav_tgt else NULL

    selectizeInput(
      "binding_target",
      "Target(s):",
      choices = c("All", all_targets),
      selected = preselected,
      multiple = TRUE,
      options = list(placeholder = "Select one or more targets, or 'All'")
    )
  })
  
  # Warn when many targets selected — heatmap will be slow / unreadable
  output$binding_target_warning <- renderUI({
    tgts <- input$binding_target
    if (is.null(tgts) || length(tgts) == 0) return(NULL)
    n <- if ("All" %in% tgts) {
      data <- binding_data()
      rbp  <- input$binding_search
      length(names(data[[rbp]]))
    } else {
      length(tgts)
    }
    if (n > 10) {
      div(
        style = "background:#fff3cd; border:1px solid #ffc107; border-radius:6px;
                 padding:7px 10px; margin-top:5px; font-size:11px; color:#856404;",
        tags$b("\u26a0 Large selection:"),
        sprintf(" %d targets selected. Generating the heatmap may take
                 30\u201360 seconds. Consider selecting fewer targets for faster results.", n)
      )
    }
  })
  
  # =========================
  # 4. dPSI selector (corrected)
  # =========================
  
  
  output$binding_dpsi_ui <- renderUI({
    req(input$binding_search, input$binding_target)
    
    data <- binding_data()
    rbp <- input$binding_search
    targets <- input$binding_target
    
    # Cleanly handle "All"
    hide_labels <- FALSE
    
    if ("All" %in% targets) {
      targets <- names(data[[rbp]])
    }
    # Single-target mode
    if (length(targets) == 1) {
      dpsi_names <- names(data[[rbp]][[targets]])
      dpsi_nums  <- suppressWarnings(as.numeric(sub(" .*", "", dpsi_names)))
      selected   <- dpsi_names[which.min(abs(dpsi_nums))]
      
      return(selectInput(
        "binding_dpsi",
        "dPSI:",
        choices  = setNames(dpsi_names, dpsi_names),
        selected = selected
      ))
    }
    
    # Multi-target mode
    options_multi <- c(
      "0.05", "0.1",
      "maximised oddsratinc", "maximised oddsratdec",
      "maximised pvalinc", "maximised pvaldec"
    )
    
    selectInput(
      "binding_dpsi",
      "dPSI:",
      choices = setNames(options_multi, options_multi),
      selected = "0.05"
    )
  })
  
  # =========================
  # 5. Metric mapping
  # =========================
  binding_metric_name <- reactive({
    if (input$binding_metric == "FDR")        return("pval")
    if (input$binding_metric == "EffectSize") return("oddsrat")
  })
  
  
  # =========================
  # 6. RUN HEAVY FUNCTION WHEN SEARCH IS PRESSED (updated with metric fix)
  # =========================
  
  schematic_exon_skipping <- function() {
    list(
      rects = data.frame(
        xmin = c(0, 450, 950),
        xmax = c(50, 550, 1000),
        ymin = -0.8,
        ymax = -0.2,
        fill = c("#1B4F72", "#BA3B46", "#1B4F72")
      ),
      introns = data.frame(
        x = c(50, 550),
        xend = c(450, 950),
        y = -0.5,
        yend = -0.5
      ),
      arrow = data.frame(
        x = 50, xend = 950, y = -0.1, yend = -0.1
      ),max_x = 1000 
    )
  }
  
  
  schematic_intron_retention <- function() {
    list(
      rects = data.frame(
        xmin = c(0, 450),
        xmax = c(50, 500),
        ymin = -0.8,
        ymax = -0.2,
        fill = c("#1B4F72", "#1B4F72")
      ),
      intron_rect = data.frame(
        xmin = 50,
        xmax = 450,
        ymin = -0.7,
        ymax = -0.3
      ),max_x = 500
    )
  }
  
  
  # -------------------------------
  # Helper: pick dPSI per target
  # -------------------------------
  resolve_dpsi_key <- function(dpsi_names, choice){
    if (is.null(dpsi_names) || length(dpsi_names) == 0 || is.null(choice)) return(NULL)
    
    # ---- 1. Exact match (this covers single-target mode, where `choice` IS
    #         one of the real keys, e.g. "0.02 (maximised pvaldec)") ----
    if (choice %in% dpsi_names) return(list(key = choice, fallback = FALSE))
    
    choice_lc <- tolower(choice)
    
    # ---- 2. Plain numeric threshold choice, e.g. "0.05" or "0.1" ----
    # Keys look like "0.05 (oddsratinc)" -- match the number exactly (escaped),
    # followed by a space or open-paren, and make sure it is NOT a "maximised" key.
    if (grepl("^[0-9.]+$", choice)) {
      pattern <- paste0("^", gsub("\\.", "\\\\.", choice), "[ (]")
      hits <- dpsi_names[grepl(pattern, dpsi_names) & !grepl("maximised", tolower(dpsi_names))]
      if (length(hits) > 0) return(list(key = hits[1], fallback = FALSE))
      return(NULL)  # that exact threshold isn't available for this target
    }
    
    # ---- 3. "maximised <metric><direction>" choice, e.g. "maximised pvaldec" ----
    # The metric+direction suffix (oddsratinc/oddsratdec/pvalinc/pvaldec) must
    # match exactly -- the leading number is target-specific and not part of
    # the match (that's the whole point of "maximised").
    if (grepl("^maximised", choice_lc)) {
      suffix <- trimws(sub("^maximised", "", choice_lc))  # e.g. "pvaldec", "oddsratinc"
      # Match keys containing "maximised <suffix>" specifically (parenthesised or not)
      pattern <- paste0("maximised\\s+", suffix, "\\b")
      hits <- dpsi_names[grepl(pattern, tolower(dpsi_names))]
      if (length(hits) > 0) return(list(key = hits[1], fallback = FALSE))
      
      # Fallback: this target has no "maximised <suffix>" key -- use its 0.1 entry instead
      hits_01 <- dpsi_names[grepl("^0\\.1[ (]", dpsi_names) & !grepl("maximised", tolower(dpsi_names))]
      if (length(hits_01) > 0) return(list(key = hits_01[1], fallback = TRUE))
      return(NULL)  # neither the maximised key nor a 0.1 fallback is available
    }
    
    # ---- 4. Unrecognised choice: no safe fallback: 4 maximised options and the two
    #         numeric ones are mutually exclusive categories, so silently picking a
    #         different one would be misleading. ----
    NULL
  }
  
  # -------------------------------
  # build_binding_heatmap (refactored)
  # -------------------------------
  build_binding_heatmap <- function(data, rbp, targets, dpsi_choice, metric, hide_row_labels = FALSE) {
    
    # -----------------------------
    # 1. Helper: pad list to matrix
    # -----------------------------
    pad_and_rbind <- function(lst) {
      if (length(lst) == 0) return(NULL)
      lengths <- sapply(lst, length)
      maxlen <- max(lengths)
      mat <- t(sapply(lst, function(v) {
        if (is.null(v)) return(rep(NA_real_, maxlen))
        if (length(v) < maxlen) return(c(v, rep(NA_real_, maxlen - length(v))))
        v
      }))
      rownames(mat) <- names(lst)
      
      # Drop rows that are entirely NA
      mat <- mat[rowSums(!is.na(mat)) > 0, , drop = FALSE]
      mat
    }
    
    # -----------------------------
    # 2. Build numeric vectors
    # -----------------------------
    inc_list <- list()
    dec_list <- list()
    fallback_targets <- character(0)  # targets that fell back to 0.1 (no maximised key)
    
    # Totals depend only on the shRNA (rbp) + dPSI threshold, not on target,
    # so we grab them once from the first target that has valid data (see the
    # "maximised" handling below for why that threshold isn't always `key`).
    inc_total <- NA
    maint_total <- NA
    dec_total <- NA
    totals_captured <- FALSE
    
    get_value <- function(df) {
      if (is.null(df) || nrow(df) == 0) return(NA)
      if ("Value" %in% names(df)) return(df$Value[1])
      as.numeric(df[1,1])
    }
    
    # For "maximised" dPSI choices, the threshold (and therefore the totals)
    # genuinely differs per target, so there is no single shared total. In
    # that case we always report the totals from the 0.1 threshold instead,
    # which is target-independent and stable.
    is_maximised_choice <- grepl("^maximised", tolower(dpsi_choice))
    
    for (tgt in targets) {
      if (is.null(data[[rbp]][[tgt]])) next
      dpsi_names <- names(data[[rbp]][[tgt]])
      resolved <- resolve_dpsi_key(dpsi_names, dpsi_choice)
      if (is.null(resolved)) next
      key <- resolved$key
      if (isTRUE(resolved$fallback)) fallback_targets <- c(fallback_targets, tgt)
      
      if (metric == "FDR") {
        inc_raw <- tryCatch(data[[rbp]][[tgt]][[key]]$pval$pvalinc, error = function(e) NULL)
        dec_raw <- tryCatch(data[[rbp]][[tgt]][[key]]$pval$pvaldec, error = function(e) NULL)
      } else {
        inc_raw <- tryCatch(data[[rbp]][[tgt]][[key]]$oddsrat$oddsratinc, error = function(e) NULL)
        dec_raw <- tryCatch(data[[rbp]][[tgt]][[key]]$oddsrat$oddsratdec, error = function(e) NULL)
      }
      
      inc_vec <- if (is.data.frame(inc_raw) && "Value" %in% names(inc_raw)) as.numeric(inc_raw$Value) else suppressWarnings(as.numeric(inc_raw))
      dec_vec <- if (is.data.frame(dec_raw) && "Value" %in% names(dec_raw)) as.numeric(dec_raw$Value) else suppressWarnings(as.numeric(dec_raw))
      
      if (!is.null(inc_vec) && !all(is.na(inc_vec))) inc_list[[tgt]] <- inc_vec
      if (!is.null(dec_vec) && !all(is.na(dec_vec))) dec_list[[tgt]] <- dec_vec
      
      # Capture the shRNA-level totals once (same for every target/dPSI key).
      # For "maximised" choices, always pull the totals from the 0.1 key
      # specifically, since the actual maximised threshold is target-specific
      # and there's no single shared total across targets in that case.
      if (!totals_captured) {
        totals_key <- key
        if (is_maximised_choice) {
          totals_resolved <- resolve_dpsi_key(dpsi_names, "0.1")
          if (!is.null(totals_resolved)) totals_key <- totals_resolved$key else totals_key <- NULL
        }
        if (!is.null(totals_key)) {
          this_data_tot <- tryCatch(data[[rbp]][[tgt]][[totals_key]], error = function(e) NULL)
          if (!is.null(this_data_tot)) {
            inc_total_tmp   <- get_value(this_data_tot$raw$TotalIncreasedEvents)
            maint_total_tmp <- get_value(this_data_tot$raw$TotalMaintainedEvents)
            dec_total_tmp   <- get_value(this_data_tot$raw$TotalDecreasedEvents)
            if (!is.na(inc_total_tmp) || !is.na(maint_total_tmp) || !is.na(dec_total_tmp)) {
              inc_total   <- inc_total_tmp
              maint_total <- maint_total_tmp
              dec_total   <- dec_total_tmp
              totals_captured <- TRUE
            }
          }
        }
      }
    }
    
    if (length(inc_list) == 0 && length(dec_list) == 0) {
      return(ggplot() + theme_void() + ggtitle("No heatmap data"))
    }
    
    inc_mat <- pad_and_rbind(inc_list)
    dec_mat <- pad_and_rbind(dec_list)
    
    # Mark fallback targets in row names with an asterisk so it's visible on the
    # heatmap y-axis (no "maximised" key was available; used 0.1 instead)
    mark_fallback <- function(mat) {
      if (is.null(mat)) return(NULL)
      rn <- rownames(mat)
      rownames(mat) <- ifelse(rn %in% fallback_targets, paste0(rn, " *"), rn)
      mat
    }
    inc_mat <- mark_fallback(inc_mat)
    dec_mat <- mark_fallback(dec_mat)
    
    # -----------------------------
    # 3. Determine color limits
    # -----------------------------
    all_vals <- c(as.numeric(inc_mat), as.numeric(dec_mat))
    all_vals <- all_vals[is.finite(all_vals)]
    lim <- if (length(all_vals) == 0) c(-1,1) else { mx <- max(abs(all_vals)); if(mx==0) c(-1,1) else c(-mx,mx) }
    
    fill_label <- if(metric=="FDR") "FDR" else "Chi-Squared\nStatistic"
    
    # -----------------------------
    # 4. Determine event type & schematic
    # -----------------------------
    etype <- NULL
    try({
      if (exists("input", envir = parent.frame())) {
        pf_input <- get("input", envir = parent.frame())
        if (!is.null(pf_input$binding_eventtype)) etype <- pf_input$binding_eventtype
      }
    }, silent = TRUE)
    if (is.null(etype)) etype <- "Exon Skipping"
    
    schem <- if (grepl("Intron Retention", etype, ignore.case = TRUE)) {
      schematic_intron_retention()
    } else {
      schematic_exon_skipping()
    }
    
    # -----------------------------
    # 5. Plot matrix helper
    # -----------------------------
    plot_matrix <- function(mat, title_text, schematic) {
      if (is.null(mat)) return(NULL)
      
      mdf <- reshape2::melt(mat)
      colnames(mdf) <- c("Target","Position","Value")
      mdf$Position <- as.numeric(mdf$Position)
      
      # Preserve matrix row order for y-axis
      mdf$Target <- factor(mdf$Target, levels = rownames(mat))
      
      max_pos <- max(mdf$Position, na.rm = TRUE)
      n_rows <- nrow(mat)
      
      # Dashed vertical lines depending on event type
      vlines_raw <- if (grepl("Intron Retention", etype, ignore.case = TRUE)) {
        c(50, 250, 450)
      } else {
        c(50, 250, 450, 500, 550, 750, 950)
      }
      vlines_scaled <- vlines_raw * (max_pos / schematic$max_x)
      
      # Use geom_segment to restrict vertical lines to actual rows
      vline_df <- data.frame(x = vlines_scaled)
      
      p <- ggplot(mdf, aes(x=Position, y=Target, fill=Value)) +
        geom_tile() +
        scale_fill_gradient2(low="#A6B1E1", mid="#EEEEEE", high="#B07156", limits=lim, oob=scales::squish) +
        labs(x="Genomic position", y="", fill=fill_label, title=title_text) +
        theme_minimal(base_size=14) +
        theme(
          text = element_text(size = 16, face = "bold"),
          axis.text.x = element_text(size = 14, face = "bold"),
          axis.text.y = if(hide_row_labels) element_blank() else element_text(size = 14, face = "bold"),
          axis.title = element_text(size = 16, face = "bold"),
          panel.grid = element_blank()
        ) +
        geom_segment(data=vline_df,
                     aes(x=x, xend=x, y=1, yend=n_rows),
                     linetype="dashed", color="black", linewidth=0.4,
                     inherit.aes=FALSE) +
        scale_y_discrete(expand = expansion(mult=c(0, 0)))+
        coord_cartesian(clip="off")
      
      # Add schematic under heatmap
      if (!is.null(schematic$rects)) {
        rects <- schematic$rects
        rects$xmin <- rects$xmin * (max_pos / schematic$max_x)
        rects$xmax <- rects$xmax * (max_pos / schematic$max_x)
        p <- p + geom_rect(data=rects,
                           aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                           fill=rects$fill, color="black", linewidth=0.3,
                           inherit.aes=FALSE)
      }
      if (!is.null(schematic$intron_rect)) {
        ir <- schematic$intron_rect
        ir$xmin <- ir$xmin * (max_pos / schematic$max_x)
        ir$xmax <- ir$xmax * (max_pos / schematic$max_x)
        p <- p + geom_rect(data=ir,
                           aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                           fill="#E5E7E9", color="black", linewidth=0.3,
                           inherit.aes=FALSE)
      }
      if (!is.null(schematic$introns)) {
        intr <- schematic$introns
        intr$x   <- intr$x   * (max_pos / schematic$max_x)
        intr$xend <- intr$xend * (max_pos / schematic$max_x)
        p <- p + geom_segment(data=intr,
                              aes(x=x, xend=xend, y=y, yend=yend),
                              linewidth=0.8, inherit.aes=FALSE)
      }
      
      return(p)
    }
    
    # -----------------------------
    # 6. Build plots
    # -----------------------------
    heat_inc_plot <- plot_matrix(inc_mat, "Increased Events", schem)
    heat_dec_plot <- plot_matrix(dec_mat, "Decreased Events", schem)
    
    # -----------------------------
    # 6b. Labels row (shRNA-level totals, shared across all targets)
    # -----------------------------
    inc_color   <- "#de425b"
    dec_color   <- "#769fca"
    maint_color <- "black"
    
    labels_plot <- ggplot() +
      annotate("text", x = 1, y = 0,
               label = paste0("Increased: ", inc_total, "\n"),
               color = inc_color, size = 5, fontface = "bold", hjust = -0.75) +
      annotate("text", x = 1, y = 0,
               label = paste0("Maintained: ", maint_total, "\n"),
               color = maint_color, size = 5, fontface = "bold", hjust = 0.5) +
      annotate("text", x = 1, y = 0,
               label = paste0("Decreased: ", dec_total, "\n"),
               color = dec_color, size = 5, fontface = "bold", hjust = 1.75) +
      theme_void() +
      theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0))
    
    # -----------------------------
    # 7. Combine side-by-side
    # -----------------------------
    heat_row <- NULL
    if (!is.null(heat_dec_plot) && !is.null(heat_inc_plot)) {
      heat_row <- ggpubr::ggarrange(heat_dec_plot, heat_inc_plot, ncol=2, widths=c(0.5,0.5))
    } else if (!is.null(heat_dec_plot)) {
      heat_row <- heat_dec_plot
    } else if (!is.null(heat_inc_plot)) {
      heat_row <- heat_inc_plot
    }
    
    if (!is.null(heat_row)) {
      footnote_lines <- character(0)
      if (is_maximised_choice) {
        footnote_lines <- c(footnote_lines,
          "Increased/Maintained/Decreased totals above use the 0.1 threshold (no single shared total exists across targets for 'maximised' choices).")
      }
      if (length(fallback_targets) > 0) {
        footnote_lines <- c(footnote_lines,
          "* No 'maximised' dPSI value available for this target — showing its 0.1 threshold instead.")
      }
      
      if (length(footnote_lines) > 0) {
        footnote <- ggplot() +
          annotate("text", x = 0.5, y = 0.5,
                   label = paste(footnote_lines, collapse = "\n"),
                   size = 4, fontface = "italic", hjust = 0.5) +
          theme_void()
        return(ggpubr::ggarrange(labels_plot, heat_row, footnote, ncol = 1, heights = c(0.05, 0.88, 0.07)))
      }
      return(ggpubr::ggarrange(labels_plot, heat_row, ncol = 1, heights = c(0.06, 0.94)))
    }
    
    ggplot() + theme_void() + ggtitle("No heatmap data")
  }
  
  
  # -------------------------------
  # results reactive
  # -------------------------------
  results <- eventReactive(input$search_btn, {
    withProgress(message = "Generating Plot...", value = 0, {
      incProgress(0.05)
      
      # Use binding_nav_target if set (navigation case), otherwise read from UI
      nav_tgt <- binding_nav_target()
      targets <- if (!is.null(nav_tgt)) {
        binding_nav_target(NULL)   # consume and clear
        nav_tgt
      } else {
        input$binding_target
      }
      
      req(input$binding_search, targets, input$binding_dpsi, input$binding_metric)
      
      data <- binding_data()
      rbp  <- input$binding_search
      
      # NEW: define hide_labels
      hide_labels <- FALSE
      
      if ("All" %in% targets) {
        targets <- names(data[[rbp]])
      }
      
      metric <- input$binding_metric
      if (is.null(metric) || !(metric %in% c("FDR","EffectSize"))) metric <- "FDR"
      
      incProgress(0.2)
      
      if (length(targets) == 1) {
        # Single target: regular eCLIPSE plot
        dpsi_names <- names(data[[rbp]][[targets]])
        dpsi_resolved <- resolve_dpsi_key(dpsi_names, input$binding_dpsi)
        if (is.null(dpsi_resolved)) {
          return(ggplot() +
                   annotate("text", x = 0.5, y = 0.5,
                            label = paste0("dPSI option '", input$binding_dpsi,
                                           "' is not available for ", rbp, " - ", targets, "."),
                            size = 6, hjust = 0.5) +
                   theme_void())
        }
        dpsi_key   <- dpsi_resolved$key
        dpsi_value <- suppressWarnings(as.numeric(sub(" .*", "", dpsi_key)))
        plot_title <- paste(rbp, targets)
        if (isTRUE(dpsi_resolved$fallback)) {
          plot_title <- paste0(plot_title, " [* fell back to 0.1, no maximised value available]")
        }
        p <- eCLIPSE_full(
          bindingvalues_nested = data,
          rnaBP = rbp,
          target = targets,
          dPSI = dpsi_value,
          metric = metric,
          title = plot_title,
          schematic_exon_skipping = schematic_exon_skipping,
          schematic_intron_retention = schematic_intron_retention
        )
        return(p)
      } else {
        # Multiple targets: heatmap
        return(build_binding_heatmap(
          data,
          rbp,
          targets,
          input$binding_dpsi,
          metric,
          hide_row_labels = hide_labels   # FIXED
        ))
      }
      
      incProgress(1)
    })
  })
  
  # ---- Helper reactive: pick the correct binding similarity object list ----
  binding_sim_objects <- reactive({
    req(input$binding_eventtype)
    is_IR  <- grepl("Intron Retention", input$binding_eventtype, ignore.case = TRUE)
    prefix <- if (is_IR) "IR" else "ES"
    cache  <- binding_sim_cache()   # triggers lazy load on first call
    
    list(
      "Both Cells" = list(
        inc = cache[[paste0("similar_binding_", prefix, "_both_inc")]],
        dec = cache[[paste0("similar_binding_", prefix, "_both_dec")]]
      ),
      "K562" = list(
        inc = cache[[paste0("similar_binding_", prefix, "_K562_inc")]],
        dec = cache[[paste0("similar_binding_", prefix, "_K562_dec")]]
      ),
      "HEPG2" = list(
        inc = cache[[paste0("similar_binding_", prefix, "_HEPG2_inc")]],
        dec = cache[[paste0("similar_binding_", prefix, "_HEPG2_dec")]]
      )
    )
  })
  
  # ─────────────────────────────────────────────────────────────────────────────
  # EXPLORE MODE — Similar Profiles (binding_dataset == "Similar Profiles")
  # ─────────────────────────────────────────────────────────────────────────────
  
  # ---- Target selector: targets available for the chosen RBP, current event type ----
  output$binding_explore_sim_target_ui <- renderUI({
    req(input$binding_search, input$binding_dataset == "Similar Profiles")
    
    sim_objs <- binding_sim_objects()
    # Use "Both Cells" / "inc" as the universe of profile IDs (RBP__Target pairs
    # are the same set across inc/dec for a given event type + cell line).
    all_ids <- rownames(sim_objs[["Both Cells"]]$inc$mean_profiles)
    if (is.null(all_ids)) all_ids <- character(0)
    
    rbp <- input$binding_search
    rbp_targets <- sort(unique(sub("^[^_]*__", "", all_ids[startsWith(all_ids, paste0(rbp, "__"))])))
    
    selectizeInput(
      "binding_explore_sim_target",
      "Target:",
      choices  = rbp_targets,
      multiple = FALSE,
      options  = list(placeholder = "Select a target")
    )
  })
  
  # ---- Plot button: render the 3 cell-line panels for the chosen pair/direction ----
  observeEvent(input$binding_explore_sim_plot_btn, {
    req(input$binding_dataset == "Similar Profiles",
        input$binding_search,
        input$binding_explore_sim_target,
        input$binding_eventtype)
    
    rbp      <- input$binding_search
    target   <- input$binding_explore_sim_target
    query_id <- paste0(rbp, "__", target)
    dir_sel  <- input$binding_explore_sim_direction   # "inc" or "dec"
    
    topN    <- input$binding_explore_sim_topN
    n_close <- input$binding_explore_sim_n_close
    n_far   <- input$binding_explore_sim_n_far
    
    # If the user left every filter blank, default to the top 10 closest
    # profiles instead of plotting every profile in the dataset.
    if (is.na(topN) && is.na(n_close) && is.na(n_far)) n_close <- 10
    
    sim_objs      <- binding_sim_objects()
    dataset_names <- c("Both Cells", "K562", "HEPG2")
    
    output$binding_explore_similar_plot_ui <- renderUI({
      tagList(
        tags$h4(paste0("Similar Profiles to ", query_id,
                       " (", if (dir_sel == "inc") "Increased" else "Decreased",
                       " Events)"),
                style = "text-align:center;"),
        lapply(dataset_names, function(ds_name) {
          fluidRow(
            column(width = 12,
                   tags$h5(ds_name, style = "text-align:center;"),
                   shinycssloaders::withSpinner(
                     plotOutput(paste0("binding_explore_sim_plot_", gsub(" ", "", ds_name)),
                                height = "450px"),
                     type = 6
                   )
            )
          )
        }) %>% tagList()
      )
    })
    
    for (ds in dataset_names) {
      local({
        ds_local <- ds
        plot_id  <- paste0("binding_explore_sim_plot_", gsub(" ", "", ds_local))
        
        output[[plot_id]] <- renderPlot({
          sim_obj <- sim_objs[[ds_local]][[dir_sel]]
          res <- binding_profile_correl(
            sim_obj       = sim_obj,
            query_id      = query_id,
            top_n         = if (!is.na(topN))    topN    else NULL,
            n_close       = if (!is.na(n_close)) n_close else NULL,
            n_far         = if (!is.na(n_far))   n_far   else NULL,
            other_ids     = NULL,
            direction     = dir_sel,
            dataset_label = ds_local
          )
          res$heatmap
        })
      })
    }
  })
  
  # ---- Reset button ----
  observeEvent(input$binding_explore_sim_reset_btn, {
    output$binding_explore_similar_plot_ui <- renderUI(NULL)
    updateSelectizeInput(session, "binding_explore_sim_target", selected = "")
    updateNumericInput(session,   "binding_explore_sim_topN",    value = NA)
    updateNumericInput(session,   "binding_explore_sim_n_close", value = NA)
    updateNumericInput(session,   "binding_explore_sim_n_far",   value = NA)
  })
  
  # -------------------------------
  # Render plot
  # -------------------------------
  output$eclipse_plot <- renderPlot({
    req(results())
    p <- results()
    if (inherits(p, "ggplot") || inherits(p, "ggarrange")) {
      dl_store$binding(p)
      hide_status()
      print(p)
    }
  })
  
  # Track plot visibility
  show_plot <- reactiveVal(FALSE)
  
  # When search button is pressed, show plot
  observeEvent(input$search_btn, {
    show_status("Generating binding map \u2014 this may take a moment\u2026")
    show_plot(TRUE)
  })
  
  # When reset button is pressed, hide plot
  observeEvent(input$reset_btn, {
    show_plot(FALSE)
    dl_store$binding(NULL)
    binding_nav_target(NULL)
  })
  output$binding_plot_ui <- renderUI({
    req(show_plot())
    
    targets <- input$binding_target
    data <- binding_data()
    rbp <- input$binding_search
    
    if ("All" %in% targets) {
      n_targets <- length(names(data[[rbp]]))
    } else {
      n_targets <- length(targets)
    }
    
    height_px <- min(900 + (n_targets - 1) * 10, 2000)
    
    tagList(
      plotOutput("eclipse_plot", height = paste0(height_px, "px")),
      uiOutput("dl_ui_binding")
    )
  })
  
  
  output$binding_similar_profiles_options <- renderUI({
    if (!upload_ok_binding()) return(NULL)
    
    # Determine available profile IDs from the inc object (Both Cells, current event type)
    sim_objs <- binding_sim_objects()
    all_profile_ids <- rownames(sim_objs[["Both Cells"]]$inc$mean_profiles)
    if (is.null(all_profile_ids)) all_profile_ids <- character(0)
    
    tagList(
      hr(),
      tags$h5("Similar Profiles"),
      tags$p(
        "Find RBP–Target binding profiles that are most similar to a chosen query profile.",
        style = "font-size: 12px; color: #666; margin-top: -5px;"
      ),
      
      # Query profile selector
      selectizeInput(
        "binding_sim_query_id",
        "Query profile (RBP__Target):",
        choices  = all_profile_ids,
        multiple = FALSE,
        options  = list(placeholder = "e.g. HNRNPC__HsaEX0001234")
      ),
      
      # Compare with specific profiles (optional)
      selectizeInput(
        "binding_sim_compare",
        "Compare with specific profiles (optional):",
        choices  = all_profile_ids,
        multiple = TRUE,
        options  = list(placeholder = "Select one or more profiles")
      ),
      
      numericInput(
        "binding_sim_topN",
        "Show top N most similar profiles (optional):",
        value = NA, min = 1
      ),
      helpText("Tip: selecting a number here takes precedence over the two below."),
      numericInput(
        "binding_sim_n_close",
        "Show top N closest profiles (optional):",
        value = NA, min = 1
      ),
      numericInput(
        "binding_sim_n_far",
        "Show top N most dissimilar profiles (optional):",
        value = NA, min = 1
      ),
      
      div(
        style = "display: flex; align-items: center; margin-top: 15px;",
        actionButton(
          "binding_sim_plot_btn",
          tagList(fa("chart-line"), " Plot"),
          class = "btn btn-primary",
          style = "margin-right: 10px; border-radius: 20px;"
        ),
        actionButton(
          "binding_sim_reset_btn",
          tagList(fa("redo"), " Reset"),
          class = "btn btn-secondary",
          style = "border-radius: 20px;"
        )
      )
    )
  })
  
  
  # ─────────────────────────────────────────────────────────────────────────────
  # PATCH 8 — SERVER: plot button for Similar Profiles (6 plots)
  # ─────────────────────────────────────────────────────────────────────────────
  
  observeEvent(input$binding_sim_plot_btn, {
    req(upload_ok_binding(),
        input$binding_sim_query_id,
        input$binding_eventtype)
    
    query_id    <- input$binding_sim_query_id
    compare_ids <- if (length(input$binding_sim_compare) > 0) input$binding_sim_compare else NULL
    topN        <- input$binding_sim_topN
    n_close     <- input$binding_sim_n_close
    n_far       <- input$binding_sim_n_far
    
    sim_objs <- binding_sim_objects()
    dataset_names <- c("Both Cells", "K562", "HEPG2")
    directions    <- c("inc", "dec")
    
    output$binding_similar_profiles_plot_ui <- renderUI({
      tagList(
        hr(),
        tags$h4("Similar Profiles Results",
                style = "font-weight:bold; color:#A10702; margin-bottom:15px;"),
        
        # Row 1 — Both Cells
        fluidRow(
          column(6,
                 tags$h5("Both Cells — Increased Events", style = "text-align:center;"),
                 shinycssloaders::withSpinner(
                   plotOutput("binding_sim_plot_both_inc", height = "400px"), type = 6)
          ),
          column(6,
                 tags$h5("Both Cells — Decreased Events", style = "text-align:center;"),
                 shinycssloaders::withSpinner(
                   plotOutput("binding_sim_plot_both_dec", height = "400px"), type = 6)
          )
        ),
        
        # Row 2 — K562
        fluidRow(
          column(6,
                 tags$h5("K562 — Increased Events", style = "text-align:center;"),
                 shinycssloaders::withSpinner(
                   plotOutput("binding_sim_plot_K562_inc", height = "400px"), type = 6)
          ),
          column(6,
                 tags$h5("K562 — Decreased Events", style = "text-align:center;"),
                 shinycssloaders::withSpinner(
                   plotOutput("binding_sim_plot_K562_dec", height = "400px"), type = 6)
          )
        ),
        
        # Row 3 — HEPG2
        fluidRow(
          column(6,
                 tags$h5("HEPG2 — Increased Events", style = "text-align:center;"),
                 shinycssloaders::withSpinner(
                   plotOutput("binding_sim_plot_HEPG2_inc", height = "400px"), type = 6)
          ),
          column(6,
                 tags$h5("HEPG2 — Decreased Events", style = "text-align:center;"),
                 shinycssloaders::withSpinner(
                   plotOutput("binding_sim_plot_HEPG2_dec", height = "400px"), type = 6)
          )
        )
      )
    })
    
    # Render each of the 6 plots
    for (ds in dataset_names) {
      for (dir in directions) {
        local({
          ds_local  <- ds
          dir_local <- dir
          plot_id   <- paste0("binding_sim_plot_", gsub(" ", "", ds_local), "_", dir_local)
          
          output[[plot_id]] <- renderPlot({
            sim_obj <- sim_objs[[ds_local]][[dir_local]]
            res <- binding_profile_correl(
              sim_obj       = sim_obj,
              query_id      = query_id,
              top_n         = if (!is.na(topN))    topN    else NULL,
              n_close       = if (!is.na(n_close)) n_close else NULL,
              n_far         = if (!is.na(n_far))   n_far   else NULL,
              other_ids     = compare_ids,
              direction     = dir_local,
              dataset_label = ds_local
            )
            res$heatmap
          })
        })
      }
    }
  })
  
  
  # ─────────────────────────────────────────────────────────────────────────────
  # PATCH 9 — SERVER: reset button for Similar Profiles
  # ─────────────────────────────────────────────────────────────────────────────
  
  observeEvent(input$binding_sim_reset_btn, {
    output$binding_similar_profiles_plot_ui <- renderUI(NULL)
    updateSelectizeInput(session, "binding_sim_query_id", selected = "")
    updateSelectizeInput(session, "binding_sim_compare",  selected = "")
    updateNumericInput(session,   "binding_sim_topN",     value    = NA)
    updateNumericInput(session,   "binding_sim_n_close",  value    = NA)
    updateNumericInput(session,   "binding_sim_n_far",    value    = NA)
  })

  # ─────────────────────────────────────────────────────────────────────────────
  # NETWORK TAB — MDS of RBP binding-profile similarity
  # ─────────────────────────────────────────────────────────────────────────────

  output$network_mds_plot <- renderPlot({
    input$network_plot_btn

    isolate({
      show_status("Building RBP similarity network \u2014 please wait\u2026")
      layers    <- input$network_layers
      cell      <- input$network_cellline
      btypes    <- input$network_binding_type
      dirs      <- input$network_binding_dir
      plot_type <- input$network_plot_type
      hl_rbps   <- input$network_highlight_rbps

      if (is.null(layers) || length(layers) == 0) {
        return(ggplot() +
                 annotate("text", x = 0.5, y = 0.5,
                          label = "Please select at least one data layer.",
                          size = 6, hjust = 0.5) + theme_void())
      }

      withProgress(message = "Building RBP similarity network\u2026", value = 0, {

      # ── Cell line mapping ──────────────────────────────────────────────────
      cell_key_binding <- switch(cell, "Both" = "both", "K562" = "K562", "HEPG2" = "HEPG2")
      charm_obj <- get_charm_object(cell)
      expr_obj  <- switch(cell,
                          "Both"  = similar_expression_all,
                          "K562"  = similar_expression_K562,
                          "HEPG2" = similar_expression_HEPG2)

      # ── Weights ───────────────────────────────────────────────────────────
      raw_weights <- c()
      if ("Expression" %in% layers)
        raw_weights["Expression"] <- max(0, if (is.null(input$network_weight_expr))    1 else input$network_weight_expr)
      if ("Splicing"   %in% layers)
        raw_weights["Splicing"]   <- max(0, if (is.null(input$network_weight_splice))  1 else input$network_weight_splice)
      if ("Binding"    %in% layers)
        raw_weights["Binding"]    <- max(0, if (is.null(input$network_weight_binding)) 1 else input$network_weight_binding)
      total_w <- sum(raw_weights)
      if (total_w == 0) raw_weights[] <- 1 / length(raw_weights)
      weights <- raw_weights / total_w

      # ── Build correlation matrices per layer ───────────────────────────────
      cor_matrices <- list()
      n_layers <- length(layers)
      step_size <- 0.6 / max(n_layers, 1)

      if ("Expression" %in% layers) {
        incProgress(step_size, message = "Computing expression correlations\u2026")
        cm <- network_cor_matrix(type = "expression", rbp_list = expr_obj)
        if (!is.null(cm) && nrow(cm) >= 3) cor_matrices[["Expression"]] <- cm
      }

      if ("Splicing" %in% layers) {
        incProgress(step_size, message = "Computing splicing correlations\u2026")
        cm <- network_cor_matrix(type = "splicing", rbp_list = charm_obj)
        if (!is.null(cm) && nrow(cm) >= 3) cor_matrices[["Splicing"]] <- cm
      }

      if ("Binding" %in% layers) {
        incProgress(step_size, message = "Computing binding correlations\u2026")
        if (!is.null(btypes) && length(btypes) > 0 &&
            !is.null(dirs)   && length(dirs)   > 0) {

          type_keys <- c()
          if ("Exon Skipping"    %in% btypes) type_keys <- c(type_keys, "ES")
          if ("Intron Retention" %in% btypes) type_keys <- c(type_keys, "IR")
          dir_keys <- c()
          if ("Increased" %in% dirs) dir_keys <- c(dir_keys, "inc")
          if ("Decreased" %in% dirs) dir_keys <- c(dir_keys, "dec")

          binding_cor_list <- list()
          cache <- binding_sim_cache()
          for (tk in type_keys) {
            for (dk in dir_keys) {
              obj_name <- paste0("similar_binding_", tk, "_", cell_key_binding, "_", dk)
              obj <- cache[[obj_name]]
              if (is.null(obj) || is.null(obj$mean_profiles)) next
              mat       <- obj$mean_profiles
              rbp_names <- sub("__.*", "", rownames(mat))
              unique_rbps <- unique(rbp_names)
              agg <- do.call(rbind, lapply(unique_rbps, function(r) {
                colMeans(mat[which(rbp_names == r), , drop = FALSE], na.rm = TRUE)
              }))
              rownames(agg) <- unique_rbps
              cm <- network_cor_matrix(type = "binding", rbp_mat = agg)
              if (!is.null(cm) && nrow(cm) >= 3)
                binding_cor_list[[paste0(tk, "_", dk)]] <- cm
            }
          }

          if (length(binding_cor_list) > 0) {
            common_b <- Reduce(intersect, lapply(binding_cor_list, rownames))
            avg_b <- Reduce("+", lapply(binding_cor_list, function(m) m[common_b, common_b])) /
                     length(binding_cor_list)
            cor_matrices[["Binding"]] <- avg_b
          }
        }
      }

      if (length(cor_matrices) == 0) {
        return(ggplot() +
                 annotate("text", x = 0.5, y = 0.5,
                          label = "No data available for the selected combination.",
                          size = 6, hjust = 0.5) + theme_void())
      }

      # ── Common RBPs ───────────────────────────────────────────────────────
      common_rbps <- Reduce(intersect, lapply(cor_matrices, rownames))
      if (length(common_rbps) < 3) {
        return(ggplot() +
                 annotate("text", x = 0.5, y = 0.5,
                          label = paste0("Only ", length(common_rbps),
                                         " RBP(s) shared across all layers. Need \u2265 3."),
                          size = 5, hjust = 0.5) + theme_void())
      }

      # ── Weighted average correlation → dissimilarity ───────────────────────
      active_layers <- names(cor_matrices)
      w <- weights[active_layers]; w <- w / sum(w)

      avg_cor <- Reduce("+", mapply(
        function(cm, wi) cm[common_rbps, common_rbps] * wi,
        cor_matrices[active_layers], w, SIMPLIFY = FALSE
      ))
      diss_mat <- 1 - avg_cor   # dissimilarity in [0, 2]

      # ── Highlight setup ───────────────────────────────────────────────────
      hl_valid  <- intersect(hl_rbps, common_rbps)
      highlight <- common_rbps %in% hl_valid

      weight_str <- paste(
        sapply(active_layers, function(l) sprintf("%s (%.0f%%)", l, w[l] * 100)),
        collapse = " | "
      )
      caption_base <- paste0("Dissimilarity = 1 \u2212 weighted Spearman r  |  n = ",
                             length(common_rbps), " RBPs")

      # ── Helper: base scatter ggplot ────────────────────────────────────────
      make_scatter <- function(df, xlab, ylab, title, subtitle) {
        df$Highlight <- ifelse(highlight, "Highlighted", "Other")
        has_hl       <- length(hl_valid) > 0

        p <- ggplot(df, aes(x = Dim1, y = Dim2)) +
          # Points: colour and size driven by Highlight column
          geom_point(
            aes(color = Highlight, size = Highlight),
            alpha = 0.8
          ) +
          scale_color_manual(values = c("Highlighted" = "#BA3B46", "Other" = "#2c3e50")) +
          scale_size_manual( values = c("Highlighted" = 4,         "Other" = 2.5)) +
          labs(title = title, subtitle = subtitle,
               x = xlab, y = ylab, caption = caption_base) +
          theme_bw() +
          theme(
            plot.title    = element_text(hjust = 0.5, face = "bold", size = 15),
            plot.subtitle = element_text(hjust = 0.5, size = 11, color = "#555555"),
            plot.caption  = element_text(hjust = 0.5, size = 9,  color = "#888888"),
            text          = element_text(size = 12),
            panel.grid    = element_blank(),
            axis.line     = element_line(colour = "black"),
            panel.border  = element_blank(),
            legend.position = if (has_hl) "right" else "none"
          )

        if (has_hl) {
          # Label only the highlighted RBPs, in red
          p <- p + ggrepel::geom_text_repel(
            data         = df[df$Highlight == "Highlighted", ],
            aes(label    = RBP),
            color        = "#BA3B46",
            fontface     = "bold",
            size         = 4,
            max.overlaps = 30
          )
        } else {
          # No highlights — label every RBP in dark
          p <- p + ggrepel::geom_text_repel(
            aes(label    = RBP),
            color        = "#2c3e50",
            fontface     = "bold",
            size         = 3.5,
            max.overlaps = 20
          )
        }
        p
      }

      # ── Plot type dispatch ─────────────────────────────────────────────────

      incProgress(0.2, message = "Running ordination\u2026")

      if (plot_type == "MDS") {
        fit     <- cmdscale(as.dist(diss_mat), k = 2, eig = TRUE)
        eig_pos <- pmax(fit$eig, 0)
        ve      <- round(eig_pos[1:2] / sum(eig_pos) * 100, 1)
        df <- data.frame(RBP = common_rbps,
                         Dim1 = fit$points[,1], Dim2 = fit$points[,2])
        p <- make_scatter(df,
          xlab     = paste0("MDS Dim 1 (", ve[1], "% var.)"),
          ylab     = paste0("MDS Dim 2 (", ve[2], "% var.)"),
          title    = paste0("RBP Similarity — MDS (", cell, ")"),
          subtitle = weight_str)
        dl_store$network(p)
        hide_status()
        return(p)
      }

      if (plot_type == "PCA") {
        pca <- prcomp(avg_cor[common_rbps, common_rbps], scale. = TRUE, center = TRUE)
        ve  <- round(summary(pca)$importance[2, 1:2] * 100, 1)
        df  <- data.frame(RBP  = common_rbps,
                          Dim1 = pca$x[, 1],
                          Dim2 = pca$x[, 2])
        p <- make_scatter(df,
          xlab     = paste0("PC1 (", ve[1], "% var.)"),
          ylab     = paste0("PC2 (", ve[2], "% var.)"),
          title    = paste0("RBP Similarity — PCA (", cell, ")"),
          subtitle = weight_str)
        dl_store$network(p)
        hide_status()
        return(p)
      }

      if (plot_type == "t-SNE") {
        if (!requireNamespace("Rtsne", quietly = TRUE)) {
          return(ggplot() +
                   annotate("text", x = 0.5, y = 0.5, size = 5, hjust = 0.5,
                            label = "Package 'Rtsne' is required for t-SNE.\nInstall it with: install.packages('Rtsne')") +
                   theme_void())
        }
        set.seed(42)
        tsne <- Rtsne::Rtsne(as.dist(diss_mat), is_distance = TRUE,
                              perplexity = min(30, floor((length(common_rbps) - 1) / 3)))
        df   <- data.frame(RBP  = common_rbps,
                           Dim1 = tsne$Y[, 1],
                           Dim2 = tsne$Y[, 2])
        p <- make_scatter(df,
          xlab     = "t-SNE 1",
          ylab     = "t-SNE 2",
          title    = paste0("RBP Similarity — t-SNE (", cell, ")"),
          subtitle = weight_str)
        dl_store$network(p)
        hide_status()
        return(p)
      }

      if (plot_type == "Dendrogram") {
        hc  <- hclust(as.dist(diss_mat), method = "ward.D2")
        dend <- as.dendrogram(hc)

        # Colour highlighted leaves
        if (length(hl_valid) > 0) {
          color_leaf <- function(node) {
            if (is.leaf(node)) {
              lbl <- attr(node, "label")
              if (lbl %in% hl_valid) {
                attr(node, "nodePar") <- list(lab.col = "#BA3B46",
                                              lab.font = 2, pch = NA)
              }
            }
            node
          }
          dend <- dendrapply(dend, color_leaf)
        }

        par(mar = c(12, 4, 4, 2))
        plot(dend, main = paste0("RBP Similarity — Dendrogram (", cell, ")"),
             sub  = weight_str,
             ylab = "Height (1 \u2212 Spearman r)",
             cex.main = 1.3, font.main = 2,
             cex.sub  = 0.9, col.sub  = "#555555",
             cex.lab  = 1.0, cex.axis = 0.9)
        if (length(hl_valid) > 0)
          legend("topright", legend = "Highlighted RBPs",
                 text.col = "#BA3B46", text.font = 2, bty = "n", cex = 0.9)
        hide_status()
        return(invisible(NULL))
      }

      }) # end withProgress
    })
  })

}


# end server


shinyApp(ui, server)

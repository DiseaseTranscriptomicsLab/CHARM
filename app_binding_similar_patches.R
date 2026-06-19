# =============================================================================
#  PATCH GUIDE — Binding Tab "Similar Profiles" in Discovery Mode
#  Apply each block in order. Each block shows WHERE to find the text and
#  what to REPLACE it with (or INSERT after it).
# =============================================================================

# ─────────────────────────────────────────────────────────────────────────────
# PATCH 1 — Load the 12 pre-computed binding similarity .qs files at app start
#           (add after the existing eCLIPSE / similar_gsea loading block)
# ─────────────────────────────────────────────────────────────────────────────
# FIND (last line of the existing similarity loading block):
#   similar_gsea_HEPG2 <- qs::qread("data/QS_Files/RBPs.gsea_HEPG2.qs")
#
# INSERT AFTER:

similar_binding_ES_both_inc  <- qs::qread("data/QS_Files/similar_binding_ES_both_inc.qs")
similar_binding_ES_both_dec  <- qs::qread("data/QS_Files/similar_binding_ES_both_dec.qs")
similar_binding_ES_K562_inc  <- qs::qread("data/QS_Files/similar_binding_ES_K562_inc.qs")
similar_binding_ES_K562_dec  <- qs::qread("data/QS_Files/similar_binding_ES_K562_dec.qs")
similar_binding_ES_HEPG2_inc <- qs::qread("data/QS_Files/similar_binding_ES_HEPG2_inc.qs")
similar_binding_ES_HEPG2_dec <- qs::qread("data/QS_Files/similar_binding_ES_HEPG2_dec.qs")
similar_binding_IR_both_inc  <- qs::qread("data/QS_Files/similar_binding_IR_both_inc.qs")
similar_binding_IR_both_dec  <- qs::qread("data/QS_Files/similar_binding_IR_both_dec.qs")
similar_binding_IR_K562_inc  <- qs::qread("data/QS_Files/similar_binding_IR_K562_inc.qs")
similar_binding_IR_K562_dec  <- qs::qread("data/QS_Files/similar_binding_IR_K562_dec.qs")
similar_binding_IR_HEPG2_inc <- qs::qread("data/QS_Files/similar_binding_IR_HEPG2_inc.qs")
similar_binding_IR_HEPG2_dec <- qs::qread("data/QS_Files/similar_binding_IR_HEPG2_dec.qs")


# ─────────────────────────────────────────────────────────────────────────────
# PATCH 2 — Rename "Similar RBPs" to "Similar Profiles" in Explore Mode
#           selectInput for the Binding tab
# ─────────────────────────────────────────────────────────────────────────────
# FIND (inside the EXPLORE MODE conditionalPanel of the Binding tab sidebar):
#   selectInput("binding_dataset", "Select option:",
#               choices = c("Both Cells", "K562", "HEPG2", "Similar RBPs")),
#
# REPLACE WITH:

selectInput("binding_dataset", "Select option:",
            choices = c("Both Cells", "K562", "HEPG2", "Similar Profiles")),


# ─────────────────────────────────────────────────────────────────────────────
# PATCH 3 — Update the conditionalPanel condition that hides eCLIPSE options
#           when "Similar Profiles" is selected in Explore Mode
# ─────────────────────────────────────────────────────────────────────────────
# FIND:
#   condition = "input.binding_dataset != 'Similar RBPs'",
#
# REPLACE WITH:

condition = "input.binding_dataset != 'Similar Profiles'",


# ─────────────────────────────────────────────────────────────────────────────
# PATCH 4 — Add "Similar Profiles" UI controls in Discovery Mode sidebar
#           (add inside the Discovery Mode conditionalPanel, after the
#            existing uiOutput("binding_discovery_options"))
# ─────────────────────────────────────────────────────────────────────────────
# FIND (end of the Discovery Mode conditionalPanel in binding sidebar):
#   uiOutput("binding_discovery_options")
#
# REPLACE WITH:

uiOutput("binding_discovery_options"),
uiOutput("binding_similar_profiles_options")


# ─────────────────────────────────────────────────────────────────────────────
# PATCH 5 — Add Similar Profiles plot area in the right column, inside
#           the Discovery Mode conditionalPanel
# ─────────────────────────────────────────────────────────────────────────────
# FIND (inside the right column, Discovery Mode conditionalPanel):
#   uiOutput("binding_discovery_plot_ui")
#
# REPLACE WITH:

uiOutput("binding_discovery_plot_ui"),
uiOutput("binding_similar_profiles_plot_ui")


# ─────────────────────────────────────────────────────────────────────────────
# PATCH 6 — SERVER: reactive to pick the 6 correct sim objects based on
#           event type + direction.  Add this near the other binding reactives.
# ─────────────────────────────────────────────────────────────────────────────

# ---- Helper reactive: pick the correct binding similarity object list ----
binding_sim_objects <- reactive({
  req(input$binding_eventtype)
  is_IR <- grepl("Intron Retention", input$binding_eventtype, ignore.case = TRUE)
  prefix <- if (is_IR) "IR" else "ES"

  list(
    "Both Cells" = list(
      inc = get(paste0("similar_binding_", prefix, "_both_inc")),
      dec = get(paste0("similar_binding_", prefix, "_both_dec"))
    ),
    "K562" = list(
      inc = get(paste0("similar_binding_", prefix, "_K562_inc")),
      dec = get(paste0("similar_binding_", prefix, "_K562_dec"))
    ),
    "HEPG2" = list(
      inc = get(paste0("similar_binding_", prefix, "_HEPG2_inc")),
      dec = get(paste0("similar_binding_", prefix, "_HEPG2_dec"))
    )
  )
})


# ─────────────────────────────────────────────────────────────────────────────
# PATCH 7 — SERVER: renderUI for the Similar Profiles sidebar options
#           (shown after a successful file upload + eCLIPSE plot has been run)
# ─────────────────────────────────────────────────────────────────────────────

output$binding_similar_profiles_options <- renderUI({
  if (!upload_ok_binding()) return(NULL)

  # Determine available profile IDs from the inc object (Both Cells, current event type)
  sim_objs <- binding_sim_objects()
  all_profile_ids <- rownames(sim_objs[["Both Cells"]]$inc$cor_matrix)
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

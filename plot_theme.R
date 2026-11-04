# =============================================================================
# plot_theme.R -- ONE visual style for every CHARM plot
# -----------------------------------------------------------------------------
# Every static (ggplot2) and interactive (Highcharts) figure in the app uses the
# font and sizes defined here, so all plots look like they come from one source.
#
#   Font family : "sans" -- the one family name every R graphics device knows
#                 (Helvetica/Arial-like on Mac, Liberation Sans in Docker). A named
#                 family such as "Arial" triggers "font family not found in PostScript
#                 font database" warnings and broken plots on devices that lack it.
#   Plot title  : 16 pt bold
#   Axis titles : 14 pt bold
#   Everything else (tick labels, legends, in-plot labels, subtitles,
#   annotations, tooltips) : 12 pt regular
#
# Style: white background, black axis lines and ticks, no grid lines.
# shiny::renderPlot() draws at 72 dpi, so 1 pt in ggplot2 == 1 px in Highcharts.
# To change the look of the whole app, edit only the constants below.
# =============================================================================

CHARM_FONT       <- "sans"
# NOTE: no quotes inside this string -- highcharter builds a Google-Fonts <link> selector
# from it and a quote character breaks that selector in the browser (chart never draws).
CHARM_FONT_CSS   <- "Helvetica, Arial, sans-serif"
CHARM_SIZE_TITLE <- 16
CHARM_SIZE_AXIS  <- 14
CHARM_SIZE_TEXT  <- 12
# geom_text()/annotate()/geom_text_repel() take 'size' in mm, not points:
CHARM_GEOM_TEXT  <- CHARM_SIZE_TEXT / ggplot2::.pt

# ---- ggplot2 ----------------------------------------------------------------
theme_charm <- function() {
  ggplot2::theme_bw(base_size = CHARM_SIZE_TEXT, base_family = CHARM_FONT) +
    ggplot2::theme(
      text             = ggplot2::element_text(family = CHARM_FONT, size = CHARM_SIZE_TEXT,
                                               face = "plain", colour = "black"),
      plot.title       = ggplot2::element_text(size = CHARM_SIZE_TITLE, face = "bold", hjust = 0.5),
      plot.subtitle    = ggplot2::element_text(size = CHARM_SIZE_TEXT, hjust = 0.5),
      plot.caption     = ggplot2::element_text(size = CHARM_SIZE_TEXT, hjust = 0.5),
      axis.title       = ggplot2::element_text(size = CHARM_SIZE_AXIS, face = "bold"),
      axis.text        = ggplot2::element_text(size = CHARM_SIZE_TEXT, colour = "black"),
      legend.title     = ggplot2::element_text(size = CHARM_SIZE_TEXT, face = "bold"),
      legend.text      = ggplot2::element_text(size = CHARM_SIZE_TEXT),
      strip.text       = ggplot2::element_text(size = CHARM_SIZE_TEXT, face = "bold"),
      strip.background = ggplot2::element_blank(),
      axis.line        = ggplot2::element_line(colour = "black"),
      axis.ticks       = ggplot2::element_line(colour = "black"),
      panel.grid       = ggplot2::element_blank(),
      panel.border     = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      plot.background  = ggplot2::element_rect(fill = "white", colour = NA),
      legend.key       = ggplot2::element_blank()
    )
}

# Same font/sizes, but without axes or panel (heatmap strips, empty-state notes).
theme_charm_void <- function() {
  ggplot2::theme_void(base_size = CHARM_SIZE_TEXT, base_family = CHARM_FONT) +
    ggplot2::theme(
      text          = ggplot2::element_text(family = CHARM_FONT, size = CHARM_SIZE_TEXT,
                                            face = "plain", colour = "black"),
      plot.title    = ggplot2::element_text(size = CHARM_SIZE_TITLE, face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = CHARM_SIZE_TEXT, hjust = 0.5),
      legend.title  = ggplot2::element_text(size = CHARM_SIZE_TEXT, face = "bold"),
      legend.text   = ggplot2::element_text(size = CHARM_SIZE_TEXT)
    )
  # NOTE: axis.text/axis.title are deliberately NOT set here. theme_void() hides axes,
  # and setting them would switch the axis numbers back on. A plot that needs an axis
  # in a void theme must add it explicitly (see eCLIPSE_full).
}

# ---- Highcharts -------------------------------------------------------------
charm_hc_style <- function(hc) {
  txt <- function(size, weight = "normal")
    list(fontFamily = CHARM_FONT_CSS, fontSize = paste0(size, "px"),
         fontWeight = weight, color = "#000000")
  axis_opts <- function(hc_axis_fn) {
    hc <<- hc_axis_fn(hc,
      lineColor = "#000000", lineWidth = 1, tickColor = "#000000", gridLineWidth = 0,
      labels = list(style = txt(CHARM_SIZE_TEXT)),
      title  = list(style = txt(CHARM_SIZE_AXIS, "bold")))
  }
  axis_opts(highcharter::hc_xAxis)
  axis_opts(highcharter::hc_yAxis)
  hc <- highcharter::hc_chart(hc, backgroundColor = "#FFFFFF", style = txt(CHARM_SIZE_TEXT))
  hc <- highcharter::hc_title(hc, style = txt(CHARM_SIZE_TITLE, "bold"), align = "center")
  hc <- highcharter::hc_subtitle(hc, style = txt(CHARM_SIZE_TEXT), align = "center")
  hc <- highcharter::hc_legend(hc, itemStyle = txt(CHARM_SIZE_TEXT))
  hc <- highcharter::hc_tooltip(hc, style = txt(CHARM_SIZE_TEXT))
  highcharter::hc_credits(hc, enabled = FALSE)
}

# Make the ggplot theme the default for every ggplot in the session.
ggplot2::theme_set(theme_charm())

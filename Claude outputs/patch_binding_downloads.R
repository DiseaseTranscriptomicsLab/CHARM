# Binding > Similar Profiles: separate bar-plot / density-plot downloads with figure-friendly sizes.
# Layout: 2x2 grid (bar plots on top, density plots below, each 8 in wide) + pooled density (16 in wide) underneath.
# Run from the CHARM folder (R console or: Rscript patch_binding_downloads.R). Edits app.R and
# helper_functions.R in place (backups: *.bak). Safe to run more than once.

rd  <- function(p) { x <- rawToChar(readBin(p, "raw", file.size(p))); Encoding(x) <- "UTF-8"; x }
wr  <- function(p, x) { file.copy(p, paste0(p, ".bak"), overwrite = TRUE); writeBin(charToRaw(x), p) }
cnt <- function(x, pat) { m <- gregexpr(pat, x, fixed = TRUE)[[1]]; if (m[1] == -1) 0L else length(m) }
sub1 <- function(x, old, new, what) {
  if (cnt(x, old) != 1) stop("Could not find expected code for: ", what, " (found ", cnt(x, old), "x). Nothing was changed in this file.")
  sub(old, gsub("\\", "\\\\", new, fixed = TRUE), x, fixed = TRUE)
}
stopifnot(file.exists("app.R"), file.exists("helper_functions.R"))

# ---------------- helper_functions.R ----------------
h <- rd("helper_functions.R")
if (grepl("bar_plot_only", h, fixed = TRUE)) {
  cat("helper_functions.R: already patched\n")
} else {
  h <- sub1(h,
    '  tbl_cols <- c("Profile", "RBP", "Target", "Correlation", if (has_ref) c("Percentile", "SameTarget"))\n  if (has_ref) {',
    paste0('  tbl_cols <- c("Profile", "RBP", "Target", "Correlation", if (has_ref) c("Percentile", "SameTarget"))\n',
           '  bar_plot_only <- heatmap_plot   # bar plot alone (separate download)\n',
           '  dens_only     <- NULL           # density plot alone (separate download)\n',
           '  if (has_ref) {'), "helper anchor 1")
  h <- sub1(h,
    '    heatmap_plot <- ggpubr::ggarrange(heatmap_plot, dens, ncol = 1, heights = c(1, 0.9))\n  }\n  list(\n    correlation_table = cor_df[, tbl_cols],\n    heatmap           = heatmap_plot\n  )',
    paste0('    bar_plot_only <- heatmap_plot\n    dens_only     <- dens\n',
           '    heatmap_plot <- ggpubr::ggarrange(heatmap_plot, dens, ncol = 1, heights = c(1, 0.9))\n  }\n',
           '  list(\n    correlation_table = cor_df[, tbl_cols],\n    heatmap           = heatmap_plot,\n',
           '    bar_plot          = bar_plot_only,\n    density_plot      = dens_only\n  )'), "helper anchor 2")
  wr("helper_functions.R", h); cat("helper_functions.R: patched\n")
}

# ---------------- app.R ----------------
a <- rd("app.R")
v1 <- '"_", dir_local, "_density"),\n            width = 16, height = 4'
if (grepl(v1, a, fixed = TRUE)) {            # an earlier (Python) version was applied: only the width changes
  a <- sub(v1, '"_", dir_local, "_density"),\n            width = 8, height = 4', a, fixed = TRUE)
  wr("app.R", a); cat("app.R: density width corrected (16 -> 8)\n")
} else if (grepl('dir_local, "_barplot"', a, fixed = TRUE)) {
  cat("app.R: already patched\n")
} else {
  # 1. Expression/Splicing density downloads: full figure width
  a <- sub('paste0(prefix, "_density"), width = 8, height = 4)', 'paste0(prefix, "_density"), width = 16, height = 4)', a, fixed = TRUE)
  # 2. store bar and density separately
  a <- sub1(a, "            dyn_dl_store[[plot_id]] <- res$heatmap\n            res$heatmap",
               "            corr_store_plots(plot_id, res)\n            res$heatmap", "store plots")
  # 3. handlers for the six bar + density panels
  a <- sub1(a,
'          dl_id <- paste0("dl_", plot_id)
          output[[dl_id]] <- make_dyn_dl_handler(
            plot_id,
            paste0("binding_similar_profile_", ds_key, "_", dir_local)
          )',
'          # bar plot and density plot are downloaded as two independent images
          output[[paste0("dl_", plot_id, "_bar")]] <- make_dyn_dl_handler(
            paste0(plot_id, "_bar"),
            paste0("binding_similar_profile_", ds_key, "_", dir_local, "_barplot"),
            width = 8, height = 6
          )
          output[[paste0("dl_", plot_id, "_density")]] <- make_dyn_dl_handler(
            paste0(plot_id, "_density"),
            paste0("binding_similar_profile_", ds_key, "_", dir_local, "_density"),
            width = 8, height = 4
          )', "per-panel handlers")
  # 4. pooled (increased + decreased) density: full width
  a <- sub1(a, 'plot_id, paste0("binding_similar_profile_", ds_key, "_together")\n          )',
               'plot_id, paste0("binding_similar_profile_", ds_key, "_together"),\n            width = 16, height = 4\n          )', "pooled density handler")
  # 5. the six buttons
  pat <- 'downloadButton\\("dl_binding_sim_plot_(both|K562|HEPG2)_(inc|dec)", "[^"]*Download plot",\n\\s*class = "btn btn-sm btn-default"\\)'
  if (length(gregexpr(pat, a, perl = TRUE)[[1]]) != 6) stop("Expected 6 download buttons; nothing was changed in app.R.")
  a <- gsub(pat,
    'tagList(downloadButton("dl_binding_sim_plot_\\1_\\2_bar", "\\\\u2193 Download bar plot", class = "btn btn-sm btn-default", style = "margin-right:6px;"), downloadButton("dl_binding_sim_plot_\\1_\\2_density", "\\\\u2193 Download density plot", class = "btn btn-sm btn-default"))',
    a, perl = TRUE)
  wr("app.R", a); cat("app.R: patched\n")
}
invisible(parse("app.R")); invisible(parse("helper_functions.R")); cat("Both files parse correctly. Restart the app.\n")

# =============================================================================
# explore_binding_targets.R -- look at binding-profile correlations for chosen targets
#   source("explore_binding_targets.R")          # from the CHARM folder, in R
# Plots, per cell line (columns) and direction (rows):
#   grey   : ALL profile pairs (from Binding_Correlation_Reference.qs2)
#   colour : pairs of profiles sharing the SAME eCLIP target (different KD RBPs),
#            one curve per target in TARGETS
# Also prints/saves a summary table. Edit the settings below, or set them before sourcing.
# =============================================================================
if (!exists("TARGETS"))   TARGETS   <- c("PTBP1", "RBFOX2")
if (!exists("PREFIX"))    PREFIX    <- "IR"                         # "ES" or "IR"
if (!exists("DATASETS"))  DATASETS  <- c("Both Cells", "K562", "HEPG2")
if (!exists("DIRS"))      DIRS      <- c("inc", "dec")
if (!exists("DATA_DIR"))  DATA_DIR  <- "data/QS2_Files"
if (!exists("OUT_PREFIX")) OUT_PREFIX <- "binding_target_distributions"

suppressMessages({ library(ggplot2); library(qs2) })
source("plot_theme.R")
source("correlation_reference.R")     # bind_key(), bind_hist_density(), CORR_BIN_WIDTH

# correlations among profiles of ONE target (rows "KD__target"), one value per pair of KD RBPs
target_pair_r <- function(M, target) {
  M  <- M[!duplicated(rownames(M)), , drop = FALSE]
  tg <- sub("^[^_]*__", "", rownames(M))
  hit <- which(tg == target)
  if (length(hit) == 0) hit <- which(toupper(tg) == toupper(target))
  if (length(hit) < 3) return(list(r = numeric(0), n_profiles = length(hit)))
  C <- suppressWarnings(stats::cor(t(M[hit, , drop = FALSE]), use = "pairwise.complete.obs"))
  r <- C[upper.tri(C)]
  list(r = r[is.finite(r)], n_profiles = length(hit))
}

# value r such that a fraction q of all pairs lie below it (from the histogram)
hist_quantile <- function(counts, q) {
  cs <- cumsum(counts) / sum(counts)
  -1 + CORR_BIN_WIDTH * (which(cs >= q)[1])
}

# density curves + summary for one cell line / direction
target_panel <- function(M, ref_entry, targets, panel_label) {
  dens <- list(); summ <- list()
  d_all <- bind_hist_density(ref_entry$counts_all, ref_entry$n_all)
  dens[[1]] <- data.frame(x = d_all$x, y = d_all$y, Group = "All profile pairs", Panel = panel_label)
  q95 <- hist_quantile(ref_entry$counts_all, 0.95); q99 <- hist_quantile(ref_entry$counts_all, 0.99)
  for (tg in targets) {
    res <- target_pair_r(M, tg)
    if (length(res$r) < 5) {
      summ[[tg]] <- data.frame(Panel = panel_label, Target = tg, Profiles = res$n_profiles, Pairs = length(res$r),
                               Median_r = NA, Pct_above_all_p95 = NA, Pct_above_all_p99 = NA)
      next
    }
    d <- stats::density(res$r, from = -1, to = 1, n = 512, bw = max(0.02, stats::bw.nrd0(res$r)))
    dens[[length(dens) + 1]] <- data.frame(x = d$x, y = d$y, Panel = panel_label,
                                           Group = sprintf("%s (%s pairs)", tg, format(length(res$r), big.mark = ",")))
    summ[[tg]] <- data.frame(Panel = panel_label, Target = tg, Profiles = res$n_profiles, Pairs = length(res$r),
                             Median_r = round(stats::median(res$r), 3),
                             Pct_above_all_p95 = round(100 * mean(res$r > q95), 1),
                             Pct_above_all_p99 = round(100 * mean(res$r > q99), 1))
  }
  list(dens = do.call(rbind, dens), summary = do.call(rbind, summ))
}

if (!exists("EXPLORE_NO_RUN")) {
  ref <- qs2::qs_read(file.path(DATA_DIR, "Binding_Correlation_Reference.qs2"))
  all_d <- list(); all_s <- list()
  for (ds in DATASETS) for (dr in DIRS) {
    key <- bind_key(PREFIX, ds, dr)
    if (is.null(ref[[key]])) { cat("no reference entry for", key, "- skipped\n"); next }
    f <- sprintf("similar_binding_%s_%s_%s", PREFIX, ifelse(ds == "Both Cells", "both", ds), dr)
    cat("reading", f, "...\n")
    M <- qs2::qs_read(file.path(DATA_DIR, paste0(f, ".qs2")))$mean_profiles
    if (!any(sub("^[^_]*__", "", rownames(M)) %in% TARGETS) &&
        !any(toupper(sub("^[^_]*__", "", rownames(M))) %in% toupper(TARGETS))) {
      tgs <- sort(unique(sub("^[^_]*__", "", rownames(M))))
      cat("None of", paste(TARGETS, collapse = ", "), "found as a target. Available targets include:\n",
          paste(head(tgs, 60), collapse = ", "), "\n")
    }
    pn <- target_panel(M, ref[[key]], TARGETS, paste(ds, dr, sep = "|"))
    pn$dens$Dataset   <- ds
    pn$dens$Direction <- ifelse(dr == "inc", "Increased events", "Decreased events")
    pn$summary$Dataset <- ds; pn$summary$Direction <- dr
    all_d[[key]] <- pn$dens; all_s[[key]] <- pn$summary
    rm(M); invisible(gc())
  }
  dd <- do.call(rbind, all_d); ss <- do.call(rbind, all_s); rownames(ss) <- NULL
  dd$Dataset <- factor(dd$Dataset, levels = DATASETS)
  grp <- unique(dd$Group); grp_t <- setdiff(grp, "All profile pairs")
  # same colour per target across panels (the pair count in the label differs between panels)
  dd$Target <- ifelse(dd$Group == "All profile pairs", "All profile pairs", sub(" \\(.*$", "", dd$Group))
  cols <- c("All profile pairs" = "#8C8C8C", stats::setNames(c("#BA3B46", "#2E7FB8", "#E0A030", "#3E9B6E", "#7272AB", "#B8639B")[seq_along(TARGETS)], TARGETS))
  p <- ggplot(dd, aes(x, y, colour = Target, fill = Target)) +
    geom_area(data = dd[dd$Target == "All profile pairs", ], alpha = 0.35, colour = NA, position = "identity") +
    geom_line(linewidth = 0.9) +
    geom_area(data = dd[dd$Target != "All profile pairs", ], alpha = 0.18, colour = NA, position = "identity") +
    facet_grid(Direction ~ Dataset, scales = "free_y") +
    scale_colour_manual(values = cols, name = NULL, breaks = names(cols)) +
    scale_fill_manual(values = cols, guide = "none") +
    coord_cartesian(xlim = c(-1, 1)) +
    labs(x = "Pearson's correlation (r)", y = "Density",
         title = sprintf("Binding-profile correlations: all pairs vs same eCLIP target (%s)", if (PREFIX == "ES") "exon skipping" else "intron retention"),
         caption = "Colour curves: pairs of profiles with the same eCLIP target and different KD RBPs. Densities are area-normalised.") +
    theme_charm() + theme(legend.position = "bottom", strip.background = element_blank())
  out_png <- sprintf("%s_%s.png", OUT_PREFIX, PREFIX)
  ggsave(out_png, p, width = 13, height = 7, dpi = 150, bg = "white")
  write.csv(ss, sprintf("%s_%s_summary.csv", OUT_PREFIX, PREFIX), row.names = FALSE)
  print(ss)
  cat("\nSaved", out_png, "and the summary CSV in", getwd(), "\n")
  print(p)
}

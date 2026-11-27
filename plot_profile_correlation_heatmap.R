# =============================================================================
# plot_profile_correlation_heatmap.R -- pairwise Pearson's r between binding profiles
#   source("plot_profile_correlation_heatmap.R")          # from the CHARM folder, in R
#
# For a chosen set of binding profiles ("KDRBP__Target", e.g. "U2AF1__U2AF2"), correlates every pair
# of mean binding profiles (Pearson's r across metagene positions, as in Binding > Similar Profiles)
# and draws one heatmap per direction: Increased Events and Decreased Events, side by side.
# Profiles are shown as "KDRBP - Target". Writes one PNG per event type.
#
# Settings -- set before sourcing, or edit the defaults:
# =============================================================================
if (!exists("PROFILES")) PROFILES <- c("RBM39__U2AF1", "RBM39__U2AF2", "U2AF1__U2AF1",
                                       "U2AF1__U2AF2", "U2AF2__U2AF1", "U2AF2__U2AF2")
if (!exists("PREFIXES")) PREFIXES <- c("ES")          # "ES" (exon skipping) and/or "IR" (intron retention)
if (!exists("DATASET"))  DATASET  <- "Both Cells"     # "Both Cells", "K562" or "HEPG2"
if (!exists("DATA_DIR")) DATA_DIR <- "data/QS2_Files"
if (!exists("OUT_DIR"))  OUT_DIR  <- "."
if (!exists("WIDTH"))    WIDTH    <- 11               # inches, both panels together
if (!exists("HEIGHT"))   HEIGHT   <- 4.5

suppressMessages({ library(ggplot2); library(qs2); library(ggpubr) })
source("plot_theme.R")

disp <- function(x) gsub("__", " - ", x, fixed = TRUE)

load_profiles <- function(prefix, dataset, direction) {
  f <- file.path(DATA_DIR, sprintf("similar_binding_%s_%s_%s.qs2", prefix,
                                   if (dataset == "Both Cells") "both" else dataset, direction))
  if (!file.exists(f)) stop("Missing ", f)
  mp <- qs2::qs_read(f)$mean_profiles
  mp[!duplicated(rownames(mp)), , drop = FALSE]
}

cor_heatmap <- function(mp, profiles, title) {
  miss <- setdiff(profiles, rownames(mp))
  if (length(miss)) warning(title, ": not in this dataset, left out: ", paste(miss, collapse = ", "))
  ids <- intersect(profiles, rownames(mp))
  if (length(ids) < 2) stop(title, ": fewer than 2 of the requested profiles are available.")
  C <- stats::cor(t(mp[ids, , drop = FALSE]), use = "pairwise.complete.obs")
  lev <- sort(disp(ids))
  d <- data.frame(x = disp(rep(ids, times = length(ids))), y = disp(rep(ids, each = length(ids))),
                  r = as.vector(C))
  d$x <- factor(d$x, levels = lev); d$y <- factor(d$y, levels = lev)
  ggplot(d, aes(x, y, fill = r)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%.2f", r)), size = CHARM_GEOM_TEXT, family = CHARM_FONT) +
    scale_fill_gradient2(low = "#53A2BE", mid = "white", high = "#BA3B46", midpoint = 0,
                         limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1), name = "Pearson's r") +
    coord_fixed() +
    labs(title = title, x = NULL, y = NULL) +
    theme_charm() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.line = element_blank(),
          axis.ticks = element_blank())
}

for (prefix in PREFIXES) {
  p_inc <- cor_heatmap(load_profiles(prefix, DATASET, "inc"), PROFILES, "Increased Events")
  p_dec <- cor_heatmap(load_profiles(prefix, DATASET, "dec"), PROFILES, "Decreased Events")
  p <- ggpubr::ggarrange(p_inc, p_dec, ncol = 2, common.legend = TRUE, legend = "right")
  out <- file.path(OUT_DIR, sprintf("profile_correlation_heatmap_%s_%s.png", prefix, gsub(" ", "", DATASET)))
  ggsave(out, p, width = WIDTH, height = HEIGHT, dpi = 300, bg = "white")
  cat("Saved", out, "\n")
}

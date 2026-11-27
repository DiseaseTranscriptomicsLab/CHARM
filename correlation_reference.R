# =============================================================================
# correlation_reference.R -- context for "how unusual is this correlation?"
# -----------------------------------------------------------------------------
# Correlating many things guarantees some "good" correlations by chance. For each
# analysis (Expression t-statistics, Splicing dPSI) and cell line, a reference file
# built OFFLINE by build_correlation_reference.R holds:
#   * observed : every RBP-vs-RBP correlation in that dataset (sorted)
#   * null     : the same correlations after shuffling the values WITHIN each gene
#                (or event) across RBPs, B times. This keeps each gene's own value
#                distribution and the missing-data pattern, but removes any real
#                relationship between RBPs. Stored as a histogram (bin width 0.001).
# From these the app reports, for any correlation r:
#   * Percentile : % of observed RBP pairs at least as extreme, in the same
#                  direction ("top 1%" of positive / of negative correlations)
#   * NullP      : share of null pairs at least as extreme (an empirical false-
#                  positive rate); never exactly 0 -- floor is 1/(1 + n_null).
# Caveat: each RBP takes part in many pairs, so pairs are not independent. These
# are descriptive measures of how unusual a correlation is, not formal p-values.
# If Correlation_Reference.qs2 is missing, everything here silently switches off.
# =============================================================================

CORR_BIN_WIDTH <- 0.001

# keys look like "expr|Both Cells", "splice|K562"
corr_key <- function(type, dataset) paste0(type, "|", dataset)

corr_has_ref <- function(key, ref = get0("CORR_REF", inherits = TRUE)) {
  !is.null(key) && !is.null(ref) && !is.null(ref[[key]])
}

corr_context <- function(r, key, ref = get0("CORR_REF", inherits = TRUE)) {
  n   <- length(r)
  out <- data.frame(Percentile = rep(NA_real_, n), NullP = rep(NA_real_, n))
  if (n == 0 || !corr_has_ref(key, ref)) return(out)
  e <- ref[[key]]

  # empirical percentile among observed RBP pairs
  obs   <- e$observed                       # stored sorted
  n_obs <- length(obs)
  n_ge  <- n_obs - findInterval(r, obs, left.open = TRUE)   # observed >= r
  n_le  <- findInterval(r, obs)                             # observed <= r
  out$Percentile <- 100 * pmax(ifelse(r >= 0, n_ge, n_le), 1) / n_obs

  # tail probability under the permutation null (histogram, conservative: the bin
  # holding r counts as "at least as extreme")
  counts <- e$null_counts
  nb     <- length(counts)
  cs     <- cumsum(counts)
  idx    <- pmin(pmax(floor((r + 1) / CORR_BIN_WIDTH) + 1, 1), nb)
  k      <- ifelse(r >= 0, e$null_n - c(0, cs)[idx], cs[idx])
  out$NullP <- (1 + k) / (1 + e$null_n)
  out
}

corr_fmt_pct <- function(p) {
  ifelse(is.na(p), "", ifelse(p >= 10, sprintf("%.0f%%", p),
                       ifelse(p >= 1,  sprintf("%.1f%%", p), sprintf("%.2f%%", p))))
}
corr_fmt_p <- function(p, floor_p = NA_real_) {
  txt <- ifelse(p < 0.001, formatC(p, format = "e", digits = 1), sprintf("%.3f", p))
  ifelse(is.na(p), "",
         ifelse(!is.na(floor_p) & p <= floor_p * 1.0000001, paste0("< ", formatC(floor_p, format = "e", digits = 0)), txt))
}

# Bar-plot label: standard percentile among all observed pairs (% of pairs with a correlation <= this one).
# `pct` is the stored "at least as extreme" share (top x% for r >= 0, bottom x% for r < 0), converted here.
corr_std_pctl <- function(r, pct) ifelse(r >= 0, 100 - pct, pct)   # -> % of pairs with a correlation <= r
corr_fmt_pctl <- function(p) format(signif(p, 3), scientific = FALSE, trim = TRUE, drop0trailing = TRUE)
corr_pctl_label <- function(r, pct) {
  ifelse(is.na(pct), "", paste0("Percentile ", corr_fmt_pctl(corr_std_pctl(r, pct))))
}

# Short label for a bar: "top 0.8%" (positive) / "bottom 0.8%" (negative)
corr_bar_label <- function(r, pct) {          # lower-case form for use inside longer labels
  ifelse(is.na(pct), "", paste0("percentile ", corr_fmt_pctl(corr_std_pctl(r, pct))))
}

# One-line context for a single scatter plot (one RBP vs one RBP / uploaded data): where the Spearman rho falls among
# all observed RBP pairs, and its tail probability under the permutation null. NULL when no reference is available.
corr_scatter_subtitle <- function(rho, ref_key, ref = get0("CORR_REF", inherits = TRUE)) {
  if (is.na(rho) || !corr_has_ref(ref_key, ref)) return(NULL)
  ctx     <- corr_context(rho, ref_key, ref)
  floor_p <- 1 / (1 + ref[[ref_key]]$null_n)
  p_txt   <- corr_fmt_p(ctx$NullP, floor_p)
  paste0("Percentile ", corr_fmt_pctl(corr_std_pctl(rho, ctx$Percentile)), " among all CHARM RBP pairs  |  FPR ",
         ifelse(startsWith(p_txt, "<"), "", "= "), p_txt, " (permutation null)")
}

# Density panel: observed (all RBP pairs) vs null (shuffled), a vertical line per
# displayed correlation.
corr_density_panel <- function(top_cor, key, up_color = "#BA3B46", down_color = "#53A2BE",
                               ref = get0("CORR_REF", inherits = TRUE)) {
  e <- ref[[key]]
  dens_obs <- stats::density(e$observed, from = -1, to = 1, n = 512,
                             bw = max(0.004, stats::bw.nrd0(e$observed)))
  centers  <- seq(-1 + CORR_BIN_WIDTH / 2, 1 - CORR_BIN_WIDTH / 2, by = CORR_BIN_WIDTH)
  w        <- e$null_counts / sum(e$null_counts)
  mu       <- sum(w * centers)
  sd_w     <- sqrt(sum(w * (centers - mu)^2))
  dens_nul <- stats::density(centers, weights = w, from = -1, to = 1, n = 512,
                             bw = max(0.002, 0.9 * sd_w * e$null_n^(-1/5)))
  n_pairs <- format(length(e$observed), big.mark = ",")
  lab_obs <- sprintf("Observed (%s RBP pairs)", n_pairs)
  is_gsea <- grepl("^gsea", key)
  lab_nul <- sprintf("Null (%s RBP pairs \u00d7 %d %s)", n_pairs, e$B,
                     if (is_gsea) "random gene-set draws" else "permutations")
  # Each curve is scaled to its own peak: the null is far narrower (so far taller)
  # than the observed distribution, and on a shared density axis the observed one
  # would look flat. What matters here is where the lines fall relative to each shape.
  dd <- rbind(data.frame(x = dens_obs$x, y = dens_obs$y / max(dens_obs$y), Distribution = lab_obs),
              data.frame(x = dens_nul$x, y = dens_nul$y / max(dens_nul$y), Distribution = lab_nul))
  dd$Distribution <- factor(dd$Distribution, levels = c(lab_obs, lab_nul))
  ymax <- 1
  rng  <- range(c(e$observed, top_cor$Correlation, 0))
  pad  <- max(0.05, 0.06 * diff(rng))
  xlim <- c(max(-1, rng[1] - pad), min(1, rng[2] + pad))

  # Displayed RBPs: one vertical line each, red (rho > 0) / blue (rho < 0), with a legend key next to the fills
  lines_df <- data.frame(x = top_cor$Correlation,
                         Status = ifelse(top_cor$Correlation > 0, "Positive", "Negative"))
  st <- intersect(c("Positive", "Negative"), unique(lines_df$Status))
  line_labs <- if (length(st) == 1) c(Positive = "Displayed RBPs", Negative = "Displayed RBPs")
               else c(Positive = "Displayed RBPs (\u03c1 > 0)", Negative = "Displayed RBPs (\u03c1 < 0)")
  ggplot2::ggplot() +
    ggplot2::geom_area(data = dd, ggplot2::aes(x = x, y = y, fill = Distribution),
                       alpha = 0.55, position = "identity", colour = NA) +
    ggplot2::geom_vline(data = lines_df, ggplot2::aes(xintercept = x, colour = Status),
                        linewidth = 0.5, alpha = 0.8, key_glyph = "vline") +
    ggplot2::scale_fill_manual(values = stats::setNames(c("#7272AB", "#8C8C8C"), c(lab_obs, lab_nul)),
                               name = NULL, guide = ggplot2::guide_legend(order = 1)) +
    ggplot2::scale_colour_manual(values = c("Positive" = up_color, "Negative" = down_color), breaks = st,
                                 labels = line_labs[st], name = NULL, guide = ggplot2::guide_legend(order = 2)) +
    ggplot2::scale_x_continuous(limits = c(-1, 1)) +
    ggplot2::coord_cartesian(xlim = xlim, ylim = c(0, 1.12)) +
    ggplot2::labs(x = "Spearman's correlation (\u03c1)", y = "Relative density") +
    theme_charm() +
    ggplot2::theme(legend.position = "bottom", legend.box = "horizontal",
                   legend.spacing.x = ggplot2::unit(2, "pt"),
                   legend.text = ggplot2::element_text(margin = ggplot2::margin(l = 2, r = 14)))
}

# Add Percentile / NullP columns and the density panel to a correlation result.
# Returns NULL (-> caller keeps its plain bar plot) when no reference is available.
corr_decorate <- function(cor_results, top_cor, bar_plot, ref_key,
                          up_color = "#BA3B46", down_color = "#53A2BE") {
  if (!corr_has_ref(ref_key) || nrow(top_cor) == 0) return(NULL)
  ctx_all <- corr_context(cor_results$Correlation, ref_key)
  cor_results$Percentile <- ctx_all$Percentile
  cor_results$NullP      <- ctx_all$NullP
  ctx <- corr_context(top_cor$Correlation, ref_key)
  top_cor$Percentile <- ctx$Percentile
  top_cor$NullP      <- ctx$NullP

  floor_p <- 1 / (1 + CORR_REF[[ref_key]]$null_n)
  lab_df  <- top_cor
  p_txt <- corr_fmt_p(lab_df$NullP, floor_p)
  lab_df$ctx_label <- paste0("\u03c1 = ", sprintf("%.2f", lab_df$Correlation), "  |  ",
                             corr_pctl_label(lab_df$Correlation, lab_df$Percentile), "  |  FPR ",
                             ifelse(startsWith(p_txt, "<"), "", "= "), p_txt)

  bar <- bar_plot +
    ggplot2::geom_text(data = lab_df,
                       ggplot2::aes(x = RBP, y = pmax(Correlation, 0), label = ctx_label),
                       inherit.aes = FALSE, hjust = -0.08,      # starts at the bar end (positive) or at 0 (negative)
                       colour = "black", size = CHARM_GEOM_TEXT, family = CHARM_FONT) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 1.15)))
  dens <- corr_density_panel(top_cor, ref_key, up_color, down_color)
  combined <- ggpubr::ggarrange(bar, dens, ncol = 1, heights = c(1, 0.7))
  # tables: Percentile as a standard percentile (% of RBP pairs with a correlation <= this one), as in the plot labels
  cor_results$Percentile <- corr_std_pctl(cor_results$Correlation, cor_results$Percentile)
  top_cor$Percentile     <- corr_std_pctl(top_cor$Correlation, top_cor$Percentile)
  list(correlation_table = cor_results, top_table = top_cor, heatmap = combined,
       bar_plot = bar, density_plot = dens)
}


# =============================================================================
# BINDING PROFILES -- Similar Profiles context
# -----------------------------------------------------------------------------
# Binding profiles are rows "KDRBP__Target" (mean binding signal of eCLIP Target
# around the events changed by knocking down KDRBP) over metagene positions.
# build_binding_reference.R stores, per event type x cell line x direction, the
# histogram (bin 0.001) of the Pearson correlation of EVERY pair of profiles:
#   counts_all  : all pairs
#   counts_same : only pairs that share the same eCLIP Target (different KD RBP)
# Here positions are the "rows"; profiles are smooth, so neighbouring positions are
# not independent and no permutation null is provided for binding (see notes).
# File: Binding_Correlation_Reference.qs2  ->  BIND_REF (NULL if missing = feature off)
# =============================================================================

bind_key <- function(prefix, dataset, direction) paste("binding", prefix, dataset, direction, sep = "|")

bind_has_ref <- function(key, ref = get0("BIND_REF", inherits = TRUE)) {
  !is.null(key) && !is.null(ref) && !is.null(ref[[key]])
}

# % of all profile pairs at least as extreme as r, in the same direction
bind_percentile <- function(r, key, ref = get0("BIND_REF", inherits = TRUE)) {
  e  <- ref[[key]]
  cs <- cumsum(e$counts_all)
  nb <- length(cs)
  idx <- pmin(pmax(floor((r + 1) / CORR_BIN_WIDTH) + 1, 1), nb)
  k   <- ifelse(r >= 0, e$n_all - c(0, cs)[idx], cs[idx])
  100 * pmax(k, 1) / e$n_all
}

bind_plot_height <- function() if (!is.null(get0("BIND_REF", inherits = TRUE))) "760px" else "400px"

# weighted density estimate from a histogram of correlations (area = 1)
bind_hist_density <- function(counts, n) {
  centers <- seq(-1 + CORR_BIN_WIDTH / 2, 1 - CORR_BIN_WIDTH / 2, by = CORR_BIN_WIDTH)
  w  <- counts / sum(counts)
  mu <- sum(w * centers)
  s  <- sqrt(sum(w * (centers - mu)^2))
  stats::density(centers, weights = w, from = -1, to = 1, n = 512,
                 bw = max(0.01, 0.9 * s * n^(-1/5)))
}

bind_density_panel_basic <- function(cor_df, key, ref = get0("BIND_REF", inherits = TRUE),
                               up_color = "#BA3B46", down_color = "#53A2BE") {
  e  <- ref[[key]]
  fn <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
  d_all  <- bind_hist_density(e$counts_all, e$n_all)
  has_same <- !is.null(e$n_same) && e$n_same > 20
  lab_all  <- sprintf("All profile pairs (n = %s)", fn(e$n_all))
  lab_same <- sprintf("Same eCLIP target, different KD RBP (n = %s)", fn(e$n_same))
  dd <- data.frame(x = d_all$x, y = d_all$y, Distribution = lab_all)
  if (has_same) {
    d_same <- bind_hist_density(e$counts_same, e$n_same)
    dd <- rbind(dd, data.frame(x = d_same$x, y = d_same$y, Distribution = lab_same))
  }
  levs <- c(lab_all, if (has_same) lab_same)
  dd$Distribution <- factor(dd$Distribution, levels = levs)

  rng  <- range(c(cor_df$Correlation, dd$x[dd$y > 0.01 * max(dd$y)], 0))
  pad  <- max(0.05, 0.06 * diff(rng))
  xlim <- c(max(-1, rng[1] - pad), min(1, rng[2] + pad))
  ymax <- max(dd$y[dd$x >= xlim[1] & dd$x <= xlim[2]])

  lines_df <- data.frame(x = cor_df$Correlation,
                         Status = ifelse(cor_df$Correlation > 0, "Positive", "Negative"))
  lab_df <- utils::head(cor_df[order(-abs(cor_df$Correlation)), , drop = FALSE], 3)
  lab_df <- data.frame(x = lab_df$Correlation, y = ymax, label = as.character(lab_df$Profile))

  ggplot2::ggplot() +
    ggplot2::geom_area(data = dd, ggplot2::aes(x = x, y = y, fill = Distribution),
                       alpha = 0.5, position = "identity", colour = NA) +
    ggplot2::geom_vline(data = lines_df, ggplot2::aes(xintercept = x, colour = Status),
                        linewidth = 0.5, alpha = 0.8) +
    ggrepel::geom_text_repel(data = lab_df, ggplot2::aes(x = x, y = y, label = label),
                             direction = "x", size = CHARM_GEOM_TEXT, family = CHARM_FONT,
                             fontface = "italic", segment.size = 0.2, min.segment.length = 0,
                             box.padding = 0.2, max.overlaps = Inf) +
    ggplot2::scale_fill_manual(values = stats::setNames(c("#7272AB", "#E0A030")[seq_along(levs)], levs), name = NULL) +
    ggplot2::scale_colour_manual(values = c("Positive" = up_color, "Negative" = down_color), guide = "none") +
    ggplot2::scale_x_continuous(limits = c(-1, 1)) +
    ggplot2::coord_cartesian(xlim = xlim, ylim = c(0, ymax * 1.15)) +
    ggplot2::labs(x = "Pearson's correlation (r)", y = "Density",
                  title = "Where these profiles fall among all profile pairs",
                  caption = "Lines: displayed profiles. Densities are area-normalised, so the two groups are directly comparable.") +
    theme_charm() +
    ggplot2::theme(legend.position = "bottom", legend.direction = "vertical")
}



# -----------------------------------------------------------------------------
# Three-density view for ONE query profile "KD__Target" (e.g. PTBP1__RBFOX2):
#   1. all profile pairs
#   2. profiles of the same eCLIP target (RBFOX2 profiles) against each other
#   3. the "self" profile of that target (RBFOX2__RBFOX2) against the other RBFOX2 profiles
# Vertical line = correlation of the query with the self profile (RBFOX2__RBFOX2).
# A combined view pools the increased and decreased events (two lines, one per direction).
# -----------------------------------------------------------------------------
bind_has_target_ref <- function(key, ref = get0("BIND_REF", inherits = TRUE)) {
  bind_has_ref(key, ref) && !is.null(ref[[key]]$by_target)
}

# % of pairs in a histogram at least as extreme as r (same direction)
bind_pct_counts <- function(r, counts) {
  cs <- cumsum(counts); n <- cs[length(cs)]; nb <- length(cs)
  idx <- pmin(pmax(floor((r + 1) / CORR_BIN_WIDTH) + 1, 1), nb)
  k   <- ifelse(r >= 0, n - c(0, cs)[idx], cs[idx])
  100 * pmax(k, 1) / n
}

# histogram counts (bin CORR_BIN_WIDTH) from a vector of correlations
bind_r_to_counts <- function(r) {
  nb <- as.integer(round(2 / CORR_BIN_WIDTH))
  as.numeric(tabulate(pmin(pmax(floor((pmin(pmax(r, -1), 1) + 1) / CORR_BIN_WIDTH) + 1, 1), nb), nbins = nb))
}

# legend names of the three groups (used by the app panel and the exploration script)
bind_group_names <- function(target) {
  c(all    = "All profile pairs",
    target = sprintf("%s profiles vs each other (different KD RBPs)", target),
    self   = sprintf("%s__%s vs the other %s profiles", target, target, target))
}

# r below which a fraction q of the pairs in a histogram lie
bind_hist_quantile <- function(counts, q) {
  cs <- cumsum(counts) / sum(counts)
  -1 + CORR_BIN_WIDTH * (which(cs >= q)[1] - 0.5)
}

# Zoomed x range: trimmed (0.5%-99.5%) range of the same-target and self curves plus the
# query line(s); falls back to the all-pairs curve if neither exists.
bind_zoom_window <- function(win, x_lines, zoom = TRUE) {
  if (!zoom) return(c(-1, 1))
  rng <- range(c(win, x_lines[is.finite(x_lines)]))
  pad <- max(0.03, 0.08 * diff(rng))
  c(max(-1, rng[1] - pad), min(1, rng[2] + pad))
}

BIND_GROUP_COLOURS <- c(all = "#8C8C8C", target = "#E0A030", self = "#7272AB")

# density curves + % position of each line, for a list of groups(counts, n, name)
bind_make_curves <- function(grp, r_lines) {
  dd <- list(); info <- list()
  for (g in names(grp)) {
    x <- grp[[g]]
    if (is.null(x$n) || is.na(x$n) || x$n < 5) next            # too few pairs for a density
    d <- bind_hist_density(x$counts, x$n)
    dd[[g]]   <- data.frame(x = d$x, y = d$y, Group = g, stringsAsFactors = FALSE)
    row <- data.frame(Group = g, Name = x$name, n = x$n, stringsAsFactors = FALSE)
    for (nm in names(r_lines)) {
      row[[paste0("pct_", nm)]] <- if (is.finite(r_lines[[nm]])) bind_pct_counts(r_lines[[nm]], x$counts) else NA_real_
    }
    info[[g]] <- row
  }
  if (length(dd) == 0) return(NULL)
  focus <- intersect(c("target", "self"), names(dd))
  if (length(focus) == 0) focus <- names(dd)
  win <- range(unlist(lapply(focus, function(g) c(bind_hist_quantile(grp[[g]]$counts, 0.005),
                                                  bind_hist_quantile(grp[[g]]$counts, 0.995)))))
  info <- do.call(rbind, info); rownames(info) <- NULL
  list(dens = do.call(rbind, dd), info = info, win = win)
}

bind_density_data <- function(query_id, key, mean_profiles, ref = get0("BIND_REF", inherits = TRUE)) {
  e       <- ref[[key]]
  target  <- sub("^[^_]*__", "", query_id)
  self_id <- paste0(target, "__", target)
  r_line  <- NA_real_
  if (!identical(self_id, query_id) && all(c(query_id, self_id) %in% rownames(mean_profiles))) {
    r_line <- suppressWarnings(stats::cor(as.numeric(mean_profiles[query_id, ]),
                                          as.numeric(mean_profiles[self_id, ]),
                                          use = "pairwise.complete.obs"))
    if (!is.finite(r_line)) r_line <- NA_real_
  }
  nm  <- bind_group_names(target)
  grp <- list(all = list(counts = e$counts_all, n = e$n_all, name = nm[["all"]]))
  # same-target pairs: profiles of this eCLIP target under different KD RBPs (computed here,
  # a few thousand correlations at most)
  tg_ids <- rownames(mean_profiles)[sub("^[^_]*__", "", rownames(mean_profiles)) == target]
  if (length(tg_ids) >= 4) {
    Ct <- suppressWarnings(stats::cor(t(mean_profiles[tg_ids, , drop = FALSE]), use = "pairwise.complete.obs"))
    rt <- Ct[upper.tri(Ct)]; rt <- rt[is.finite(rt)]
    if (length(rt) >= 5) grp$target <- list(counts = bind_r_to_counts(rt), n = length(rt), name = nm[["target"]])
    # self profile (e.g. RBFOX2__RBFOX2) against the other profiles of the same target
    if (self_id %in% tg_ids) {
      rs <- Ct[self_id, setdiff(tg_ids, self_id)]; rs <- rs[is.finite(rs)]
      if (length(rs) >= 5) grp$self <- list(counts = bind_r_to_counts(rs), n = length(rs), name = nm[["self"]])
    }
  }
  cv <- bind_make_curves(grp, c(line = r_line))
  if (is.null(cv)) return(NULL)
  cv$info$Percentile <- cv$info$pct_line
  c(cv, list(r_line = r_line, grp = grp, query = query_id, self_id = self_id, target = target))
}

# --- inclusion + exclusion pooled: ALL profiles of both directions in one pool -------------
# Rows are tagged "inc|KD__Target" / "dec|KD__Target"; every correlation among the pooled
# profiles counts, including inc-vs-dec pairs. The reference (key from bind_pool_key()) holds
# the all-pairs histogram of the pool; the same-target and self curves are computed on the fly.
bind_pool_key <- function(prefix, dataset) paste("binding", prefix, dataset, "pool", sep = "|")

bind_pool_matrix <- function(mp_inc, mp_dec) {
  if (is.null(mp_inc) || is.null(mp_dec) || ncol(mp_inc) != ncol(mp_dec)) return(NULL)
  mp_inc <- mp_inc[!duplicated(rownames(mp_inc)), , drop = FALSE]
  mp_dec <- mp_dec[!duplicated(rownames(mp_dec)), , drop = FALSE]
  P <- rbind(mp_inc, mp_dec)
  rownames(P) <- c(paste0("inc|", rownames(mp_inc)), paste0("dec|", rownames(mp_dec)))
  P
}

# Same result structure as bind_density_data(), but r_lines = c(inc = , dec = )
bind_density_data_pool <- function(query_id, key, mp_inc, mp_dec, ref = get0("BIND_REF", inherits = TRUE)) {
  e <- ref[[key]]
  P <- bind_pool_matrix(mp_inc, mp_dec)
  if (is.null(e) || is.null(P)) return(NULL)
  target  <- sub("^[^_]*__", "", query_id)
  self_id <- paste0(target, "__", target)
  line_of <- function(mp) {
    if (is.null(mp) || identical(self_id, query_id) || !all(c(query_id, self_id) %in% rownames(mp))) return(NA_real_)
    r <- suppressWarnings(stats::cor(as.numeric(mp[query_id, ]), as.numeric(mp[self_id, ]), use = "pairwise.complete.obs"))
    if (is.finite(r)) r else NA_real_
  }
  r_lines <- c(inc = line_of(mp_inc), dec = line_of(mp_dec))
  nm  <- bind_group_names(target)
  grp <- list(all = list(counts = e$counts_all, n = e$n_all, name = nm[["all"]]))
  tg_ids <- rownames(P)[sub("^[^_]*__", "", rownames(P)) == target]
  if (length(tg_ids) >= 4) {
    Ct <- suppressWarnings(stats::cor(t(P[tg_ids, , drop = FALSE]), use = "pairwise.complete.obs"))
    rt <- Ct[upper.tri(Ct)]; rt <- rt[is.finite(rt)]
    if (length(rt) >= 5) grp$target <- list(counts = bind_r_to_counts(rt), n = length(rt), name = nm[["target"]])
    selfs <- tg_ids[sub("^(inc|dec)\\|", "", tg_ids) == self_id]       # self profile of each direction
    if (length(selfs) > 0) {
      rs <- unlist(lapply(seq_along(selfs), function(i) Ct[selfs[i], setdiff(tg_ids, selfs[seq_len(i)])]))
      rs <- rs[is.finite(rs)]
      if (length(rs) >= 5) grp$self <- list(counts = bind_r_to_counts(rs), n = length(rs), name = nm[["self"]])
    }
  }
  cv <- bind_make_curves(grp, r_lines)
  if (is.null(cv)) return(NULL)
  c(cv, list(r_lines = r_lines, query = query_id, self_id = self_id, target = target))
}

bind_fmt_n <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)

# One direction: three densities + the line of the query
bind_density_panel <- function(query_id, key, mean_profiles, ref = get0("BIND_REF", inherits = TRUE),
                               facet = FALSE, zoom = TRUE) {
  D <- bind_density_data(query_id, key, mean_profiles, ref)
  if (is.null(D)) return(ggplot2::ggplot() + theme_charm_void() + ggplot2::ggtitle("No reference data"))
  has_line <- is.finite(D$r_line)
  D$info$Label <- sprintf("%s (n = %s)%s", D$info$Name, bind_fmt_n(D$info$n),
                          ifelse(is.na(D$info$Percentile), "",
                                 paste0("  |  line: ", corr_bar_label(D$r_line, D$info$Percentile))))
  dd <- D$dens
  dd$Label <- factor(D$info$Label[match(dd$Group, D$info$Group)], levels = D$info$Label)
  cols <- stats::setNames(BIND_GROUP_COLOURS[D$info$Group], D$info$Label)
  xlim <- bind_zoom_window(D$win, D$r_line, zoom)
  ymax <- max(dd$y[dd$x >= xlim[1] & dd$x <= xlim[2]])

  caption <- if (has_line) {
    sprintf("Line: %s vs %s, r = %.2f. Densities are area-normalised.", D$query, D$self_id, D$r_line)
  } else if (identical(D$query, D$self_id)) {
    "This is the self profile itself, so there is no line to draw."
  } else {
    sprintf("%s is not available in this dataset, so no line is drawn.", D$self_id)
  }

  p <- ggplot2::ggplot(dd, ggplot2::aes(x = x, y = y, fill = Label)) +
    ggplot2::geom_area(alpha = 0.5, position = "identity", colour = NA)
  if (has_line) {
    p <- p +
      ggplot2::geom_vline(xintercept = D$r_line, colour = "#BA3B46", linewidth = 0.8) +
      ggplot2::annotate("text", x = D$r_line, y = ymax * 1.12,
                        label = sprintf("%s\nr = %.2f", D$query, D$r_line),
                        hjust = if (D$r_line > mean(xlim)) 1.05 else -0.05, vjust = 1,
                        colour = "#BA3B46", size = CHARM_GEOM_TEXT, family = CHARM_FONT)
  }
  p <- p +
    ggplot2::scale_fill_manual(values = cols, name = NULL) +
    ggplot2::scale_x_continuous(limits = c(-1, 1)) +
    ggplot2::coord_cartesian(xlim = xlim, ylim = c(0, ymax * 1.25)) +
    ggplot2::labs(x = "Pearson's correlation (r)", y = "Density",
                  title = "Where this profile falls", caption = caption) +
    theme_charm() +
    ggplot2::theme(legend.position = "bottom", legend.direction = "vertical")
  if (facet) {
    p <- p + ggplot2::facet_wrap(~ Label, ncol = 1, scales = "free_y") +
      ggplot2::theme(legend.position = "none", strip.background = ggplot2::element_blank())
  }
  p
}

# Increased and decreased events together: all profiles of both directions in one pool,
# one line per direction (the profile's inc r and dec r with its self profile)
bind_density_panel_combined <- function(query_id, key_pool, mp_inc, mp_dec,
                                        ref = get0("BIND_REF", inherits = TRUE), zoom = TRUE) {
  if (!bind_has_ref(key_pool, ref))
    return(ggplot2::ggplot() + theme_charm_void() +
             ggplot2::ggtitle("Pooled reference not built yet (run build_binding_reference.R)"))
  C <- bind_density_data_pool(query_id, key_pool, mp_inc, mp_dec, ref)
  if (is.null(C)) return(ggplot2::ggplot() + theme_charm_void() + ggplot2::ggtitle("No reference data"))
  pc <- function(r, p) if (is.finite(r) && is.finite(p)) corr_bar_label(r, p) else "n/a"
  C$info$Label <- sprintf("%s (n = %s)  |  inc: %s, dec: %s", C$info$Name, bind_fmt_n(C$info$n),
                          mapply(pc, C$r_lines[["inc"]], C$info$pct_inc), mapply(pc, C$r_lines[["dec"]], C$info$pct_dec))
  dd <- C$dens
  dd$Label <- factor(C$info$Label[match(dd$Group, C$info$Group)], levels = C$info$Label)
  cols <- stats::setNames(BIND_GROUP_COLOURS[C$info$Group], C$info$Label)
  xlim <- bind_zoom_window(C$win, C$r_lines, zoom)
  ymax <- max(dd$y[dd$x >= xlim[1] & dd$x <= xlim[2]])

  lines_df <- data.frame(Direction = c("Increased events", "Decreased events"), x = unname(C$r_lines),
                         short = c("inc", "dec"), stringsAsFactors = FALSE)
  lines_df <- lines_df[is.finite(lines_df$x), , drop = FALSE]
  caption <- if (nrow(lines_df) > 0) {
    sprintf("Lines: %s vs %s. Curves use all increased and decreased profiles pooled (including inc-vs-dec pairs). Densities are area-normalised.", C$query, C$self_id)
  } else if (identical(C$query, C$self_id)) {
    "This is the self profile itself, so there is no line to draw."
  } else {
    sprintf("%s is not available in these datasets, so no lines are drawn.", C$self_id)
  }

  p <- ggplot2::ggplot() +
    ggplot2::geom_area(data = dd, ggplot2::aes(x = x, y = y, fill = Label),
                       alpha = 0.5, position = "identity", colour = NA)
  if (nrow(lines_df) > 0) {
    lines_df$Direction <- factor(lines_df$Direction, levels = c("Increased events", "Decreased events"))
    lines_df$label <- sprintf("%s r = %.2f", lines_df$short, lines_df$x)
    lines_df$y <- ymax * 1.15
    p <- p +
      ggplot2::geom_vline(data = lines_df, ggplot2::aes(xintercept = x, linetype = Direction),
                          colour = "#BA3B46", linewidth = 0.8) +
      ggrepel::geom_text_repel(data = lines_df, ggplot2::aes(x = x, y = y, label = label),
                               direction = "x", colour = "#BA3B46", size = CHARM_GEOM_TEXT, family = CHARM_FONT,
                               segment.size = 0.2, min.segment.length = 0, box.padding = 0.2, max.overlaps = Inf) +
      ggplot2::scale_linetype_manual(values = c("Increased events" = "solid", "Decreased events" = "dashed"),
                                     name = NULL, drop = FALSE)
  }
  p +
    ggplot2::scale_fill_manual(values = cols, name = NULL) +
    ggplot2::scale_x_continuous(limits = c(-1, 1)) +
    ggplot2::coord_cartesian(xlim = xlim, ylim = c(0, ymax * 1.3)) +
    ggplot2::labs(x = "Pearson's correlation (r)", y = "Density",
                  title = "Increased and decreased events together (all profiles pooled)", caption = caption) +
    theme_charm() +
    ggplot2::theme(legend.position = "bottom", legend.direction = "vertical", legend.box = "vertical")
}


# =============================================================================
# GSEA gene-set window (Expression > GSEA plot: click a gene set)
# -----------------------------------------------------------------------------
# build_gsea_reference.R stores, per cell line ("gsea|K562" ... in CORR_REF), a "gsea_sets" element:
#   obs         : NES of every RBP for every hallmark set (set x RBP)
#   null_counts : per hallmark set, histogram (bin 0.02, -6..6) of the NES that random gene sets of the same size
#                 get against ALL RBP rankings. It does not depend on the RBP, so it is stored once.
# The window shows (1) where the selected RBP's NES sits among all RBPs for that set and (2) its false-positive
# rate: the share of size-matched random gene sets with an NES at least as extreme in the same direction.
# Without the reference only (1) is shown (computed from the loaded per-RBP GSEA tables).
# =============================================================================
gsea_pw_id <- function(x) sub("^HALLMARK_", "", gsub(" ", "_", as.character(x)))

gsea_pw_has_ref <- function(key, ref = get0("CORR_REF", inherits = TRUE)) {
  !is.null(key) && !is.null(ref) && !is.null(ref[[key]]) && !is.null(ref[[key]]$gsea_sets)
}

# obs_fallback: named list (one data.frame per RBP with columns pathway, NES), used when no reference is built
gsea_pw_context <- function(key, pathway, nes, ref = get0("CORR_REF", inherits = TRUE), obs_fallback = NULL) {
  id <- gsea_pw_id(pathway)
  out <- list(id = id, nes = nes, obs = NULL, pct = NA_real_, null = NULL, null_n = NA_real_, size = NA_integer_,
              fpr = NA_real_, floor_p = NA_real_)
  out$why <- if (is.null(ref)) "Correlation_Reference.qs2 was not found when the app started."
             else if (is.null(ref[[key]])) sprintf("The reference has no entry '%s' (available: %s).", key, paste(names(ref), collapse = ", "))
             else if (is.null(ref[[key]]$gsea_sets)) sprintf("Entry '%s' has no gene-set backgrounds: rerun build_gsea_reference.R, then restart the app.", key)
             else NULL
  if (gsea_pw_has_ref(key, ref)) {
    g <- ref[[key]]$gsea_sets
    i <- match(id, g$ids)
    if (is.na(i)) out$why <- sprintf("Gene set '%s' not found in the reference (it has: %s ...).", id, paste(utils::head(g$ids, 3), collapse = ", "))
    if (!is.na(i)) {
      out$obs    <- stats::setNames(as.numeric(g$obs[i, ]), colnames(g$obs))
      out$size   <- g$size[i]
      out$B      <- ref[[key]]$B          # random gene sets drawn per hallmark set
      out$n_rank <- ref[[key]]$n_rbps     # RBP rankings each random set is scored on
      cnt        <- g$null_counts[i, ]
      nb         <- length(cnt)
      out$null_n <- g$null_n[i]
      out$null   <- data.frame(x = g$lo + (seq_len(nb) - 0.5) * g$bin, w = cnt / sum(cnt))
      idx        <- min(max(floor((nes - g$lo) / g$bin) + 1, 1), nb)
      cs         <- cumsum(cnt)
      k          <- if (nes >= 0) out$null_n - c(0, cs)[idx] else cs[idx]   # bin holding nes counts as "at least as extreme"
      out$fpr    <- (1 + k) / (1 + out$null_n)
      out$floor_p <- 1 / (1 + out$null_n)
    }
  }
  if (is.null(out$obs) && !is.null(obs_fallback)) {
    out$obs <- vapply(obs_fallback, function(d) {
      if (is.null(d) || !all(c("pathway", "NES") %in% colnames(d))) return(NA_real_)
      v <- d$NES[match(id, gsea_pw_id(d$pathway))]
      if (length(v)) as.numeric(v[1]) else NA_real_
    }, numeric(1))
  }
  if (!is.null(out$obs)) {
    o <- out$obs[is.finite(out$obs)]
    out$obs <- o
    if (length(o)) out$pct <- 100 * max(if (nes >= 0) sum(o >= nes) else sum(o <= nes), 1) / length(o)
  }
  out
}

gsea_pw_plot <- function(ctx, rbp, dataset = NULL, up_color = "#BA3B46", down_color = "#53A2BE",
                         obs_color = "#7272AB", null_color = "#8C8C8C") {
  # One panel: observed NES of this gene set across all RBP knockdowns (purple) and NES of size-matched random
  # gene sets (grey; same colours as the similar-RBP density plots). Each curve is scaled to a maximum of 1 so
  # both stay readable (the random-set curve is much narrower). Definitions are in the figure legend.
  line_col <- if (ctx$nes >= 0) up_color else down_color
  pw_label <- gsub("_", " ", ctx$id)
  has_obs  <- length(ctx$obs) >= 3
  has_null <- !is.null(ctx$null)
  if (!has_obs && !has_null)
    return(ggplot2::ggplot() + ggplot2::annotate("text", x = 0, y = 0, label = "No NES values available for this gene set", family = CHARM_FONT) + theme_charm_void())
  xr <- range(c(ctx$obs, ctx$nes, if (has_null) ctx$null$x[ctx$null$w > 1e-4]), na.rm = TRUE)
  pad <- 0.08 * diff(xr); xlim <- c(xr[1] - pad, xr[2] + pad)
  fmt3 <- function(x) format(signif(x, 3), scientific = FALSE, trim = TRUE, drop0trailing = TRUE)

  lab_obs  <- sprintf("RBP KDs (n = %d)", length(ctx$obs))
  lab_null <- if (has_null) {
    if (!is.null(ctx$B) && !is.null(ctx$n_rank)) sprintf("Random gene sets (%s sets \u00d7 %s RBP KDs)", ctx$B, ctx$n_rank)
    else "Random gene sets"
  } else NULL
  curves <- list()
  if (has_obs) {
    d <- stats::density(ctx$obs, adjust = 1.2, from = xlim[1], to = xlim[2], n = 512)
    curves$obs <- data.frame(x = d$x, y = d$y / max(d$y), grp = lab_obs)
  }
  if (has_null) {
    d <- stats::density(ctx$null$x, weights = ctx$null$w, bw = 0.12, from = xlim[1], to = xlim[2], n = 512)
    curves$null <- data.frame(x = d$x, y = d$y / max(d$y), grp = lab_null)
  }
  cd <- do.call(rbind, curves); cd$grp <- factor(cd$grp, levels = c(lab_obs, lab_null))

  stat_lines <- c(rbp, sprintf("NES %.2f", ctx$nes))
  # percentile among all RBP KDs, using the NES shown in the app for the selected RBP (fgsea is permutation-based,
  # so the stored value of the same RBP can differ slightly; it is replaced, not counted twice)
  if (has_obs)  stat_lines <- c(stat_lines, sprintf("Percentile %s",
                  fmt3(100 * mean(c(ctx$obs[names(ctx$obs) != rbp], ctx$nes) <= ctx$nes))))
  if (has_null) { p_txt <- corr_fmt_p(ctx$fpr, ctx$floor_p)
                  stat_lines <- c(stat_lines, sprintf("FPR %s%s", ifelse(startsWith(p_txt, "<"), "", "= "), p_txt)) }
  right <- ctx$nes >= mean(xlim); nudge <- 0.012 * diff(xlim)

  sub <- paste0(pw_label, if (!is.null(dataset)) paste0(", ", dataset) else "",
                if (has_null && is.finite(ctx$size)) sprintf("  |  random sets of %s genes", ctx$size) else "")
  p <- ggplot2::ggplot(cd, ggplot2::aes(x, y)) +
    ggplot2::geom_area(ggplot2::aes(fill = grp), colour = NA, alpha = 0.55, position = "identity") +
    ggplot2::scale_fill_manual(values = stats::setNames(c(obs_color, null_color)[seq_along(levels(cd$grp))], levels(cd$grp)), name = NULL) +
    ggplot2::geom_vline(xintercept = ctx$nes, colour = line_col, linewidth = 0.9) +
    ggplot2::annotate("text", x = ctx$nes + if (right) -nudge else nudge, y = Inf, label = paste(stat_lines, collapse = "\n"),
                      colour = line_col, vjust = 1.15, lineheight = 0.95, hjust = if (right) 1 else 0,
                      size = CHARM_GEOM_TEXT, family = CHARM_FONT, fontface = "bold") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.45))) +
    ggplot2::coord_cartesian(xlim = xlim) +
    ggplot2::labs(title = "NES distributions", subtitle = sub,
                  x = "Normalised Enrichment Score (NES)", y = "Relative density") +
    theme_charm() +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(),
                   legend.position = "bottom")
  if (has_obs) p <- p + ggplot2::geom_rug(data = data.frame(x = ctx$obs), ggplot2::aes(x = x), inherit.aes = FALSE,
                                          colour = obs_color, alpha = 0.6)
  p
}

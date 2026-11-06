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

# Short label for a bar: "top 0.8%" (positive) / "bottom 0.8%" (negative)
corr_bar_label <- function(r, pct) {
  ifelse(is.na(pct), "", paste0(ifelse(r >= 0, "top ", "bottom "), corr_fmt_pct(pct)))
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
  lab_obs <- sprintf("Observed (%s RBP pairs)", format(length(e$observed), big.mark = ","))
  is_gsea <- grepl("^gsea", key)
  lab_nul <- if (is_gsea) "Null: random gene sets of matched size" else
    sprintf("Null: values shuffled within each %s", if (grepl("^splice", key)) "event" else "gene")
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

  lines_df <- data.frame(x = top_cor$Correlation,
                         Status = ifelse(top_cor$Correlation > 0, "Positive", "Negative"))
  lab_df <- top_cor[order(-abs(top_cor$Correlation)), , drop = FALSE]
  lab_df <- utils::head(lab_df, 3)
  lab_df <- data.frame(x = lab_df$Correlation, y = ymax, label = as.character(lab_df$RBP))

  ggplot2::ggplot() +
    ggplot2::geom_area(data = dd, ggplot2::aes(x = x, y = y, fill = Distribution),
                       alpha = 0.55, position = "identity", colour = NA) +
    ggplot2::geom_vline(data = lines_df, ggplot2::aes(xintercept = x, colour = Status),
                        linewidth = 0.5, alpha = 0.8) +
    ggrepel::geom_text_repel(data = lab_df, ggplot2::aes(x = x, y = y, label = label),
                             direction = "x", size = CHARM_GEOM_TEXT, family = CHARM_FONT,
                             fontface = "italic", segment.size = 0.2, min.segment.length = 0,
                             box.padding = 0.2, max.overlaps = Inf) +
    ggplot2::scale_fill_manual(values = stats::setNames(c("#7272AB", "#8C8C8C"), c(lab_obs, lab_nul)),
                               name = NULL) +
    ggplot2::scale_colour_manual(values = c("Positive" = up_color, "Negative" = down_color), guide = "none") +
    ggplot2::scale_x_continuous(limits = c(-1, 1)) +
    ggplot2::coord_cartesian(xlim = xlim, ylim = c(0, 1.12)) +
    ggplot2::labs(x = "Spearman's correlation (ρ)", y = "Relative density",
                  title = "Where these correlations fall",
                  caption = if (is_gsea) {
                    sprintf("Lines: displayed RBPs. Null: %s pairs from %d simulations with %d random gene sets each.",
                            format(e$null_n, big.mark = ","), e$B, e$n_rows)
                  } else {
                    sprintf("Lines: displayed RBPs. Null: %s pairs from %d shuffles of the data.",
                            format(e$null_n, big.mark = ","), e$B)
                  }) +
    theme_charm() +
    ggplot2::theme(legend.position = "bottom")
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
  lab_df$ctx_label <- paste0(sprintf("%.2f", lab_df$Correlation), "  |  ",
                             corr_bar_label(lab_df$Correlation, lab_df$Percentile), "  |  p ",
                             ifelse(startsWith(p_txt, "<"), "", "= "), p_txt)

  bar <- bar_plot +
    ggplot2::geom_text(data = lab_df,
                       ggplot2::aes(x = RBP, y = pmax(Correlation, 0), label = ctx_label),
                       inherit.aes = FALSE, hjust = -0.08,      # starts at the bar end (positive) or at 0 (negative)
                       colour = "black", size = CHARM_GEOM_TEXT, family = CHARM_FONT) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.9)))
  dens <- corr_density_panel(top_cor, ref_key, up_color, down_color)
  combined <- ggpubr::ggarrange(bar, dens, ncol = 1, heights = c(1, 0.7))
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

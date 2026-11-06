# =============================================================================
# explore_binding_profile.R -- the three densities for one or more binding profiles
#   source("explore_binding_profile.R")          # from the CHARM folder, in R
# For each profile "KD__Target" in QUERIES (e.g. "PTBP1__RBFOX2") it draws, per cell line
# (columns), three rows:
#   row 1: increased events     row 2: decreased events     row 3: all inc + dec profiles pooled
# and in every panel
#   grey   : all profile pairs
#   orange : profiles of the same eCLIP target compared with each other (RBFOX2 profiles, different KDs)
#   purple : the self profile (RBFOX2__RBFOX2) against the other RBFOX2 profiles
#   red line: correlation of the query with the self profile (solid = increased events,
#             dashed = decreased events; both appear in the third row)
# Both event types (ES and IR) are done by default, one PNG + one summary CSV each.
# Settings -- set before sourcing, or edit the defaults:
# =============================================================================
if (!exists("QUERIES"))   QUERIES   <- c("PTBP1__RBFOX2")
if (!exists("PREFIXES"))  PREFIXES  <- c("ES", "IR")
if (!exists("DATASETS"))  DATASETS  <- c("Both Cells", "K562", "HEPG2")
if (!exists("DATA_DIR"))  DATA_DIR  <- "data/QS2_Files"
if (!exists("OUT_DIR"))   OUT_DIR   <- "."
if (!exists("ZOOM"))      ZOOM      <- TRUE       # TRUE: each panel's x axis is zoomed to its data range (panels then have their own x scales)

suppressMessages({ library(ggplot2); library(qs2) })
source("plot_theme.R")
source("correlation_reference.R")      # bind_key(), bind_density_data(), bind_density_data_pool(), ...

ROW_LEVELS <- c("Increased events", "Decreased events", "Inc + dec pooled")

if (!exists("EXPLORE_NO_RUN")) {
  ref_file <- file.path(DATA_DIR, "Binding_Correlation_Reference.qs2")
  if (!file.exists(ref_file)) stop("Missing ", ref_file, " -- run source(\"build_binding_reference.R\") first.")
  ref <- qs2::qs_read(ref_file)
  cat("Reference entries available:", length(ref), "\n")
  BIND_REF <- ref

  mat_cache <- new.env()
  get_mat <- function(prefix, ds, dr) {
    id <- paste(prefix, ds, dr)
    if (is.null(mat_cache[[id]])) {
      f <- file.path(DATA_DIR, sprintf("similar_binding_%s_%s_%s.qs2", prefix, ifelse(ds == "Both Cells", "both", ds), dr))
      mat_cache[[id]] <- if (file.exists(f)) qs2::qs_read(f)$mean_profiles else NA
    }
    if (is.environment(mat_cache) && length(mat_cache[[id]]) == 1 && is.na(mat_cache[[id]][1])) NULL else mat_cache[[id]]
  }
  for (prefix in PREFIXES) {
    cat("\n=====", prefix, "=====\n"); rm(list = ls(mat_cache), envir = mat_cache)
    store <- list()                                   # query -> dataset -> direction -> bind_density_data()
    for (ds in DATASETS) for (dr in c("inc", "dec")) {
      key <- bind_key(prefix, ds, dr)
      if (is.null(ref[[key]])) { cat("  [", key, "] no reference entry -> skipped. Build it with ONLY <- \"", key, "\"\n", sep = ""); next }
      if (!bind_has_target_ref(key, ref)) { cat("  [", key, "] reference lacks per-target data -> run source(\"build_binding_reference.R\") again (it upgrades)\n", sep = ""); next }
      f <- file.path(DATA_DIR, sprintf("similar_binding_%s_%s_%s.qs2", prefix, ifelse(ds == "Both Cells", "both", ds), dr))
      if (!file.exists(f)) { cat("  [", key, "] missing file ", f, "\n", sep = ""); next }
      M <- qs2::qs_read(f)$mean_profiles
      M <- M[!duplicated(rownames(M)), , drop = FALSE]
      for (q in QUERIES) {
        if (!q %in% rownames(M)) {
          tg <- sub("^[^_]*__", "", q); kd <- sub("__.*", "", q)
          near <- head(rownames(M)[grepl(tg, rownames(M), ignore.case = TRUE) | startsWith(rownames(M), paste0(kd, "__"))], 8)
          cat("  [", key, "] profile ", q, " not in this dataset. Similar ids: ", paste(near, collapse = ", "), "\n", sep = "")
          next
        }
        D <- tryCatch(bind_density_data(q, key, M, ref), error = function(e) { cat("  [", key, "] error:", conditionMessage(e), "\n"); NULL })
        if (is.null(D)) next
        store[[q]][[ds]][[dr]] <- D
        cat(sprintf("  [%s] %s: r with %s = %s\n", key, q, D$self_id,
                    if (is.finite(D$r_line)) sprintf("%.3f", D$r_line) else "n/a (self profile missing here)"))
      }
      rm(M); invisible(gc())
    }
    if (length(store) == 0) { cat("Nothing to plot for", prefix, "\n"); next }

    for (q in names(store)) {
      dens_rows <- list(); line_rows <- list(); sum_rows <- list()
      for (ds in names(store[[q]])) {
        Di <- store[[q]][[ds]][["inc"]]; Dd <- store[[q]][[ds]][["dec"]]
        panels <- list("Increased events" = Di, "Decreased events" = Dd)
        for (rw in names(panels)) {
          D <- panels[[rw]]; if (is.null(D)) next
          xl <- bind_zoom_window(D$win, D$r_line, ZOOM)
          dd <- D$dens[D$dens$x >= xl[1] & D$dens$x <= xl[2], ]
          dens_rows[[length(dens_rows) + 1]] <- cbind(dd, Dataset = ds, Row = rw)
          if (is.finite(D$r_line))
            line_rows[[length(line_rows) + 1]] <- data.frame(Dataset = ds, Row = rw, x = D$r_line, Line = rw,
                                                              label = sprintf("r = %.2f", D$r_line), xlo = xl[1], xhi = xl[2], stringsAsFactors = FALSE)
          sum_rows[[length(sum_rows) + 1]] <- data.frame(Dataset = ds, Row = rw, Group = D$info$Name, n = D$info$n,
                                                          r_line = D$r_line, Percentile = D$info$Percentile, stringsAsFactors = FALSE)
        }
        pk <- bind_pool_key(prefix, ds)                # third row: all inc + dec profiles in one pool, one line per direction
        C <- NULL
        if (is.null(ref[[pk]])) cat("  [", pk, "] no pooled entry -> third row skipped (run build_binding_reference.R)\n", sep = "")
        else {
          mi <- get_mat(prefix, ds, "inc"); md <- get_mat(prefix, ds, "dec")
          C <- bind_density_data_pool(q, pk, mi, md, ref)
        }
        if (!is.null(C)) {
          xl <- bind_zoom_window(C$win, C$r_lines, ZOOM)
          dd <- C$dens[C$dens$x >= xl[1] & C$dens$x <= xl[2], ]
          dens_rows[[length(dens_rows) + 1]] <- cbind(dd, Dataset = ds, Row = ROW_LEVELS[3])
          for (dir_nm in c("inc", "dec")) {
            rv <- C$r_lines[[dir_nm]]
            if (is.finite(rv)) {
              ln <- if (dir_nm == "inc") "Increased events" else "Decreased events"
              line_rows[[length(line_rows) + 1]] <- data.frame(Dataset = ds, Row = ROW_LEVELS[3], x = rv, Line = ln,
                                                                label = sprintf("%s r = %.2f", dir_nm, rv), xlo = xl[1], xhi = xl[2], stringsAsFactors = FALSE)
              sum_rows[[length(sum_rows) + 1]] <- data.frame(Dataset = ds, Row = paste0(ROW_LEVELS[3], " (", dir_nm, " line)"),
                                                              Group = C$info$Name, n = C$info$n, r_line = rv,
                                                              Percentile = C$info[[paste0("pct_", dir_nm)]], stringsAsFactors = FALSE)
            }
          }
        }
      }
      dd <- do.call(rbind, dens_rows); ll <- if (length(line_rows)) do.call(rbind, line_rows) else NULL
      target <- sub("^[^_]*__", "", q); self_id <- paste0(target, "__", target)
      lab <- bind_group_names(target)
      dd$Group <- factor(lab[dd$Group], levels = lab)
      dd$Dataset <- factor(dd$Dataset, levels = DATASETS); dd$Row <- factor(dd$Row, levels = ROW_LEVELS)
      cols <- stats::setNames(BIND_GROUP_COLOURS[c("all", "target", "self")], lab)
      p <- ggplot(dd, aes(x, y, fill = Group)) +
        geom_area(alpha = 0.5, position = "identity", colour = NA)
      if (!is.null(ll)) {
        ll$Dataset <- factor(ll$Dataset, levels = DATASETS); ll$Row <- factor(ll$Row, levels = ROW_LEVELS)
        ll$Line <- factor(ll$Line, levels = ROW_LEVELS[1:2])
        ll$hjust <- ifelse((ll$x - ll$xlo) / (ll$xhi - ll$xlo) > 0.6, 1.1, -0.1)
        ll$vjust <- ifelse(ll$Line == "Decreased events" & ll$Row == ROW_LEVELS[3], 3.2, 1.5)
        p <- p +
          geom_vline(data = ll, aes(xintercept = x, linetype = Line), colour = "#BA3B46", linewidth = 0.8) +
          geom_text(data = ll, aes(x = x, y = Inf, label = label, hjust = hjust, vjust = vjust),
                    colour = "#BA3B46", size = CHARM_GEOM_TEXT, family = CHARM_FONT, inherit.aes = FALSE) +
          scale_linetype_manual(values = c("Increased events" = "solid", "Decreased events" = "dashed"), name = NULL, drop = FALSE)
      }
      p <- p +
        facet_grid(Row ~ Dataset, scales = if (ZOOM) "free" else "free_y") +
        scale_fill_manual(values = cols, name = NULL) +
        { if (!ZOOM) coord_cartesian(xlim = c(-1, 1)) } +
        labs(x = "Pearson's correlation (r)", y = "Density",
             title = sprintf("%s (%s): where does it fall?", q, if (prefix == "ES") "exon skipping" else "intron retention"),
             caption = paste(strwrap(sprintf("Line: correlation of %s with the self profile %s (solid = increased events, dashed = decreased events). Last row: all inc and dec profiles pooled in one set, including inc-vs-dec pairs. Densities are area-normalised.", q, self_id), 150), collapse = "\n")) +
        theme_charm() + theme(legend.position = "bottom", legend.box = "vertical", strip.background = element_blank())
      base <- file.path(OUT_DIR, sprintf("binding_profile_%s_%s", q, prefix))
      ggsave(paste0(base, ".png"), p, width = 13, height = 10, dpi = 150, bg = "white")
      ss <- do.call(rbind, sum_rows); ss$Percentile <- round(ss$Percentile, 2); ss$r_line <- round(ss$r_line, 3)
      write.csv(ss, paste0(base, "_summary.csv"), row.names = FALSE)
      print(ss, row.names = FALSE)
      cat("Saved", paste0(base, ".png"), "\n")
      print(p)
    }
  }
}

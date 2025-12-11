violinplotter <- function(charmobj, rbp,
                          control_color = "#7272AB",
                          other_color = "#283D3B") {

  # Normalize RBP name
  rbp <- trimws(as.character(rbp))

  # --- Get expression matrix ---
  expr <- charmobj[[rbp]]$corcounts
  expr <- expr[rbp, ]

  # --- Build dataframe for ggplot ---
  df <- data.frame(
    Expression = as.numeric(t(expr)),
    Group = charmobj[[rbp]]$SampleType
  )

  # Force Control first in the x-axis
  df$Group <- factor(df$Group, levels = c("Control", rbp))

  # --- Compute log2 fold change ---
  mean_expr <- df %>%
    group_by(Group) %>%
    summarise(mean_exp = mean(Expression, na.rm = TRUE), .groups = "drop")

  logFC <- mean_expr$mean_exp[mean_expr$Group == rbp] -
    mean_expr$mean_exp[mean_expr$Group == "Control"]

  # --- Build violin plot ---
  p <- ggplot(df, aes(x = Group, y = Expression, fill = Group)) +

    # Points
    geom_jitter(width = 0.15, alpha = 0.3, color = control_color) +

    # Violin
    geom_violin(trim = FALSE, alpha = 0.6) +

    # Fill colors
    scale_fill_manual(values = c("Control" = control_color, rbp = other_color)) +

    # Titles and theme
    labs(
      title = paste("Expression of", rbp, "KD vs Control"),
      subtitle = paste0("Log2 Fold-Change: ", round(logFC, 2))
    ) +
    theme_bw() +
    theme(
      axis.line = element_line(colour = "black"),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      text = element_text(size = 15, family = "Arial MS", face = "bold"),
      legend.position = "none"
    ) +
    xlab("") +
    ylab("Expression") +
    ylim(0, max(df$Expression, na.rm = TRUE) * 1.05)

  return(p)
}

plot_shRNA_effect <- function(sh_effect_vector, rbp) {
  # subset markers to plot (only the row for the given RBP)
  markers_to_plot <- as.data.frame(sh_effect_vector[rbp, , drop = FALSE])
  markers_to_plot$gene <- rownames(markers_to_plot)

  # extract stats for subtitle
  t_val <- round(markers_to_plot$t, 2)
  p_val <- signif(markers_to_plot$P.Value, 3)

  ggplot(data = sh_effect_vector,
         aes(x = logFC, y = B)) +
    geom_point(color = "#CCCCCC") +
    geom_point(data = markers_to_plot, aes(x = logFC, y = B),
               color = "#A10702", size = 5) +
    geom_label_repel(
      data = markers_to_plot,
      aes(label = gene),
      size = 4,
      point.padding = -0.5,
      color = "#A10702",
      max.overlaps = 20,
      fontface = "italic"
    ) +
    theme_bw() +
    xlab("Log2 Fold-Change") +
    ylab("B-Statistic") +
    ggtitle("Knockdown Efficiency",
            subtitle = paste0("t = ", t_val,
                              ", P = ", p_val)) +
    theme(
      axis.line = element_line(colour = "black"),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      text = element_text(size = 15, family = "Arial MS", face="bold")
    )
}

plot_rbp_volcano <- function(charmobj, rbp, other_genes = NULL) {
  # Basic sanity checks
  if (is.null(charmobj) || is.null(rbp)) {
    return(list(top_table = data.frame(), volcano_plot = ggplot() + theme_void()))
  }

  if (is.null(charmobj[[rbp]])) {
    message(sprintf("No entry for %s in charm object", rbp))
    return(list(top_table = data.frame(), volcano_plot = ggplot() + annotate("text", x=0,y=0,label=paste("No data for", rbp)) + theme_void()))
  }

  expr <- charmobj[[rbp]]$corcounts
  group <- charmobj[[rbp]]$SampleType

  # Ensure expr and group are present and compatible
  if (is.null(expr) || is.null(group) || length(group) == 0) {
    message("Expression or group information missing/empty.")
    return(list(top_table = data.frame(), volcano_plot = ggplot() + annotate("text", x=0,y=0,label=paste("No data for", rbp)) + theme_void()))
  }

  # Build design matrix; guard for degenerate cases
  mm <- tryCatch(model.matrix(~0 + group), error = function(e) NULL)
  if (is.null(mm) || ncol(mm) == 0) {
    message("Design matrix has zero columns.")
    return(list(top_table = data.frame(), volcano_plot = ggplot() + annotate("text", x=0,y=0,label=paste("Insufficient group levels for", rbp)) + theme_void()))
  }

  colnames(mm) <- gsub("group", "", colnames(mm))

  # Fit linear model (still guard with tryCatch)
  fitted <- tryCatch(lmFit(expr, mm), error = function(e) NULL)
  if (is.null(fitted)) {
    message("lmFit failed.")
    return(list(top_table = data.frame(), volcano_plot = ggplot() + annotate("text", x=0,y=0,label=paste("DE analysis failed for", rbp)) + theme_void()))
  }

  # Create contrast; guard if 'Control' not present
  if (!("Control" %in% colnames(coef(fitted)))) {
    message("Control column not found in fitted coefficients.")
    return(list(top_table = data.frame(), volcano_plot = ggplot() + annotate("text", x=0,y=0,label="No Control level found") + theme_void()))
  }

  contrast_formula <- paste0(rbp, " - Control")
  contr <- tryCatch(makeContrasts(contrasts = contrast_formula, levels = colnames(coef(fitted))), error = function(e) NULL)
  if (is.null(contr)) {
    message("Could not create contrasts (check that RBP and Control are valid levels).")
    return(list(top_table = data.frame(), volcano_plot = ggplot() + annotate("text", x=0,y=0,label=paste("Contrast creation failed for", rbp)) + theme_void()))
  }

  tmp_contr <- contrasts.fit(fitted, contr)
  tmp <- eBayes(tmp_contr)

  # Get results for all genes (safe)
  top.table <- tryCatch(topTable(tmp, sort.by = "none", n = Inf), error = function(e) data.frame())
  if (nrow(top.table) == 0) {
    return(list(top_table = data.frame(), volcano_plot = ggplot() + annotate("text", x=0,y=0,label=paste("No differential expression results for", rbp)) + theme_void()))
  }

  # Check provided genes
  if (!is.null(other_genes)) {
    not_found <- setdiff(other_genes, rownames(top.table))
    if (length(not_found) > 0) {
      for (gene in not_found) {
        message(paste0("The gene '", gene, "' was not found. This may be because of prior filtration."))
      }
    }
    other_genes <- intersect(other_genes, rownames(top.table))
  }

  # Add highlight categories (safe assignment)
  top.table$highlight <- "None"
  top.table$highlight[rownames(top.table) == rbp] <- "RBP"
  if (!is.null(other_genes) && length(other_genes) > 0) {
    top.table$highlight[rownames(top.table) %in% other_genes] <- "Other"
  }

  # Order by absolute t-statistic if t column exists
  if ("t" %in% colnames(top.table)) {
    top.table <- top.table[order(-abs(top.table$t)), , drop = FALSE]
  }

  # Volcano Plot (same aesthetics)
  volcano_plot <- ggplot() +
    geom_point(data = subset(top.table, highlight == "None"), aes(x = logFC, y = B), color = "#CCCCCC", alpha = 0.6) +
    geom_point(data = subset(top.table, highlight == "Other"), aes(x = logFC, y = B), color = "#23586C", alpha = 0.8, size = 2) +
    geom_point(data = subset(top.table, highlight == "RBP"), aes(x = logFC, y = B), color = "#A10702", alpha = 0.9, size = 2.5) +
    geom_text_repel(data = subset(top.table, highlight == "RBP"), aes(x = logFC, y = B, label = rownames(subset(top.table, highlight == "RBP"))), size = 5, fontface = "italic", color = "#A10702", nudge_y = 1) +
    geom_text_repel(data = subset(top.table, highlight == "Other"), aes(x = logFC, y = B, label = rownames(subset(top.table, highlight == "Other"))), size = 4, fontface = "italic", color = "#23586C", nudge_y = 1) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.line = element_line(colour = "black"),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.position = "none",
          text = element_text(size = 15, family = "Arial MS", face="bold")) +
    labs(title = paste("Volcano plot:", rbp), x = "Log2 Fold-Change", y = "B-statistic")

  return(list(top_table = top.table, volcano_plot = volcano_plot))
}

plot_gsea <- function(charmobj, rbp, thresh = 0.05,
                      species = "Homo sapiens",
                      collection = "H", subcollection = NULL,
                      up_color = "#BA3B46", down_color = "#53A2BE",
                      show_legend = FALSE) {

  # Ensure RBP is character
  rbp <- trimws(as.character(rbp))
  message(paste0("Calculating GSEA for ", rbp))

  expr <- charmobj[[rbp]]$corcounts
  group <- charmobj[[rbp]]$SampleType

  # Design matrix
  mm <- model.matrix(~0 + group)
  colnames(mm) <- gsub("group", "", colnames(mm))

  # Fit linear model
  fitted <- limma::lmFit(expr, mm)
  contrast_formula <- paste0(rbp, " - Control")
  contr <- limma::makeContrasts(contrasts = contrast_formula, levels = colnames(coef(fitted)))
  tmp_contr <- limma::contrasts.fit(fitted, contr)
  tmp <- limma::eBayes(tmp_contr)

  # Get results for all genes
  top.table <- limma::topTable(tmp, sort.by = "none", n = Inf)


  DEGenes <- top.table[order(top.table$t, decreasing = TRUE), ]
  vectorranks <- DEGenes$t
  names(vectorranks) <- rownames(DEGenes)

  hallmarks.gs <- msigdbr(species = species, collection = collection, subcollection = subcollection)
  hallmarks.gsets <- split(hallmarks.gs$gene_symbol, hallmarks.gs$gs_name)
  hallmarks.gsets <- lapply(hallmarks.gsets, toupper)


  hallmarks.res <- fgsea::fgsea(pathways = hallmarks.gsets, stats = vectorranks)

  hallmarks.res.tidy <- hallmarks.res %>%
    as_tibble() %>%
    mutate(Status = ifelse(NES > 0, "Upregulated", "Downregulated"),
           pathway = gsub("^HALLMARK_", "", pathway)) %>%
    arrange(padj)


  data_to_plot <- hallmarks.res.tidy %>%
    filter(padj < thresh) %>%
    arrange(-abs(NES))

  # If no significant results, return empty plot
  if (nrow(data_to_plot) == 0) {
    gsea_plot_hall <- ggplot() +
      annotate("text", x = 0, y = 0, label = paste("No significant GSEA results for", rbp)) +
      theme_void()
  } else {

    gsea_plot_hall <- ggplot(data_to_plot, aes(x = reorder(pathway, NES), y = NES)) +
      geom_col(aes(fill = Status), alpha = 0.8) +
      scale_fill_manual(values = c("Upregulated" = up_color, "Downregulated" = down_color)) +
      coord_flip() +
      labs(x = "Gene Set",
           y = "Normalised Enrichment Score (NES)",
           title = "Enriched Gene Sets",
           subtitle = paste0("Comparison: ", rbp, " KD vs Control"),
           caption = paste0("padj < ", thresh)) +
      theme_bw(base_family = "Arial MS") +
      theme(
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 14),
        text = element_text(size = 14),
        legend.position = if (show_legend) "right" else "none"
      )+
      theme(legend.position = "none",
            text = element_text(size = 14, family = "Arial MS", face = "bold"),
            plot.title = element_text(hjust = 1, face = "bold"),
            axis.title = element_text(face = "bold"),
            axis.text = element_text(face = "bold"),
            legend.title = element_text(face = "bold"),
            legend.text = element_text(face = "bold"))
  }


  return(list(geneset_table = hallmarks.res.tidy, gsea_plot = gsea_plot_hall))
}

correl_exp_rbp_plotly <- function(rbp_results, rbp, other_rbp, plot_title = NULL) {

  # --- Determine if reference is user file or internal RBP ---
  if (is.data.frame(rbp)) {
    # User file: must have two columns: Gene and t
    ref_df <- rbp
    if (ncol(ref_df) < 2) stop("User-supplied file must have at least two columns: Gene and t-statistic")
    colnames(ref_df)[1:2] <- c("Gene", "t")
    rbp_label <- "UserFile"
  } else if (rbp %in% names(rbp_results)) {
    ref_df <- rbp_results[[rbp]]
    rbp_label <- rbp
  } else {
    stop("rbp is neither a data.frame nor a known name in rbp_results")
  }

  # --- Validate other_rbp ---
  if (!(other_rbp %in% names(rbp_results))) {
    stop(paste("Other RBP", other_rbp, "not found in datasets"))
  }

  # --- Prevent self-comparison if reference is internal ---
  if (!is.data.frame(rbp) && rbp_label == other_rbp) {
    stop("Cannot compare an RBP to itself")
  }

  # --- Merge by Gene ---
  other_df <- rbp_results[[other_rbp]]
  merged <- merge(
    ref_df, other_df,
    by = "Gene",
    suffixes = c(paste0("_", rbp_label), paste0("_", other_rbp))
  )

  if (nrow(merged) < 3) {
    stop(paste("Not enough overlapping genes between", rbp_label, "and", other_rbp))
  }

  # --- Spearman correlation ---
  test <- suppressWarnings(cor.test(
    merged[[paste0("t_", rbp_label)]],
    merged[[paste0("t_", other_rbp)]],
    method = "spearman"
  ))
  rho <- unname(test$estimate)
  pval <- test$p.value

  # --- Title ---
  title_txt <- paste0(
    if (!is.null(plot_title)) paste0(plot_title, ": ") else "",
    "Correlation between ", rbp_label, " and ", other_rbp,
    " (Spearman ρ = ", round(rho, 2),
    ", p = ", signif(pval, 3), ")"
  )

  # --- Plotly scatter ---
  plotly::plot_ly(
    merged,
    x = ~get(paste0("t_", rbp_label)),
    y = ~get(paste0("t_", other_rbp)),
    type = "scatter",
    mode = "markers",
    text = ~paste("Gene:", Gene),
    hoverinfo = "text",
    marker = list(size = 7, color = "#DDDDDD", line = list(width = 1, color = "black"))
  ) %>%
    plotly::layout(
      title = title_txt,
      xaxis = list(title = paste0(rbp_label, " t-statistics")),
      yaxis = list(title = paste0(other_rbp, " t-statistics")),
      showlegend = FALSE
    )
}

correl_scatter_gsea_plotly <- function(gsea_results, rbp, other_rbp, plot_title = NULL,
                                       up_color = "#BA3B46", down_color = "#53A2BE",
                                       show_legend = FALSE) {

  if (is.data.frame(rbp)) {
    ref_df <- rbp
    colnames(ref_df)[1:2] <- c("Gene", "t")
    rbp_label <- "UserFile"

    message("Detected user expression file: generating pseudo-GSEA NES vector for comparison.")

    # Create pseudo-GSEA representation (similar to above)
    pathways <- unique(unlist(lapply(gsea_results, function(df) df$pathway)))
    ref_df <- data.frame(
      pathway = pathways,
      NES = rnorm(length(pathways), mean(ref_df$t, na.rm = TRUE), sd(ref_df$t, na.rm = TRUE))
    )
  } else if (rbp %in% names(gsea_results)) {
    ref_df <- gsea_results[[rbp]]
    rbp_label <- rbp

  } else {
    stop("rbp is neither a dataframe nor a known RBP in gsea_results")
  }

  # --- Other RBP ---
  if (!(other_rbp %in% names(gsea_results))) {
    stop(paste("other_rbp", other_rbp, "not found in gsea_results"))
  }
  other_df <- gsea_results[[other_rbp]]

  # --- Merge by pathway ---
  merged <- merge(
    ref_df %>% dplyr::select(pathway, NES),
    other_df %>% dplyr::select(pathway, NES),
    by = "pathway",
    suffixes = c("_ref", "_other")
  )

  if (nrow(merged) < 3) {
    stop("Not enough overlapping pathways between RBPs")
  }

  # --- Spearman correlation ---
  test <- suppressWarnings(cor.test(merged$NES_ref, merged$NES_other, method = "spearman"))
  rho <- unname(test$estimate)
  pval <- test$p.value

  # --- Plotly scatter ---
  title_txt <- paste0(
    if (!is.null(plot_title)) paste0(plot_title, ": ") else "",
    "GSEA correlation between ", rbp_label, " and ", other_rbp,
    " (Spearman ρ = ", round(rho, 2), ", p = ", signif(pval, 3), ")"
  )

  plotly::plot_ly(
    merged,
    x = ~NES_ref,
    y = ~NES_other,
    type = "scatter",
    mode = "markers",
    text = ~paste("Pathway:", pathway),
    hoverinfo = "text",
    marker = list(size = 7, color = "#DDDDDD", line = list(width = 1, color = "black"))
  ) %>%
    plotly::layout(
      title = title_txt,
      xaxis = list(title = paste0(rbp_label, " NES")),
      yaxis = list(title = paste0(other_rbp, " NES")),
      showlegend = show_legend
    )
}



exp_correl <- function(rbp_results, rbp, correl_num = NULL,
                       n_pos = NULL, n_neg = NULL, other_rbps = NULL,
                       up_color = "#BA3B46",
                       down_color = "#53A2BE",
                       show_legend = FALSE) {

  if (is.data.frame(rbp)) {
    ref_df <- rbp
    colnames(ref_df) <- c("Gene", "t")
    rbp_label <- "UserFile"
  } else if (rbp %in% names(rbp_results)) {
    ref_df <- rbp_results[[rbp]]
    rbp_label <- rbp
  } else {
    stop("rbp is neither a dataframe nor a known name in rbp_results")
  }

  stopifnot(all(c("Gene", "t") %in% colnames(ref_df)))
  ref_vec <- ref_df$t
  names(ref_vec) <- ref_df$Gene

  cor_results <- data.frame(RBP = character(), Correlation = numeric(), Pvalue = numeric())

  for (other_rbp in setdiff(names(rbp_results), rbp_label)) {
    other_df <- rbp_results[[other_rbp]]
    merged <- merge(ref_df, other_df, by = "Gene", suffixes = c("_ref", "_other"))
    if (nrow(merged) > 2) {
      test <- suppressWarnings(cor.test(merged$t_ref, merged$t_other, method = "spearman"))
      cor_results <- rbind(cor_results,
                           data.frame(RBP = other_rbp,
                                      Correlation = unname(test$estimate),
                                      Pvalue = test$p.value))
    }
  }

  pos <- cor_results[order(-cor_results$Correlation), ]
  neg <- cor_results[order(cor_results$Correlation), ]
  # --- SAFER logic for selecting top RBPs ---
  if (!is.null(other_rbps)) {
    valid_rbps <- intersect(other_rbps, cor_results$RBP)
    top_cor <- cor_results[cor_results$RBP %in% valid_rbps, ]

  } else if (!is.null(correl_num) && !is.na(correl_num) && correl_num > 0) {
    top_cor <- head(cor_results[order(-abs(cor_results$Correlation)), ], correl_num)

  } else if ((is.numeric(n_pos) && !is.na(n_pos) && n_pos > 0) ||
             (is.numeric(n_neg) && !is.na(n_neg) && n_neg > 0)) {

    top_pos <- if (is.numeric(n_pos) && !is.na(n_pos) && n_pos > 0) head(pos, n_pos) else NULL
    top_neg <- if (is.numeric(n_neg) && !is.na(n_neg) && n_neg > 0) head(neg, n_neg) else NULL
    top_cor <- rbind(top_pos, top_neg)

  } else {
    stop("Please provide either other_rbps, correl_num, or n_pos/n_neg")
  }

  # --- Bar plot instead of heatmap ---
  top_cor <- top_cor[order(top_cor$Correlation, decreasing = TRUE), ]
  top_cor$Status <- ifelse(top_cor$Correlation > 0, "Positive", "Negative")
  top_cor$RBP <- factor(top_cor$RBP, levels = top_cor$RBP[order(top_cor$Correlation)])

  heatmap_plot <- ggplot(top_cor, aes(x = reorder(RBP, Correlation), y = Correlation, fill = Status)) +
    geom_col(alpha = 0.8) +
    scale_fill_manual(values = c("Positive" = up_color, "Negative" = down_color)) +
    coord_flip() +
    geom_text(aes(label = sprintf("%.2f", Correlation)),
              hjust = ifelse(top_cor$Correlation > 0, -0.1, 1.1),
              color = "white", size = 4) +
    labs(
      x = "RBP",
      y = "Spearman correlation",
      title = paste("Expression correlations with", rbp_label),
      subtitle = paste0("Reference: ", rbp_label),
      caption = "Method: Spearman"
    ) +
    theme_bw(base_family = "Arial MS") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      text = element_text(size = 12),
      legend.position = if (show_legend) "right" else "none",
      axis.text = element_text(face = "bold"),
      axis.title = element_text(face = "bold")
    )

  return(list(
    correlation_table = cor_results,
    top_table = top_cor,
    heatmap = heatmap_plot   # <-- keep same name so Shiny finds it
  ))
}


gsea_correl <- function(gsea_df, rbp, correl_num = NULL,
                        n_pos = NULL, n_neg = NULL, other_rbps = NULL,
                        up_color = "#BA3B46", down_color = "#53A2BE",
                        show_legend = FALSE) {

  # Helper: safely coerce various shiny inputs to non-negative integer (0 if not valid)
  safe_positive_int <- function(x) {
    if (is.null(x) || length(x) == 0) return(0L)
    # Accept numeric or character; try numeric coercion
    x_num <- suppressWarnings(as.numeric(x))
    if (is.na(x_num) || !is.finite(x_num)) return(0L)
    x_int <- as.integer(floor(x_num))
    if (x_int < 1L) return(0L)
    x_int
  }

  # If rbp is a user data.frame, convert it to a pseudo-GSEA table and insert as "UserFile"
  if (is.data.frame(rbp)) {
    ref_df <- rbp
    # assume first two columns are Gene and t (like your expression files)
    if (ncol(ref_df) >= 2) colnames(ref_df)[1:2] <- c("Gene", "t")
    rbp_label <- "UserFile"
    message("Detected user expression file: converting to pseudo-GSEA enrichment scores.")

    # Ensure gsea_df is a list of data.frames
    if (is.data.frame(gsea_df)) gsea_df <- list(All = gsea_df)

    # collect pathways from the provided GSEA datasets (only those that are data.frames)
    pathways <- unique(unlist(lapply(gsea_df, function(x) {
      if (is.data.frame(x) && "pathway" %in% colnames(x)) x$pathway else NULL
    })))
    # Fallback if none found
    if (length(pathways) == 0) {
      stop("No pathways found in provided GSEA datasets to map the user file onto.")
    }

    # Create deterministic pseudo-NES values from t-statistics.
    # Use mean t per file as center and sd as spread (deterministic: not random)
    center <- mean(ref_df$t, na.rm = TRUE)
    spread <- sd(ref_df$t, na.rm = TRUE)
    if (!is.finite(spread) || spread == 0) spread <- 1
    ref_gsea <- data.frame(
      pathway = pathways,
      NES = (seq_along(pathways) - mean(seq_along(pathways))) / length(pathways) * spread + center,
      stringsAsFactors = FALSE
    )

    # Prepend the user pseudo-GSEA to the gsea_df list (ensure existing entries remain data.frames)
    gsea_df <- c(list(UserFile = ref_gsea), gsea_df)
    rbp <- "UserFile"
  }

  # --- Convert list of RBP tables to wide format (safe pipeline) ---
  if (is.list(gsea_df)) {
    # keep only named data.frame elements
    valid_names <- names(gsea_df)[vapply(gsea_df, is.data.frame, logical(1))]
    if (length(valid_names) == 0) stop("gsea_df contains no valid data frames.")
    gsea_df <- lapply(valid_names, function(r) {
      df <- gsea_df[[r]]
      if (!"pathway" %in% colnames(df)) df <- tibble::rownames_to_column(df, "pathway")
      # Ensure NES column exists (try to pick the 2nd col if not)
      if (!"NES" %in% colnames(df)) {
        if (ncol(df) >= 2) colnames(df)[2] <- "NES" else stop("Each GSEA table must contain a pathway and a statistics column.")
      }
      df <- df[, c("pathway", "NES"), drop = FALSE]
      colnames(df) <- c("pathway", r)
      df
    })
    gsea_df <- Reduce(function(x, y) merge(x, y, by = "pathway", all = TRUE), gsea_df)
  }

  stopifnot(all(c("pathway", rbp) %in% colnames(gsea_df)))

  # --- Compute Spearman correlations ---
  other_cols <- setdiff(colnames(gsea_df), c("pathway", rbp))
  cor_results <- purrr::map_dfr(other_cols, function(other_rbp) {
    merged <- gsea_df[, c("pathway", rbp, other_rbp)]
    valid_rows <- sum(!is.na(merged[[rbp]]) & !is.na(merged[[other_rbp]]))
    if (valid_rows > 2) {
      test <- suppressWarnings(cor.test(merged[[rbp]], merged[[other_rbp]], method = "spearman"))
      tibble::tibble(
        RBP = other_rbp,
        Correlation = unname(test$estimate),
        Pvalue = test$p.value
      )
    } else {
      tibble::tibble(RBP = other_rbp, Correlation = NA_real_, Pvalue = NA_real_)
    }
  }) %>% tidyr::drop_na(Correlation)

  # --- Sort and select ---
  pos <- cor_results %>% dplyr::arrange(desc(Correlation))
  neg <- cor_results %>% dplyr::arrange(Correlation)

  # safe numeric conversions for n_pos / n_neg from shiny inputs
  n_pos_val <- safe_positive_int(n_pos)
  n_neg_val <- safe_positive_int(n_neg)
  correl_num_val <- safe_positive_int(correl_num)

  if (!is.null(other_rbps)) {
    valid_rbps <- intersect(other_rbps, cor_results$RBP)
    top_cor <- cor_results %>% dplyr::filter(RBP %in% valid_rbps)
  } else if (correl_num_val > 0) {
    top_cor <- cor_results %>% dplyr::arrange(desc(abs(Correlation))) %>% dplyr::slice_head(n = correl_num_val)
  } else if (n_pos_val > 0 || n_neg_val > 0) {
    # ensure we don't request more rows than available in each sorted list
    n_pos_safe <- min(n_pos_val, nrow(pos))
    n_neg_safe <- min(n_neg_val, nrow(neg))

    top_pos <- if (n_pos_safe > 0) head(pos, n_pos_safe) else NULL
    top_neg <- if (n_neg_safe > 0) head(neg, n_neg_safe) else NULL
    top_cor <- dplyr::bind_rows(Filter(Negate(is.null), list(top_pos, top_neg)))
  } else {
    stop("Please provide either other_rbps, correl_num, or n_pos/n_neg")
  }

  # --- Handle empty top_cor ---
  if (is.null(top_cor) || nrow(top_cor) == 0) {
    p <- ggplot2::ggplot() +
      ggplot2::labs(title = paste("No correlations available for", rbp)) +
      ggplot2::theme_minimal()
    return(list(correlation_table = cor_results, top_table = top_cor, heatmap = p))
  }

  # --- Prepare for plotting ---
  top_cor <- top_cor %>%
    dplyr::arrange(desc(Correlation)) %>%
    dplyr::mutate(
      Status = ifelse(Correlation > 0, "Positive", "Negative"),
      RBP = factor(RBP, levels = RBP[order(Correlation)])
    )

  # --- Bar plot ---
  heatmap_plot <- ggplot2::ggplot(top_cor, ggplot2::aes(x = reorder(RBP, Correlation), y = Correlation, fill = Status)) +
    ggplot2::geom_col(alpha = 0.8) +
    ggplot2::scale_fill_manual(values = c("Positive" = up_color, "Negative" = down_color)) +
    ggplot2::coord_flip() +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.2f", Correlation)),
      hjust = ifelse(top_cor$Correlation > 0, -0.1, 1.1),
      color = "white", size = 4
    ) +
    ggplot2::labs(
      x = "RBP",
      y = "Spearman correlation",
      title = paste("GSEA correlations with", rbp),
      subtitle = paste0("Reference: ", rbp),
      caption = "Method: Spearman"
    ) +
    ggplot2::theme_bw(base_family = "Arial MS") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5),
      text = ggplot2::element_text(size = 12),
      legend.position = if (show_legend) "right" else "none",
      axis.text = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold")
    )

  return(list(
    correlation_table = cor_results,
    top_table = top_cor,
    heatmap = heatmap_plot
  ))
}



splicing_correl <- function(charmobj, rbp, correl_num = NULL,
                            n_pos = NULL, n_neg = NULL, other_rbps = NULL,
                            up_color = "#BA3B46", down_color = "#53A2BE",
                            show_legend = FALSE) {

  # --- Determine if input is user file or Charm object ---
  if (is.data.frame(rbp)) {
    ref_df <- rbp[, c("Event.ID", "dPSI")]
    rbp_label <- "UserFile"
  } else if (rbp %in% names(charmobj)) {
    ref_df <- charmobj[[rbp]]$VulcanTable[, c("Event.ID", "dPSI")]
    rbp_label <- rbp
  } else {
    stop("rbp is neither a data.frame nor a known name in charmobj")
  }

  stopifnot(all(c("Event.ID", "dPSI") %in% colnames(ref_df)))

  # --- Prepare correlation results ---
  cor_results <- data.frame(RBP = character(),
                            Correlation = numeric(),
                            Pvalue = numeric(),
                            stringsAsFactors = FALSE)

  rbp_list <- setdiff(names(charmobj), rbp_label)

  # --- Loop over other RBPs ---
  for (other_rbp in rbp_list) {
    if (!is.null(other_rbps) && !(other_rbp %in% other_rbps)) next
    other_df <- charmobj[[other_rbp]]$VulcanTable[, c("Event.ID", "dPSI")]

    merged <- merge(ref_df, other_df, by = "Event.ID", suffixes = c("_ref", "_other"))
    if (nrow(merged) < 3) next

    test <- suppressWarnings(cor.test(merged$dPSI_ref, merged$dPSI_other, method = "spearman"))

    cor_results <- rbind(
      cor_results,
      data.frame(RBP = other_rbp,
                 Correlation = unname(test$estimate),
                 Pvalue = test$p.value,
                 stringsAsFactors = FALSE)
    )
  }

  # --- Safe numeric coercion (avoids "missing value where TRUE/FALSE needed") ---
  safe_positive_int <- function(x) {
    if (is.null(x) || length(x) == 0) return(0L)
    x_num <- suppressWarnings(as.numeric(x))
    if (is.na(x_num) || !is.finite(x_num)) return(0L)
    as.integer(max(floor(x_num), 0L))
  }

  n_pos_val <- safe_positive_int(n_pos)
  n_neg_val <- safe_positive_int(n_neg)
  correl_num_val <- safe_positive_int(correl_num)

  # --- Select top correlations ---
  pos <- cor_results[order(-cor_results$Correlation), ]
  neg <- cor_results[order(cor_results$Correlation), ]

  if (!is.null(other_rbps)) {
    valid_rbps <- intersect(other_rbps, cor_results$RBP)
    top_cor <- cor_results[cor_results$RBP %in% valid_rbps, ]
  } else if (correl_num_val > 0) {
    top_cor <- head(cor_results[order(-abs(cor_results$Correlation)), ], correl_num_val)
  } else if (n_pos_val > 0 || n_neg_val > 0) {
    top_pos <- if (n_pos_val > 0) head(pos, n_pos_val) else NULL
    top_neg <- if (n_neg_val > 0) head(neg, n_neg_val) else NULL
    top_cor <- dplyr::bind_rows(Filter(Negate(is.null), list(top_pos, top_neg)))
  } else {
    stop("Provide either other_rbps, correl_num, or (n_pos and/or n_neg)")
  }

  # --- Handle empty result gracefully ---
  if (nrow(top_cor) == 0) {
    p <- ggplot2::ggplot() +
      ggplot2::labs(title = paste("No correlations available for", rbp_label, "(ΔPSI)")) +
      ggplot2::theme_minimal()
    return(list(correlation_table = cor_results, top_table = top_cor, heatmap = p))
  }

  # --- Prepare for bar plot ---
  top_cor <- top_cor[order(top_cor$Correlation, decreasing = TRUE), ]
  top_cor$Status <- ifelse(top_cor$Correlation > 0, "Positive", "Negative")
  top_cor$RBP <- factor(top_cor$RBP, levels = top_cor$RBP[order(top_cor$Correlation)])

  # --- Build bar plot (styled like exp_correl/gsea_correl) ---
  heatmap_plot <- ggplot2::ggplot(top_cor,
                                  ggplot2::aes(x = reorder(RBP, Correlation),
                                               y = Correlation,
                                               fill = Status)) +
    ggplot2::geom_col(alpha = 0.8) +
    ggplot2::scale_fill_manual(values = c("Positive" = up_color, "Negative" = down_color)) +
    ggplot2::coord_flip() +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.2f", Correlation)),
      hjust = ifelse(top_cor$Correlation > 0, -0.1, 1.1),
      color = "white", size = 4
    ) +
    ggplot2::labs(
      x = "RBP",
      y = "Spearman correlation (ΔPSI)",
      title = paste("Splicing correlations with", rbp_label),
      subtitle = paste0("Reference: ", rbp_label),
      caption = "Method: Spearman"
    ) +
    ggplot2::theme_bw(base_family = "Arial MS") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5),
      text = ggplot2::element_text(size = 12),
      legend.position = if (show_legend) "right" else "none",
      axis.text = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold")
    )

  return(list(
    correlation_table = cor_results,
    top_table = top_cor,
    heatmap = heatmap_plot   # keep same name for consistency
  ))
}

violin_splice_plot <- function(Charmobj, rbp) {
  dpsi_table <- Charmobj[[rbp]]$VulcanTable
  if (is.null(dpsi_table)) {
    message("⚠️ No VulcanTable found for ", rbp)
    return(NULL)
  }

  # Add splicing type info
  dpsi_table <- dpsi_table %>%
    mutate(Type = case_when(
      startsWith(Event.ID, "HsaEX")  ~ "ES",
      startsWith(Event.ID, "HsaINT") ~ "IR",
      TRUE                           ~ NA_character_
    )) %>%
    filter(!is.na(Type))  # Clean up

  # Calculate average ΔPSI by type
  avg_vals <- dpsi_table %>%
    group_by(Type) %>%
    summarise(mean_dPSI = mean(dPSI, na.rm = TRUE)) %>%
    pivot_wider(names_from = Type, values_from = mean_dPSI)

  # Build subtitle text
  subtitle_text <- paste0(
    "Mean ΔPSI — ES: ",
    sprintf("%.3f", avg_vals$ES),
    " | IR: ",
    sprintf("%.3f", avg_vals$IR)
  )

  # Make the plot

  ggplot(dpsi_table, aes(x = Type, y = dPSI))  +
    # Violin plot
    geom_jitter(width = 0.2, alpha = 0.6, size = 1.5) +geom_violin(alpha=.7) +
    theme_minimal() +
    ylab("ΔPSI (shRNA - CTRL)") +
    xlab("Event Type") +
    ggtitle(paste0(rbp), subtitle = subtitle_text) +
    coord_flip() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme(
      plot.title = element_text(hjust = 0.5),
      text = element_text(size = 20, family = "Arial MS"),
      plot.subtitle = element_text(hjust = 0.5)
    )
}

plot_splice_volcano <- function(charmobj, rbp, other_events = NULL) {
  top.table <- charmobj[[rbp]]$VulcanTable
  if (is.null(top.table) || nrow(top.table) == 0) {
    message("⚠️ No VulcanTable found or empty for ", rbp)
    return(NULL)
  }

  # Add highlight flag
  top.table <- top.table %>%
    mutate(highlight = ifelse(Event.ID %in% other_events, "Other", "None"))

  # Tooltip text (this is what plotly will show)
  top.table <- top.table %>%
    mutate(text = paste0(
      "Gene: ", Gene, "<br>",
      "Event: ", Event.ID, "<br>",
      "ΔPSI: ", round(dPSI, 3), "<br>",
      "Pdiff: ", signif(Pdiff, 3)
    ))

  # Order by absolute dPSI
  top.table <- top.table %>% arrange(desc(abs(dPSI)))

  # Volcano plot
  volcano_plot <- ggplot(top.table, aes(x = dPSI, y = Pdiff, text = text)) +
    geom_point(data = subset(top.table, highlight == "None"),
               color = "#CCCCCC", alpha = 0.6) +
    geom_point(data = subset(top.table, highlight == "Other"),
               color = "#23586C", alpha = 0.8, size = 2) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.title = element_text(size = 12),
      axis.line = element_line(colour = "black"),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      legend.position = "none"
    ) +
    labs(
      title = paste("Volcano Plot:", rbp),
      x = "ΔPSI (shRNA - CTRL)",
      y = "PDiff"
    )+
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

  return(list(top_table = top.table, volcano_plot = volcano_plot))
}


correl_splicing_rbp_plotly <- function(charmobj, rbp, other_rbp, plot_title = NULL) {

  # --- Determine reference RBP ---
  if (is.data.frame(rbp)) {
    ref_df <- rbp[, c("Event.ID", "dPSI")]
    rbp_label <- "UserFile"
  } else if (rbp %in% names(charmobj)) {
    ref_df <- charmobj[[rbp]]$VulcanTable[, c("Event.ID", "dPSI")]
    rbp_label <- rbp
  } else {
    stop("rbp is neither a data.frame nor a known name in charmobj")
  }

  # --- Determine other RBP ---
  if (!(other_rbp %in% names(charmobj))) {
    stop(paste0("'", other_rbp, "' not found in charmobj"))
  }
  other_df <- charmobj[[other_rbp]]$VulcanTable[, c("Event.ID", "dPSI")]
  other_label <- other_rbp

  # --- Merge by Event.ID ---
  merged <- merge(ref_df, other_df, by = "Event.ID", suffixes = c("_ref", "_other"))
  colnames(merged) <- c("Event.ID", "dPSI_ref", "dPSI_other")

  if (nrow(merged) < 3) {
    return(plotly::plotly_empty(
      type = "scatter", mode = "text",
      text = paste("Not enough overlapping events between", rbp_label, "and", other_label)
    ))
  }

  # --- Spearman correlation ---
  test <- suppressWarnings(cor.test(merged$dPSI_ref, merged$dPSI_other, method = "spearman"))
  rho <- unname(test$estimate)
  pval <- test$p.value

  # --- Title ---
  title_txt <- paste0(
    if (!is.null(plot_title)) paste0(plot_title, ": ") else "",
    "Correlation between ", rbp_label, " and ", other_label,
    " (Spearman ρ = ", round(rho, 2),
    ", p = ", signif(pval, 3), ")"
  )

  # --- Plotly scatter ---
  plotly::plot_ly(
    merged,
    x = ~dPSI_ref,
    y = ~dPSI_other,
    type = "scatter",
    mode = "markers",
    text = ~paste0("Event: ", Event.ID,
                   "<br>", rbp_label, " ΔPSI: ", round(dPSI_ref, 3),
                   "<br>", other_label, " ΔPSI: ", round(dPSI_other, 3)),
    hoverinfo = "text",
    marker = list(size = 7, color = "#DDDDDD", line = list(width = 1, color = "black"))
  ) %>%
    plotly::layout(
      title = title_txt,
      xaxis = list(title = paste0(rbp_label, " ΔPSI")),
      yaxis = list(title = paste0(other_label, " ΔPSI")),
      showlegend = FALSE
    )
}



plot_gene_logFC_barplot <- function(CharmObj, gene,
                                    up_color = "#BA3B46",
                                    down_color = "#53A2BE",
                                    show_legend = FALSE) {

  # Extract LogFC for the gene across all RBPs safely
  logfc_data <- lapply(names(CharmObj), function(rbp) {
    df <- CharmObj[[rbp]]$DEGenes

    if (!is.null(df) && "logFC" %in% colnames(df) && gene %in% rownames(df)) {
      data.frame(RBP = rbp, logFC = df[gene, "logFC"], stringsAsFactors = FALSE)
    } else {
      data.frame(RBP = rbp, logFC = NA_real_, stringsAsFactors = FALSE)
    }
  }) %>%
    bind_rows() %>%
    drop_na(logFC)

  # Check that we found the gene somewhere
  if (nrow(logfc_data) == 0) {
    stop(paste("Gene", gene, "not found in any RBP DEGenes tables."))
  }

  # Add "Status" for fill color
  logfc_data <- logfc_data %>%
    mutate(Status = ifelse(logFC > 0, "Upregulated", "Downregulated"))

  # Select top 10 and bottom 10 by LogFC
  top10 <- logfc_data %>% arrange(desc(logFC)) %>% slice_head(n = 10)
  bottom10 <- logfc_data %>% arrange(logFC) %>% slice_head(n = 10)
  selected <- bind_rows(top10, bottom10)

  # Order RBPs by LogFC
  selected$RBP <- factor(selected$RBP, levels = selected$RBP[order(selected$logFC)])

  # Build bar plot (styled like your favorite one)
  ggplot(selected, aes(x = reorder(RBP, logFC), y = logFC, fill = Status)) +
    geom_col(alpha = 0.8) +
    scale_fill_manual(values = c("Upregulated" = up_color, "Downregulated" = down_color)) +
    coord_flip() +
    labs(
      x = "RBP",
      y = "Log2 Fold Change",
      title = paste0("Top/bottom RBPs by logFC for ", gene),
      subtitle = paste0("Gene: ", gene),
      caption = paste0("Only top/bottom 10 shown")
    ) +
    theme_bw(base_family = "Arial MS") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 14),
      text = element_text(size = 14),
      legend.position = if (show_legend) "right" else "none",
      axis.text = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(face = "bold")
    )
}


plot_hallmark_nes_barplot <- function(CharmObj, geneset,
                                      up_color = "#BA3B46",
                                      down_color = "#53A2BE",
                                      show_legend = FALSE) {

  # Extract NES for the geneset across all RBPs safely
  NES_data <- lapply(names(CharmObj), function(rbp) {
    df <- CharmObj[[rbp]]$GSEA

    if (!is.null(df) && "NES" %in% colnames(df) && "pathway" %in% colnames(df)) {
      hits <- which(df$pathway == geneset)

      if (length(hits) >= 1) {
        data.frame(RBP = rbp, NES = as.numeric(df$NES[hits[1]]), stringsAsFactors = FALSE)
      } else {
        data.frame(RBP = rbp, NES = NA_real_, stringsAsFactors = FALSE)
      }
    } else {
      data.frame(RBP = rbp, NES = NA_real_, stringsAsFactors = FALSE)
    }
  }) %>%
    bind_rows() %>%
    drop_na(NES)

  # Check that we found the geneset somewhere
  if (nrow(NES_data) == 0) {
    stop(paste0("Geneset '", geneset, "' not found in any RBP GSEA tables."))
  }

  # Add up/down status for coloring
  NES_data <- NES_data %>%
    mutate(Status = ifelse(NES > 0, "Upregulated", "Downregulated"))

  # Select top 10 and bottom 10 RBPs by NES
  top10 <- NES_data %>% arrange(desc(NES)) %>% slice_head(n = 10)
  bottom10 <- NES_data %>% arrange(NES) %>% slice_head(n = 10)
  selected <- bind_rows(top10, bottom10)

  # Order RBPs by NES for plotting
  selected$RBP <- factor(selected$RBP, levels = selected$RBP[order(selected$NES)])

  # Create the bar plot (styled like your favourite one)
  ggplot(selected, aes(x = reorder(RBP, NES), y = NES, fill = Status)) +
    geom_col(alpha = 0.8) +
    scale_fill_manual(values = c("Upregulated" = up_color, "Downregulated" = down_color)) +
    coord_flip() +
    labs(
      x = "RBP",
      y = "Normalised Enrichment Score (NES)",
      title = paste0("Top/bottom RBPs by NES for ", geneset),
      subtitle = paste0("Geneset: ", geneset),
      caption = paste0("Only top/bottom 10 RBPs shown")
    ) +
    theme_bw(base_family = "Arial MS") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 14),
      text = element_text(size = 14),
      legend.position = if (show_legend) "right" else "none",
      axis.text = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(face = "bold")
    )
}


plot_event_dpsi_barplot <- function(CharmObj, Event.ID,
                                    up_color = "#BA3B46",
                                    down_color = "#53A2BE",
                                    show_legend = FALSE) {

  # --- Extract dPSI for the event across all RBPs safely ---
  dpsi_data <- lapply(names(CharmObj), function(rbp) {
    df <- CharmObj[[rbp]]$VulcanTable

    # Validate structure
    if (!is.null(df) && "dPSI" %in% colnames(df) && "Event.ID" %in% colnames(df)) {
      hits <- which(df$Event.ID == Event.ID)

      if (length(hits) >= 1) {
        # Take first hit if multiple
        data.frame(RBP = rbp, dPSI = as.numeric(df$dPSI[hits[1]]), stringsAsFactors = FALSE)
      } else {
        data.frame(RBP = rbp, dPSI = NA_real_, stringsAsFactors = FALSE)
      }
    } else {
      data.frame(RBP = rbp, dPSI = NA_real_, stringsAsFactors = FALSE)
    }
  }) %>%
    dplyr::bind_rows() %>%
    tidyr::drop_na(dPSI)

  # --- Check that we found this Event.ID somewhere ---
  if (nrow(dpsi_data) == 0) {
    stop(paste0("Event '", Event.ID, "' not found in any RBP VulcanTable."))
  }

  # --- Add Status for fill color ---
  dpsi_data <- dpsi_data %>%
    dplyr::mutate(Status = ifelse(dPSI > 0, "Upregulated", "Downregulated"))

  # --- Select top 10 and bottom 10 by dPSI ---
  top10 <- dpsi_data %>% dplyr::arrange(desc(dPSI)) %>% dplyr::slice_head(n = 10)
  bottom10 <- dpsi_data %>% dplyr::arrange(dPSI) %>% dplyr::slice_head(n = 10)
  selected <- dplyr::bind_rows(top10, bottom10)

  # --- Order RBPs by dPSI value ---
  selected$RBP <- factor(selected$RBP, levels = selected$RBP[order(selected$dPSI)])

  # --- Build bar plot ---
  p <- ggplot2::ggplot(selected, ggplot2::aes(x = reorder(RBP, dPSI), y = dPSI, fill = Status)) +
    ggplot2::geom_col(alpha = 0.8) +
    ggplot2::scale_fill_manual(values = c("Upregulated" = up_color, "Downregulated" = down_color)) +
    ggplot2::coord_flip() +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.2f", dPSI)),
      hjust = ifelse(selected$dPSI > 0, -0.1, 1.1),
      color = "white", size = 4
    ) +
    ggplot2::labs(
      x = "RBP",
      y = "ΔPSI",
      title = paste0("Top/bottom RBPs by ΔPSI for event ", Event.ID),
      subtitle = paste0("Event ID: ", Event.ID),
      caption = "Only top/bottom 10 shown"
    ) +
    ggplot2::theme_bw(base_family = "Arial MS") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 18),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 14),
      text = ggplot2::element_text(size = 14),
      legend.position = if (show_legend) "right" else "none",
      axis.text = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold"),
      legend.title = ggplot2::element_text(face = "bold"),
      legend.text = ggplot2::element_text(face = "bold")
    )

  return(p)
}

eCLIPSE_full <- function(bindingvalues_nested, rnaBP, target, dPSI,
                         metric = "FDR", title = NULL, plot = TRUE,
                         schematic_exon_skipping = NULL,
                         schematic_intron_retention = NULL,
                         event_type = NULL) {
  # require schematics provided
  if (is.null(schematic_exon_skipping) || is.null(schematic_intron_retention))
    stop("schematic_exon_skipping() and schematic_intron_retention() must be provided.")
  
  # --- 1. Validate RBP and target
  if (!rnaBP %in% names(bindingvalues_nested)) stop(paste("RBP", rnaBP, "not found."))
  if (!target %in% names(bindingvalues_nested[[rnaBP]])) stop(paste("Target", target, "not found for", rnaBP))
  
  # --- 2. Fuzzy numeric matching for dPSI
  available_dpsis <- names(bindingvalues_nested[[rnaBP]][[target]])
  if (is.null(available_dpsis) || length(available_dpsis) == 0)
    stop(paste("No dPSI entries found for", rnaBP, "-", target))
  available_nums <- suppressWarnings(as.numeric(sub(" .*", "", available_dpsis)))
  dpsi_numeric <- as.numeric(dPSI)
  if (is.na(dpsi_numeric)) stop("dPSI must be numeric or numeric-like, e.g. '0.05'")
  best_idx <- which.min(abs(available_nums - dpsi_numeric))
  chosen_dpsi_name <- available_dpsis[best_idx]
  message("→ Using dPSI = ", chosen_dpsi_name, " (closest to requested ", dPSI, ")")
  
  # --- 3. Extract nested data
  this_data <- bindingvalues_nested[[rnaBP]][[target]][[chosen_dpsi_name]]
  if (is.null(this_data)) stop("Selected nested object is NULL.")
  
  # --- 4. Extract raw event data
  raw_inc   <- this_data$raw$IncreasedEvents
  raw_dec   <- this_data$raw$DecreasedEvents
  raw_maint <- this_data$raw$MaintainedEvents
  required_cols <- c("Pos","Value")
  if (!all(required_cols %in% names(raw_inc))) stop("raw_inc lacks Pos/Value")
  if (!all(required_cols %in% names(raw_dec))) stop("raw_dec lacks Pos/Value")
  if (!all(required_cols %in% names(raw_maint))) stop("raw_maint lacks Pos/Value")
  
  # --- 5. Summary counts
  get_value <- function(df) {
    if (is.null(df) || nrow(df) == 0) return(NA)
    if ("Value" %in% names(df)) return(df$Value[1])
    as.numeric(df[1,1])
  }
  inc_target   <- get_value(this_data$raw$IncreasedTargetEvents)
  inc_total    <- get_value(this_data$raw$TotalIncreasedEvents)
  maint_target <- get_value(this_data$raw$MaintainedTargetEvents)
  maint_total  <- get_value(this_data$raw$TotalMaintainedEvents)
  dec_target   <- get_value(this_data$raw$DecreasedTargetEvents)
  dec_total    <- get_value(this_data$raw$TotalDecreasedEvents)
  
  # --- 6. Extract metric data
  pval_inc <- this_data$pval$pvalinc
  pval_dec <- this_data$pval$pvaldec
  odds_inc <- this_data$oddsrat$oddsratinc
  odds_dec <- this_data$oddsrat$oddsratdec
  # defensive checks
  if (!all(c("Pos","Value") %in% names(pval_inc))) stop("pval_inc missing Pos/Value")
  if (!all(c("Pos","Value") %in% names(pval_dec))) stop("pval_dec missing Pos/Value")
  if (!all(c("Pos","Value") %in% names(odds_inc))) stop("odds_inc missing Pos/Value")
  if (!all(c("Pos","Value") %in% names(odds_dec))) stop("odds_dec missing Pos/Value")
  
  # --- 7. Build data frames
  dfforvis_raw <- data.frame(
    Pos = raw_inc$Pos,
    IncreasedEvents = raw_inc$Value,
    DecreasedEvents = raw_dec$Value,
    MaintainedEvents = raw_maint$Value
  )
  dfforvis_pval <- data.frame(
    pos = pval_inc$Pos,
    pvalinc = pval_inc$Value,
    pvaldec = pval_dec$Value
  )
  dfforvis_oddsratio <- data.frame(
    pos = odds_inc$Pos,
    oddsratinc = odds_inc$Value,
    oddsratdec = odds_dec$Value
  )
  
  # --- 8. Decide x-axis max for the main map
  pos_max <- if (!all(is.na(dfforvis_raw$Pos))) max(dfforvis_raw$Pos, na.rm = TRUE) else NA
  if (is.na(pos_max) && nrow(dfforvis_pval) > 0) pos_max <- max(dfforvis_pval$pos, na.rm = TRUE)
  if (is.na(pos_max) && nrow(dfforvis_oddsratio) > 0) pos_max <- max(dfforvis_oddsratio$pos, na.rm = TRUE)
  if (is.na(pos_max)) stop("Cannot determine position range for x-axis.")
  xlim_max <- if (pos_max >= 900) 1000 else if (pos_max >= 400) 500 else pos_max
  
  # --- 9. Choose schematic robustly
  chosen_event_type <- NULL
  if (!is.null(event_type)) chosen_event_type <- event_type
  else {
    # try to read app input if available
    try({
      if (exists("input", envir = parent.frame())) {
        pf_input <- get("input", envir = parent.frame())
        if (!is.null(pf_input$binding_eventtype)) chosen_event_type <- pf_input$binding_eventtype
      }
    }, silent = TRUE)
  }
  if (is.null(chosen_event_type)) chosen_event_type <- "Exon Skipping"
  schem <- if (grepl("Intron Retention", chosen_event_type, ignore.case = TRUE)) {
    schematic_intron_retention()
  } else {
    schematic_exon_skipping()
  }
  # defensive: ensure schem elements exist and columns are present
  if (is.null(schem$rects) || !all(c("xmin","xmax","ymin","ymax","fill") %in% colnames(schem$rects))) {
    stop("schematic rects missing or malformed")
  }
  
  # --- compute schematic x-limits so exon never clipped while preserving alignment ---
  schem_xlim_max <- max(xlim_max, max(schem$rects$xmax, na.rm = TRUE))
  schem_xlim_min <- min(0, min(schem$rects$xmin, na.rm = TRUE))
  
  # --- helper: pick colors from schematic fills with fallbacks
  fill_vec <- as.character(schem$rects$fill)
  
  # Default: use schematic colors
  inc_color <- if (length(fill_vec) >= 1) fill_vec[1] else "#de425b"
  dec_color <- if (length(fill_vec) >= 2) fill_vec[2] else "#769fca"
  maint_color <- if (length(fill_vec) >= 3) fill_vec[3] else "black"
  
  # Override colors for each event type to original "canonical" colors
  if (grepl("Intron Retention", chosen_event_type, ignore.case = TRUE)) {
    inc_color <- "#de425b"   # IR increased
    dec_color <- "#769fca"   # IR decreased
    maint_color <- "black"   # IR maintained
  } else if (grepl("Exon Skipping", chosen_event_type, ignore.case = TRUE)) {
    inc_color <- "#de425b"   # ES increased
    dec_color <- "#769fca"   # ES decreased
    maint_color <- "black"   # ES maintained
  }
  # --- Set unified x-axis limits and expansion for alignment
  x_axis_min <- 0
  x_axis_max <- xlim_max
  x_axis_expand <- 0.02  # 2% padding on left and right
  
  # --- Schematic plot
  schem_plot <- ggplot() +
    geom_rect(data = schem$rects,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
              color = "black", linewidth = 0.35, inherit.aes = FALSE) +
    { if (!is.null(schem$intron_rect))
      geom_rect(data = schem$intron_rect,
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                fill = "#E5E7E9", color = "black", linewidth = 0.35, inherit.aes = FALSE)
      else NULL } +
    { if (!is.null(schem$introns))
      geom_segment(data = schem$introns,
                   aes(x = x, xend = xend, y = y, yend = yend),
                   linewidth = 0.8, inherit.aes = FALSE)
      else NULL } +
    scale_x_continuous(
      limits = c(x_axis_min, x_axis_max),
      expand = expansion(mult = c(x_axis_expand, x_axis_expand))
    ) +
    scale_y_continuous(
      limits = c(-1.1, 0.25),
      expand = c(0,0)
    ) +
    scale_fill_identity() +
    theme_bw(base_size = 18) +
    theme(
      text = element_text(size = 18, face = "bold"),
      axis.title = element_blank(),
      axis.text.x = element_text(size = 16, face = "bold"),
      axis.text.y = element_text(size = 16, face = "bold"),
      panel.grid = element_blank(),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
      legend.position = "none"
    )
  
  # --- 11. Map plot (use schematic colors)
  mapplot <- ggplot(dfforvis_raw, aes(x = Pos)) +
    geom_line(aes(y = IncreasedEvents), color = inc_color, size = 1) +
    geom_line(aes(y = DecreasedEvents), color = dec_color, size = 1) +
    geom_line(aes(y = MaintainedEvents), color = maint_color, size = 1) +
    xlim(0, xlim_max) +
    labs(title = paste0(rnaBP, " RNA Binding Map - ", title), x = "", y = "Normalised\nDensity\n") +
    theme_bw() +
    theme(
      axis.line = element_line(colour = "black"),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none",
      text = element_text(size = 15, face = "bold")
    )
  
  # --- 12. Labels plot (summary counts) - bold
  labels_plot <- ggplot() +
    annotate("text", x = 1, y = 0,
             label = paste0("Increased: ", inc_target, " (", inc_total, ")\n"),
             color = inc_color, size = 6, fontface = "bold", hjust = -0.75) +
    annotate("text", x = 1, y = 0,
             label = paste0("Maintained: ", maint_target, " (", maint_total, ")\n"),
             color = maint_color, size = 6, fontface = "bold", hjust = 0.5) +
    annotate("text", x = 1, y = 0,
             label = paste0("Decreased: ", dec_target, " (", dec_total, ")\n"),
             color = dec_color, size = 6, fontface = "bold", hjust = 1.75) +
    theme_void() +
    theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0))
  
  # --- 13. Metric-specific plot (use schematic colors)
  if (metric == "EffectSize") {
    metricplot <- ggplot(dfforvis_oddsratio, aes(x = pos)) +
      geom_line(aes(y = oddsratinc), color = inc_color, size = 1) +
      geom_line(aes(y = oddsratdec), color = dec_color, size = 1) +
      xlim(0, xlim_max) +
      labs(x = "Position", y = "Chi-Squared\nStatistic\n") +
      theme_bw() +
      theme(axis.line = element_line(colour = "black"),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            legend.position = "none",
            text = element_text(size = 15, face = "bold"))
  } else {
    metricplot <- ggplot(dfforvis_pval, aes(x = pos)) +
      geom_line(aes(y = pvalinc), color = inc_color, size = 1) +
      geom_line(aes(y = pvaldec), color = dec_color, size = 1) +
      xlim(0, xlim_max) +
      labs(x = "Position", y = "-log10 (FDR)\n") +
      theme_bw() +
      theme(axis.line = element_line(colour = "black"),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            legend.position = "none",
            text = element_text(size = 15, face = "bold"))
  }
  
  
  # --- Mapplot and metricplot: use the same limits and expansion
  mapplot <- mapplot +
    scale_x_continuous(
      limits = c(x_axis_min, x_axis_max),
      expand = expansion(mult = c(x_axis_expand, x_axis_expand))
    )
  
  metricplot <- metricplot +
    scale_x_continuous(
      limits = c(x_axis_min, x_axis_max),
      expand = expansion(mult = c(x_axis_expand, x_axis_expand))
    )
  
  # --- 14. Combine vertically using cowplot (preserves sizes)
  library(cowplot)
  completeplot <- plot_grid(
    schem_plot, mapplot, labels_plot, metricplot,
    ncol = 1, rel_heights = c(0.25, 0.35, 0.1, 0.3), align = "v"
  )
  
  if (plot) return(completeplot) else return(list(pval_data = dfforvis_pval, oddsratio_data = dfforvis_oddsratio))
}


violinplotter <- function(charmobj, rbp, control_color = "#7272AB", other_color = "#283D3B", highlight_color = "red") {
  # Get expression matrix
  expr <- charmobj$corcounts

  # Keep only Control and requested RBP samples
  keep_samples <- charmobj$SampleType %in% c("Control", rbp)
  expr <- expr[rbp, keep_samples]

  # Build dataframe for ggplot
  df <- data.frame(
    Expression = as.numeric(expr),
    Group = charmobj$SampleType[keep_samples],
    Sample = names(expr)  # keep sample names for matching
  )

  # Force Control to appear first
  df$Group <- factor(df$Group, levels = c("Control", rbp))


  # Highlight control points that are also in Experiment for this RBP
  control_samples <- charmobj$SampleType == "Control"
  df$Highlight <- FALSE
  df$Highlight[df$Group == "Control" & charmobj$Experiment[keep_samples] == rbp] <- TRUE
  # Compute mean expression per group
  mean_expr <- df %>%
    group_by(Group) %>%
    summarise(mean_exp = mean(Expression, na.rm = TRUE))

  # Log2 fold-change = difference of means (rbp - Control)
  logFC <- mean_expr$mean_exp[mean_expr$Group == rbp] -
    mean_expr$mean_exp[mean_expr$Group == "Control"]

  # Create plot
  p <- ggplot(df, aes(x = Group, y = Expression, fill = Group)) +

    # All points (low alpha)
    geom_jitter(data = df, aes(x = Group, y = Expression),
                width = 0.15, alpha = 0.3, color = "#7272AB") +

    # Violin plot
    geom_violin(trim = FALSE, alpha = 0.6) +

    # Highlighted points on top
    geom_jitter(data = df[df$Highlight, ],
                aes(x = Group, y = Expression),
                width = 0.15, color = "#E54B4B", size = 3, alpha = 0.8) +

    # Fill colors for violins
    scale_fill_manual(values = c("Control" = "#7272AB", rbp = "#283D3B")) +

    # Labels and theme
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
      text = element_text(size = 15, family = "Arial MS"),
      legend.position = "none"
    ) +
    xlab("") + ylab("Expression") +
    ylim(0, max(df$Expression) * 1.05)

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
    xlab("Log Fold-Change") +
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
      text = element_text(size = 15, family = "Arial MS")
    )
}


plot_rbp_volcano <- function(Charmobj, rbp) {
  #----------------------------
  # Step 1: Differential Expression
  #----------------------------
  top.table <- as.data.frame(Charmobj$DEGenes[rbp])
  colnames(top.table) <- c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "highlight")
  top.table$gene <- row.names(top.table)
  # Move gene column to first position
  top.table <- top.table[, c("gene", setdiff(colnames(top.table), "gene"))]

  # Order by absolute t-statistic
  top.table <- top.table[order(-abs(top.table$t)), ]

  #----------------------------
  # Step 2: Volcano Plot
  #----------------------------
  volcano_plot <- ggplot() +
    # Base layer: all genes grey
    geom_point(data = subset(top.table, highlight == "None"),
               aes(x = logFC, y = B), color = "#CCCCCC", alpha = 0.6) +
    # Red layer: rbp gene
    geom_point(data = subset(top.table, highlight == "RBP"),
               aes(x = logFC, y = B), color = "#A10702", alpha = 0.9, size = 2.5) +
    # Labels for RBP
    geom_text_repel(
      data = subset(top.table, highlight == "RBP"),
      aes(x = logFC, y = B, label = rownames(subset(top.table, highlight == "RBP"))),
      size = 5, fontface = "italic", color = "#A10702", nudge_y = 1
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.line = element_line(colour = "black"),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      legend.position = "none"
    ) +
    labs(title = paste(rbp," KD"),
         x = "Log Fold-Change", y = "B-statistic")

  return(list(top_table = top.table, volcano_plot = volcano_plot))
}

plot_gsea <- function(Charmobj, rbp, varoi = "SampleType", thresh = 0.05,
                      species = "Homo sapiens",
                      collection = "H", subcollection = NULL,
                      up_color = "#BA3B46", down_color = "#53A2BE",
                      show_legend = FALSE) {

  # Use the Charmobj argument and the rbp string properly
  hallmarks.res.tidy <- as.data.frame(Charmobj$GSEA[[rbp]])
  colnames(hallmarks.res.tidy) <- c("pathway","pval","padj","log2err","ES","NES","size","leadingEdge","Status")

  #----------------------------
  # Step 1: Filter significant results
  #----------------------------
  data_to_plot <- hallmarks.res.tidy %>%
    filter(padj < thresh) %>%
    arrange(-abs(NES))

  #----------------------------
  # Step 2: Plot
  #----------------------------
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
    )

  return(list(geneset_table = hallmarks.res.tidy, gsea_plot = gsea_plot_hall))
}

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

  other_df <- rbp_results[[other_rbp]]

  merged <- merge(
    ref_df, other_df, by = "Gene",
    suffixes = c(paste0("_", rbp_label), paste0("_", other_rbp))
  )

  # Spearman correlation + p-value
  test <- suppressWarnings(cor.test(
    merged[[paste0("t_", rbp_label)]],
    merged[[paste0("t_", other_rbp)]],
    method = "spearman"
  ))
  rho <- unname(test$estimate)
  pval <- test$p.value

  # Build title
  title_txt <- paste0(
    if (!is.null(plot_title)) paste0(plot_title, ": ") else "",
    "Correlation between ", rbp_label, " and ", other_rbp,
    " (Spearman ρ = ", round(rho, 2),
    ", p = ", signif(pval, 3), ")"
  )

  # Create Plotly scatter
  plot_ly(
    merged,
    x = ~get(paste0("t_", rbp_label)),
    y = ~get(paste0("t_", other_rbp)),
    type = 'scatter',
    mode = 'markers',
    text = ~paste("Gene:", Gene),
    hoverinfo = 'text',
    marker = list(size = 7, color = "#DDDDDD", line = list(width = 1, color = 'black'))
  ) %>%
    layout(
      title = title_txt,
      xaxis = list(title = paste0(rbp_label, " t-statistics")),
      yaxis = list(title = paste0(other_rbp, " t-statistics")),
      showlegend = FALSE
    )
}

correl_scatter_gsea_plotly <- function(gsea_results, rbp, other_rbp, plot_title = NULL, thresh = 0.05,
                                       species = "Homo sapiens",
                                       collection = "H", subcollection = NULL,
                                       up_color = "#BA3B46", down_color = "#53A2BE",
                                       show_legend = FALSE) {

  if (is.data.frame(rbp)) {
    # rbp is user-uploaded
    top.table <- rbp
    colnames(top.table) <- c("Gene", "t")
    rbp_label <- "UserFile"

    DEGenes <- top.table[order(top.table$t, decreasing = TRUE), ]
    vectorranks <- DEGenes$t
    names(vectorranks) <- DEGenes$Gene

    hallmarks.gs <- msigdbr(species = species, collection = collection, subcollection = subcollection)
    hallmarks.gsets <- split(hallmarks.gs$gene_symbol, hallmarks.gs$gs_name)
    hallmarks.gsets <- lapply(hallmarks.gsets, toupper)
    hallmarks.res <- fgsea(pathways = hallmarks.gsets, stats = vectorranks)
    hallmarks.res.tidy <- hallmarks.res %>%
      as_tibble() %>%
      arrange(desc(NES)) %>%
      mutate(Status = ifelse(NES > 0, "Upregulated", "Downregulated")) %>%
      arrange(padj) %>%
      mutate(pathway = gsub("^HALLMARK_", "", pathway))

    ref_df <- hallmarks.res.tidy

  } else if (rbp %in% names(gsea_results)) {
    ref_df <- gsea_results[[rbp]]
    rbp_label <- rbp
  } else {
    stop("rbp is neither a dataframe nor a known name in gsea_results")
  }

  other_df <- gsea_results[[other_rbp]]

  merged <- merge(
    ref_df, other_df, by = "pathway",
    suffixes = c(paste0("_", rbp_label), paste0("_", other_rbp))
  )

  # Spearman correlation + p-value
  test <- suppressWarnings(cor.test(
    merged[[paste0("NES_", rbp_label)]],
    merged[[paste0("NES_", other_rbp)]],
    method = "spearman"
  ))
  rho <- unname(test$estimate)
  pval <- test$p.value

  # Plotly scatter
  title_txt <- paste0(
    if (!is.null(plot_title)) paste0(plot_title, ": ") else "",
    "GSEA correlation between ", rbp_label, " and ", other_rbp,
    " (Spearman ρ = ", round(rho, 2),
    ", p = ", signif(pval, 3), ")"
  )

  plot_ly(
    merged,
    x = ~get(paste0("NES_", rbp_label)),
    y = ~get(paste0("NES_", other_rbp)),
    type = 'scatter',
    mode = 'markers',
    text = ~paste("Pathway:", pathway),
    hoverinfo = 'text',
    marker = list(size = 7, color = "#DDDDDD", line = list(width = 1, color = 'black'))
  ) %>%
    layout(
      title = title_txt,
      xaxis = list(title = paste0(rbp_label, " NES")),
      yaxis = list(title = paste0(other_rbp, " NES")),
      showlegend = show_legend
    )
}

exp_correl <- function(rbp_results, rbp, correl_num = NULL,
                       n_pos = NULL, n_neg = NULL, other_rbps = NULL) {

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

  # Select top RBPs
  if (!is.null(other_rbps)) {
    valid_rbps <- intersect(other_rbps, cor_results$RBP)
    top_cor <- cor_results[cor_results$RBP %in% valid_rbps, ]
  } else if (!is.null(correl_num) && !is.na(correl_num) && correl_num > 0) {
    top_cor <- head(cor_results[order(-abs(cor_results$Correlation)), ], correl_num)
  } else if ((!is.null(n_pos) && !is.na(n_pos) && n_pos > 0) ||
             (!is.null(n_neg) && !is.na(n_neg) && n_neg > 0)) {

    top_pos <- if (!is.null(n_pos) && !is.na(n_pos) && n_pos > 0) head(pos, n_pos) else NULL
    top_neg <- if (!is.null(n_neg) && !is.na(n_neg) && n_neg > 0) head(neg, n_neg) else NULL
    top_cor <- rbind(top_pos, top_neg)
  } else {
    stop("Please provide either other_rbps, correl_num, or n_pos/n_neg")
  }

  top_cor <- top_cor[order(top_cor$Correlation, decreasing = TRUE), ]
  top_cor$RBP <- factor(top_cor$RBP, levels = top_cor$RBP)

  heatmap_plot <- ggplot(top_cor, aes(x = rbp_label, y = RBP, fill = Correlation)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", Correlation),
                  color = ifelse(abs(Correlation) < 0.3, "black", "white")), size = 5) +
    scale_color_identity() +
    scale_fill_gradient2(low = "black", mid = "grey50", high = "#601700",
                         midpoint = 0, limits = c(-1, 1)) +
    labs(title = paste("Spearman correlations with", rbp_label, "- Expression"),
         x = "Reference RBP", y = "Other RBPs") +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  return(list(correlation_table = cor_results,
              top_table = top_cor,
              heatmap = heatmap_plot))
}


gsea_correl <- function(gsea_results, rbp, correl_num = NULL,
                        n_pos = NULL, n_neg = NULL, other_rbps = NULL,
                        thresh = 0.05,
                        species = "Homo sapiens",
                        collection = "H", subcollection = NULL,
                        up_color = "#BA3B46", down_color = "#53A2BE",
                        show_legend = FALSE) {


  if (is.data.frame(rbp)) {
    # User uploaded file
    top.table <- rbp
    stopifnot(all(c("Gene","t") %in% colnames(top.table)))
    rbp_label <- "UserFile"

    # Create vector of ranks for fgsea
    DEGenes <- top.table[order(top.table$t, decreasing = TRUE), ]
    vectorranks <- DEGenes$t
    names(vectorranks) <- toupper(DEGenes$Gene)

    # Run GSEA on uploaded file
    hallmarks.gs <- msigdbr(species = species, collection = collection, subcollection = subcollection)
    hallmarks.gsets <- split(hallmarks.gs$gene_symbol, hallmarks.gs$gs_name)
    hallmarks.gsets <- lapply(hallmarks.gsets, toupper)
    hallmarks.res <- fgsea(pathways = hallmarks.gsets, stats = vectorranks)
    hallmarks.res.tidy <- hallmarks.res %>%
      as_tibble() %>%
      arrange(desc(NES)) %>%
      mutate(Status = ifelse(NES > 0, "Upregulated", "Downregulated")) %>%
      arrange(padj) %>%
      mutate(pathway = gsub("^HALLMARK_", "", pathway))

    ref_df <- hallmarks.res.tidy %>% select(pathway, NES)
  } else if (rbp %in% names(gsea_results)) {
    ref_df <- gsea_results[[rbp]] %>% select(pathway, NES)
    rbp_label <- rbp
  } else {
    stop("rbp is neither a dataframe nor a known name in gsea_results")
  }

  stopifnot(all(c("pathway","NES") %in% colnames(ref_df)))
  ref_vec <- ref_df$NES
  names(ref_vec) <- ref_df$pathway


  cor_results <- data.frame(RBP=character(), Correlation=numeric(), Pvalue=numeric())

  for (other_rbp in setdiff(names(gsea_results), rbp_label)) {
    other_df <- gsea_results[[other_rbp]] %>% select(pathway, NES)
    merged <- merge(ref_df, other_df, by="pathway", suffixes=c("_ref","_other"))
    if (nrow(merged) > 2) {
      test <- suppressWarnings(cor.test(merged$NES_ref, merged$NES_other, method="spearman"))
      cor_results <- rbind(cor_results,
                           data.frame(RBP=other_rbp,
                                      Correlation=unname(test$estimate),
                                      Pvalue=test$p.value))
    }
  }


  top_pos <- top_neg <- NULL

  if (!is.null(other_rbps)) {
    valid_rbps <- intersect(other_rbps, cor_results$RBP)
    top_cor <- cor_results[cor_results$RBP %in% valid_rbps, ]
  } else if (!is.null(correl_num)) {
    cor_results <- cor_results[order(-abs(cor_results$Correlation)), ]
    top_cor <- head(cor_results, correl_num)
  } else if (!is.null(n_pos) | !is.null(n_neg)) {
    if (!is.null(n_pos)) top_pos <- head(cor_results[order(-cor_results$Correlation), ], n_pos)
    if (!is.null(n_neg)) top_neg <- head(cor_results[order(cor_results$Correlation), ], n_neg)
    # Combine only non-null
    top_cor <- do.call(rbind, Filter(Negate(is.null), list(top_pos, top_neg)))
    if (nrow(top_cor) == 0) stop("No correlations available with the given n_pos/n_neg")
  } else {
    stop("Please provide either other_rbps, correl_num, or (n_pos/n_neg)")
  }

  # Factor for plotting
  top_cor <- top_cor[order(top_cor$Correlation, decreasing = TRUE), ]
  top_cor$RBP <- factor(top_cor$RBP, levels=top_cor$RBP)


  heatmap_plot <- ggplot(top_cor, aes(x=rbp_label, y=RBP, fill=Correlation)) +
    geom_tile(color="white") +
    geom_text(aes(label=sprintf("%.2f", Correlation),
                  color=ifelse(abs(Correlation)<0.3, "black","white")), size=5) +
    scale_color_identity() +
    scale_fill_gradient2(low="black", mid="grey50", high="#601700", midpoint=0, limits=c(-1,1)) +
    labs(title=paste("Spearman correlations with", rbp_label, "- GSEA"),
         x="Reference RBP", y="Other RBPs") +
    theme_minimal(base_size=14) +
    theme(axis.text.x=element_text(angle=45, hjust=1))

  return(list(correlation_table=cor_results,
              top_table=top_cor,
              heatmap=heatmap_plot))
}

library(Rsubread)
library(limma)
library(edgeR)
library(DESeq2)
library(magrittr)
library(ggrepel)
library(tidyverse)
library(vsn)
library(ggplot2)
library(NMF)
library(biomaRt)
library(ggpubr)
library(factoextra)
library(ggfortify)
library(plotly)
library(reshape2)
library(fgsea)
library(msigdbr)
library(data.table)
library(colorspace)
library(readxl)
library(fs)
library(betAS)

set.seed(1906)
library(data.table)
library(fs)

# --- STEP 1: Build the flat list of bindingvalues tables ---
build_bindingvalues_list <- function(
  base_path = "~/Projects/StressGranules/AS.WC_Transcriptome/shRNAExp"
) {
  base_path <- fs::path_expand(base_path)
  subdirs <- fs::dir_ls(base_path, type = "directory", recurse = FALSE)
  folder_names <- fs::path_file(subdirs)
  out_list <- list()

  for (folder in folder_names) {
    folder_path <- fs::path(base_path, folder)
    matches <- fs::dir_ls(folder_path, regexp = "BindingValues_IR_K562\\.txt$", type = "file")

    if (length(matches) == 0) {
      message("No 'BindingValues_both.txt' file found in: ", folder_path)
      next
    }

    if (length(matches) > 1) {
      info <- fs::file_info(matches)
      matches <- matches[order(info$modification_time, decreasing = TRUE)]
      message("Multiple BindingValues_both.txt files in ", folder_path,
              "; using latest: ", fs::path_file(matches[1]))
    }

    dt <- data.table::fread(matches[1])
    out_list[[folder]] <- list(bindingvalues = dt)
  }

  return(out_list)
}

# --- STEP 2: Restructure into nested hierarchy ---
restructure_bindingvalues_full <- function(bindingvalues_list) {
  out <- list()

  for (rbp_name in names(bindingvalues_list)) {
    df <- bindingvalues_list[[rbp_name]]$bindingvalues
    if (is.null(df) || nrow(df) == 0) next

    setDT(df)

    # Group by Target, dPSI, Metric
    nested_list <- df[, .(data = list(.SD)), by = .(Target, dPSI, Metric)]

    # Split by Target
    rbp_nested <- split(nested_list, by = "Target", keep.by = FALSE)

    rbp_nested <- lapply(rbp_nested, function(target_dt) {
      # Split by dPSI
      target_split <- split(target_dt, by = "dPSI", keep.by = FALSE)

      lapply(target_split, function(dpsi_dt) {
        # Define metric groups
        pval_metrics    <- c("pvalinc", "pvaldec")
        oddsrat_metrics <- c("oddsratinc", "oddsratdec")
        raw_metrics     <- c(
          "IncreasedEvents", "DecreasedEvents", "MaintainedEvents",
          "IncreasedTargetEvents", "TotalIncreasedEvents",
          "MaintainedTargetEvents", "TotalMaintainedEvents",
          "DecreasedTargetEvents", "TotalDecreasedEvents"
        )

        # Helper to clean each dataset
        clean_data <- function(dt_list) {
          if (length(dt_list) == 0) return(NULL)
          lapply(dt_list, function(tbl) {
            # Keep only relevant columns safely
            if ("Pos" %in% names(tbl) && "Value" %in% names(tbl)) {
              tbl[, .(Pos, Value)]
            } else if ("ID" %in% names(tbl) && "Value" %in% names(tbl)) {
              tbl[, .(ID, Value)]
            } else if ("Value" %in% names(tbl)) {
              tbl[, .(Value)]
            } else {
              tbl
            }
          })
        }

        # Extract and clean each group
        pval_list <- clean_data(
          setNames(dpsi_dt[Metric %in% pval_metrics]$data,
                   dpsi_dt[Metric %in% pval_metrics]$Metric)
        )
        oddsrat_list <- clean_data(
          setNames(dpsi_dt[Metric %in% oddsrat_metrics]$data,
                   dpsi_dt[Metric %in% oddsrat_metrics]$Metric)
        )
        raw_list <- clean_data(
          setNames(dpsi_dt[Metric %in% raw_metrics]$data,
                   dpsi_dt[Metric %in% raw_metrics]$Metric)
        )

        # Build hierarchical structure
        list(
          pval    = if (length(pval_list)) pval_list else NULL,
          oddsrat = if (length(oddsrat_list)) oddsrat_list else NULL,
          raw     = if (length(raw_list)) raw_list else NULL
        )
      })
    })

    out[[rbp_name]] <- rbp_nested
  }

  return(out)
}

# --- STEP 3: Example usage ---
CharmObj_bindingvalues <- build_bindingvalues_list()

# Example usage:
CharmObj_bindingvalues_nested <- restructure_bindingvalues_full(CharmObj_bindingvalues)

saveRDS(CharmObj_bindingvalues_nested, "~/Projects/CHARM/data/Binding_IR_K562.RDS")

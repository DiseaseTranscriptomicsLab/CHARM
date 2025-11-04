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

build_bindingvalues_list <- function(
  base_path = "~/Projects/StressGranules/AS.WC_Transcriptome/shRNAExp"
) {
  # Expand base path
  base_path <- fs::path_expand(base_path)

  # List all immediate subdirectories (each RBP)
  subdirs <- fs::dir_ls(base_path, type = "directory", recurse = FALSE)
  folder_names <- fs::path_file(subdirs)

  # Initialize list for results
  out_list <- list()

  # Iterate through each folder and look for BindingValues_both.txt
  for (folder in folder_names) {
    folder_path <- fs::path(base_path, folder)

    # Find BindingValues_both.txt file
    matches <- fs::dir_ls(folder_path, regexp = "BindingValues_IR_K562\\.txt$", type = "file")

    if (length(matches) == 0) {
      message("No 'BindingValues_both.txt' file found in: ", folder_path)
      next  # skip this folder entirely
    }

    if (length(matches) > 1) {
      info <- fs::file_info(matches)
      matches <- matches[order(info$modification_time, decreasing = TRUE)]
      message("Multiple BindingValues_both.txt files in ", folder_path,
              "; using latest: ", fs::path_file(matches[1]))
    }

    # Read the selected file
    dt <- data.table::fread(matches[1])

    # Add to output list, using folder name as key
    out_list[[folder]] <- list(bindingvalues = dt)
  }

  return(out_list)
}

# Example usage
CharmObj_bindingvalues <- build_bindingvalues_list()




restructure_bindingvalues_full <- function(bindingvalues_list) {
  out <- list()

  for (rbp_name in names(bindingvalues_list)) {
    df <- bindingvalues_list[[rbp_name]]$bindingvalues
    if (is.null(df) || nrow(df) == 0) next

    setDT(df)  # ensure it's a data.table

    # Group by Target, dPSI, Metric → store only Pos & Value
    nested_list <- df[, .(data = list(.SD[, .(Pos, Value)])),
                      by = .(Target, dPSI, Metric)]

    # Now structure as: RBP → Target → dPSI → Metric → data.frame
    rbp_nested <- split(nested_list, by = "Target", keep.by = FALSE)

    rbp_nested <- lapply(rbp_nested, function(target_dt) {
      target_split <- split(target_dt, by = "dPSI", keep.by = FALSE)

      lapply(target_split, function(dpsi_dt) {
        setNames(dpsi_dt$data, dpsi_dt$Metric)
      })
    })

    out[[rbp_name]] <- rbp_nested
  }

  return(out)
}

# Example usage:
CharmObj_bindingvalues_nested <- restructure_bindingvalues_full(CharmObj_bindingvalues)

saveRDS(CharmObj_bindingvalues_nested, "~/Projects/CHARM/data/Binding_IR_K562.RDS")

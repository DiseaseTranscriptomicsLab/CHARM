library(data.table)
library(tidyverse)
library(betAS)
set.seed(1906)

CharmObj_rbp_All <- readRDS("~/Projects/CHARM/data/CharmObj_rbp_All.rds")
CharmObj_rbp_HEPG2 <- readRDS("~/Projects/CHARM/data/CharmObj_rbp_HEPG2.rds")
CharmObj_rbp_K562 <- readRDS("~/Projects/CHARM/data/CharmObj_rbp_K562.rds")
Charm.object <- readRDS("data/Charm.object.RDS")
Charm.object_K562 <- readRDS("data/Charm.object_K562.RDS")
Charm.object_HEPG2 <- readRDS("data/Charm.object_HEPG2.RDS")

#Both Cells
for (rbp in names(CharmObj_rbp_All)) {
  dataset_All_filtered <- CharmObj_rbp_All[[rbp]]
  psi <- dataset_All_filtered$PSI
  qual <- dataset_All_filtered$Qual

  # --- Automatically detect control and shRNA columns ---
  cols_CTRL_names <- grep("control", colnames(psi), value = TRUE, ignore.case = TRUE)
  cols_shRNA_names <- grep("shrna", colnames(psi), value = TRUE, ignore.case = TRUE)

  # Skip safely if nothing is found
  if (length(cols_CTRL_names) == 0 || length(cols_shRNA_names) == 0) {
    message("⚠️  Skipping ", rbp, ": no control or shRNA columns found.")
    next
  }

  # Convert to column indices for betAS
  cols_CTRL <- convertCols(psi, cols_CTRL_names)
  cols_shRNA <- convertCols(psi, cols_shRNA_names)

  # --- Build volcano table ---
  volcanoTable_Pdiff_All <- tryCatch({
    prepareTableVolcano(
      psitable = psi,
      qualtable = qual,
      npoints = 500,
      colsA = cols_CTRL,
      colsB = cols_shRNA,
      labA = "CTRL",
      labB = "shRNA",
      basalColor = "#89C0AE",
      interestColor = "#E69A9C",
      maxDevTable = maxDevSimulationN100,
      seed = TRUE,
      CoverageWeight = FALSE
    )
  }, error = function(e) {
    message("⚠️  Skipping ", rbp, " due to error: ", e$message)
    return(NULL)
  })

  # Skip if failed
  if (is.null(volcanoTable_Pdiff_All)) next

  # --- Format output ---
  volcanoTable_Pdiff_All <- volcanoTable_Pdiff_All[, c("EVENT", "GENE", "deltapsi", "Pdiff")]
  colnames(volcanoTable_Pdiff_All) <- c("Event.ID", "Gene", "dPSI", "Pdiff")

  # --- Store result ---
  Charm.object[[rbp]]$VulcanTable <- volcanoTable_Pdiff_All
}

# K562
for (rbp in names(CharmObj_rbp_K562)) {
  dataset_All_filtered <- CharmObj_rbp_K562[[rbp]]
  psi <- dataset_All_filtered$PSI
  qual <- dataset_All_filtered$Qual

  # --- Automatically detect control and shRNA columns ---
  cols_CTRL_names <- grep("control", colnames(psi), value = TRUE, ignore.case = TRUE)
  cols_shRNA_names <- grep("shrna", colnames(psi), value = TRUE, ignore.case = TRUE)

  # Skip safely if nothing is found
  if (length(cols_CTRL_names) == 0 || length(cols_shRNA_names) == 0) {
    message("⚠️  Skipping ", rbp, ": no control or shRNA columns found.")
    next
  }

  # Convert to column indices for betAS
  cols_CTRL <- convertCols(psi, cols_CTRL_names)
  cols_shRNA <- convertCols(psi, cols_shRNA_names)

  # --- Build volcano table ---
  volcanoTable_Pdiff_All <- tryCatch({
    prepareTableVolcano(
      psitable = psi,
      qualtable = qual,
      npoints = 500,
      colsA = cols_CTRL,
      colsB = cols_shRNA,
      labA = "CTRL",
      labB = "shRNA",
      basalColor = "#89C0AE",
      interestColor = "#E69A9C",
      maxDevTable = maxDevSimulationN100,
      seed = TRUE,
      CoverageWeight = FALSE
    )
  }, error = function(e) {
    message("⚠️  Skipping ", rbp, " due to error: ", e$message)
    return(NULL)
  })

  # Skip if failed
  if (is.null(volcanoTable_Pdiff_All)) next

  # --- Format output ---
  volcanoTable_Pdiff_All <- volcanoTable_Pdiff_All[, c("EVENT", "GENE", "deltapsi", "Pdiff")]
  colnames(volcanoTable_Pdiff_All) <- c("Event.ID", "Gene", "dPSI", "Pdiff")

  # --- Store result ---
  Charm.object_K562[[rbp]]$VulcanTable <- volcanoTable_Pdiff_All
}


# HEPG2
for (rbp in names(CharmObj_rbp_HEPG2)) {
  dataset_All_filtered <- CharmObj_rbp_HEPG2[[rbp]]
  psi <- dataset_All_filtered$PSI
  qual <- dataset_All_filtered$Qual

  # --- Automatically detect control and shRNA columns ---
  cols_CTRL_names <- grep("control", colnames(psi), value = TRUE, ignore.case = TRUE)
  cols_shRNA_names <- grep("shrna", colnames(psi), value = TRUE, ignore.case = TRUE)

  # Skip safely if nothing is found
  if (length(cols_CTRL_names) == 0 || length(cols_shRNA_names) == 0) {
    message("⚠️  Skipping ", rbp, ": no control or shRNA columns found.")
    next
  }

  # Convert to column indices for betAS
  cols_CTRL <- convertCols(psi, cols_CTRL_names)
  cols_shRNA <- convertCols(psi, cols_shRNA_names)

  # --- Build volcano table ---
  volcanoTable_Pdiff_All <- tryCatch({
    prepareTableVolcano(
      psitable = psi,
      qualtable = qual,
      npoints = 500,
      colsA = cols_CTRL,
      colsB = cols_shRNA,
      labA = "CTRL",
      labB = "shRNA",
      basalColor = "#89C0AE",
      interestColor = "#E69A9C",
      maxDevTable = maxDevSimulationN100,
      seed = TRUE,
      CoverageWeight = FALSE
    )
  }, error = function(e) {
    message("⚠️  Skipping ", rbp, " due to error: ", e$message)
    return(NULL)
  })

  # Skip if failed
  if (is.null(volcanoTable_Pdiff_All)) next

  # --- Format output ---
  volcanoTable_Pdiff_All <- volcanoTable_Pdiff_All[, c("EVENT", "GENE", "deltapsi", "Pdiff")]
  colnames(volcanoTable_Pdiff_All) <- c("Event.ID", "Gene", "dPSI", "Pdiff")

  # --- Store result ---
  Charm.object_HEPG2[[rbp]]$VulcanTable <- volcanoTable_Pdiff_All
}


# Save separately
saveRDS(Charm.object,   file = "~/Projects/CHARM/data/CharmObj_ALL_DEDS.rds")
saveRDS(Charm.object_HEPG2, file = "~/Projects/CHARM/data/CharmObj_HEPG2_DEDS.rds")
saveRDS(Charm.object_K562,  file = "~/Projects/CHARM/data/CharmObj_K562_DEDS.rds")

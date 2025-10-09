library(data.table)
library(tidyverse)
library(betAS)
set.seed(1906)

# --- Load main event list ---
FullTrial <- fread("~/Projects/StressGranules/AS.WC_Transcriptome/eCLIPSE/FullCassetteEvents.txt")

# --- Helper: filter dataset ---
NewFilterFunction <- function(filteredList, eventlist) {
  psiTable  <- filteredList$PSI %>% filter(EVENT %in% eventlist$EVENT)
  qualTable <- filteredList$Qual %>% filter(EVENT %in% psiTable$EVENT)
  list(
    PSI           = psiTable,
    Qual          = qualTable,
    EventsPerType = table(psiTable$COMPLEX),
    Samples       = colnames(psiTable)[-c(1:6)]
  )
}

# --- Helper: subset dataset by keyword (HEPG2/K562) ---
make_subset <- function(dataset, keyword) {
  psi_cols  <- c(1:6, grep(keyword, colnames(dataset$PSI)))
  qual_cols <- c(1:6, grep(keyword, colnames(dataset$Qual)))
  list(
    PSI     = dataset$PSI[, psi_cols, drop = FALSE],
    Qual    = dataset$Qual[, qual_cols, drop = FALSE],
    Samples = grep(keyword, dataset$Samples, value = TRUE)
  )
}

# --- Process each RBP folder ---
process_shRNAExp <- function(base_path = "~/Projects/StressGranules/AS.WC_Transcriptome/shRNAExp") {
  rbp_dirs <- fs::dir_ls(base_path, type = "directory")

  map(rbp_dirs, function(rbp_dir) {
    rbp_name <- basename(rbp_dir)
    vast_file <- list.files(file.path(rbp_dir, "TrimmedSamples/vast_out"),
                            pattern = "^INCLUSION_LEVELS", full.names = TRUE)[1]
    if (is.na(vast_file)) return(NULL)

    message("Processing ", rbp_name)

    dataset <- getDataset(pathTables = vast_file, tool = "vast-tools") %>%
      getEvents(tool = "vast-tools") %>%
      filterEvents(types = c("S", "IR"), N = 0)

    dataset$PSI  <- na.omit(dataset$PSI)
    dataset$Qual <- dataset$Qual[dataset$Qual$EVENT %in% dataset$PSI$EVENT, ]

    dataset_All   <- dataset
    dataset_HEPG2 <- make_subset(dataset, "HEPG2")
    dataset_K562  <- make_subset(dataset, "K562")

    list(All = dataset_All, HEPG2 = dataset_HEPG2, K562 = dataset_K562)
  }) %>% set_names(basename(rbp_dirs))
}

CharmObj_rbp <- process_shRNAExp()

# --- Apply event filter to all ---
CharmObj_rbp <- map(CharmObj_rbp, ~ map(.x, NewFilterFunction, eventlist = FullTrial))

# --- Safer Volcano builder ---
make_volcano <- function(psitable, qualtable, ctrl, shrna) {
  # --- Safety checks ---
  if (length(ctrl) == 0 || length(shrna) == 0) {
    message("⚠️  Skipping volcano: no control or shRNA columns found.")
    return(NULL)
  }
  if (nrow(psitable) == 0 || nrow(qualtable) == 0) {
    message("⚠️  Skipping volcano: empty PSI or Qual table.")
    return(NULL)
  }

  # Wrap the volcano creation in tryCatch for robustness
  safe_tab <- tryCatch({
    tab <- prepareTableVolcano(
      psitable = psitable,
      qualtable = qualtable,
      npoints = 500,
      colsA = convertCols(psitable, ctrl),
      colsB = convertCols(psitable, shrna),
      labA = "CTRL", labB = "shRNA",
      basalColor = "#89C0AE", interestColor = "#E69A9C",
      maxDevTable = maxDevSimulationN100,
      seed = TRUE, CoverageWeight = FALSE
    )
    tab <- tab[, c("EVENT", "GENE", "deltapsi", "Pdiff")]
    colnames(tab) <- c("Event.ID", "Gene", "dPSI", "Pdiff")
    tab
  }, error = function(e) {
    message("⚠️  Skipping due to error in prepareTableVolcano: ", e$message)
    return(NULL)
  })

  return(safe_tab)
}

# --- Build volcano tables for each RBP/subset safely ---
CharmObj_volcano <- map(CharmObj_rbp, function(rbp_data) {
  psitable <- rbp_data$All$PSI
  volcano_cfg <- make_volcano_cfg(psitable)

  imap(volcano_cfg, function(cfg, name) {
    make_volcano(
      psitable = rbp_data[[name]]$PSI,
      qualtable = rbp_data[[name]]$Qual,
      ctrl = cfg$ctrl,
      shrna = cfg$shrna
    )
  })
})

CharmObj_volcano <- map(CharmObj_volcano, purrr::compact)

# --- Save results ---
saveRDS(CharmObj_volcano %>% map("All"),   "~/Projects/CHARM/data/CharmObj_rbp_All.rds")
saveRDS(CharmObj_volcano %>% map("HEPG2"), "~/Projects/CHARM/data/CharmObj_rbp_HEPG2.rds")
saveRDS(CharmObj_volcano %>% map("K562"),  "~/Projects/CHARM/data/CharmObj_rbp_K562.rds")

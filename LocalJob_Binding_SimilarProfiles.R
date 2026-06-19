# =============================================================================
#  Pre-computation: Binding Profile Similarity Matrices
#  ─────────────────────────────────────────────────────────────────────────────
#  For each of the 6 binding datasets (ES/IR × Both/K562/HEPG2):
#    1. Smooth positional profiles (sliding-window, same as LocalJob_File_Binding_New.R)
#    2. Build per (RBP, Target) mean profile matrix per Direction (inc / dec)
#    3. Save the mean_profiles matrix — Euclidean distances are computed at
#       query time in the app (one row vs all others), so no N×N matrix is stored
#
#  Distance metric: Euclidean (same as LocalJob_File_Binding_New.R MDS + heatmaps)
#  Profiles are ranked by SMALLEST distance (= most similar shape + magnitude).
#
#  Output files (save to data/QS_Files/):
#    similar_binding_ES_both_inc.qs   – Exon Skipping, Both cells, oddsratinc
#    similar_binding_ES_both_dec.qs   – Exon Skipping, Both cells, oddsratdec
#    similar_binding_ES_K562_inc.qs
#    similar_binding_ES_K562_dec.qs
#    similar_binding_ES_HEPG2_inc.qs
#    similar_binding_ES_HEPG2_dec.qs
#    similar_binding_IR_both_inc.qs
#    similar_binding_IR_both_dec.qs
#    similar_binding_IR_K562_inc.qs
#    similar_binding_IR_K562_dec.qs
#    similar_binding_IR_HEPG2_inc.qs
#    similar_binding_IR_HEPG2_dec.qs
#
#  Each saved object is a named list with a single element:
#    $mean_profiles : numeric matrix [n_profiles × n_positions], smoothed mean profiles
#                     rownames = "RBP__Target"
#
# =============================================================================

# ── 0. Dependencies ───────────────────────────────────────────────────────────
required_packages <- c("data.table", "dplyr", "qs")
invisible(lapply(required_packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  library(pkg, character.only = TRUE)
}))

# ── 1. Paths — adjust as needed ───────────────────────────────────────────────
DATA_DIR   <- "~/Projects/CHARM/intermediary_files"   # where raw binding CSVs live
OUTPUT_DIR <- "data/QS_Files"                         # where .qs files will be saved

# Raw binding CSV files (same structure as LocalJob_File_Binding_New.R)
# Columns 1-4: RBP | Target | dPSI | Direction
# Columns 5-1004: pos_1 … pos_1000
input_files <- list(
  ES_both  = file.path(DATA_DIR, "Charm.binding.filtered.both_new_full.csv"),
  ES_K562  = file.path(DATA_DIR, "Charm.binding.filtered.both_new_full_K562.csv"),
  ES_HEPG2 = file.path(DATA_DIR, "Charm.binding.filtered.both_new_full_HEPG2.csv"),
  IR_both  = file.path(DATA_DIR, "Charm.binding.filtered.both_new_full_IR.csv"),
  IR_K562  = file.path(DATA_DIR, "Charm.binding.filtered.both_new_full_IR_K562.csv"),
  IR_HEPG2 = file.path(DATA_DIR, "Charm.binding.filtered.both_new_full_IR_HEPG2.csv")
)

# ── 2. Sliding-window smoother (identical to LocalJob_File_Binding_New.R) ────
smooth_profiles <- function(df, half_win = 4L, min_win = 5L) {
  df <- as.data.frame(df)
  pos_cols <- grep("^pos_\\d+$", colnames(df), value = TRUE)
  pos_cols <- pos_cols[order(as.integer(sub("pos_", "", pos_cols)))]
  n        <- length(pos_cols)
  meta_cols <- c("RBP", "Target", "dPSI", "Direction")
  mat       <- as.matrix(df[, pos_cols])
  smoothed  <- matrix(NA_real_, nrow = nrow(mat), ncol = n,
                      dimnames = list(NULL, pos_cols))
  for (i in seq_len(n)) {
    left  <- max(1L, i - half_win)
    right <- min(n,  i + half_win)
    win_size <- right - left + 1L
    if (win_size < min_win) {
      if (left == 1L) right <- min(n, left + min_win - 1L)
      else             left  <- max(1L, right - min_win + 1L)
    }
    smoothed[, i] <- rowMeans(mat[, left:right, drop = FALSE], na.rm = TRUE)
  }
  cbind(df[, meta_cols], as.data.frame(smoothed))
}

# ── 3. Build similarity object for one dataset + direction ───────────────────
#
#  Strategy:
#    • For each (RBP, Target) pair, average the smoothed profiles over all rows
#      that share the same RBP, Target, and Direction.  (Typically one row each,
#      but averaging is safe if duplicates exist.)
#    • This gives a matrix of shape [n_profiles × n_positions].
#    • Only the mean_profiles matrix is saved.  Euclidean distances between the
#      query and all other profiles are computed at query time in the app — this
#      is a single dist(rbind(query, others)) call and is fast even for large
#      profile sets, while avoiding storing an N×N distance matrix on disk.
#
build_similarity_object <- function(df_smooth, direction) {
  pos_cols <- grep("^pos_\\d+$", colnames(df_smooth), value = TRUE)
  pos_cols <- pos_cols[order(as.integer(sub("pos_", "", pos_cols)))]

  sub <- df_smooth[df_smooth$Direction == direction, , drop = FALSE]
  if (nrow(sub) == 0) {
    warning("No rows for Direction = '", direction, "'. Skipping.")
    return(NULL)
  }

  # Mean profile per (RBP, Target)
  sub$profile_id <- paste(sub$RBP, sub$Target, sep = "__")
  ids <- unique(sub$profile_id)

  mean_mat <- do.call(rbind, lapply(ids, function(pid) {
    rows <- sub[sub$profile_id == pid, pos_cols, drop = FALSE]
    colMeans(as.matrix(rows), na.rm = TRUE)
  }))
  rownames(mean_mat) <- ids

  list(mean_profiles = mean_mat)
}

# ── 4. Main loop ──────────────────────────────────────────────────────────────
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

for (dataset_name in names(input_files)) {
  fpath <- input_files[[dataset_name]]
  if (!file.exists(fpath)) {
    warning("File not found, skipping: ", fpath)
    next
  }

  message("\n▶  Processing: ", dataset_name)
  message("   Reading … ", fpath)
  raw_df <- data.table::fread(fpath)

  message("   Smoothing profiles …")
  df_smooth <- smooth_profiles(raw_df, half_win = 4L, min_win = 5L)
  message("   Done. ", nrow(df_smooth), " rows.")

  for (direction in c("oddsratinc", "oddsratdec")) {
    dir_suffix <- if (direction == "oddsratinc") "inc" else "dec"
    out_name   <- paste0("similar_binding_", dataset_name, "_", dir_suffix, ".qs")
    out_path   <- file.path(OUTPUT_DIR, out_name)

    message("   Building similarity object for direction: ", direction)
    sim_obj <- build_similarity_object(df_smooth, direction)

    if (!is.null(sim_obj)) {
      qs::qsave(sim_obj, out_path)
      message("   Saved → ", out_path,
              "  (", nrow(sim_obj$mean_profiles), " profiles)")
    } else {
      message("   Skipped (no data).")
    }
  }
}

message("\n✔  All similarity objects saved to: ", OUTPUT_DIR)
message("Files created:")
for (dataset_name in names(input_files)) {
  for (dir_suffix in c("inc", "dec")) {
    message("  similar_binding_", dataset_name, "_", dir_suffix, ".qs")
  }
}

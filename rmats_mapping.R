# =============================================================================
# rmats_mapping.R -- rMATS event IDs anchored on VAST-TOOLS event IDs (Splicing tab)
#
# CHARM's splicing data are quantified with VAST-TOOLS, so every event is identified by a VastDB ID
# (HsaEX..., HsaINT...). The rMATS IDs shown in the app are NOT measured: they come from the
# coordinate cross-reference built by eCLIPSE/map_vast_rmats_events.R and are attached to each VAST
# ID. This file holds the (small) helpers shared by the app and that builder:
#   rmats_prepare_map()   mapping table (one row per VAST/rMATS pair)  ->  compact table, one row per VAST ID
#   rmats_add_columns()   adds "rMATS ID" / "rMATS match" columns to any table that has an Event.ID column
#   rmats_label()         one-line label for a single event (plot subtitles)
# =============================================================================

RMATS_DISCLAIMER_SHORT <- paste(
  "rMATS IDs are not measured in CHARM: they are matched to each VAST-TOOLS event by genomic coordinates",
  "(\u00b11 nt). All dPSI values and statistics come from VAST-TOOLS.")

RMATS_DISCLAIMER_LONG <- paste(
  "The rMATS IDs shown are anchored on the VAST-TOOLS event IDs. They were assigned by matching splice-site",
  "coordinates (hg38, \u00b11 nt) between the VastDB event table and the rMATS annotation (GENCODE v42; exon",
  "skipping and intron retention only). dPSI values and p-values always come from VAST-TOOLS, not from rMATS.",
  "One VAST-TOOLS event can match several rMATS events, and some have no match. \"Same exon, different flanks\"",
  "means the alternative exon is identical but rMATS pairs it with different neighbouring exons.")

# labels for the match classes written by map_vast_rmats_events.R
RMATS_MATCH_LABELS <- c(identical                    = "identical exon and flanks",
                        exon_identical_flanks_differ = "same exon, different flanks",
                        within_tolerance             = "within \u00b11 nt")

# mapping: data.frame with vast_id, rmats_id, event_type ("ES"/"IR"), match_class (as written by the builder)
# Returns one row per VAST event: vast_id, rmats_ids ("SE:123; SE:456"), rmats_match, n_rmats
rmats_prepare_map <- function(mapping, max_ids = 5L) {
  m <- as.data.frame(mapping)[, c("vast_id", "rmats_id", "event_type", "match_class")]
  m$rank <- match(m$match_class, names(RMATS_MATCH_LABELS))
  m$tag  <- sprintf("%s:%s", ifelse(m$event_type == "ES", "SE", "RI"), m$rmats_id)
  m <- m[order(m$vast_id, m$rank, suppressWarnings(as.numeric(m$rmats_id))), ]
  pieces <- split(m, m$vast_id)
  out <- lapply(pieces, function(d) {
    ids <- unique(d$tag)
    shown <- if (length(ids) > max_ids) c(ids[seq_len(max_ids)], sprintf("+%d more", length(ids) - max_ids)) else ids
    data.frame(vast_id = d$vast_id[1], rmats_ids = paste(shown, collapse = "; "),
               rmats_match = unname(RMATS_MATCH_LABELS[names(RMATS_MATCH_LABELS)[min(d$rank)]]),
               n_rmats = length(ids), stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, out); rownames(out) <- NULL
  out
}

# Adds "rMATS ID" and "rMATS match" right after the Event.ID column (row order is untouched).
# note = TRUE also adds a column carrying the disclaimer (used for downloads).
rmats_add_columns <- function(tbl, map, id_col = "Event.ID", note = FALSE) {
  if (is.null(tbl) || is.null(map) || !id_col %in% names(tbl)) return(tbl)
  i <- match(tbl[[id_col]], map$vast_id)
  tbl <- as.data.frame(tbl, check.names = FALSE)
  add <- data.frame(`rMATS ID`    = ifelse(is.na(i), "no match", map$rmats_ids[i]),
                    `rMATS match` = ifelse(is.na(i), "\u2013",   map$rmats_match[i]),
                    check.names = FALSE, stringsAsFactors = FALSE)
  if (note) add[["rMATS note"]] <- "rMATS IDs mapped from VAST-TOOLS IDs by genomic coordinates (+/-1 nt); dPSI and statistics are from VAST-TOOLS"
  pos <- match(id_col, names(tbl))
  cbind(tbl[seq_len(pos)], add, tbl[setdiff(seq_along(tbl), seq_len(pos))], stringsAsFactors = FALSE)
}

# One-line label for a single VAST event, e.g. "SE:37285; SE:37286 (identical exon and flanks)"
rmats_label <- function(event_id, map) {
  if (is.null(map)) return(NULL)
  i <- match(event_id, map$vast_id)
  if (is.na(i)) return("no rMATS match")
  sprintf("%s (%s)", map$rmats_ids[i], map$rmats_match[i])
}

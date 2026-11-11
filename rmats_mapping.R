# =============================================================================
# rmats_mapping.R -- rMATS event IDs anchored on VAST-TOOLS event IDs (Splicing tab)
#
# CHARM's splicing data are quantified with VAST-TOOLS, so every event is identified by a VastDB ID
# (HsaEX..., HsaINT...). The rMATS IDs shown in the app are NOT measured: they come from the
# coordinate cross-reference built by eCLIPSE/map_vast_rmats_events.R and are attached to each VAST
# ID. This file holds the (small) helpers shared by the app and that builder:
#   rmats_prepare_map()    mapping table (one row per VAST/rMATS pair) -> list(forward, reverse)
#   rmats_unpack()         accepts that list (or an older forward-only table) and returns both parts
#   rmats_add_columns()    adds an "rMATS ID" column to any table that has an Event.ID column
#   rmats_label()          one-line label for a single event (plot subtitles)
#   rmats_event_choices()  search-box choices: VAST IDs plus the rMATS IDs that map onto them
#   rmats_resolve()        rMATS ID typed by the user -> the VAST-TOOLS event it is anchored on
# =============================================================================

RMATS_DISCLAIMER_LONG <- paste(
  "rMATS IDs are not measured in CHARM: they are anchored on the VAST-TOOLS event IDs and were assigned by",
  "matching splice-site coordinates (hg38, \u00b11 nt) between the VastDB event table and the rMATS annotation",
  "(GENCODE v42; exon skipping and intron retention only). The match is based on the alternative exon (or intron);",
  "rMATS may pair it with different neighbouring exons. dPSI values and p-values always come from VAST-TOOLS.",
  "One VAST-TOOLS event can match several rMATS events, and some have no match.")

RMATS_DISCLAIMER_SHORT <- "rMATS IDs are mapped from the VAST-TOOLS IDs by genomic coordinates (\u00b11 nt); dPSI and statistics are from VAST-TOOLS"

# ranking of the match classes written by map_vast_rmats_events.R (best first)
RMATS_MATCH_ORDER <- c("identical", "exon_identical_flanks_differ", "within_tolerance")
RMATS_MATCH_LABELS <- c(identical = "identical exon and flanks", exon_identical_flanks_differ = "same exon, different flanks",
                        within_tolerance = "within \u00b11 nt")

# mapping: data.frame with vast_id, rmats_id, event_type ("ES"/"IR"), match_class (as written by the builder)
# Returns list(forward = one row per VAST event [vast_id, rmats_ids "SE:123; SE:456", rmats_match, n_rmats],
#              reverse = one row per VAST/rMATS pair [rmats_tag, vast_id, rank, best]) ; best = the VAST event an
#              rMATS ID is anchored on when it matches several.
rmats_prepare_map <- function(mapping, max_ids = 5L) {
  m <- as.data.frame(mapping)[, c("vast_id", "rmats_id", "event_type", "match_class")]
  m$rank <- match(m$match_class, RMATS_MATCH_ORDER)
  m$tag  <- sprintf("%s:%s", ifelse(m$event_type == "ES", "SE", "RI"), m$rmats_id)
  m <- m[order(m$vast_id, m$rank, suppressWarnings(as.numeric(m$rmats_id))), ]
  fw <- lapply(split(m, m$vast_id), function(d) {
    ids <- unique(d$tag)
    shown <- if (length(ids) > max_ids) c(ids[seq_len(max_ids)], sprintf("+%d more", length(ids) - max_ids)) else ids
    data.frame(vast_id = d$vast_id[1], rmats_ids = paste(shown, collapse = "; "),
               rmats_match = unname(RMATS_MATCH_LABELS[RMATS_MATCH_ORDER[min(d$rank)]]),
               n_rmats = length(ids), stringsAsFactors = FALSE)
  })
  fw <- do.call(rbind, fw); rownames(fw) <- NULL
  rv <- data.frame(rmats_tag = m$tag, vast_id = m$vast_id, rank = m$rank, stringsAsFactors = FALSE)
  rv <- rv[order(rv$rmats_tag, rv$rank, rv$vast_id), ]
  rv$best <- !duplicated(rv$rmats_tag); rownames(rv) <- NULL
  list(forward = fw, reverse = rv)
}

# Older forward-only tables: rebuild the reverse lookup from the "SE:123; SE:456" strings
rmats_build_reverse <- function(forward) {
  tags <- strsplit(forward$rmats_ids, "; ", fixed = TRUE)
  n <- lengths(tags)
  rv <- data.frame(rmats_tag = unlist(tags), vast_id = rep(forward$vast_id, n),
                   rank = rep(match(forward$rmats_match, RMATS_MATCH_LABELS), n), stringsAsFactors = FALSE)
  rv <- rv[grepl("^(SE|RI):", rv$rmats_tag), ]                      # drops the "+N more" placeholders
  rv$rank[is.na(rv$rank)] <- length(RMATS_MATCH_ORDER) + 1L
  rv <- rv[order(rv$rmats_tag, rv$rank, rv$vast_id), ]
  rv$best <- !duplicated(rv$rmats_tag); rownames(rv) <- NULL
  rv
}

rmats_unpack <- function(x) {
  if (is.null(x)) return(list(forward = NULL, reverse = NULL))
  if (is.data.frame(x)) return(list(forward = x, reverse = rmats_build_reverse(x)))
  list(forward = x$forward, reverse = x$reverse)
}

# Adds an "rMATS ID" column right after the Event.ID column (row order is untouched).
# note = TRUE also adds a column carrying the disclaimer (used for downloads).
rmats_add_columns <- function(tbl, map, id_col = "Event.ID", note = FALSE) {
  if (is.null(tbl) || is.null(map) || !id_col %in% names(tbl)) return(tbl)
  i <- match(tbl[[id_col]], map$vast_id)
  tbl <- as.data.frame(tbl, check.names = FALSE)
  add <- data.frame(`rMATS ID` = ifelse(is.na(i), "no match", map$rmats_ids[i]), check.names = FALSE, stringsAsFactors = FALSE)
  if (note) add[["rMATS note"]] <- "rMATS IDs mapped from VAST-TOOLS IDs by genomic coordinates (+/-1 nt); dPSI and statistics are from VAST-TOOLS"
  pos <- match(id_col, names(tbl))
  cbind(tbl[seq_len(pos)], add, tbl[setdiff(seq_along(tbl), seq_len(pos))], stringsAsFactors = FALSE)
}

# One-line label for a single VAST event, e.g. "SE:37285; SE:37286" (or "no rMATS match")
rmats_label <- function(event_id, map) {
  if (is.null(map)) return(NULL)
  i <- match(event_id, map$vast_id)
  if (is.na(i)) "no rMATS match" else map$rmats_ids[i]
}

# Choices for the Event ID search box: the dataset's VAST IDs plus every rMATS ID anchored on one of them
# (labelled "SE:37285  (HsaEX0000025)" so it is clear which VAST event it opens).
rmats_event_choices <- function(event_ids, reverse) {
  event_ids <- as.character(event_ids)
  base <- stats::setNames(event_ids, event_ids)
  if (is.null(reverse)) return(base)
  r <- reverse[reverse$vast_id %in% event_ids, ]
  r <- r[order(r$rmats_tag, r$rank, r$vast_id), ]
  r <- r[!duplicated(r$rmats_tag), ]
  r <- r[order(r$rmats_tag), ]
  c(base, stats::setNames(r$rmats_tag, sprintf("%s  (%s)", r$rmats_tag, r$vast_id)))
}

# rMATS ID -> list(vast_id, tag, others) ; others = further VAST events the same rMATS ID matches. NULL if unknown.
rmats_resolve <- function(x, reverse, event_ids = NULL) {
  if (is.null(reverse) || length(x) != 1 || is.na(x)) return(NULL)
  r <- reverse[reverse$rmats_tag == x, ]
  if (!is.null(event_ids)) r <- r[r$vast_id %in% event_ids, ]
  if (nrow(r) == 0) return(NULL)
  r <- r[order(r$rank, r$vast_id), ]
  list(vast_id = r$vast_id[1], tag = x, others = setdiff(r$vast_id, r$vast_id[1]))
}

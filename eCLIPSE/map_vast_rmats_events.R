# =============================================================================
# map_vast_rmats_events.R
#   Cross-reference VAST-TOOLS and rMATS event IDs by genomic coordinates (SE and RI).
#
#   Rscript map_vast_rmats_events.R          (or source() it from the CHARM folder in R)
#
# WHY COORDINATES: VAST-TOOLS IDs (HsaEX0067681, HsaINT0140508) are stable VastDB identifiers,
# rMATS IDs (1, 2, 3 ...) are row numbers that restart in every run and every event type.
# The only thing the two share is the genome, so events are matched on splice-site coordinates.
#
# CONVENTIONS (checked on the eCLIPSE data, see the README block at the end of this file)
#   * Common format everywhere: 1-based, closed [start, end], genomic order (start <= end),
#     strand "+" / "-".
#   * VAST-TOOLS is already 1-based closed. rMATS is 0-based half-open: every *_0base START is
#     shifted by +1, every END is left unchanged.
#   * Both tools list the flanking exons in genomic order (upstream = lower coordinate, also
#     on the "-" strand).
#
# WHAT IS MATCHED
#   ES (exon skipping)    VAST HsaEX*  <->  rMATS SE.MATS.JC.txt
#       VAST gives, for every event, a SET of possible sites (e.g. "100+250" = two alternative
#       donors). An rMATS exon matches when its start is one of the VAST exon starts AND its end
#       is one of the VAST exon ends. The flanking exons are compared too, but they only
#       change the match class, they never prevent a match (rMATS lists one SE per flank pair).
#   IR (intron retention) VAST HsaINT* <->  rMATS RI.MATS.JC.txt
#       retained intron = [upstreamEE + 1, downstreamES(0-based)]  vs  VAST COORD.
#
# TOLERANCE: every coordinate may differ by up to MAP_TOL nt (default 1, to absorb 0/1-based slips).
#   Each pair records the exact offset used, so tolerant matches can be filtered out later.
#
# OUTPUT (MAP_OUT_DIR): Event_ID_Mapping.csv (+ .rds), Event_ID_Mapping_summary.csv,
#   Event_ID_Mapping_by_vast_class.csv, Event_ID_Mapping_offsets.csv, Event_ID_Unmapped_rmats.csv
#   Key columns are kept as-is: vast_id, rmats_id (+ rmats_key = "<type>_<id>@<run>").
# =============================================================================

# ---- settings (set MAP_* variables before sourcing to override; otherwise the defaults below are used) ----
# Defaults chosen by a previous run are discarded first, so re-sourcing always picks up fresh defaults.
if (exists(".MAP_AUTO", envir = globalenv())) rm(list = intersect(get(".MAP_AUTO", envir = globalenv()), ls(globalenv(), all.names = TRUE)), envir = globalenv())
.MAP_AUTO <- character()
.set <- function(name, value) {                       # value is evaluated lazily, after earlier defaults exist
  if (!exists(name, envir = globalenv(), inherits = FALSE)) { assign(name, value, envir = globalenv()); .MAP_AUTO <<- c(.MAP_AUTO, name) }
}
.set("MAP_DIR", "eCLIPSE")
.set("MAP_VAST_FILE", file.path(MAP_DIR, "INCLUSION_LEVELS_FULL-hg38-25-v251.tab.gz"))
# rMATS event universe: fromGTF.<type>.txt = EVERY annotated event (independent of the reads); falls back
# to <type>.MATS.JC.txt (only events that passed the read filters) when the fromGTF file is absent.
.pick <- function(a, b) if (file.exists(a)) a else b
.set("MAP_SE_FILE", .pick(file.path(MAP_DIR, "fromGTF.SE.txt"), file.path(MAP_DIR, "SE.MATS.JC.txt")))
.set("MAP_RI_FILE", .pick(file.path(MAP_DIR, "fromGTF.RI.txt"), file.path(MAP_DIR, "RI.MATS.JC.txt")))   # missing file = skipped
# optional: the quantified output of the same run (adds FDR / IncLevelDifference and flags quantified events)
.set("MAP_SE_QUANT_FILE", file.path(MAP_DIR, "SE.MATS.JC.txt"))
.set("MAP_RI_QUANT_FILE", file.path(MAP_DIR, "RI.MATS.JC.txt"))
.set("MAP_OUT_DIR", MAP_DIR)
.set("MAP_APP_DIR", "data/QS2_Files")   # where the Shiny app reads RMATS_Event_Mapping.qs2 (NULL = do not write)
.set("MAP_TOL", 1L)                 # max offset (nt) per coordinate
.set("MAP_MIN_LEN", 0L)             # drop SE events shorter than this (the eCLIPSE tables use 10)
.set("MAP_RUN_LABEL", "rMATS")      # tags rMATS IDs, which are only unique within one run
.set("MAP_MICROEXON_MAX", 27L)      # VAST's microexon range (3-27 nt)
.MAP_AUTO <- get(".MAP_AUTO")       # (kept in the global environment for the next run)
cat("Mapping settings: SE =", MAP_SE_FILE, "| RI =", MAP_RI_FILE, "| tolerance =", MAP_TOL, "nt\n")

suppressPackageStartupMessages(library(data.table))

# ---- small helpers -------------------------------------------------------------------------
msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "|", sprintf(...), "\n")

# "100+250" -> long table (id, pos): one row per alternative site
expand_sites <- function(id, txt) {
  l <- strsplit(txt, "+", fixed = TRUE)
  data.table(id = rep(id, lengths(l)), pos = as.integer(unlist(l)))
}

# smallest |pos - value| over the sites of each pair's VAST event
min_abs_diff <- function(pairs, sites, value_col) {
  x <- merge(pairs[, .(rmats_key, vast_id, v = get(value_col))], sites, by.x = "vast_id", by.y = "id",
             allow.cartesian = TRUE)
  x[, .(d = min(abs(pos - v))), by = .(rmats_key, vast_id)]
}

# ---- readers -------------------------------------------------------------------------------
read_vast <- function(file) {
  if (!file.exists(file)) stop("VAST file not found: ", file, "\n  working directory: ", getwd(), "\n  set MAP_DIR (folder) or MAP_VAST_FILE (full path) before sourcing")
  msg("reading VAST-TOOLS table (first 6 columns only) ...")
  v <- tryCatch(fread(cmd = paste("gzip -dc", shQuote(file), "| cut -f1-6"), showProgress = FALSE),
                error = function(e) fread(file, select = 1:6, showProgress = FALSE))
  need <- c("GENE", "EVENT", "COORD", "LENGTH", "FullCO", "COMPLEX")
  if (!all(need %in% names(v))) stop("VAST table lacks columns: ", paste(setdiff(need, names(v)), collapse = ", "))
  v[, .SD, .SDcols = need]
}

# Exon skipping: FullCO = "chr:<C1 donors>,<C2 acceptors>-<C2 donors>,<C3 acceptors>", sites joined by "+",
# written in TRANSCRIPT order. Strand = "-" when the first flank lies at higher coordinates than the last.
# Genomic left flank  = C1 donors (+ strand) or C3 acceptors (- strand); right flank = the other one.
parse_vast_ex <- function(v) {
  ex <- v[grepl("^HsaEX", EVENT)]
  p  <- tstrsplit(ex$FullCO, "[:,]")
  if (length(p) != 4) stop("unexpected FullCO layout for HsaEX events")
  bc <- tstrsplit(p[[3]], "-", fixed = TRUE)
  ex[, `:=`(chr = p[[1]], A = p[[2]], B = bc[[1]], C = bc[[2]], D = p[[4]])]
  ex[, strand := ifelse(as.integer(sub("\\+.*", "", A)) > as.integer(sub("\\+.*", "", D)), "-", "+")]
  ev <- ex[, .(vast_id = EVENT, gene_vast = GENE, chr, strand, vast_length = LENGTH, vast_complex = COMPLEX)]
  ev[, is_microexon := vast_complex == "MIC" | vast_length <= MAP_MICROEXON_MAX]
  starts <- expand_sites(ex$EVENT, ex$B); setnames(starts, "pos", "st")
  ends   <- expand_sites(ex$EVENT, ex$C); setnames(ends,   "pos", "en")
  alts <- merge(starts, ends, by = "id", allow.cartesian = TRUE)         # every start x end the event allows
  alts <- merge(alts, ev[, .(id = vast_id, chr, strand)], by = "id")
  setnames(alts, "id", "vast_id")
  a <- expand_sites(ex$EVENT, ex$A); d <- expand_sites(ex$EVENT, ex$D)
  neg <- ex[strand == "-", EVENT]
  left  <- rbind(a[!id %in% neg], d[id %in% neg])
  right <- rbind(d[!id %in% neg], a[id %in% neg])
  list(events = ev, alts = alts, left = left, right = right)
}

# Intron retention: COORD = the intron itself ("chr:start-end"), strand is the last field of FullCO.
parse_vast_ir <- function(v) {
  ir <- v[grepl("^HsaINT", EVENT)]
  co <- tstrsplit(ir$COORD, "[:-]")
  ir[, `:=`(chr = co[[1]], st = as.integer(co[[2]]), en = as.integer(co[[3]]),
            strand = sub(".*:", "", FullCO))]
  ir <- ir[strand %in% c("+", "-")]
  ev <- ir[, .(vast_id = EVENT, gene_vast = GENE, chr, strand, vast_length = LENGTH, vast_complex = COMPLEX,
               st, en)]
  ev[, is_microexon := FALSE]
  ev
}

read_rmats <- function(file, type, quant_file = NULL) {
  if (is.null(file) || !file.exists(file)) { msg("%s: file not found, skipped (%s)", type, file); return(NULL) }
  r <- fread(file, showProgress = FALSE)
  setnames(r, make.unique(names(r)))                                    # rMATS repeats the "ID" column
  base <- c("ID", "GeneID", "geneSymbol", "chr", "strand", "upstreamEE", "downstreamES")
  need <- c(base, if (type == "SE") c("exonStart_0base", "exonEnd") else c("riExonStart_0base", "riExonEnd"))
  if (!all(need %in% names(r))) stop(type, " file lacks columns: ", paste(setdiff(need, names(r)), collapse = ", "))
  out <- r[, .(rmats_id = ID, gene_id = GeneID, gene_rmats = geneSymbol, chr, strand,
               FDR = if ("FDR" %in% names(r)) FDR else NA_real_,
               IncLevelDifference = if ("IncLevelDifference" %in% names(r)) IncLevelDifference else NA_real_,
               up_end = upstreamEE, dn_start = downstreamES + 1L)]      # 1-based: last base / first base of the flanking exons
  if (type == "SE") {
    out[, `:=`(st = r$exonStart_0base + 1L, en = r$exonEnd)]            # 0-based start -> 1-based; end unchanged
  } else {
    out[, `:=`(st = r$upstreamEE + 1L, en = r$downstreamES)]            # retained intron: between the two flanking exons
  }
  out[, `:=`(event_type = ifelse(type == "SE", "ES", "IR"),
             rmats_key  = sprintf("%s_%s@%s", type, rmats_id, MAP_RUN_LABEL),
             rmats_len  = en - st + 1L)]
  # same run's quantified table (IDs are identical to the fromGTF IDs): FDR, dPSI and a "quantified" flag
  qf <- !is.null(quant_file) && file.exists(quant_file)
  if (qf && normalizePath(quant_file) != normalizePath(file)) {
    q <- fread(quant_file, showProgress = FALSE); setnames(q, make.unique(names(q)))
    i <- match(out$rmats_id, q$ID)
    out[, `:=`(FDR = q$FDR[i], IncLevelDifference = q$IncLevelDifference[i], quantified_in_run = !is.na(i))]
  } else out[, quantified_in_run := TRUE]
  if (type == "SE") out <- out[rmats_len >= MAP_MIN_LEN]
  out
}

# ---- matching ------------------------------------------------------------------------------
# alts: (vast_id, chr, strand, st, en) one row per start x end the VAST event allows.
# Returns one row per (rmats_key, vast_id) with the smallest offset found within +-MAP_TOL.
match_coords <- function(rm, alts, tol) {
  shifts <- CJ(ds = -tol:tol, de = -tol:tol)
  hits <- rbindlist(lapply(seq_len(nrow(shifts)), function(i) {
    x <- rm[, .(rmats_key, chr, strand, st = st + shifts$ds[i], en = en + shifts$de[i])]
    m <- merge(x, alts, by = c("chr", "strand", "st", "en"), allow.cartesian = TRUE)
    if (nrow(m) == 0) return(NULL)
    m[, `:=`(offset_start = -shifts$ds[i], offset_end = -shifts$de[i])]  # rMATS - VAST
    m
  }))
  if (is.null(hits) || nrow(hits) == 0) return(NULL)
  setnames(hits, c("st", "en"), c("vast_start", "vast_end"))
  hits[, off_sum := abs(offset_start) + abs(offset_end)]
  setorder(hits, rmats_key, vast_id, off_sum)
  hits[, .SD[1], by = .(rmats_key, vast_id)][, off_sum := NULL][]
}

CLASS_LEVELS <- c("identical", "exon_identical_flanks_differ", "within_tolerance")

map_se <- function(rm, vx, tol) {
  hits <- match_coords(rm, vx$alts, tol)
  if (is.null(hits)) return(NULL)
  hits <- merge(hits, rm[, .(rmats_key, rmats_up_end = up_end, rmats_down_start = dn_start)], by = "rmats_key")
  fu <- min_abs_diff(hits, vx$left,  "rmats_up_end");   setnames(fu, "d", "flank_up_diff")
  fd <- min_abs_diff(hits, vx$right, "rmats_down_start"); setnames(fd, "d", "flank_down_diff")
  hits <- merge(merge(hits, fu, by = c("rmats_key", "vast_id"), all.x = TRUE), fd, by = c("rmats_key", "vast_id"), all.x = TRUE)
  hits[, exon_off := abs(offset_start) + abs(offset_end)]
  hits[, match_class := fifelse(exon_off == 0 & flank_up_diff == 0 & flank_down_diff == 0, CLASS_LEVELS[1],
                        fifelse(exon_off == 0, CLASS_LEVELS[2], CLASS_LEVELS[3]))]
  hits[, exon_off := NULL]
  hits
}

map_ir <- function(rm, ir, tol) {
  alts <- ir[, .(vast_id, chr, strand, st, en)]
  hits <- match_coords(rm, alts, tol)
  if (is.null(hits)) return(NULL)
  hits <- merge(hits, rm[, .(rmats_key, rmats_up_end = up_end, rmats_down_start = dn_start)], by = "rmats_key")
  hits[, `:=`(flank_up_diff = NA_integer_, flank_down_diff = NA_integer_)]
  hits[, match_class := fifelse(offset_start == 0 & offset_end == 0, CLASS_LEVELS[1], CLASS_LEVELS[3])]
  hits
}

# one best VAST partner per rMATS event and vice versa (class, then offset, then flank agreement)
rank_pairs <- function(map) {
  map[, score := match(match_class, CLASS_LEVELS) * 1e6 + (abs(offset_start) + abs(offset_end)) * 1e3 +
                 pmin(fifelse(is.na(flank_up_diff), 0L, flank_up_diff) + fifelse(is.na(flank_down_diff), 0L, flank_down_diff), 999L)]
  setorder(map, score, vast_id, rmats_key)
  map[, rank_for_rmats := seq_len(.N), by = rmats_key]
  map[, rank_for_vast  := seq_len(.N), by = vast_id]
  map[, `:=`(n_vast_for_rmats = .N), by = rmats_key]
  map[, `:=`(n_rmats_for_vast = .N), by = vast_id]
  map[, `:=`(best_for_rmats = rank_for_rmats == 1L, best_for_vast = rank_for_vast == 1L)]
  map[, reciprocal_best := best_for_rmats & best_for_vast]
  map[, c("score", "rank_for_rmats", "rank_for_vast") := NULL]
  map[]
}

# unmatched rMATS events: does any VAST event at least overlap them on the same strand?
overlap_diagnostic <- function(unm, ranges) {
  if (nrow(unm) == 0) return(character(0))
  q <- unm[, .(rmats_key, chr, strand, st, en)]
  setkey(ranges, chr, strand, st, en)
  o <- foverlaps(q, ranges, by.x = c("chr", "strand", "st", "en"), type = "any", nomatch = NULL)
  fifelse(unm$rmats_key %in% o$rmats_key, "overlaps_other_VAST_event", "no_VAST_event_overlap")
}

# ---- run ----------------------------------------------------------------------------------
run_mapping <- function() {
  vast <- read_vast(MAP_VAST_FILE)
  vx   <- parse_vast_ex(vast)
  vir  <- parse_vast_ir(vast)
  msg("VAST: %s HsaEX events, %s HsaINT events", format(nrow(vx$events), big.mark = ","), format(nrow(vir), big.mark = ","))
  if (MAP_MIN_LEN > 0) { vx$events <- vx$events[vast_length >= MAP_MIN_LEN]; keep <- vx$events$vast_id
    vx$alts <- vx$alts[vast_id %in% keep]; vir <- vir[vast_length >= MAP_MIN_LEN] }

  maps <- list(); unmapped <- list(); n_rm <- list(); n_vast_total <- c(ES = nrow(vx$events), IR = nrow(vir))

  # ---- ES
  rs <- read_rmats(MAP_SE_FILE, "SE", MAP_SE_QUANT_FILE)
  if (!is.null(rs)) {
    n_rm$ES <- nrow(rs); msg("ES: %s rMATS SE events, matching (tolerance %d nt) ...", format(nrow(rs), big.mark = ","), MAP_TOL)
    h <- map_se(rs, vx, MAP_TOL)
    h <- merge(h, vx$events[, .(vast_id, gene_vast, vast_complex, vast_length, is_microexon)], by = "vast_id")
    maps$ES <- merge(h, rs[, .(rmats_key, rmats_id, gene_rmats, gene_id, rmats_start = st, rmats_end = en, FDR, IncLevelDifference, quantified_in_run)], by = "rmats_key")
    maps$ES[, event_type := "ES"]
    unm <- rs[!rmats_key %in% maps$ES$rmats_key]
    rng <- vx$alts[, .(st = min(st), en = max(en)), by = .(vast_id, chr, strand)]
    unm[, reason := overlap_diagnostic(unm, rng)]
    unmapped$ES <- unm
  }
  # ---- IR
  ri <- read_rmats(MAP_RI_FILE, "RI", MAP_RI_QUANT_FILE)
  if (!is.null(ri)) {
    n_rm$IR <- nrow(ri); msg("IR: %s rMATS RI events, matching (tolerance %d nt) ...", format(nrow(ri), big.mark = ","), MAP_TOL)
    h <- map_ir(ri, vir, MAP_TOL)
    if (!is.null(h)) {
      h <- merge(h, vir[, .(vast_id, gene_vast, vast_complex, vast_length, is_microexon)], by = "vast_id")
      maps$IR <- merge(h, ri[, .(rmats_key, rmats_id, gene_rmats, gene_id, rmats_start = st, rmats_end = en, FDR, IncLevelDifference, quantified_in_run)], by = "rmats_key")
      maps$IR[, event_type := "IR"]
    }
    unm <- ri[!rmats_key %in% (if (is.null(maps$IR)) character(0) else maps$IR$rmats_key)]
    rng <- vir[, .(vast_id, chr, strand, st, en)]
    unm[, reason := overlap_diagnostic(unm, rng)]
    unmapped$IR <- unm
  }
  if (length(maps) == 0) stop("nothing was matched - check the input files")

  # ---- tidy mapping table (one row per matched pair)
  map <- rbindlist(maps, use.names = TRUE, fill = TRUE)
  map[, `:=`(rmats_run = MAP_RUN_LABEL, pair_id = paste(vast_id, rmats_key, sep = "|"))]
  map <- rbindlist(lapply(split(map, map$event_type), rank_pairs))
  setcolorder(map, c("event_type", "pair_id", "vast_id", "rmats_id", "rmats_key", "rmats_run", "match_class",
                     "best_for_rmats", "best_for_vast", "reciprocal_best", "n_vast_for_rmats", "n_rmats_for_vast",
                     "gene_vast", "gene_rmats", "chr", "strand", "rmats_start", "rmats_end", "vast_start", "vast_end",
                     "offset_start", "offset_end", "flank_up_diff", "flank_down_diff", "vast_complex", "vast_length",
                     "is_microexon", "gene_id", "rmats_up_end", "rmats_down_start", "FDR", "IncLevelDifference", "quantified_in_run"))
  setorder(map, event_type, vast_id, -best_for_vast, rmats_key)
  map[, gene_match := fifelse(is.na(gene_vast) | is.na(gene_rmats), NA, toupper(gene_vast) == toupper(gene_rmats))]

  # ---- summaries
  summ <- rbindlist(lapply(names(maps), function(t) {
    m <- map[event_type == t]; b <- m[best_for_rmats == TRUE]; u <- unmapped[[t]]
    nq <- if (t == "ES") sum(rs$quantified_in_run) else sum(ri$quantified_in_run)
    mq <- uniqueN(m[quantified_in_run == TRUE, rmats_key])
    data.table(event_type = t, rmats_events = n_rm[[t]], rmats_quantified = nq, rmats_quantified_mapped = mq,
               rmats_mapped = uniqueN(m$rmats_key),
               pct_rmats_mapped = round(100 * uniqueN(m$rmats_key) / n_rm[[t]], 1),
               best_identical = sum(b$match_class == CLASS_LEVELS[1]),
               best_exon_identical_flanks_differ = sum(b$match_class == CLASS_LEVELS[2]),
               best_within_tolerance = sum(b$match_class == CLASS_LEVELS[3]),
               rmats_unmapped_overlapping_VAST = sum(u$reason == "overlaps_other_VAST_event"),
               rmats_unmapped_no_VAST_overlap = sum(u$reason == "no_VAST_event_overlap"),
               vast_events_total = n_vast_total[[t]], vast_events_mapped = uniqueN(m$vast_id),
               pairs = nrow(m), reciprocal_best_pairs = sum(m$reciprocal_best),
               rmats_with_multiple_vast = sum(b$n_vast_for_rmats > 1),
               vast_with_multiple_rmats = uniqueN(m[n_rmats_for_vast > 1, vast_id]),
               gene_symbol_agrees_pct = round(100 * mean(b$gene_match, na.rm = TRUE), 1))
  }))
  by_class <- map[best_for_rmats == TRUE, .(rmats_events = .N, identical = sum(match_class == CLASS_LEVELS[1])),
                  by = .(event_type, vast_complex)][order(event_type, -rmats_events)]
  offs <- map[, .N, by = .(event_type, offset_start, offset_end)][order(event_type, -N)]
  setnames(offs, "N", "pairs")
  unm_all <- rbindlist(lapply(names(unmapped), function(t) unmapped[[t]][, .(event_type, rmats_id, rmats_key, gene_rmats, chr, strand, start = st, end = en, reason)]))

  dir.create(MAP_OUT_DIR, showWarnings = FALSE, recursive = TRUE)
  fwrite(map,      file.path(MAP_OUT_DIR, "Event_ID_Mapping.csv"))
  fwrite(summ,     file.path(MAP_OUT_DIR, "Event_ID_Mapping_summary.csv"))
  fwrite(by_class, file.path(MAP_OUT_DIR, "Event_ID_Mapping_by_vast_class.csv"))
  fwrite(offs,     file.path(MAP_OUT_DIR, "Event_ID_Mapping_offsets.csv"))
  fwrite(unm_all,  file.path(MAP_OUT_DIR, "Event_ID_Unmapped_rmats.csv"))
  res <- list(mapping = map, summary = summ, by_vast_class = by_class, offsets = offs, unmapped_rmats = unm_all,
              settings = list(MAP_TOL = MAP_TOL, MAP_MIN_LEN = MAP_MIN_LEN, MAP_RUN_LABEL = MAP_RUN_LABEL, built = Sys.time()))
  saveRDS(res, file.path(MAP_OUT_DIR, "Event_ID_Mapping.rds"))
  if (requireNamespace("qs2", quietly = TRUE)) qs2::qs_save(res, file.path(MAP_OUT_DIR, "Event_ID_Mapping.qs2"))
  msg("written to %s", normalizePath(MAP_OUT_DIR))

  # compact table for the app (one row per VAST event, "SE:123; SE:456") -> data/QS2_Files/RMATS_Event_Mapping.qs2
  if (!is.null(MAP_APP_DIR) && dir.exists(MAP_APP_DIR) && file.exists("rmats_mapping.R") && requireNamespace("qs2", quietly = TRUE)) {
    source("rmats_mapping.R", local = TRUE)
    app_map <- rmats_prepare_map(map)
    qs2::qs_save(app_map, file.path(MAP_APP_DIR, "RMATS_Event_Mapping.qs2"))
    msg("app file written: %s (%s VAST events with at least one rMATS match)",
        file.path(MAP_APP_DIR, "RMATS_Event_Mapping.qs2"), format(nrow(app_map), big.mark = ","))
  } else msg("app file NOT written (needs: run from the CHARM folder, %s present, qs2 installed)", MAP_APP_DIR)
  res
}

# ---- lookups for the Shiny app / downstream binding analysis -----------------------------------
# VAST IDs -> rMATS IDs (and back). only_best = TRUE keeps one partner per event.
lookup_rmats_from_vast <- function(vast_ids, mapping, only_best = TRUE, min_class = "within_tolerance") {
  m <- mapping[vast_id %in% vast_ids & match(match_class, CLASS_LEVELS) <= match(min_class, CLASS_LEVELS)]
  if (only_best) m <- m[best_for_vast == TRUE]
  m[, .(vast_id, rmats_id, rmats_key, match_class, event_type)]
}
lookup_vast_from_rmats <- function(rmats_ids, mapping, only_best = TRUE, min_class = "within_tolerance") {
  m <- mapping[rmats_id %in% rmats_ids & match(match_class, CLASS_LEVELS) <= match(min_class, CLASS_LEVELS)]
  if (only_best) m <- m[best_for_rmats == TRUE]
  m[, .(rmats_id, rmats_key, vast_id, match_class, event_type)]
}

if (!exists("MAPPING_NO_RUN")) {
  res <- run_mapping()
  cat("\n=== SUMMARY ===\n"); print(res$summary, width = 200)
  cat("\n=== OFFSETS USED (rMATS - VAST; 0/0 expected to dominate) ===\n"); print(res$offsets[, head(.SD, 5), by = event_type])
}

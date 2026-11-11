# =============================================================================
# build_gsea_reference.R -- run ONCE (and again whenever the expression/GSEA data change)
#   source("build_gsea_reference.R")          # from the CHARM folder, in R
# Adds the GSEA entries ("gsea|Both Cells", "gsea|K562", "gsea|HEPG2") to
# data/QS2_Files/Correlation_Reference.qs2 (the same file the expression and splicing
# entries live in; existing entries are kept).
#
# What it builds, per cell-line dataset
#   observed : every RBP-vs-RBP Spearman correlation of the hallmark NES profiles, exactly as
#              gsea_correl() computes them (from the stored RBPs.gsea_* tables)
#   null     : B times, draw one random gene set for EACH hallmark set, with the same size, from
#              the genes that are ranked in (almost) every RBP. Run fgsea on every RBP's ranked
#              t-statistics with these random sets, obtaining a (hallmark sets x RBPs) matrix of
#              null NES values, and correlate the RBP columns (Spearman) -- the same
#              calculation as for the real data, just with random signatures.
#              The null correlations of all replicates are pooled into one histogram.
#   Because the random sets are scored against each RBP's real ranking, genome-wide similarity
#   between two RBPs' rankings is kept in the null; only the hallmark biology is removed.
#
#   gene-set background (used by the "click a gene set" window of the Expression > GSEA plot): the SAME random
#   NES values are also kept per hallmark set -- a histogram (bin 0.02) of the NES that random gene sets of
#   that set's size get against ALL RBP rankings (identical for every RBP, so stored once) -- together with
#   the observed NES of every RBP for that set ("gsea_sets" element of each entry).
#
# Settings (set before sourcing, or edit): B = replicates (each = one random set per hallmark set;
#   100 as in the earlier mito script), CORES, ONLY (e.g. "gsea|K562"), SEED, MIN_SIZE.
# Needs: fgsea, msigdbr (the app already uses both), qs2.
# =============================================================================
if (!exists("B"))        B        <- 100
if (!exists("SEED"))     SEED     <- 1
if (!exists("ONLY"))     ONLY     <- NULL
if (!exists("MIN_SIZE")) MIN_SIZE <- 5
if (!exists("CORES"))    CORES    <- if (.Platform$OS.type == "windows") 1 else max(1, min(4, parallel::detectCores() - 1))
if (!exists("DATA_DIR")) DATA_DIR <- "data/QS2_Files"
UNIVERSE_FRACTION <- 0.95          # random genes are drawn from genes ranked in >= 95% of the RBPs
GSEA_PW_LO <- -6; GSEA_PW_HI <- 6; GSEA_PW_BIN <- 0.02   # NES histogram of the per-gene-set null
MIN_SHARED <- 3                    # same minimum overlap gsea_correl() requires
source("correlation_reference.R", local = TRUE)    # CORR_BIN_WIDTH

# NES of every pathway for ONE ranking; returns a named numeric vector (NA where fgsea drops the set).
# Replaceable (GSEA_NES_FUN) so the plumbing can be tested without fgsea.
if (!exists("GSEA_NES_FUN")) {
  GSEA_NES_FUN <- function(ranks, pathways) {
    res <- suppressWarnings(fgsea::fgsea(pathways = pathways, stats = ranks, minSize = MIN_SIZE, maxSize = 5000,
                                         eps = 1e-10, BPPARAM = BiocParallel::SerialParam()))
    stats::setNames(res$NES, res$pathway)
  }
}

gsea_to_counts <- function(r) {
  nb <- as.integer(round(2 / CORR_BIN_WIDTH))
  tabulate(pmin(pmax(floor((r + 1) / CORR_BIN_WIDTH) + 1, 1), nb), nbins = nb)
}

pairwise_spearman_nes <- function(M) {             # M: gene sets x RBPs
  n_shared <- crossprod(!is.na(M))
  cm <- suppressWarnings(stats::cor(M, method = "spearman", use = "pairwise.complete.obs"))
  cm[n_shared < MIN_SHARED] <- NA
  cm[upper.tri(cm)][is.finite(cm[upper.tri(cm)])]
}

fmt_dur <- function(sec) {
  sec <- round(sec)
  if (sec < 90) sprintf("%ds", sec) else if (sec < 5400) sprintf("%.0f min", sec / 60) else sprintf("%.1f h", sec / 3600)
}

# ranked statistics: list of data.frames(Gene, t) -> named numeric vectors (upper-case symbols)
prepare_ranks <- function(lst) {
  out <- lapply(lst, function(df) {
    if (is.null(df) || !all(c("Gene", "t") %in% colnames(df))) return(NULL)
    v <- stats::setNames(as.numeric(df$t), toupper(as.character(df$Gene)))
    v <- v[is.finite(v) & !duplicated(names(v))]
    if (length(v) < 100) NULL else v
  })
  out[!vapply(out, is.null, logical(1))]
}

# observed NES matrix (pathway x RBP) from the stored GSEA tables
nes_matrix <- function(lst) {
  vecs <- lapply(lst, function(df) {
    if (is.null(df) || !all(c("pathway", "NES") %in% colnames(df))) return(NULL)
    stats::setNames(as.numeric(df$NES), as.character(df$pathway))
  })
  vecs <- vecs[!vapply(vecs, is.null, logical(1))]
  paths <- unique(unlist(lapply(vecs, names), use.names = FALSE))
  M <- matrix(NA_real_, length(paths), length(vecs), dimnames = list(paths, names(vecs)))
  for (j in seq_along(vecs)) M[names(vecs[[j]]), j] <- vecs[[j]]
  M
}

build_gsea_entry <- function(label, gsea_list, rank_list, hallmark_sets) {
  t0 <- Sys.time()
  cat(sprintf("[%s] observed correlations\n", label))
  Mobs <- nes_matrix(gsea_list)
  observed <- sort(round(pairwise_spearman_nes(Mobs), 4))
  cat(sprintf("[%s] %d hallmark rows x %d RBPs -> %s observed pairs\n", label, nrow(Mobs), ncol(Mobs),
              format(length(observed), big.mark = ",")))

  ranks <- prepare_ranks(rank_list)
  ranks <- ranks[intersect(names(ranks), colnames(Mobs))]            # same RBPs as the observed matrix
  gene_freq <- table(unlist(lapply(ranks, names), use.names = FALSE))
  universe  <- names(gene_freq)[gene_freq >= UNIVERSE_FRACTION * length(ranks)]
  sizes <- vapply(hallmark_sets, function(g) sum(g %in% universe), integer(1))
  sizes <- sizes[sizes >= MIN_SIZE]
  K <- length(sizes)
  cat(sprintf("[%s] %d RBPs, universe %s genes, %d hallmark sets (sizes %d-%d)\n", label, length(ranks),
              format(length(universe), big.mark = ","), K, min(sizes), max(sizes)))

  # sanity check: our NES function reproduces the stored NES (tiny sampling noise expected)
  chk <- intersect(names(ranks), colnames(Mobs))[seq_len(min(3, length(ranks)))]
  for (r in chk) {
    nn <- GSEA_NES_FUN(ranks[[r]], lapply(hallmark_sets, toupper))
    names(nn) <- sub("^HALLMARK_", "", names(nn))
    sto <- Mobs[, r]; names(sto) <- sub("^HALLMARK_", "", names(sto))
    common <- intersect(names(nn), names(sto))
    cat(sprintf("[%s] check %s: NES recomputed vs stored, Spearman %.3f over %d sets\n", label, r,
                suppressWarnings(stats::cor(nn[common], sto[common], use = "complete.obs", method = "spearman")), length(common)))
  }

  # B x K random gene sets (same sets for every RBP)
  set.seed(SEED)
  rand_sets <- vector("list", B * K); nm <- character(B * K)
  for (b in seq_len(B)) for (k in seq_len(K)) {
    i <- (b - 1L) * K + k
    rand_sets[[i]] <- sample(universe, sizes[k]); nm[i] <- sprintf("r%d_%d", b, k)
  }
  names(rand_sets) <- nm

  one_rbp <- function(r) {
    set.seed(SEED + match(r, names(ranks)))
    nn <- GSEA_NES_FUN(ranks[[r]], rand_sets)
    as.numeric(nn[nm])
  }
  rbp_names <- names(ranks)
  chunks <- split(rbp_names, ceiling(seq_along(rbp_names) / max(1, CORES * 2)))
  nes_list <- list(); t1 <- Sys.time(); done <- 0
  for (ch in chunks) {
    res <- if (CORES > 1) parallel::mclapply(ch, one_rbp, mc.cores = CORES, mc.preschedule = FALSE) else lapply(ch, one_rbp)
    bad <- vapply(res, function(x) !is.numeric(x), logical(1))
    if (any(bad)) stop("fgsea failed for some RBPs in ", label, " (try CORES <- 1 to see the error)")
    names(res) <- ch; nes_list <- c(nes_list, res); done <- done + length(ch)
    el <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
    cat(sprintf("[%s]   null NES: %d / %d RBPs | elapsed %s | est. remaining %s (done ~%s)\n", label, done, length(rbp_names),
                fmt_dur(el), fmt_dur(el / done * (length(rbp_names) - done)),
                format(Sys.time() + el / done * (length(rbp_names) - done), "%H:%M")))
  }
  NES <- do.call(cbind, nes_list)                                   # (B*K) x RBPs
  counts <- numeric(as.integer(round(2 / CORR_BIN_WIDTH)))
  for (b in seq_len(B)) {
    Mb <- NES[(b - 1L) * K + seq_len(K), , drop = FALSE]
    counts <- counts + gsea_to_counts(pairwise_spearman_nes(Mb))
  }
  cat(sprintf("[%s] null: %d replicates x %d sets, %s null pairs (%s total)\n", label, B, K,
              format(sum(counts), big.mark = ","), fmt_dur(as.numeric(difftime(Sys.time(), t0, units = "secs")))))
  # ---- per-gene-set background: observed NES of every RBP + null NES histogram of size-matched random sets ----
  pw_ids  <- names(sizes)                                            # K sets, same order as the random sets
  nbin    <- as.integer(round((GSEA_PW_HI - GSEA_PW_LO) / GSEA_PW_BIN))
  null_pw <- matrix(0L, K, nbin, dimnames = list(sub("^HALLMARK_", "", pw_ids), NULL))
  for (k in seq_len(K)) {
    v <- as.numeric(NES[(seq_len(B) - 1L) * K + k, , drop = FALSE])  # B replicates x all RBPs
    v <- v[is.finite(v)]
    null_pw[k, ] <- tabulate(pmin(pmax(floor((v - GSEA_PW_LO) / GSEA_PW_BIN) + 1L, 1L), nbin), nbins = nbin)
  }
  obs_pw <- Mobs[match(sub("^HALLMARK_", "", pw_ids), sub("^HALLMARK_", "", rownames(Mobs))), , drop = FALSE]
  rownames(obs_pw) <- sub("^HALLMARK_", "", pw_ids)
  gsea_sets <- list(ids = rownames(null_pw), size = as.integer(sizes), obs = round(obs_pw, 4),
                    null_counts = null_pw, null_n = as.numeric(rowSums(null_pw)),
                    lo = GSEA_PW_LO, hi = GSEA_PW_HI, bin = GSEA_PW_BIN)
  list(observed = observed, null_counts = as.integer(counts), null_n = sum(counts), B = B,
       n_rows = K, n_rbps = length(ranks), null_type = "gsea", gsea_sets = gsea_sets, built = Sys.time())
}

if (!exists("CORR_BUILD_NO_RUN")) {
  suppressMessages({ library(qs2); library(fgsea); library(msigdbr) })
  rd <- function(name) qs2::qs_read(file.path(DATA_DIR, paste0(name, ".qs2")))
  hm <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
  hallmark_sets <- lapply(split(hm$gene_symbol, hm$gs_name), function(g) unique(toupper(g)))
  jobs <- list(list(key = "gsea|Both Cells", gsea = "RBPs.gsea_All",   ranks = "RBPs.t_All"),
               list(key = "gsea|K562",       gsea = "RBPs.gsea_K562",  ranks = "RBPs.t_K562"),
               list(key = "gsea|HEPG2",      gsea = "RBPs.gsea_HEPG2", ranks = "RBPs.t_HEPG2"))
  if (!is.null(ONLY)) jobs <- Filter(function(j) j$key %in% ONLY, jobs)
  out_file <- file.path(DATA_DIR, "Correlation_Reference.qs2")
  ref <- if (file.exists(out_file)) qs2::qs_read(out_file) else list()
  for (j in jobs) {
    ref[[j$key]] <- build_gsea_entry(j$key, rd(j$gsea), rd(j$ranks), hallmark_sets)
    qs2::qs_save(ref, out_file)
    cat("saved ->", out_file, "\n")
  }
  cat("\nDone. Entries:", paste(names(ref), collapse = ", "), "\n")
}

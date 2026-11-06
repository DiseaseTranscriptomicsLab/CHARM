# =============================================================================
# build_correlation_reference.R  -- run ONCE (and again whenever the data change)
#   source("build_correlation_reference.R")      # from the CHARM folder, in R
# Writes data/QS2_Files/Correlation_Reference.qs2, which the app reads at startup
# (see correlation_reference.R for what it contains and how it is used).
#
# For each analysis x cell line it builds a (genes|events) x RBP matrix, then:
#   1. computes every RBP-vs-RBP Spearman correlation  -> "observed"
#   2. B times: shuffles the values WITHIN each row (gene/event) across the RBPs
#      that have a value for it, and recomputes all correlations -> "null"
# Shuffling only among observed cells keeps each pair's number of shared rows the
# same as in the real data, so null and observed are directly comparable.
#
# Settings (edit, or set before sourcing):  B = number of shuffles (every shuffle
# contributes ~12,000 null pairs, pooled; 20-100 are all fine); CORES = maximum
# parallel workers; RAM_BUDGET_GB = memory the workers may use together.
# Runtime: minutes per shuffle for the largest splicing set -- run it while you
# do something else. Set ONLY <- "expr|K562" (etc.) to build a single entry.
# =============================================================================

if (!exists("B"))      B      <- 100
if (!exists("SEED"))   SEED   <- 1906
if (!exists("ONLY"))   ONLY   <- NULL
if (!exists("CORES"))  CORES  <- if (.Platform$OS.type == "windows") 1 else max(1, min(4, parallel::detectCores() - 1))
# Every parallel worker holds ~3 copies of the matrix, so the number of workers is
# capped per entry so that workers x 3 x matrix size stays below this budget (GB).
# Lower it (e.g. 3) if you still run out of memory; raise it if you have plenty of RAM.
if (!exists("RAM_BUDGET_GB")) RAM_BUDGET_GB <- 6
if (!exists("DATA_DIR")) DATA_DIR <- "data/QS2_Files"
MIN_SHARED <- 3     # same minimum overlap the app requires for a correlation
source("correlation_reference.R", local = TRUE)   # CORR_BIN_WIDTH

# list of data.frames (one per RBP)  ->  rows x RBP matrix, NA where not measured
list_to_matrix <- function(lst, id_col, val_col) {
  vecs <- lapply(lst, function(df) {
    if (is.null(df) || nrow(df) == 0 || !all(c(id_col, val_col) %in% colnames(df))) return(NULL)
    df <- df[!duplicated(df[[id_col]]), , drop = FALSE]          # first occurrence, as in the app
    v  <- stats::setNames(as.numeric(df[[val_col]]), as.character(df[[id_col]]))
    v[is.finite(v)]
  })
  vecs <- vecs[vapply(vecs, function(v) !is.null(v) && length(v) >= MIN_SHARED, logical(1))]
  ids  <- unique(unlist(lapply(vecs, names), use.names = FALSE))
  M <- matrix(NA_real_, nrow = length(ids), ncol = length(vecs), dimnames = list(ids, names(vecs)))
  for (j in seq_along(vecs)) M[names(vecs[[j]]), j] <- vecs[[j]]
  M[rowSums(!is.na(M)) >= 2, , drop = FALSE]    # a row seen in <2 RBPs never forms a pair
}

pairwise_spearman <- function(M, n_shared) {
  cm <- suppressWarnings(stats::cor(M, method = "spearman", use = "pairwise.complete.obs"))
  cm[n_shared < MIN_SHARED] <- NA
  cm[upper.tri(cm)][is.finite(cm[upper.tri(cm)])]
}

# shuffle values within each row, among the RBPs that have a value for that row
shuffle_rows <- function(Mt) {            # Mt = t(M): RBPs x rows, so a row is a contiguous column
  for (i in seq_len(ncol(Mt))) {
    idx <- which(!is.na(Mt[, i]))
    if (length(idx) > 1L) Mt[idx, i] <- Mt[idx[sample.int(length(idx))], i]
  }
  Mt
}

to_counts <- function(r) {
  nb <- as.integer(round(2 / CORR_BIN_WIDTH))
  tabulate(pmin(pmax(floor((r + 1) / CORR_BIN_WIDTH) + 1, 1), nb), nbins = nb)
}

fmt_dur <- function(sec) {
  sec <- round(sec)
  if (sec < 90) sprintf("%ds", sec) else if (sec < 5400) sprintf("%.0f min", sec / 60) else sprintf("%.1f h", sec / 3600)
}

build_entry <- function(M, label) {
  t0 <- Sys.time()
  cat(sprintf("[%s] %s rows x %d RBPs\n", label, format(nrow(M), big.mark = ","), ncol(M)))
  n_shared <- crossprod(!is.na(M))
  observed <- sort(round(pairwise_spearman(M, n_shared), 4))
  cat(sprintf("[%s] observed: %s pairs (%.0fs)\n", label, format(length(observed), big.mark = ","),
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  Mt <- t(M)
  mat_gb  <- as.numeric(object.size(M)) / 1e9
  workers <- max(1L, min(CORES, floor(RAM_BUDGET_GB / (3 * mat_gb))))
  cat(sprintf("[%s] matrix %.2f GB -> %d parallel worker(s)\n", label, mat_gb, workers))
  one_shuffle <- function(b) {
    set.seed(SEED + b)
    P <- t(shuffle_rows(Mt))
    to_counts(pairwise_spearman(P, n_shared))
  }
  # run in rounds of `workers` shuffles and free memory between rounds
  counts <- 0L; done <- 0L; t_shuf <- Sys.time()
  for (idx in split(seq_len(B), ceiling(seq_len(B) / workers))) {
    res <- if (workers > 1) parallel::mclapply(idx, one_shuffle, mc.cores = workers, mc.preschedule = FALSE)
           else lapply(idx, one_shuffle)
    bad <- vapply(res, function(x) !is.numeric(x), logical(1))
    if (any(bad)) stop("a shuffle failed in ", label, ": try CORES <- 1 to see the error (often = out of memory)")
    counts <- counts + Reduce(`+`, res)
    done <- done + length(idx); rm(res); invisible(gc())
    el_total <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    el_shuf  <- as.numeric(difftime(Sys.time(), t_shuf, units = "secs"))
    eta      <- el_shuf / done * (B - done)
    cat(sprintf("[%s]   %d / %d shuffles | elapsed %s | est. remaining for this entry %s (done ~%s)\n",
                label, done, B, fmt_dur(el_total), fmt_dur(eta), format(Sys.time() + eta, "%H:%M")))
  }
  cat(sprintf("[%s] null: %d shuffles, %s pairs (%.0fs total)\n", label, B,
              format(sum(counts), big.mark = ","), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  list(observed = observed, null_counts = as.integer(counts), null_n = sum(counts), B = B,
       n_rows = nrow(M), n_rbps = ncol(M), built = Sys.time())
}

if (!exists("CORR_BUILD_NO_RUN")) {
  suppressMessages(library(qs2))
  rd <- function(name) qs2::qs_read(file.path(DATA_DIR, paste0(name, ".qs2")))
  jobs <- list(
    list(key = "expr|Both Cells",   file = "RBPs.t_All",        type = "expr"),
    list(key = "expr|K562",         file = "RBPs.t_K562",       type = "expr"),
    list(key = "expr|HEPG2",        file = "RBPs.t_HEPG2",      type = "expr"),
    list(key = "splice|Both Cells", file = "Charm.object",       type = "splice"),
    list(key = "splice|K562",       file = "Charm.object_K562",  type = "splice"),
    list(key = "splice|HEPG2",      file = "Charm.object_HEPG2", type = "splice")
  )
  if (!is.null(ONLY)) jobs <- Filter(function(j) j$key %in% ONLY, jobs)
  out_file <- file.path(DATA_DIR, "Correlation_Reference.qs2")
  ref <- if (!is.null(ONLY) && file.exists(out_file)) qs2::qs_read(out_file) else list()
  for (j in jobs) {
    obj <- rd(j$file)
    M <- if (j$type == "expr") list_to_matrix(obj, "Gene", "t")
         else list_to_matrix(lapply(obj, function(x) x$VulcanTable), "Event.ID", "dPSI")
    rm(obj); invisible(gc())
    ref[[j$key]] <- build_entry(M, j$key)
    qs2::qs_save(ref, out_file)                    # save after each entry: nothing is lost if interrupted
    cat("saved ->", out_file, "\n")
  }
  cat("\nDone. Entries:", paste(names(ref), collapse = ", "), "\n")
}

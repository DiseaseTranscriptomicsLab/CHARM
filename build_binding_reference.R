# =============================================================================
# build_binding_reference.R -- run ONCE (again when the binding data change)
#   source("build_binding_reference.R")        # from the CHARM folder, in R
# Writes data/QS2_Files/Binding_Correlation_Reference.qs2 (see correlation_reference.R).
#
# For each of the 12 similar_binding_* objects (ES/IR x Both/K562/HEPG2 x inc/dec) it
# computes the Pearson correlation between binding profiles (rows of $mean_profiles,
# over positions) and stores only histograms (bin 0.001):
#   counts_all  : every pair of profiles
#   counts_same : pairs sharing the same eCLIP Target (different KD RBP)
#   by_target   : for each eCLIP Target T, every pair that involves a T profile
#   self        : for each Target T, the self profile T__T against all other profiles
# Also builds ONE pooled entry per event type x dataset (key binding|<ES|IR>|<dataset>|pool): the inc and
# dec profiles in a single pool, every pair counted (including inc-vs-dec) -- used by the "together" plot.
# Entries built by an older version of this script are UPGRADED (only by_target/self
# are added, the all-pairs part is kept). Set REBUILD <- TRUE to redo everything, or
# ONLY <- c("binding|IR|K562|inc", ...) to build specific entries.
# Memory: blocks of CHUNK profiles (~ CHUNK x N doubles); lower CHUNK if memory is tight.
# =============================================================================
if (!exists("DATA_DIR")) DATA_DIR <- "data/QS2_Files"
if (!exists("CHUNK"))    CHUNK    <- 500
if (!exists("ONLY"))     ONLY     <- NULL
if (!exists("REBUILD"))  REBUILD  <- FALSE
source("correlation_reference.R", local = TRUE)    # CORR_BIN_WIDTH, bind_key()

bind_to_counts <- function(r) {
  nb <- as.integer(round(2 / CORR_BIN_WIDTH))
  as.numeric(tabulate(pmin(pmax(floor((pmin(pmax(r, -1), 1) + 1) / CORR_BIN_WIDTH) + 1, 1), nb), nbins = nb))
}

fmt_dur <- function(sec) {
  sec <- round(sec)
  if (sec < 90) sprintf("%ds", sec) else if (sec < 5400) sprintf("%.0f min", sec / 60) else sprintf("%.1f h", sec / 3600)
}

# de-duplicate (as the app does), drop constant profiles, standardise rows
bind_prepare <- function(M) {
  M <- M[!duplicated(rownames(M)), , drop = FALSE]
  has_na <- anyNA(M); Z <- NULL
  if (!has_na) {
    Zc <- M - rowMeans(M); s <- sqrt(rowSums(Zc^2))
    keep <- is.finite(s) & s > 0
    M <- M[keep, , drop = FALSE]; Z <- Zc[keep, , drop = FALSE] / s[keep]; rm(Zc)
  }
  list(M = M, Z = Z, has_na = has_na)
}

# every pair once: all pairs, and same-target pairs
bind_allpairs <- function(prep, label) {
  t0 <- Sys.time(); M <- prep$M; Z <- prep$Z; has_na <- prep$has_na
  if (has_na) cat(sprintf("[%s] NAs present -> slower pairwise-complete path\n", label))
  N   <- nrow(M)
  tgt <- sub("^[^_]*__", "", rownames(M))
  tg  <- match(tgt, unique(tgt))
  cat(sprintf("[%s] %s profiles x %d positions, %s targets\n", label, format(N, big.mark = ","), ncol(M),
              format(length(unique(tgt)), big.mark = ",")))
  nb <- as.integer(round(2 / CORR_BIN_WIDTH))
  cnt_all <- numeric(nb); cnt_same <- numeric(nb)
  starts <- seq(1L, N, by = CHUNK)
  for (k in seq_along(starts)) {
    I <- starts[k]:min(N, starts[k] + CHUNK - 1L)
    J <- starts[k]:N
    C <- if (!has_na) tcrossprod(Z[I, , drop = FALSE], Z[J, , drop = FALSE])
         else suppressWarnings(stats::cor(t(M[I, , drop = FALSE]), t(M[J, , drop = FALSE]), use = "pairwise.complete.obs"))
    up   <- outer(I, J, "<")
    ok   <- up & is.finite(C)
    cnt_all  <- cnt_all + bind_to_counts(C[ok])
    cnt_same <- cnt_same + bind_to_counts(C[ok & outer(tg[I], tg[J], "==")])
    rm(C, up, ok)
    if (k %% 10 == 0 || k == length(starts)) {
      el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      cat(sprintf("[%s]   all pairs: block %d / %d | elapsed %s\n", label, k, length(starts), fmt_dur(el)))
    }
  }
  list(counts_all = cnt_all, counts_same = cnt_same, n_all = sum(cnt_all), n_same = sum(cnt_same),
       n_profiles = N, n_targets = length(unique(tgt)), n_targets_multi = sum(tabulate(tg) >= 2))
}

# per target T: pairs involving a T profile (each pair once) and the self profile T__T vs all others
bind_target_hists <- function(prep, label) {
  t0 <- Sys.time(); M <- prep$M; Z <- prep$Z; has_na <- prep$has_na
  rn  <- rownames(M); tgt <- sub("^[^_]*__", "", rn)
  targets <- sort(unique(tgt)); nT <- length(targets)
  by_target <- vector("list", nT); names(by_target) <- targets
  n_by <- stats::setNames(numeric(nT), targets)
  self <- list(); n_self <- numeric(0)
  for (k in seq_len(nT)) {
    tn  <- targets[k]
    idx <- which(tgt == tn)
    C <- if (!has_na) tcrossprod(Z[idx, , drop = FALSE], Z)
         else suppressWarnings(stats::cor(t(M[idx, , drop = FALSE]), t(M), use = "pairwise.complete.obs"))
    Cin <- C[, idx, drop = FALSE]
    r   <- c(Cin[upper.tri(Cin)], if (length(idx) < ncol(C)) as.vector(C[, -idx, drop = FALSE]))
    r   <- r[is.finite(r)]
    by_target[[k]] <- as.integer(bind_to_counts(r)); n_by[k] <- length(r)
    s <- which(rn == paste0(tn, "__", tn))
    if (length(s) == 1L) {
      rs <- C[match(s, idx), -s]; rs <- rs[is.finite(rs)]
      self[[tn]] <- as.integer(bind_to_counts(rs)); n_self[tn] <- length(rs)
    }
    rm(C, Cin)
    if (k %% 20 == 0 || k == nT) {
      el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      cat(sprintf("[%s]   per-target: %d / %d targets | elapsed %s | est. remaining %s\n",
                  label, k, nT, fmt_dur(el), fmt_dur(el / k * (nT - k))))
    }
  }
  cat(sprintf("[%s] self profiles found for %d of %d targets\n", label, length(self), nT))
  list(by_target = by_target, n_by_target = n_by, self = self, n_self = n_self)
}

bind_build_entry <- function(M, label) {
  prep <- bind_prepare(M)
  e <- bind_allpairs(prep, label)
  c(e, bind_target_hists(prep, label), list(built = Sys.time()))
}

if (!exists("CORR_BUILD_NO_RUN")) {
  suppressMessages(library(qs2))
  out_file <- file.path(DATA_DIR, "Binding_Correlation_Reference.qs2")
  ref <- if (file.exists(out_file)) qs2::qs_read(out_file) else list()
  jobs <- expand.grid(prefix = c("ES", "IR"), ds = c("Both Cells", "K562", "HEPG2"), dir = c("inc", "dec"),
                      stringsAsFactors = FALSE)
  jobs$key  <- bind_key(jobs$prefix, jobs$ds, jobs$dir)
  jobs$file <- sprintf("similar_binding_%s_%s_%s", jobs$prefix, ifelse(jobs$ds == "Both Cells", "both", jobs$ds), jobs$dir)
  if (!is.null(ONLY)) jobs <- jobs[jobs$key %in% ONLY, ]
  need <- vapply(jobs$key, function(k) REBUILD || is.null(ref[[k]]) || is.null(ref[[k]]$by_target), logical(1))
  jobs <- jobs[need, ]
  cat(sprintf("%d entr%s to build or upgrade\n", nrow(jobs), if (nrow(jobs) == 1) "y" else "ies"))
  for (i in seq_len(nrow(jobs))) {
    key <- jobs$key[i]
    f   <- file.path(DATA_DIR, paste0(jobs$file[i], ".qs2"))
    if (!file.exists(f)) { cat("MISSING FILE:", f, "- skipped\n"); next }
    obj <- qs2::qs_read(f)
    M <- obj$mean_profiles; rm(obj); invisible(gc())
    if (REBUILD || is.null(ref[[key]])) {
      ref[[key]] <- bind_build_entry(M, key)
    } else {
      cat(sprintf("[%s] upgrading: adding per-target histograms\n", key))
      ref[[key]] <- utils::modifyList(ref[[key]], bind_target_hists(bind_prepare(M), key))
    }
    rm(M); invisible(gc())
    qs2::qs_save(ref, out_file)
    cat(sprintf("saved -> %s  (%d of %d)\n", out_file, i, nrow(jobs)))
  }
  # pooled inc + dec profiles (all pairs, including inc-vs-dec): key "binding|<ES|IR>|<dataset>|pool"
  pj <- expand.grid(prefix = c("ES", "IR"), ds = c("Both Cells", "K562", "HEPG2"), stringsAsFactors = FALSE)
  pj$key <- bind_pool_key(pj$prefix, pj$ds)
  if (!is.null(ONLY)) pj <- pj[pj$key %in% ONLY, ]
  pj <- pj[vapply(pj$key, function(k) REBUILD || is.null(ref[[k]]), logical(1)), ]
  cat(sprintf("%d pooled (inc + dec) entr%s to build\n", nrow(pj), if (nrow(pj) == 1) "y" else "ies"))
  for (i in seq_len(nrow(pj))) {
    fs <- sprintf("similar_binding_%s_%s_%s.qs2", pj$prefix[i], ifelse(pj$ds[i] == "Both Cells", "both", pj$ds[i]), c("inc", "dec"))
    fs <- file.path(DATA_DIR, fs)
    if (!all(file.exists(fs))) { cat("MISSING FILE for", pj$key[i], "- skipped\n"); next }
    Mi <- qs2::qs_read(fs[1])$mean_profiles; Md <- qs2::qs_read(fs[2])$mean_profiles
    P <- bind_pool_matrix(Mi, Md); rm(Mi, Md); invisible(gc())
    if (is.null(P)) { cat("inc/dec profiles have different positions for", pj$key[i], "- skipped\n"); next }
    e <- bind_allpairs(bind_prepare(P), pj$key[i]); rm(P); invisible(gc())
    ref[[pj$key[i]]] <- c(e, list(pooled = TRUE, built = Sys.time()))
    qs2::qs_save(ref, out_file)
    cat(sprintf("saved -> %s  (pooled %d of %d)\n", out_file, i, nrow(pj)))
  }
  cat("\nDone. Entries:", length(ref), "\n")
}

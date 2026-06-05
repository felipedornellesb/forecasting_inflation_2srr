# ==============================================================================
# 01_data_prep.R  --  inputs for the horse-race, in the professor's 312-window
# convention (ForecastingInflation/03_call_model.R), so the 2SRR lines up
# cell-by-cell with his benchmarks.
#
# Convention: for horizon h, window j (j=1..312) forecasts y[tau+j] from a rolling
# window, so the realized target is the same across h:
#   yout[j,h] = y[tau+j]      (inflation RATE, identical across h)
#   rw[j,h]   = y[tau+j-h]    (random walk: last rate at the origin)
# Target = RATE pi_{t+h} (advisor's instruction), not the accumulated series.
#
# Also imports the professor's 12 official benchmark forecasts (per-step RATE
# columns t+1..t+12) as the benchmark of record. 02 is an optional repro check.
# ==============================================================================
cat("== 01_data_prep.R ==\n\n")
source("00_prog/00_setup.R")

variable <- "CPIAUCSL"
nwindows <- 312
maxh     <- 12

# Professor's official run (benchmark of record)
DIR_PROF <- file.path(ROOT, "ForecastingInflation", "forecasts")
stopifnot("ForecastingInflation/forecasts not found." = dir.exists(DIR_PROF))

# ---- Data: must be byte-identical to the professor's (the 2SRR is estimated on
#      the same in-sample windows as his benchmarks, so the data must match). ----
load(file.path(DIR_DATA, "data.rda"))                      # project data.rda
prof_data_path <- file.path(ROOT, "ForecastingInflation", "data", "data.rda")
if (file.exists(prof_data_path)) {
  e <- new.env(); load(prof_data_path, envir = e)
  prof_data <- get(ls(e)[1], envir = e)
  same_cols <- identical(names(data), names(prof_data))
  num <- names(data)[vapply(names(data), function(c) is.numeric(data[[c]]), logical(1))]
  maxd <- suppressWarnings(max(vapply(intersect(num, names(prof_data)),
            function(c) max(abs(data[[c]] - prof_data[[c]]), na.rm = TRUE), numeric(1))))
  cat(sprintf("Data identity check vs professor: cols_equal=%s  max_abs_diff=%.2e\n",
              same_cols, maxd))
  if (!isTRUE(same_cols) || !is.finite(maxd) || maxd > 1e-10)
    stop("Project data.rda differs from the professor's data.rda — align them ",
         "before running (the 2SRR must see the same in-sample data as the benchmarks).")
} else {
  cat("  [note] professor's data.rda not found for the identity check (skipped).\n")
}

y   <- data[[variable]]
n   <- length(y)              # 786
tau <- n - nwindows           # 474  (= rolling window_size at h=1)
cat(sprintf("n=%d | nwindows=%d | tau=%d (window_size at h=1)\n", n, nwindows, tau))

# ---- Realized target (RATE) and Random Walk, professor's POOS convention ----
yout <- matrix(NA_real_, nwindows, maxh, dimnames = list(NULL, paste0("h", 1:maxh)))
rw   <- matrix(NA_real_, nwindows, maxh, dimnames = list(NULL, paste0("h", 1:maxh)))
for (h in 1:maxh) for (j in 1:nwindows) {
  yout[j, h] <- y[tau + j]                       # same realized rate for every h
  if (tau + j - h >= 1) rw[j, h] <- y[tau + j - h]
}
save(yout, file = file.path(DIR_FORECASTS, "yout.rda"))
save(rw,   file = file.path(DIR_FORECASTS, "rw.rda"))
cat(sprintf("Saved yout, rw: %d x %d\n", nwindows, maxh))
cat(sprintf("  yout columns identical across h: %s\n",
            all(apply(yout, 2, function(col) isTRUE(all.equal(col, yout[, 1]))))))
cat(sprintf("  yout[,1] == y[%d:%d]: %s\n", tau + 1, n,
            isTRUE(all.equal(as.numeric(yout[, 1]), y[(tau + 1):n]))))

# ---- Import the professor's 12 benchmark forecasts (per-step RATE columns) ----
benchmarks <- c("LASSO", "Ridge", "ElNET", "AdaLASSO", "AdaElNET", "RF", "Bagging",
                "Factor", "T.Factor", "CSR", "AR", "AR_BIC")
cat("\nImporting professor's benchmark forecasts (t+1..t+12 -> h1..h12):\n")
n_imported <- 0L
for (m in benchmarks) {
  src <- file.path(DIR_PROF, paste0(m, ".rda"))
  if (!file.exists(src)) { cat(sprintf("  %-10s MISSING in professor folder — skipped\n", m)); next }
  e <- new.env(); load(src, envir = e)
  M <- as.matrix(get(ls(e)[1], envir = e))
  if (ncol(M) < maxh) { cat(sprintf("  %-10s has only %d cols — skipped\n", m, ncol(M))); next }
  forecasts <- M[, 1:maxh, drop = FALSE]          # t+1..t+12 (the RATE forecasts)
  colnames(forecasts) <- paste0("h", 1:maxh)
  save(forecasts, file = file.path(DIR_FORECASTS, paste0(m, ".rda")))
  n_imported <- n_imported + 1L
  cat(sprintf("  imported %-10s -> %s.rda (%d x %d, rate)  RMSE/RW: %s\n",
              m, m, nrow(forecasts), maxh,
              paste(sprintf("h%d=%.3f", c(1, 3, 6, 12),
                sapply(c(1, 3, 6, 12), function(h) {
                  ok <- complete.cases(yout[, h], forecasts[, h], rw[, h])
                  sqrt(mean((yout[ok, h] - forecasts[ok, h])^2)) /
                  sqrt(mean((yout[ok, h] - rw[ok, h])^2))
                })), collapse = " ")))
}
cat(sprintf("\n%d/%d benchmark files imported into %s\n",
            n_imported, length(benchmarks), DIR_FORECASTS))
cat("== done ==\n")

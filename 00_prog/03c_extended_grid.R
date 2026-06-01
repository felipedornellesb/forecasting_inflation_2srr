# =============================================================================
# 03c_extended_grid.R — robustness run for Section 4.4 (the grid argument).
#
# Re-estimates ONLY the three 2SRR specifications on a deliberately OVER-EXTENDED
# penalty grid, exp(linspace(-4, 18, 25)), to compare against the bounded grid
# exp(linspace(-2, 12, 15)) used in the main run. It produces the relative-RMSE
# table that 05_article_figures.R reads for FIG6 (original vs extended grid).
#
# SAFE BY DESIGN:
#   - forecasts are written to 30_output/forecasts_extgrid/  (a SEPARATE folder;
#     the main 30_output/forecasts/ is never touched);
#   - the P1 table is written to a NEW run folder
#     40_results/run_extended_grid_overexpanded_linspace_-4_18_25/tables/.
#   - runs SERIALLY with the callr-isolated GARCH (crash-safe), no PSOCK cluster.
#
# ~4 h on a typical laptop. Run in the background.
# =============================================================================
cat("== 03c_extended_grid.R (extended-grid robustness run) ==\n\n")
source("00_prog/00_setup.R")
options(garch_isolate = TRUE)            # serial context: callr isolation is safe
set.seed(12345)

# ---- Extended grid (the only change vs the main run) ------------------------
EXT_GRID <- exp(pracma::linspace(-4, 18, n = 25))
lambda2  <- 0.1
variable <- "CPIAUCSL"; horizons <- c(1, 3, 6, 12); maxh <- 12
ly <- 2; lf <- 2; nf <- 4; kfold <- 5
cases_tvp <- c("AR", "Factor", "FAVAR")

# ---- Separate output locations (never overwrite the main run) ---------------
FC_EXT  <- file.path("30_output", "forecasts_extgrid")
RUN_EXT <- file.path("40_results", "run_extended_grid_overexpanded_linspace_-4_18_25")
dir.create(FC_EXT,  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(RUN_EXT, "tables"), recursive = TRUE, showWarnings = FALSE)

load(file.path(DIR_DATA, "data.rda"))
load(file.path(DIR_FORECASTS, "yout.rda"))
load(file.path(DIR_FORECASTS, "rw.rda"))
data$date <- NULL
bigt <- nrow(data); n_oos <- nrow(yout); tau <- bigt - n_oos
y_col <- which(names(data) == variable)

prepare_X_is <- function(X_is) {
  cv   <- apply(X_is, 2, var, na.rm = TRUE)
  good <- which(is.finite(cv) & cv > 1e-10)
  X_is <- X_is[, good, drop = FALSE]; X_is[is.na(X_is)] <- 0
  X_means <- colMeans(X_is); X_sds <- apply(X_is, 2, sd); X_sds[X_sds < 1e-10] <- 1
  scale(X_is, center = X_means, scale = X_sds)
}
rmse_fn <- function(y, f) { ok <- complete.cases(y, f)
  if (sum(ok) < 5) NA_real_ else sqrt(mean((y[ok] - f[ok])^2)) }

fc_ext <- setNames(lapply(cases_tvp, function(.) matrix(NA_real_, n_oos, maxh,
                          dimnames = list(NULL, paste0("h", 1:maxh)))), cases_tvp)

prog <- file.path(DIR_CHECKPOINTS, "run_extgrid.progress")
dir.create(dirname(prog), recursive = TRUE, showWarnings = FALSE)
cat(sprintf("ext-grid run start %s\n", Sys.time()), file = prog)
t0 <- Sys.time()

for (wi in 1:n_oos) {
  t_end    <- tau + wi - 1
  X_is_std <- prepare_X_is(as.matrix(data[1:t_end, -y_col, drop = FALSE]))
  y_is     <- as.numeric(data[1:t_end, y_col])
  for (case in cases_tvp) for (h in horizons) {
    res <- fit_2srr_window(y_is = y_is, X_is_raw = X_is_std, h = h, case = case,
                           ly = ly, lf = lf, nf = nf, kfold = kfold,
                           lambdas = EXT_GRID, lambda2 = lambda2,
                           coulombe_lambdavec = EXT_GRID, engine = "coulombe_fast")
    fc_ext[[case]][wi, h] <- res$forecast
  }
  if (wi %% 10 == 1 || wi == n_oos) {
    el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    msg <- sprintf("  win %3d/%d  %.1f min  ETA %.1f min\n", wi, n_oos, el,
                   el / wi * (n_oos - wi))
    cat(msg); cat(msg, file = prog, append = TRUE)
  }
  if (wi %% 30 == 0) for (case in cases_tvp) {
    forecasts <- fc_ext[[case]]
    save(forecasts, file = file.path(FC_EXT, paste0("2SRR_", case, ".rda")))
  }
}

# ---- Save forecasts + write the P1 (RMSE relative to RW) table for FIG6 -----
for (case in cases_tvp) {
  forecasts <- fc_ext[[case]]
  save(forecasts, file = file.path(FC_EXT, paste0("2SRR_", case, ".rda")))
}
p1 <- data.frame(model = paste0("2SRR_", cases_tvp))
for (h in 1:maxh) {
  col <- sapply(cases_tvp, function(case) {
    if (!(h %in% horizons)) return(NA_real_)
    rr <- rmse_fn(yout[, h], rw[, h])
    if (is.na(rr) || rr == 0) NA_real_ else round(rmse_fn(yout[, h], fc_ext[[case]][, h]) / rr, 4)
  })
  p1[[paste0("h", h)]] <- as.numeric(col)
}
write.csv(p1, file.path(RUN_EXT, "tables", "P1_rmsfe_relative_rw_all_h.csv"),
          row.names = FALSE)

cat(sprintf("\n  Done in %.1f min.\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cat("  Wrote:", file.path(RUN_EXT, "tables", "P1_rmsfe_relative_rw_all_h.csv"), "\n")
cat("  RMSE/RW (extended grid):\n"); print(p1[, c("model", "h1", "h3", "h6", "h12")])
cat("  -> now re-run 05_article_figures.R to generate FIG6.\n")
cat("== 03c done ==\n")

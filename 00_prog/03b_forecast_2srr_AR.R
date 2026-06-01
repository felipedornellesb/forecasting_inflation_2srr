# =============================================================================
# 03b_forecast_2srr_AR.R — runs ONLY the 2SRR-AR case.
#
# Why this exists: rugarch's solver intermittently crashes R at the C level on
# Windows. In the full parallel 03 that crash silently killed the worker that
# happened to be handling the AR case, leaving 2SRR_AR.rda unsaved (Factor and
# FAVAR completed). This script regenerates ONLY the AR case, serially, with the
# GARCH(1,1) step isolated in a callr subprocess (see tvp_functions.R): a
# subprocess crash is caught in the parent and that single window falls back to
# a rolling variance, so the run always finishes WITH GARCH on (almost) every
# window. Factor/FAVAR outputs are NOT touched.
#
# Output: 2SRR_AR.rda, Ridge_from_2SRR_AR.rda, betas_2SRR_AR.rda
# =============================================================================
cat("== 03b_forecast_2srr_AR.R ==\n\n")
source("00_prog/00_setup.R")

options(garch_isolate = TRUE,        # isolate rugarch in a callr subprocess
        garch_solver  = "hybrid")    # same solver as the Factor/FAVAR runs

# ---- Config (mirrors 03_forecast_2srr.R) ------------------------------------
variable   <- "CPIAUCSL"
horizons   <- c(1, 3, 6, 12)
maxh       <- 12
ly <- 2; lf <- 2; nf <- 4; kfold <- 5
lambda_vec <- exp(pracma::linspace(-2, 12, n = 15))
lambda2    <- 0.1
ENGINE     <- "coulombe_fast"

load(file.path(DIR_DATA, "data.rda"))
load(file.path(DIR_FORECASTS, "yout.rda"))
load(file.path(DIR_FORECASTS, "rw.rda"))
data$date <- NULL
bigt  <- nrow(data); n_oos <- nrow(yout); tau <- bigt - n_oos
y_col <- which(names(data) == variable)
cat(sprintf("  Total obs: %d | OOS windows: %d | tau: %d\n", bigt, n_oos, tau))

prepare_X_is <- function(X_is) {
  cv   <- apply(X_is, 2, var, na.rm = TRUE)
  good <- which(is.finite(cv) & cv > 1e-10)
  X_is <- X_is[, good, drop = FALSE]; X_is[is.na(X_is)] <- 0
  X_means <- colMeans(X_is); X_sds <- apply(X_is, 2, sd); X_sds[X_sds < 1e-10] <- 1
  scale(X_is, center = X_means, scale = X_sds)
}

# ---- Containers --------------------------------------------------------------
fc_AR    <- matrix(NA_real_, n_oos, maxh, dimnames = list(NULL, paste0("h", 1:maxh)))
fc_ridge <- matrix(NA_real_, n_oos, maxh, dimnames = list(NULL, paste0("h", 1:maxh)))
betas_bundle <- setNames(lapply(horizons, function(.) vector("list", n_oos)),
                         paste0("h", horizons))

prog <- file.path(DIR_CHECKPOINTS, "run_AR.progress")
dir.create(dirname(prog), recursive = TRUE, showWarnings = FALSE)
cat(sprintf("AR run start %s\n", Sys.time()), file = prog)

# ---- Serial POOS loop (AR only) ---------------------------------------------
t_total <- Sys.time(); n_fallback <- 0L
for (wi in 1:n_oos) {
  t_end    <- tau + wi - 1
  X_is_std <- prepare_X_is(as.matrix(data[1:t_end, -y_col, drop = FALSE]))
  y_is     <- as.numeric(data[1:t_end, y_col])

  for (h in horizons) {
    res <- fit_2srr_window(y_is = y_is, X_is_raw = X_is_std, h = h, case = "AR",
                           ly = ly, lf = lf, nf = nf, kfold = kfold,
                           lambdas = lambda_vec, lambda2 = lambda2,
                           coulombe_lambdavec = lambda_vec, engine = ENGINE)
    fc_AR[wi, h]    <- res$forecast
    fc_ridge[wi, h] <- res$ridge_forecast
    if (!is.null(res$status) && res$status == "fit_failed_fallback_ridge")
      n_fallback <- n_fallback + 1L
    betas_bundle[[paste0("h", h)]][[wi]] <- list(
      betas_tvp = res$betas_tvp, lambda = res$lambda,
      lambda_step1 = res$lambda_step1, omega = res$omega, sigma2 = res$sigma2,
      var_names = res$var_names, n_obs = res$n_obs, n_vars = res$n_vars,
      status = res$status)
  }

  if (wi %% 10 == 1 || wi == n_oos) {
    el <- as.numeric(difftime(Sys.time(), t_total, units = "mins"))
    eta <- el / wi * (n_oos - wi)
    msg <- sprintf("  win %3d/%d  %.1f min elapsed  ETA %.1f min\n", wi, n_oos, el, eta)
    cat(msg); cat(msg, file = prog, append = TRUE)
  }
  if (wi %% 30 == 0) {           # checkpoint (safety net)
    forecasts <- fc_AR
    save(forecasts, file = file.path(DIR_FORECASTS, "2SRR_AR.rda"))
  }
}

# ---- Save outputs ------------------------------------------------------------
forecasts <- fc_AR
save(forecasts, file = file.path(DIR_FORECASTS, "2SRR_AR.rda"))
forecasts <- fc_ridge
save(forecasts, file = file.path(DIR_FORECASTS, "Ridge_from_2SRR_AR.rda"))
save(betas_bundle, file = file.path(DIR_BETAS, "betas_2SRR_AR.rda"))

garch_session_close()

cat(sprintf("\n  Done in %.1f min. GARCH->rolling fallbacks: %d of %d fits.\n",
            as.numeric(difftime(Sys.time(), t_total, units = "mins")),
            n_fallback, n_oos * length(horizons)))
cat("  Saved: 2SRR_AR.rda, Ridge_from_2SRR_AR.rda, betas_2SRR_AR.rda\n")
cat(sprintf("  RMSE(2SRR-AR)/RMSE(RW): %s\n",
            paste(sprintf("h%d=%.3f", horizons,
              sapply(horizons, function(h) {
                ok <- complete.cases(yout[, h], fc_AR[, h], rw[, h])
                sqrt(mean((yout[ok, h] - fc_AR[ok, h])^2)) /
                sqrt(mean((yout[ok, h] - rw[ok, h])^2))
              })), collapse = "  ")))
cat("== 03b done ==\n")

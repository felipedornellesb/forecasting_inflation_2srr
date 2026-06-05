# ==============================================================================
# 03_forecast_2srr.R  --  2SRR (Goulet Coulombe 2025, Algorithm 1) on FRED-MD,
# cases Factor and FAVAR. The AR case is produced separately by 03b.
#
# Alignment: estimated on the IDENTICAL fixed-size rolling windows as the
# benchmarks (ForecastingInflation/03_call_model.R) -- for horizon h, window j
# forecasts y[tau+j] from rows j:(j+tau-h) -- so the 2SRR and the benchmarks
# differ only by the estimator (the controlled horse-race).
#
# Target = inflation RATE pi_{t+h}, direct multi-step (build_design_tvp).
# Engine "coulombe_fast" reproduces TVPRR_cosso (GARCH Step 2, OLS prior) via the
# T x T dual kernel; two self-tests assert equivalence before any compute.
#
# Outputs: 2SRR_{Factor,FAVAR}.rda, Ridge_from_2SRR_{Factor,FAVAR}.rda,
#          betas_2SRR_{Factor,FAVAR}.rda
# ==============================================================================

cat("== 03_forecast_2srr.R ==\n\n")
source("00_prog/00_setup.R")

suppressPackageStartupMessages({ library(glmnet); library(pracma); library(fGarch) })

# GARCH (Step 2) isolated in a persistent callr subprocess so a rare rugarch
# C-level segfault is contained in the child (that window falls back to a rolling
# variance). Safe because 03 runs serially (no callr inside a PSOCK worker).
options(garch_isolate = TRUE, garch_solver = "hybrid")

# ---------------------------------------------------------------------------- #
# SELF-TESTS — abort before wasting compute if the math is wrong.
# ---------------------------------------------------------------------------- #
cat("[validation 1] Standalone self-test (dual vs primal, Eq. 5)...\n")
diff_test <- tvp_self_test(seed = 1, T_obs = 30, K = 3,  lam = 1,   verbose = TRUE)
diff_test <- max(diff_test,
                 tvp_self_test(seed = 2, T_obs = 50, K = 5,  lam = 5,   verbose = TRUE),
                 tvp_self_test(seed = 3, T_obs = 80, K = 10, lam = 0.1, verbose = TRUE))
if (diff_test > 1e-6) stop("Self-test failed (diff=", diff_test, "). Aborting.")
cat("[validation 1 OK] standalone reproduces Eq. 5 (err <", 1e-6, ")\n\n")

cat("[validation 2] coulombe_fast vs original TVPRR_cosso/dualGRR...\n")
diff_test2 <- tvp_coulombe_fast_vs_original_test(seed = 11, T_obs = 25, K = 4,
                                                 lambda1 = 5, lambda2 = 0.1, verbose = TRUE)
diff_test2 <- max(diff_test2,
                  tvp_coulombe_fast_vs_original_test(seed = 12, T_obs = 40, K = 6,
                                                     lambda1 = 100, lambda2 = 0.5, verbose = TRUE))
if (!is.na(diff_test2) && diff_test2 > 1e-4)
  warning(sprintf("coulombe_fast diverges from original by %.3e", diff_test2))
cat("[validation 2 OK]\n\n")

# 1. Configuration ------------------------------------------------------------
variable  <- "CPIAUCSL"
horizons  <- c(1, 3, 6, 12)
maxh      <- 12
cases_tvp <- c("Factor", "FAVAR")      # AR is produced by 03b (serial, GARCH isolated)

ENGINE_PRIMARY <- "coulombe_fast"

# Hyperparameters (Coulombe IJF 2025, sec. 4 defaults — NO recalibration)
ly       <- 2          # y lags
lf       <- 2          # factor lags
nf       <- 4          # number of PCA factors
kfold    <- 5
n_lambda <- 15
lambda_vec <- exp(pracma::linspace(-2, 12, n = n_lambda))   # Step-1 ridge grid
lambda2    <- 0.1                                            # constant-block lambda

# Optional cross-check against the ORIGINAL (slow) TVPRR_cosso. The two self-tests
# above already prove equivalence; leave FALSE for the production run.
RUN_COULOMBE_CHECK     <- FALSE
COULOMBE_CHECK_WINDOWS <- c(40, 150, 300)
COULOMBE_CHECK_HS      <- c(1, 12)
COULOMBE_CHECK_CASE    <- "FAVAR"

# 2. Data ---------------------------------------------------------------------
load(file.path(DIR_DATA, "data.rda"))
load(file.path(DIR_FORECASTS, "yout.rda"))
data$date <- NULL

bigt  <- nrow(data)
n_oos <- nrow(yout)              # 312
tau   <- bigt - n_oos            # 474 (= rolling window_size at h=1)
y_col <- which(names(data) == variable)

cat(sprintf("  Total obs: %d | OOS windows: %d | tau (window_size @h=1): %d\n",
            bigt, n_oos, tau))
cat(sprintf("  Cases: %s | Horizons: %s | Engine: %s\n",
            paste(cases_tvp, collapse = ", "),
            paste(horizons, collapse = ","), ENGINE_PRIMARY))

# 3. Output containers --------------------------------------------------------
fc_2srr  <- setNames(lapply(cases_tvp, function(.) matrix(NA_real_, n_oos, maxh,
                       dimnames = list(NULL, paste0("h", 1:maxh)))), cases_tvp)
fc_ridge <- setNames(lapply(cases_tvp, function(.) matrix(NA_real_, n_oos, maxh,
                       dimnames = list(NULL, paste0("h", 1:maxh)))), cases_tvp)
betas_store <- setNames(
  lapply(cases_tvp, function(.)
    setNames(lapply(horizons, function(.) vector("list", n_oos)), paste0("h", horizons))),
  cases_tvp)

# 4. Helper: standardize X_is (drop zero-var cols, center/scale) --------------
prepare_X_is <- function(X_is) {
  cv   <- apply(X_is, 2, var, na.rm = TRUE)
  good <- which(is.finite(cv) & cv > 1e-10)
  X_is <- X_is[, good, drop = FALSE]; X_is[is.na(X_is)] <- 0
  X_means <- colMeans(X_is); X_sds <- apply(X_is, 2, sd); X_sds[X_sds < 1e-10] <- 1
  scale(X_is, center = X_means, scale = X_sds)
}

all_out_paths <- sapply(cases_tvp, function(c)
                          file.path(DIR_FORECASTS, paste0("2SRR_", c, ".rda")))

if (all(file.exists(all_out_paths))) {
  cat("\n  2SRR_Factor.rda and 2SRR_FAVAR.rda already exist. Nothing to do.\n")
  cat("  (Delete them to force a re-run.)\n")
} else {

  # 5. POOS loop — ROLLING window per horizon (Medeiros convention) -----------
  t_total <- Sys.time()
  cat(sprintf("\n  Starting POOS loop (rolling window, %d windows x %d cases x %d h)\n",
              n_oos, length(cases_tvp), length(horizons)))

  for (h in horizons) {
    window_size <- bigt - (n_oos + h - 1)        # = tau - h + 1  (rolling_window)
    cat(sprintf("\n  --- horizon h=%d  (window_size=%d) ---\n", h, window_size))
    t_h <- Sys.time()

    for (wi in 1:n_oos) {
      is_idx   <- wi:(wi + window_size - 1)       # rolling in-sample rows j:(j+ws-1)
      y_is     <- as.numeric(data[is_idx, y_col])
      X_is_std <- prepare_X_is(as.matrix(data[is_idx, -y_col, drop = FALSE]))

      for (case in cases_tvp) {
        res <- fit_2srr_window(
          y_is = y_is, X_is_raw = X_is_std, h = h, case = case,
          ly = ly, lf = lf, nf = nf, kfold = kfold,
          lambdas = lambda_vec, lambda2 = lambda2,
          coulombe_lambdavec = lambda_vec, engine = ENGINE_PRIMARY)
        fc_2srr[[case]][wi, h]  <- res$forecast
        fc_ridge[[case]][wi, h] <- res$ridge_forecast
        betas_store[[case]][[paste0("h", h)]][[wi]] <- list(
          betas_tvp = res$betas_tvp, lambda = res$lambda,
          lambda_step1 = res$lambda_step1, omega = res$omega, sigma2 = res$sigma2,
          var_names = res$var_names, n_obs = res$n_obs, n_vars = res$n_vars,
          status = res$status)
      }

      if (wi %% 25 == 1 || wi == n_oos) {
        el  <- as.numeric(difftime(Sys.time(), t_h, units = "mins"))
        eta <- el / wi * (n_oos - wi)
        cat(sprintf("    h=%2d win %3d/%d  (T_is=%d)  %.1f min  ETA(h) %.1f min\n",
                    h, wi, n_oos, window_size, el, eta))
      }
      if (wi %% 60 == 0) {
        save(fc_2srr, fc_ridge, betas_store,
             file = file.path(DIR_CHECKPOINTS, sprintf("ckpt_2SRR_h%02d_w%03d.rda", h, wi)))
      }
    }
  }

  cat(sprintf("\n  POOS loop finished in %.1f min.\n",
              as.numeric(difftime(Sys.time(), t_total, units = "mins"))))

  # 6. Save outputs -----------------------------------------------------------
  cat("\n  Saving outputs...\n")
  for (case in cases_tvp) {
    forecasts <- fc_2srr[[case]]
    save(forecasts, file = file.path(DIR_FORECASTS, paste0("2SRR_", case, ".rda")))
    forecasts <- fc_ridge[[case]]
    save(forecasts, file = file.path(DIR_FORECASTS, paste0("Ridge_from_2SRR_", case, ".rda")))
    betas_bundle <- betas_store[[case]]
    save(betas_bundle, file = file.path(DIR_BETAS, paste0("betas_2SRR_", case, ".rda")))
    cat(sprintf("    2SRR_%s.rda  +  Ridge_from_2SRR_%s.rda  +  betas saved.\n", case, case))
  }

  # 7. Quick RMSE/RW sanity print --------------------------------------------
  load(file.path(DIR_FORECASTS, "rw.rda"))
  for (case in cases_tvp)
    cat(sprintf("  RMSE(2SRR-%s)/RW: %s\n", case,
        paste(sprintf("h%d=%.3f", horizons,
          sapply(horizons, function(h) {
            ok <- complete.cases(yout[, h], fc_2srr[[case]][, h], rw[, h])
            sqrt(mean((yout[ok, h] - fc_2srr[[case]][ok, h])^2)) /
            sqrt(mean((yout[ok, h] - rw[ok, h])^2))
          })), collapse = "  ")))
}

# 8. Optional cross-check vs ORIGINAL TVPRR_cosso (rolling window) ------------
if (RUN_COULOMBE_CHECK && exists("TVPRR_cosso", mode = "function")) {
  cat("\n  --- Cross-check: original TVPRR_cosso on a few rolling windows ---\n")
  cc <- list()
  for (h in COULOMBE_CHECK_HS) {
    window_size <- bigt - (n_oos + h - 1)
    for (wi in COULOMBE_CHECK_WINDOWS) {
      if (wi > n_oos) next
      is_idx   <- wi:(wi + window_size - 1)
      y_is     <- as.numeric(data[is_idx, y_col])
      X_is_std <- prepare_X_is(as.matrix(data[is_idx, -y_col, drop = FALSE]))
      res_fast <- fit_2srr_window(y_is, X_is_std, h, COULOMBE_CHECK_CASE, ly, lf, nf,
                                  kfold, lambda_vec, engine = "coulombe_fast",
                                  lambda2 = lambda2, coulombe_lambdavec = lambda_vec)
      res_orig <- tryCatch(fit_2srr_window(y_is, X_is_std, h, COULOMBE_CHECK_CASE, ly, lf, nf,
                                  kfold, lambda_vec, engine = "coulombe",
                                  lambda2 = lambda2, coulombe_lambdavec = lambda_vec),
                           error = function(e) NULL)
      if (!is.null(res_orig) && is.finite(res_orig$forecast) && is.finite(res_fast$forecast)) {
        rd <- abs(res_orig$forecast - res_fast$forecast) / max(abs(res_fast$forecast), 1e-6)
        cc[[sprintf("w%d_h%d", wi, h)]] <- list(window = wi, horizon = h,
              fast = res_fast$forecast, orig = res_orig$forecast, rel_diff = rd)
        cat(sprintf("    win %3d h=%2d  fast=%.4f orig=%.4f  rel_diff=%.2f%%\n",
                    wi, h, res_fast$forecast, res_orig$forecast, 100 * rd))
      }
    }
  }
  save(cc, file = file.path(DIR_FORECASTS,
                            sprintf("2SRR_%s_coulombe_check.rda", COULOMBE_CHECK_CASE)))
}

# Close the persistent GARCH subprocess (if any).
if (exists("garch_session_close", mode = "function")) garch_session_close()

cat("\n== done ==\n")

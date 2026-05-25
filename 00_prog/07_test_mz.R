# Standalone test of the Mincer-Zarnowitz joint p-value (Table 1 / FIG8 / TABLE2)
# from 04_analysis.R PART 14b, run directly against 30_output/forecasts.
#
# Purpose: the published table showed p_joint = 0.000 for EVERY model at h=12.
# Direct h-step forecasts overlap by h-1 months, so the MZ residuals follow an
# MA(h-1) process. The current code uses plain OLS (homoskedastic) standard
# errors in the joint Wald (alpha=0, beta=1), which ignore that serial
# correlation and over-reject at long horizons. This script recomputes the
# joint p-value with a Newey-West HAC covariance (lag = h-1, the overlap length;
# at h=1 it reduces to a heteroskedasticity-robust estimator) and compares.
suppressPackageStartupMessages({ library(sandwich) })

DIR_FC   <- "30_output/forecasts"
OUT_DIR  <- "40_results/_test_mz"
RUN_FINAL <- "40_results/run_coulombe_grid_linspace(-2_12_15)"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

ld_fc <- function(f) { e <- new.env(); load(file.path(DIR_FC, f), envir = e)
                       as.matrix(get(ls(e)[1], envir = e)) }
load(file.path(DIR_FC, "yout.rda"))
horizons <- c(1, 3, 6, 12)

# --- Build the SAME model universe 04_analysis.R PART 14b uses ----------------
medeiros_models <- c("Ridge", "LASSO", "ElNET", "AdaLASSO", "AdaElNET",
                     "RF", "Bagging", "Factor", "T.Factor", "CSR", "AR", "AR_BIC")
cases_tvp <- c("AR", "Factor", "FAVAR")
all_fc <- list()
for (m in medeiros_models) {
  fp <- file.path(DIR_FC, paste0(m, ".rda"))
  if (file.exists(fp)) all_fc[[m]] <- ld_fc(paste0(m, ".rda"))
}
for (case in cases_tvp) {
  fp <- file.path(DIR_FC, paste0("2SRR_", case, ".rda"))
  if (file.exists(fp)) all_fc[[paste0("2SRR_", case)]] <- ld_fc(paste0("2SRR_", case, ".rda"))
  fp <- file.path(DIR_FC, paste0("Ridge_from_2SRR_", case, ".rda"))
  if (file.exists(fp)) all_fc[[paste0("RidgeStep1_", case)]] <- ld_fc(paste0("Ridge_from_2SRR_", case, ".rda"))
}
cat(sprintf("Loaded %d forecast matrices.\n\n", length(all_fc)))

# --- OLD joint p-value: OLS (homoskedastic) vcov ------------------------------
mz_p_old <- function(y, f) {
  ok <- complete.cases(y, f); if (sum(ok) < 30) return(NA_real_)
  y_ <- y[ok]; f_ <- f[ok]
  m <- tryCatch(lm(y_ ~ f_), error = function(e) NULL); if (is.null(m)) return(NA_real_)
  cf <- coef(m); vcv <- vcov(m); d <- cf - c(0, 1)
  w <- tryCatch(as.numeric(t(d) %*% solve(vcv) %*% d), error = function(e) NA_real_)
  if (is.na(w) || w < 0) NA_real_ else 1 - pchisq(w, df = 2)
}

# --- NEW joint p-value: Newey-West HAC, lag = h-1 (overlap length) ------------
# This is exactly the version that will replace mz_test() in 04_analysis.R.
mz_test_new <- function(y, f, h = 1) {
  ok <- complete.cases(y, f)
  if (sum(ok) < 30) return(list(alpha = NA, beta = NA, R2 = NA, p_joint = NA, n = sum(ok)))
  y_ <- y[ok]; f_ <- f[ok]
  m <- tryCatch(lm(y_ ~ f_), error = function(e) NULL)
  if (is.null(m)) return(list(alpha = NA, beta = NA, R2 = NA, p_joint = NA, n = length(y_)))
  cf <- coef(m)
  vcv <- tryCatch(
    sandwich::NeweyWest(m, lag = max(0, h - 1), prewhite = FALSE, adjust = TRUE),
    error = function(e) vcov(m))
  R <- rbind(c(1, 0), c(0, 1)); r <- c(0, 1)
  w <- tryCatch({
    d <- R %*% cf - r
    as.numeric(t(d) %*% solve(R %*% vcv %*% t(R)) %*% d)
  }, error = function(e) NA_real_)
  list(alpha = cf[1], beta = cf[2], R2 = summary(m)$r.squared,
       p_joint = if (is.na(w) || w < 0) NA_real_ else 1 - pchisq(w, df = 2),
       n = length(y_))
}

# --- Compare across every model x horizon -------------------------------------
rows <- list()
for (mn in names(all_fc)) {
  M <- all_fc[[mn]]
  for (h in horizons) {
    if (h > ncol(M)) next
    y_h <- yout[, h]; f_h <- M[, h]
    new <- mz_test_new(y_h, f_h, h)
    rows[[length(rows) + 1]] <- data.frame(
      model = mn, h = h,
      p_OLS = round(mz_p_old(y_h, f_h), 4),
      p_HAC = round(new$p_joint, 4))
  }
}
cmp <- do.call(rbind, rows)

to_wide <- function(col) {
  w <- reshape(cmp[, c("model", "h", col)], idvar = "model", timevar = "h",
               direction = "wide")
  names(w) <- c("model", paste0("h", horizons)); w
}
cat("=== OLD (OLS standard errors) MZ joint p-value ===\n")
print(to_wide("p_OLS"), row.names = FALSE)
cat("\n=== NEW (Newey-West HAC, lag = h-1) MZ joint p-value ===\n")
print(to_wide("p_HAC"), row.names = FALSE)

cat("\n--- h=12 column: how many models reject efficiency (p<0.10)? ---\n")
h12 <- cmp[cmp$h == 12, ]
cat(sprintf("  OLS: %d of %d models with p<0.10 (i.e. %d 'pass')\n",
            sum(h12$p_OLS < 0.10, na.rm = TRUE), nrow(h12),
            sum(h12$p_OLS >= 0.10, na.rm = TRUE)))
cat(sprintf("  HAC: %d of %d models with p<0.10 (i.e. %d 'pass')\n",
            sum(h12$p_HAC < 0.10, na.rm = TRUE), nrow(h12),
            sum(h12$p_HAC >= 0.10, na.rm = TRUE)))

# --- Patch RUN_FINAL P14b in place (p_joint only) into the test folder --------
# alpha/beta/R2/PT are point estimates unaffected by the vcov choice; only the
# joint p-value changes. We rewrite ONLY MZ_p_joint, preserving every other
# column verbatim, producing a drop-in replacement for the published CSV.
p14b_path <- file.path(RUN_FINAL, "tables", "P14b_MZ_PT_tests.csv")
if (file.exists(p14b_path)) {
  tab <- read.csv(p14b_path, stringsAsFactors = FALSE)
  old_pj <- tab$MZ_p_joint
  for (i in seq_len(nrow(tab))) {
    mn <- tab$model[i]; h <- tab$h[i]
    if (!is.null(all_fc[[mn]]) && h <= ncol(all_fc[[mn]])) {
      pj <- mz_test_new(yout[, h], all_fc[[mn]][, h], h)$p_joint
      tab$MZ_p_joint[i] <- round(pj, 4)
    }
  }
  changed <- which(abs(ifelse(is.na(old_pj), -1, old_pj) -
                       ifelse(is.na(tab$MZ_p_joint), -1, tab$MZ_p_joint)) > 1e-9)
  cat(sprintf("\nPatched %d of %d MZ_p_joint rows.\n", length(changed), nrow(tab)))
  out_csv <- file.path(OUT_DIR, "P14b_MZ_PT_tests_CORRECTED.csv")
  write.csv(tab, out_csv, row.names = FALSE)
  cat("Wrote corrected drop-in CSV ->", out_csv, "\n")
} else {
  cat("\n[skip] RUN_FINAL P14b not found at", p14b_path, "\n")
}
cat("done\n")

# =============================================================================
# 05_article_figures.R — single source for every article figure and the
# console diagnostics the prose draws on.
#
# Replaces the legacy 05/05b/06/07 split. All comparisons take the
# autoregression (AR) of Medeiros et al. (2021) as the benchmark and 2SRR-FAVAR
# as the championed specification (the data-rich form faithful to Goulet
# Coulombe, 2025). Forecasts are read from 30_output, with the inflation rate
# target (yout[i,h] = pi_{t+h}) produced by 01_data_prep.R.
#
# Figures
#   FIG1 — Relative RMSE vs AR with Giacomini-White conditional-test p-values
#   FIG2 — Diebold-Mariano heatmap, reference = 2SRR-FAVAR
#   FIG3 — 2SRR-FAVAR coefficient paths (h = 6, top-5 by variability)
#   FIG4 — CSFE of the three TVP specs against AR, 4-panel (h = 1, 3, 6, 12)
#   FIG5 — Cross-validation: share of windows at the grid ceiling
#   FIG6 — Original vs extended grid: RMSE/RW for the three TVP specs
#   FIG7-FIG10 — Realised inflation vs forecasts (Top-3 by RMSE + Ridge + 2SRR),
#                one panel per horizon (y-axis clipped to the realised range)
#
# Console reports (for the article prose):
#   - RMSE(model)/RMSE(AR) and GW (HAC) p-values
#   - DM p-values, 2SRR-FAVAR as reference
#   - Top-3 lowest RMSE per horizon
#   - Final cumulative CSFE vs AR
#   - Mincer-Zarnowitz table with RMSE relative to 2SRR-FAVAR
# =============================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(scales); library(sandwich)
})

# ---- Configuration ----------------------------------------------------------
OUT_FIG   <- "40_results/article_figures/figures"
DIR_FC    <- "30_output/forecasts"
DIR_BETA  <- "30_output/betas"
DIR_DATA  <- "10_data"

# Paths to the 04_analysis.R run folders. Adjust if you rename the run folders.
# If RUN_FINAL does not exist we fall back to the most recent run_final_*.
RUN_FINAL <- "40_results/run_coulombe_grid_linspace(-2_12_15)"
GRID_EXT  <- "40_results/run_extended_grid_overexpanded_linspace_-4_18_25"

dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)

LINE_W <- 0.8
MODEL_COLORS <- c(
  "Realized"          = "black",
  "Ridge"             = "#7F7F7F",
  "2SRR-AR"           = "#D62728",
  "2SRR-Factor"       = "#2CA02C",
  "2SRR-FAVAR"        = "#1F4E9C",
  "LASSO"             = "#FF7F0E",
  "Elastic Net"       = "#17BECF",
  "Adaptive LASSO"    = "#9467BD",
  "Adaptive ElNet"    = "#8C564B",
  "Random Forest"     = "#E377C2",
  "Bagging"           = "#BCBD22",
  "Factor"            = "#1A55A3",
  "Target Factor"     = "#7B4FA3",
  "CSR"               = "#FF1493",
  "AR"                = "#000080",
  "AR-BIC"            = "#4B0082",
  "Ridge-Step1-AR"    = "#FFA07A")

theme_article <- theme_minimal(base_size = 12) +
  theme(legend.position  = "bottom",
        plot.subtitle    = element_text(size = 9, colour = "grey30"),
        panel.grid.minor = element_blank())

save_fig <- function(plt, fname, width = 9, height = 5.5) {
  ggsave(file.path(OUT_FIG, paste0(fname, ".png")), plt,
         width = width, height = height, dpi = 150, bg = "white")
  cat("  [fig]", fname, "\n")
}

# ---- Helpers ----------------------------------------------------------------
ld_fc <- function(f) {
  e <- new.env(); load(file.path(DIR_FC, f), envir = e)
  as.matrix(get(ls(e)[1], envir = e))
}
rmse_fn <- function(y, f) {
  ok <- complete.cases(y, f)
  if (sum(ok) < 5) NA_real_ else sqrt(mean((y[ok] - f[ok])^2))
}
disp <- function(x) {
  m <- c("2SRR_AR" = "2SRR-AR", "2SRR_Factor" = "2SRR-Factor",
         "2SRR_FAVAR" = "2SRR-FAVAR", "AR_BIC" = "AR-BIC",
         "T.Factor" = "Target Factor", "ElNET" = "Elastic Net",
         "AdaLASSO" = "Adaptive LASSO", "AdaElNET" = "Adaptive ElNet",
         "RF" = "Random Forest", "Ridge" = "Ridge", "LASSO" = "LASSO",
         "AR" = "AR", "Bagging" = "Bagging", "Factor" = "Factor",
         "CSR" = "CSR", "RidgeStep1_AR" = "Ridge-Step1-AR")
  out <- m[x]; out[is.na(out)] <- x[is.na(out)]; unname(out)
}

# Giacomini-White (2006) CONDITIONAL predictive ability test. Test function
# h_t = (1, d_{t-1})'; statistic T_n * z_bar' Omega^-1 z_bar ~ chi^2(q = 2),
# with a Newey-West (Bartlett) HAC covariance of bandwidth h - 1. d = loss_a -
# loss_b (squared-error losses); a low p-value rejects equal CONDITIONAL
# predictive ability. This is the genuine GW test promised in Section 3.3 (and
# the one tabulated in 04_analysis.R, P13_GW_test.csv), distinct from the
# unconditional Diebold-Mariano test of Figure 2.
gw_cond <- function(loss_a, loss_b, h = 1) {
  d <- loss_a - loss_b; d <- d[complete.cases(d)]; n <- length(d)
  if (n < 20) return(NA_real_)
  z <- cbind(1, c(NA, d[-n])) * d              # h_t * d_t, n x 2
  z <- z[complete.cases(z), , drop = FALSE]; T_n <- nrow(z)
  if (T_n < 20) return(NA_real_)
  z_bar <- colMeans(z); zc <- z - rep(z_bar, each = T_n)
  Omega <- crossprod(zc) / T_n
  bw <- max(0, h - 1)
  if (bw > 0) for (l in 1:bw) {
    w_l <- 1 - l / (bw + 1)                     # Bartlett weight
    G   <- crossprod(zc[-(1:l), , drop = FALSE],
                     zc[-((T_n - l + 1):T_n), , drop = FALSE]) / T_n
    Omega <- Omega + w_l * (G + t(G))
  }
  stat <- tryCatch(as.numeric(T_n * t(z_bar) %*% solve(Omega) %*% z_bar),
                   error = function(e) NA_real_)
  if (is.na(stat)) NA_real_ else 1 - pchisq(stat, df = ncol(z))
}

# Auto-detect a run folder if the configured one is missing.
if (!dir.exists(RUN_FINAL)) {
  cand <- list.files("40_results", pattern = "^run_(coulombe|final)",
                     full.names = TRUE)
  cand <- cand[dir.exists(cand)]
  if (length(cand) > 0) {
    RUN_FINAL <- cand[which.max(file.mtime(cand))]
    cat("[info] RUN_FINAL auto-detected:", RUN_FINAL, "\n")
  }
}

# ---- Load forecasts and target series ---------------------------------------
load(file.path(DIR_FC,   "yout.rda"))
load(file.path(DIR_FC,   "rw.rda"))
load(file.path(DIR_DATA, "data.rda"))

n_oos     <- nrow(yout)
oos_dates <- tail(as.Date(data$date), n_oos)
hz        <- c(1, 3, 6, 12)

medeiros_models <- c("Ridge", "LASSO", "ElNET", "AdaLASSO", "AdaElNET",
                     "RF", "Bagging", "Factor", "T.Factor", "CSR",
                     "AR", "AR_BIC")
all_fc <- list()
for (m in medeiros_models) {
  fp <- file.path(DIR_FC, paste0(m, ".rda"))
  if (file.exists(fp)) all_fc[[m]] <- ld_fc(paste0(m, ".rda"))
}
for (case in c("AR", "Factor", "FAVAR")) {
  fp <- file.path(DIR_FC, paste0("2SRR_", case, ".rda"))
  if (file.exists(fp)) all_fc[[paste0("2SRR_", case)]] <- ld_fc(paste0("2SRR_", case, ".rda"))
}
fp <- file.path(DIR_FC, "Ridge_from_2SRR_AR.rda")
if (file.exists(fp)) all_fc[["RidgeStep1_AR"]] <- ld_fc("Ridge_from_2SRR_AR.rda")

AR_fc    <- all_fc[["AR"]]               # Medeiros autoregression (benchmark)
ref_FAVAR <- all_fc[["2SRR_FAVAR"]]      # championed spec (faithful to Coulombe)
fc_tvp <- list("2SRR-FAVAR"  = ref_FAVAR,
               "2SRR-Factor" = all_fc[["2SRR_Factor"]],
               "2SRR-AR"     = all_fc[["2SRR_AR"]])

cat(sprintf("Loaded %d forecasts; %d OOS windows; horizons %s\n",
            length(all_fc), n_oos, paste(hz, collapse = ", ")))

# Row order for the heatmaps (top = best, bottom = worst on the AR scale).
mod_levels <- rev(c(
  "Ridge", "Ridge-Step1-AR", "2SRR-FAVAR", "2SRR-Factor",
  "AR", "AR-BIC", "Bagging", "Factor", "Target Factor", "CSR",
  "Adaptive ElNet", "Adaptive LASSO", "Elastic Net", "LASSO",
  "Random Forest", "2SRR-AR"))

# ---- FIG1 — Relative RMSE vs AR + GW HAC p-value ----------------------------
cat("FIG1: GW relative-RMSE vs AR ...\n")
gw_rows <- list()
for (mn in names(all_fc)) {
  if (mn == "AR") next
  M <- all_fc[[mn]]
  for (h in hz) {
    ok <- complete.cases(yout[, h], AR_fc[, h], M[, h]); if (sum(ok) < 20) next
    rel <- rmse_fn(yout[, h], M[, h]) / rmse_fn(yout[, h], AR_fc[, h])
    pv  <- tryCatch(gw_cond((yout[ok, h] - M[ok, h])^2,
                            (yout[ok, h] - AR_fc[ok, h])^2, h),
                    error = function(e) NA_real_)
    gw_rows[[length(gw_rows) + 1]] <- data.frame(
      model = disp(mn), h = h, rel = round(rel, 3), p = round(pv, 3))
  }
}
gw_df <- do.call(rbind, gw_rows)
gw_df$stars <- ifelse(is.na(gw_df$p), "",
                ifelse(gw_df$p < 0.01, "***",
                ifelse(gw_df$p < 0.05, "**",
                ifelse(gw_df$p < 0.10, "*",  ""))))
gw_df$lab   <- sprintf("%.2f\n(%.2f)%s", gw_df$rel, gw_df$p, gw_df$stars)
gw_df$model <- factor(gw_df$model, levels = intersect(mod_levels, unique(gw_df$model)))
fig1 <- ggplot(gw_df, aes(factor(h), model,
                          fill = pmin(pmax(rel, 0.5), 1.5))) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = lab), size = 2.7, lineheight = 0.85) +
  scale_fill_gradient2(midpoint = 1, low = "#1A7F37", mid = "white",
                       high = "#B0202A", limits = c(0.5, 1.5),
                       oob = scales::squish,
                       name = "RMSE(model) / RMSE(AR)") +
  labs(x = "Horizon (h, months)", y = NULL,
       subtitle = "Cell: RMSE relative to AR (below 1, green = beats AR), with the Giacomini-White conditional-test p-value below; * p<0.10, ** p<0.05, *** p<0.01") +
  theme_article + theme(panel.grid = element_blank())
save_fig(fig1, "FIG1_gw_relrmse_vs_AR", width = 8.8, height = 7)

# ---- FIG2 — Diebold-Mariano heatmap, reference = 2SRR-FAVAR ----------------
cat("FIG2: DM heatmap (reference = 2SRR-FAVAR) ...\n")
dm_rows <- list()
for (mn in names(all_fc)) {
  if (mn == "2SRR_FAVAR") next
  M <- all_fc[[mn]]
  for (h in hz) {
    ok <- complete.cases(yout[, h], ref_FAVAR[, h], M[, h]); if (sum(ok) < 20) next
    # Pass RAW forecast errors (y - f); dm.test applies the quadratic loss
    # internally via power = 2 (MSE). Passing pre-squared errors would test an
    # e^4 loss, whose differential variance underflows to zero on the small
    # inflation-rate scale and makes dm.test abort on every cell.
    d <- tryCatch(forecast::dm.test(yout[ok, h] - ref_FAVAR[ok, h],
                                    yout[ok, h] - M[ok, h],
                                    alternative = "two.sided", h = h, power = 2),
                  error = function(e) NULL)
    if (is.null(d)) next
    dm_rows[[length(dm_rows) + 1]] <- data.frame(
      model = disp(mn), h = h, p = round(d$p.value, 4),
      ref_wins = isTRUE(rmse_fn(yout[, h], ref_FAVAR[, h]) <
                          rmse_fn(yout[, h], M[, h])))
  }
}
dm_df <- if (length(dm_rows)) do.call(rbind, dm_rows) else NULL
if (is.null(dm_df) || nrow(dm_df) == 0) {
  cat("  [skip] FIG2: no valid Diebold-Mariano cells\n")
} else {
dm_df$signed_p <- ifelse(dm_df$ref_wins, -dm_df$p, dm_df$p)
dm_df$model    <- factor(dm_df$model,
                         levels = intersect(mod_levels, unique(dm_df$model)))
fig2 <- ggplot(dm_df, aes(factor(h), model,
                          fill = pmin(pmax(signed_p, -0.5), 0.5))) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", p)), size = 3) +
  scale_fill_gradient2(midpoint = 0, low = "#1A7F37", mid = "white",
                       high = "#B0202A", limits = c(-0.5, 0.5),
                       name = "Signed p-value") +
  labs(x = "Horizon (h, months)", y = NULL,
       subtitle = "Reference: 2SRR-FAVAR. Green = 2SRR-FAVAR more accurate; red = the other model. Cell = DM p-value") +
  theme_article + theme(panel.grid = element_blank())
save_fig(fig2, "FIG2_dm_heatmap_2SRR_FAVAR", width = 8.5, height = 6.5)
}

# ---- FIG3 — 2SRR-FAVAR coefficient paths (h = 6, top-5 by variability) ------
cat("FIG3: 2SRR-FAVAR coefficient paths (h = 6) ...\n")
bp <- file.path(DIR_BETA, "betas_2SRR_FAVAR.rda")
if (file.exists(bp)) {
  e <- new.env(); load(bp, envir = e)
  betas_bundle <- get(ls(e)[1], envir = e)
  bh <- betas_bundle[["h6"]]
  n_w <- length(bh)
  fv  <- which(vapply(bh, function(b)
    !is.null(b) && !is.null(b$betas_tvp) && is.matrix(b$betas_tvp) &&
      nrow(b$betas_tvp) > 0, logical(1)))[1]
  if (!is.na(fv)) {
    K   <- ncol(bh[[fv]]$betas_tvp)
    vn  <- bh[[fv]]$var_names
    mat <- matrix(NA_real_, n_w, K)
    for (wi in seq_len(n_w)) {
      b <- bh[[wi]]
      if (!is.null(b) && !is.null(b$betas_tvp) && ncol(b$betas_tvp) == K)
        mat[wi, ] <- b$betas_tvp[nrow(b$betas_tvp), ]
    }
    colnames(mat) <- vn
    # FAVAR carries many coefficients (intercept + inflation lags + factor
    # loadings); show the five most time-varying for legibility.
    sds  <- apply(mat, 2, sd, na.rm = TRUE)
    topk <- names(sort(sds, decreasing = TRUE))[seq_len(min(5, ncol(mat)))]
    df3   <- data.frame(date = tail(oos_dates, n_w), mat[, topk, drop = FALSE],
                        check.names = FALSE)
    long3 <- tidyr::pivot_longer(df3, -date, names_to = "coef", values_to = "value")
    long3$coef <- factor(long3$coef, levels = topk)
    long3 <- long3[!is.na(long3$value), ]
    fig3 <- ggplot(long3, aes(date, value, colour = coef)) +
      geom_hline(yintercept = 0, linetype = 3, colour = "grey50") +
      geom_line(linewidth = 0.8, na.rm = TRUE) +
      scale_colour_brewer(palette = "Dark2", name = NULL) +
      scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
      labs(x = "Forecast origin", y = "Estimated coefficient",
           subtitle = "2SRR-FAVAR: five most time-varying coefficients (inflation lags + factor loadings), h = 6") +
      theme_article
    save_fig(fig3, "FIG3_coef_paths_2srrFAVAR_h6", width = 9.5, height = 5.5)
  } else cat("  [skip] no valid TVP betas in betas_2SRR_FAVAR.rda\n")
} else cat("  [skip] betas_2SRR_FAVAR.rda not found\n")

# ---- FIG4 — CSFE vs AR, 4-panel ---------------------------------------------
cat("FIG4: CSFE vs AR (h = 1, 3, 6, 12) ...\n")
# cumsum that treats an NA increment as 0, so the trailing windows where the
# realised h-step value does not yet exist (yout[,h] is NA there) do not blank
# out the whole curve — the line simply flattens over that unrealised tail.
csum0 <- function(v) cumsum(ifelse(is.na(v), 0, v))
csfe_dfs <- lapply(hz, function(hC) {
  e_b <- (yout[, hC] - AR_fc[, hC])^2
  data.frame(
    date = oos_dates,
    h    = factor(paste0("h = ", hC), levels = paste0("h = ", hz)),
    `2SRR-AR`     = csum0(e_b - (yout[, hC] - fc_tvp[["2SRR-AR"]][, hC])^2),
    `2SRR-Factor` = csum0(e_b - (yout[, hC] - fc_tvp[["2SRR-Factor"]][, hC])^2),
    `2SRR-FAVAR`  = csum0(e_b - (yout[, hC] - fc_tvp[["2SRR-FAVAR"]][, hC])^2),
    check.names = FALSE)
})
csfe_long <- tidyr::pivot_longer(do.call(rbind, csfe_dfs),
                                 cols = -c(date, h),
                                 names_to = "Model", values_to = "CSFE")
csfe_long$Model <- factor(csfe_long$Model,
                          levels = c("2SRR-AR", "2SRR-Factor", "2SRR-FAVAR"))
fig4 <- ggplot(csfe_long, aes(date, CSFE, colour = Model)) +
  annotate("rect", xmin = as.Date("2020-02-01"), xmax = as.Date("2020-09-01"),
           ymin = -Inf, ymax = Inf, alpha = 0.13, fill = "grey40") +
  annotate("rect", xmin = as.Date("2021-04-01"), xmax = as.Date("2023-01-01"),
           ymin = -Inf, ymax = Inf, alpha = 0.10, fill = "#F4A582") +
  geom_hline(yintercept = 0, linetype = 3, colour = "grey50") +
  geom_line(linewidth = LINE_W, na.rm = TRUE) +
  facet_wrap(~ h, ncol = 2, scales = "free_y") +
  scale_colour_manual(values = MODEL_COLORS, name = NULL) +
  scale_x_date(date_breaks = "3 years", date_labels = "%Y") +
  labs(x = "Forecast origin",
       y = "Cumulative squared-error difference vs. the autoregression",
       subtitle = "Rising = beats AR; falling = trails AR. Grey: COVID 2020; orange: 2021-22 surge") +
  theme_article
save_fig(fig4, "FIG4_csfe_vs_AR_4h", width = 11, height = 7)

# ---- FIG5 — Cross-validation grid saturation --------------------------------
cat("FIG5: grid saturation ...\n")
f5 <- file.path(RUN_FINAL, "tables", "P5b_lambda_saturation.csv")
if (file.exists(f5)) {
  p5 <- read.csv(f5, stringsAsFactors = FALSE)
  p5 <- p5[p5$step == "lambda", ]
  p5$case <- factor(p5$case, levels = c("AR", "Factor", "FAVAR"),
                    labels = c("2SRR-AR", "2SRR-Factor", "2SRR-FAVAR"))
  fig5 <- ggplot(p5, aes(factor(h), pct_at_top, fill = case)) +
    geom_col(position = position_dodge(0.8), width = 0.75) +
    geom_text(aes(label = sprintf("%.0f%%", pct_at_top)),
              position = position_dodge(0.8), vjust = -0.35, size = 3) +
    geom_hline(yintercept = 50, linetype = 2, colour = "grey50") +
    scale_fill_manual(values = MODEL_COLORS, name = NULL) +
    scale_y_continuous(limits = c(0, 100)) +
    labs(x = "Horizon (h, months)",
         y = "Share of windows (%) at the grid ceiling",
         subtitle = "Cross-validation picks the largest available penalty in 2/3 to 9/10 of windows") +
    theme_article
  save_fig(fig5, "FIG5_lambda_saturation", width = 9, height = 5.5)
} else cat("  [skip] P5b_lambda_saturation.csv not found\n")

# ---- FIG6 — Original vs extended grid ---------------------------------------
cat("FIG6: original vs extended grid ...\n")
f_o <- file.path(RUN_FINAL, "tables", "P1_rmsfe_relative_rw_all_h.csv")
f_e <- file.path(GRID_EXT,  "tables", "P1_rmsfe_relative_rw_all_h.csv")
if (file.exists(f_o) && file.exists(f_e)) {
  pf <- function(df, grid_lab) {
    df <- df[df$model %in% c("2SRR_AR", "2SRR_Factor", "2SRR_FAVAR"),
             c("model", "h1", "h3", "h6", "h12")]
    df$model <- factor(disp(df$model),
                       levels = c("2SRR-AR", "2SRR-Factor", "2SRR-FAVAR"))
    df$grid <- grid_lab
    tidyr::pivot_longer(df, c(-model, -grid), names_to = "h", values_to = "ratio")
  }
  d6 <- rbind(pf(read.csv(f_o, stringsAsFactors = FALSE), "Original (exp(-2, 12))"),
              pf(read.csv(f_e, stringsAsFactors = FALSE), "Extended (exp(-4, 18))"))
  d6$h <- factor(d6$h, levels = c("h1", "h3", "h6", "h12"),
                 labels = c("1", "3", "6", "12"))
  fig6 <- ggplot(d6, aes(h, ratio, fill = grid)) +
    geom_col(position = position_dodge(0.8), width = 0.75) +
    geom_hline(yintercept = 1, linetype = 2) +
    geom_text(aes(label = sprintf("%.2f", ratio)),
              position = position_dodge(0.8), vjust = -0.35, size = 3) +
    facet_wrap(~ model, nrow = 1) +
    scale_fill_manual(values = c("Original (exp(-2, 12))" = "#2CA02C",
                                 "Extended (exp(-4, 18))" = "#D62728"),
                      name = NULL) +
    labs(x = "Horizon (h, months)", y = "RMSE / RMSE(RW)",
         subtitle = "Extending the penalty grid degrades out-of-sample accuracy in 8 of 12 cells; the exceptions fall at h = 12") +
    theme_article
  save_fig(fig6, "FIG6_grid_comparison", width = 11, height = 5)
} else cat("  [skip] grid-comparison CSVs not found (need both run folders)\n")

# ---- FIG7-FIG10 — Realised inflation vs forecasts ---------------------------
cat("FIG7-FIG10: realised vs forecasts (one panel per horizon) ...\n")
top3_per_h <- function(h) {
  rmse_h <- sapply(all_fc, function(M) rmse_fn(yout[, h], M[, h]))
  rmse_h <- rmse_h[!is.na(rmse_h)]
  names(sort(rmse_h)[seq_len(min(3, length(rmse_h)))])
}
plot_real_vs_fc <- function(h) {
  must <- c("Ridge", "2SRR_FAVAR", "2SRR_Factor", "2SRR_AR")
  sel  <- unique(c(top3_per_h(h), intersect(must, names(all_fc))))
  df   <- data.frame(date = oos_dates, Realized = yout[, h])
  for (m in sel) df[[disp(m)]] <- all_fc[[m]][, h]
  long <- tidyr::pivot_longer(df, -date, names_to = "Series", values_to = "value")
  realised  <- long[long$Series == "Realized", ]
  forecasts <- long[long$Series != "Realized", ]
  forecasts$Series <- factor(forecasts$Series,
                             levels = setdiff(unique(long$Series), "Realized"))
  # Focus the y-axis on the realised range (plus margin). At long horizons the
  # high-dimensional ridge overshoots far beyond it (it chases noise); clipping
  # keeps the realised-vs-forecast comparison legible. coord_cartesian clips the
  # view without dropping data from the lines.
  yr  <- range(yout[, h], na.rm = TRUE); pad <- 0.25 * diff(yr)
  ylim <- c(yr[1] - pad, yr[2] + pad)
  ggplot(mapping = aes(date, value)) +
    geom_line(data = forecasts, aes(colour = Series), linewidth = 0.6, na.rm = TRUE) +
    geom_line(data = realised,  colour = "black",     linewidth = 0.6, na.rm = TRUE) +
    scale_colour_manual(values = MODEL_COLORS, name = NULL) +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    coord_cartesian(ylim = ylim) +
    labs(x = "Date", y = sprintf("Inflation rate (%%), h = %d", h),
         subtitle = sprintf(
           "Realised (black) vs forecasts (Top-3 by RMSE + Ridge + 2SRR); y-axis clipped to the realised range, h = %d", h)) +
    theme_article
}
for (idx in seq_along(hz)) {
  h <- hz[idx]
  save_fig(plot_real_vs_fc(h),
           sprintf("FIG%d_realized_vs_fc_h%02d", 6 + idx, h),
           width = 10, height = 4.8)
}

# ---- FIG11 — Monthly inflation vs forecasts (P7 style, natural scale) -------
# Realised MONTHLY inflation y_t (the actual observed series, the same line in
# every panel) against the h-step forecasts, in the data's natural monthly
# scale (percentage points). Unlike the legacy P7 in 04_analysis.R, the
# forecasts are NOT divided by h: with the inflation-rate target (yout = pi_t+h)
# the model already produces a monthly-rate forecast, so dividing by h would
# wrongly push the lines toward zero. 2SRR-FAVAR is the champion, shown with
# Ridge and the leading Medeiros models, in the 05 colour scheme. The realised
# line is the contemporaneous monthly inflation at the forecast origin; each
# coloured line is that model's forecast made h months ahead.
cat("FIG11: monthly inflation vs forecasts (P7 style, h = 1,3,6,12) ...\n")
y_oos_monthly <- tail(as.numeric(data$CPIAUCSL), n_oos)
mi_models <- intersect(c("2SRR_FAVAR", "2SRR_Factor", "2SRR_AR", "Ridge", "LASSO", "RF"),
                       names(all_fc))
mi_dfs <- lapply(hz, function(hC) {
  d <- data.frame(date = oos_dates,
                  h = factor(paste0("h = ", hC), levels = paste0("h = ", hz)),
                  Realized = y_oos_monthly, check.names = FALSE)
  for (m in mi_models) d[[disp(m)]] <- all_fc[[m]][, hC]
  d
})
mi_long <- tidyr::pivot_longer(do.call(rbind, mi_dfs), cols = -c(date, h),
                               names_to = "Series", values_to = "value")
mi_real <- mi_long[mi_long$Series == "Realized", ]
mi_fc   <- mi_long[mi_long$Series != "Realized", ]
# Keep all series in the legend (FAVAR first for prominence) ...
mi_fc$Series <- factor(mi_fc$Series, levels = disp(mi_models))
# ... but DRAW 2SRR-FAVAR last (on top) and thicker, so it stands out against
# the other, thinner forecast lines (the realised series sits just beneath it).
mi_fc_oth <- mi_fc[mi_fc$Series != "2SRR-FAVAR", ]
mi_fc_fav <- mi_fc[mi_fc$Series == "2SRR-FAVAR", ]
yr  <- range(y_oos_monthly, na.rm = TRUE); pad <- 0.30 * diff(yr)
fig11 <- ggplot(mapping = aes(date, value)) +
  geom_hline(yintercept = 0, linetype = 3, colour = "grey60") +
  geom_line(data = mi_fc_oth, aes(colour = Series), linewidth = 0.45,
            alpha = 0.85, na.rm = TRUE) +
  geom_line(data = mi_real,   colour = "black", linewidth = 0.45, na.rm = TRUE) +
  geom_line(data = mi_fc_fav, aes(colour = Series), linewidth = 0.45, na.rm = TRUE) +
  facet_wrap(~ h, ncol = 2) +
  scale_colour_manual(values = MODEL_COLORS, name = NULL) +
  scale_x_date(date_breaks = "3 years", date_labels = "%Y") +
  coord_cartesian(ylim = c(yr[1] - pad, yr[2] + pad)) +
  labs(x = "Forecast origin",
       y = "Monthly inflation %",
       subtitle = "Realised monthly inflation (black) vs h-step forecasts; 2SRR-FAVAR (blue). Natural monthly scale, y-axis clipped to the realised range") +
  theme_article
save_fig(fig11, "FIG11_monthly_inflation_vs_fc_4h", width = 11, height = 7)

# ---- Console diagnostics for the prose --------------------------------------
cat("\n=== RMSE(model)/RMSE(AR) and Giacomini-White CONDITIONAL p-value ===\n")
cat(sprintf("%-16s | %s\n", "spec",
            paste(sprintf("h=%-2d rel (gwp)", hz), collapse = "  ")))
for (nm in c("2SRR_FAVAR", "2SRR_Factor", "2SRR_AR")) {
  M <- all_fc[[nm]]; if (is.null(M)) next
  line <- sprintf("%-16s |", disp(nm))
  for (h in hz) {
    ok  <- complete.cases(yout[, h], AR_fc[, h], M[, h])
    rel <- rmse_fn(yout[, h], M[, h]) / rmse_fn(yout[, h], AR_fc[, h])
    gwp <- tryCatch(gw_cond((yout[ok, h] - M[ok, h])^2,
                            (yout[ok, h] - AR_fc[ok, h])^2, h),
                    error = function(e) NA_real_)
    line <- paste0(line, sprintf("  %.3f (%.3f)", rel, gwp))
  }
  cat(line, "\n")
}

cat("\n=== DM two-sided p-value with 2SRR-FAVAR as the reference ===\n")
for (lab in if (is.null(dm_df)) character(0) else
            c("AR", "2SRR-Factor", "2SRR-AR", "Random Forest",
              "LASSO", "Elastic Net", "Ridge")) {
  rr <- dm_df[dm_df$model == lab, ]
  if (nrow(rr) == 0) next
  rr <- rr[order(rr$h), ]
  cat(sprintf("  vs %-16s : %s\n", lab,
              paste(sprintf("h=%2d p=%.3f%s", rr$h, rr$p,
                            ifelse(rr$ref_wins, "(FAVAR+)", "(rival+)")),
                    collapse = "  ")))
}

cat("\n=== Top-3 lowest RMSE per horizon ===\n")
for (h in hz)
  cat(sprintf("  h=%2d: %s\n", h, paste(disp(top3_per_h(h)), collapse = ", ")))

cat("\n=== Final cumulative CSFE vs AR (positive = beats AR) ===\n")
last_finite <- function(z) { z <- z[is.finite(z)]; if (length(z)) z[length(z)] else NA_real_ }
for (hC in hz) {
  e_b  <- (yout[, hC] - AR_fc[, hC])^2
  vals <- sapply(fc_tvp, function(M)
    last_finite(cumsum(ifelse(is.na(e_b - (yout[, hC] - M[, hC])^2), 0,
                              e_b - (yout[, hC] - M[, hC])^2))))
  cat(sprintf("  h=%2d: %s\n", hC,
              paste(sprintf("%s=%+.3f", names(vals), vals), collapse = "  ")))
}

# TABLE 1 (Section 4.5) — RMSE RELATIVE TO THE AUTOREGRESSION paired with the
# Mincer-Zarnowitz forecast-efficiency p-value (Newey-West HAC). AR is the
# common benchmark, matching FIG1 and the GW analysis. Emitted both as a console
# table (the article's Table 1) and as a heatmap (FIG12) coloured by the MZ
# p-value: green (p > 0.10) = efficiency not rejected, red = rejected.
mz_test_hac <- function(y, f, h) {
  ok <- complete.cases(y, f); if (sum(ok) < 30) return(NA_real_)
  y_ <- y[ok]; f_ <- f[ok]
  m  <- tryCatch(lm(y_ ~ f_), error = function(e) NULL); if (is.null(m)) return(NA_real_)
  vcv <- tryCatch(sandwich::NeweyWest(m, lag = max(0, h - 1),
                                      prewhite = FALSE, adjust = TRUE),
                  error = function(e) vcov(m))
  d <- coef(m) - c(0, 1)
  w <- tryCatch(as.numeric(t(d) %*% solve(vcv) %*% d),
                error = function(e) NA_real_)
  if (is.na(w) || w < 0) NA_real_ else 1 - pchisq(w, df = 2)
}
mz_models <- c("AR", "2SRR_FAVAR", "2SRR_Factor", "2SRR_AR", "RF", "LASSO", "ElNET")
mz_rows <- list()
for (mn in mz_models) {
  M <- all_fc[[mn]]; if (is.null(M)) next
  for (h in hz) {
    rel <- rmse_fn(yout[, h], M[, h]) / rmse_fn(yout[, h], AR_fc[, h])
    pv  <- mz_test_hac(yout[, h], M[, h], h)
    mz_rows[[length(mz_rows) + 1]] <- data.frame(
      model = disp(mn), h = h, rel = round(rel, 3), mzp = round(pv, 3))
  }
}
mz_tab <- do.call(rbind, mz_rows)

cat("\n=== TABLE 1: RMSE relative to AR  /  Mincer-Zarnowitz HAC p-value ===\n")
cat(sprintf("%-14s | %s\n", "model",
            paste(sprintf("h=%-2d  rel / MZp", hz), collapse = "  ")))
for (mn in disp(mz_models)) {
  rr <- mz_tab[mz_tab$model == mn, ]; if (!nrow(rr)) next
  rr <- rr[match(hz, rr$h), ]
  cat(sprintf("%-14s |", mn))
  for (i in seq_along(hz)) cat(sprintf("  %5.3f / %5.3f", rr$rel[i], rr$mzp[i]))
  cat("\n")
}

# FIG12 — heatmap of Table 1 (coloured by the MZ p-value).
mz_hm <- mz_tab
mz_hm$model <- factor(mz_hm$model, levels = rev(disp(mz_models)))   # AR on top
mz_hm$lab   <- sprintf("%.2f\n(%.3f)", mz_hm$rel, mz_hm$mzp)
fig12 <- ggplot(mz_hm, aes(factor(h), model, fill = pmin(pmax(mzp, 0), 0.5))) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = lab), size = 2.8, lineheight = 0.85) +
  scale_fill_gradient2(midpoint = 0.10, low = "#B0202A", mid = "white",
                       high = "#1A7F37", limits = c(0, 0.5),
                       oob = scales::squish, name = "Mincer-Zarnowitz p-value") +
  labs(x = "Horizon (h, months)", y = NULL,
       subtitle = "RMSE relative to AR, with the MZ joint p-value below. Green (p > 0.10) = efficiency not rejected; red = rejected") +
  theme_article + theme(panel.grid = element_blank())
save_fig(fig12, "FIG12_mincer_zarnowitz_relAR", width = 8.8, height = 5.5)

cat("\n== 05_article_figures.R done ==\n")
cat("Figures:", normalizePath(OUT_FIG), "\n")

# ---- 0. Configuration: paths to the result folders --------------------------
RUN_FINAL <- "40_results/run_final_20260519_1817"
RUN_V2    <- "40_results/run_v2_20260519_1822"
GRID_ORIG <- "40_results/run_coulombe_grid_linspace(-2_12_15)"
GRID_EXT  <- "40_results/run_extended_grid_overexpanded_linspace_-4_18_25"
DIR_BETAS <- "30_output/betas"
DIR_DATA  <- "10_data"

OUT_FIG <- "40_results/article_figures/figures"
OUT_TAB <- "40_results/article_figures/tables"
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)

# ---- 1. Packages ------------------------------------------------------------
need <- c("ggplot2", "dplyr", "tidyr", "scales")
for (p in need)
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cran.r-project.org")
suppressPackageStartupMessages(invisible(lapply(need, library, character.only = TRUE)))
has_flex <- requireNamespace("flextable", quietly = TRUE)
if (!has_flex) {
  try(install.packages("flextable", repos = "https://cran.r-project.org"), silent = TRUE)
  has_flex <- requireNamespace("flextable", quietly = TRUE)
}

# ---- 2. House style: palette, theme, helpers --------------------------------
MODEL_COLORS <- c(
  "2SRR-FAVAR" = "#D62728", "2SRR-AR" = "#FF7F0E", "2SRR-Factor" = "#2CA02C",
  "LASSO" = "#1F77B4", "Elastic Net" = "#17BECF", "Adaptive LASSO" = "#9467BD",
  "Adaptive ElNet" = "#8C564B", "Random Forest" = "#E377C2", "Bagging" = "#7F7F7F",
  "Factor" = "#BCBD22", "Target Factor" = "#1A55A3", "CSR" = "#FF1493",
  "AR" = "#000080", "AR-BIC" = "#4B0082", "Ridge (117 series)" = "#8B0000",
  "Ridge-Step1-AR" = "#FFA07A", "Ridge-Step1-Factor" = "#A0522D",
  "Ridge-Step1-FAVAR" = "#DEB887")

theme_article <- theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        legend.title    = element_text(size = 10),
        plot.subtitle   = element_text(size = 9, colour = "grey30"),
        panel.grid.minor = element_blank())

# Pretty display names for the raw model codes used in the CSVs
disp <- function(x) {
  m <- c("2SRR_AR" = "2SRR-AR", "2SRR_Factor" = "2SRR-Factor",
         "2SRR_FAVAR" = "2SRR-FAVAR", "RidgeStep1_AR" = "Ridge-Step1-AR",
         "RidgeStep1_Factor" = "Ridge-Step1-Factor",
         "RidgeStep1_FAVAR" = "Ridge-Step1-FAVAR", "AR_BIC" = "AR-BIC",
         "T.Factor" = "Target Factor", "ElNET" = "Elastic Net",
         "AdaLASSO" = "Adaptive LASSO", "AdaElNET" = "Adaptive ElNet",
         "RF" = "Random Forest", "Ridge" = "Ridge (117 series)")
  out <- m[x]; out[is.na(out)] <- x[is.na(out)]; unname(out)
}

save_article_fig <- function(plt, fname, width = 9, height = 5.5) {
  ggsave(file.path(OUT_FIG, paste0(fname, ".pdf")), plt, width = width, height = height)
  ggsave(file.path(OUT_FIG, paste0(fname, ".png")), plt, width = width, height = height,
         dpi = 300, bg = "white")
  cat("  [fig]", fname, ".pdf / .png\n", sep = "")
}

rd <- function(folder, file) read.csv(file.path(folder, "tables", file),
                                       check.names = FALSE, stringsAsFactors = FALSE)

# =============================================================================
# FIGURE 1 — RMSE relative to the random walk: full field as a heatmap
#   Replaces the 18-row table. Green = beats RW (<1), red = worse (>1).
# =============================================================================
cat("\nFIGURE 1: RMSE heatmap (all models)...\n")
p1 <- rd(RUN_FINAL, "P1_rmsfe_relative_rw_all_h.csv")
p1 <- p1[, c("model", "h1", "h3", "h6", "h12")]
p1$model <- disp(p1$model)
long1 <- tidyr::pivot_longer(p1, -model, names_to = "h", values_to = "ratio")
long1$h <- factor(long1$h, levels = c("h1", "h3", "h6", "h12"),
                  labels = c("1", "3", "6", "12"))
# y-axis order (bottom -> top): worst families at the bottom, 2SRR at the top
ord <- c("Random walk", "Ridge (117 series)", "Factor", "CSR", "Target Factor",
         "Bagging", "Adaptive ElNet", "Adaptive LASSO", "Elastic Net", "LASSO",
         "Random Forest", "Ridge-Step1-FAVAR", "Ridge-Step1-Factor",
         "Ridge-Step1-AR", "AR-BIC", "AR", "2SRR-Factor", "2SRR-FAVAR", "2SRR-AR")
long1$model <- factor(ifelse(long1$model == "RW", "Random walk", long1$model),
                      levels = ord)
fig1 <- ggplot(long1, aes(h, model, fill = ratio)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", ratio)), size = 3) +
  scale_fill_gradient2(midpoint = 1, low = "#1A7F37", mid = "white",
                       high = "#B0202A", limits = c(0.5, 1.3),
                       oob = scales::squish, name = "RMSE / RMSE(RW)") +
  labs(x = "Horizon (h, months)", y = NULL,
       subtitle = "Values below 1.00 (green) beat the random walk; above 1.00 (red) are worse") +
  # ggtitle("Out-of-sample accuracy relative to the random walk") +
  theme_article + theme(panel.grid = element_blank())
save_article_fig(fig1, "FIG1_rmse_heatmap", width = 8.5, height = 7)

# =============================================================================
# FIGURE 2 — The value of parsimony + time variation: 2SRR-AR vs two ridges
#   High-dimensional Medeiros ridge (117 series) vs Step-1 ridge (2 lags, static)
#   vs full 2SRR-AR (2 lags, time-varying). RMSE relative to RW.
# =============================================================================
cat("FIGURE 2: 2SRR-AR vs static ridges...\n")
sel <- c("Ridge", "RidgeStep1_AR", "2SRR_AR")
p2 <- p1[FALSE, ]  # not used; rebuild from raw to keep raw codes
raw1 <- rd(RUN_FINAL, "P1_rmsfe_relative_rw_all_h.csv")
d2 <- raw1[raw1$model %in% sel, c("model", "h1", "h3", "h6", "h12")]
d2$model <- factor(disp(d2$model),
                   levels = c("Ridge (117 series)", "Ridge-Step1-AR", "2SRR-AR"))
long2 <- tidyr::pivot_longer(d2, -model, names_to = "h", values_to = "ratio")
long2$h <- factor(long2$h, levels = c("h1", "h3", "h6", "h12"),
                  labels = c("1", "3", "6", "12"))
fig2 <- ggplot(long2, aes(h, ratio, fill = model)) +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_text(aes(label = sprintf("%.2f", ratio)),
            position = position_dodge(0.8), vjust = -0.35, size = 3) +
  scale_fill_manual(values = MODEL_COLORS, name = NULL) +
  labs(x = "Horizon (h, months)", y = "RMSE / RMSE(RW)",
       subtitle = "Two inflation lags beat the full 117-series ridge; time variation (2SRR-AR) sharpens the edge") +
  theme_article
save_article_fig(fig2, "FIG2_2srr_vs_ridges", width = 9, height = 5.5)

# =============================================================================
# FIGURE 3 — Where the parameters move: 2SRR-AR coefficient paths (h = 6)
#   Per window, the last in-sample beta (the one used for the forecast).
# =============================================================================
cat("FIGURE 3: 2SRR-AR coefficient paths (h=6)...\n")
load(file.path(DIR_DATA, "data.rda"))            # provides `data` with $date
all_dates <- as.Date(data$date)
fig3 <- NULL
bp <- file.path(DIR_BETAS, "betas_2SRR_AR.rda")
if (file.exists(bp)) {
  e <- new.env(); load(bp, envir = e)
  betas_bundle <- get(ls(e)[1], envir = e)
  bh <- betas_bundle[["h6"]]
  n_oos <- length(bh)
  oos_dates <- tail(all_dates, n_oos)
  fv <- which(vapply(bh, function(b)
    !is.null(b) && !is.null(b$betas_tvp) && is.matrix(b$betas_tvp) &&
      nrow(b$betas_tvp) > 0, logical(1)))[1]
  K  <- ncol(bh[[fv]]$betas_tvp)
  vn <- bh[[fv]]$var_names
  mat <- matrix(NA_real_, n_oos, K)
  for (wi in seq_len(n_oos)) {
    b <- bh[[wi]]
    if (!is.null(b) && !is.null(b$betas_tvp) && ncol(b$betas_tvp) == K)
      mat[wi, ] <- b$betas_tvp[nrow(b$betas_tvp), ]
  }
  colnames(mat) <- vn
  df3 <- data.frame(date = oos_dates, mat, check.names = FALSE)
  long3 <- tidyr::pivot_longer(df3, -date, names_to = "coef", values_to = "value")
  # Friendly labels: intercept + the two own lags
  long3$coef <- dplyr::recode(long3$coef,
                              "intercept" = "Intercept (local mean)",
                              "Yh_L6" = "Inflation lag 6",
                              "Yh_L7" = "Inflation lag 7")
  # Two of the 180 windows (Nov/2021, Aug/2023) used a static-ridge fallback when
  # the time-varying optimizer did not converge under extreme volatility, so they
  # have no TVP betas. Dropping these NAs lets the line bridge them (a short,
  # near-invisible straight segment) instead of breaking. Disclose in a footnote.
  n_fallback <- sum(is.na(long3$value)) / dplyr::n_distinct(long3$coef)
  if (n_fallback > 0)
    cat(sprintf("  note: %d window(s) used ridge fallback (no TVP betas); line bridges them.\n",
                round(n_fallback)))
  long3 <- long3[!is.na(long3$value), ]
  fig3 <- ggplot(long3, aes(date, value, colour = coef)) +
    geom_hline(yintercept = 0, linetype = 3, colour = "grey50") +
    geom_line(linewidth = 0.8, na.rm = TRUE) +
    scale_colour_manual(values = c("Intercept (local mean)" = "#D62728",
                                   "Inflation lag 7" = "#1F77B4",
                                   "Inflation lag 6" = "#2CA02C"), name = NULL) +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    labs(x = "Forecast origin", y = "Estimated coefficient",
         subtitle = "The intercept re-anchors the inflation mean; the autoregressive slope stays stable") +
    theme_article
  save_article_fig(fig3, "FIG3_coef_paths_2srrAR_h6", width = 9.5, height = 5.5)
} else cat("  [skip] betas_2SRR_AR.rda not found\n")

# =============================================================================
# FIGURE 4 — Diebold-Mariano heatmap: 2SRR-AR (championed spec) vs each model
#   Computed fresh from the saved forecasts so the REFERENCE is 2SRR-AR (the
#   spec the analysis champions), with the sibling 2SRR specs shown as rows.
#   The DM call mirrors 04_v2_analysis.R exactly, so this is directly
#   comparable to the 2SRR-FAVAR heatmap already in run_v2.
# =============================================================================
cat("FIGURE 4: Diebold-Mariano heatmap (reference = 2SRR-AR)...\n")
fig4 <- NULL
DIR_FC <- "30_output/forecasts"
if (requireNamespace("forecast", quietly = TRUE) &&
    file.exists(file.path(DIR_FC, "2SRR_AR.rda"))) {
  ld_fc <- function(f) { e <- new.env(); load(file.path(DIR_FC, f), envir = e)
                         as.matrix(get(ls(e)[1], envir = e)) }
  load(file.path(DIR_FC, "yout.rda"))            # realized cumulative target
  fc_ref <- ld_fc("2SRR_AR.rda")
  bench_map <- c(
    "2SRR_FAVAR.rda" = "2SRR-FAVAR", "2SRR_Factor.rda" = "2SRR-Factor",
    "Ridge_from_2SRR_AR.rda" = "Ridge-Step1-AR", "AR.rda" = "AR",
    "AR_BIC.rda" = "AR-BIC", "RF.rda" = "Random Forest", "LASSO.rda" = "LASSO",
    "ElNET.rda" = "Elastic Net", "AdaLASSO.rda" = "Adaptive LASSO",
    "AdaElNET.rda" = "Adaptive ElNet", "Bagging.rda" = "Bagging",
    "CSR.rda" = "CSR", "Factor.rda" = "Factor", "T.Factor.rda" = "Target Factor",
    "Ridge.rda" = "Ridge (117 series)")
  hz <- c(1, 3, 6, 12)
  rmse_fn <- function(y, f) { ok <- complete.cases(y, f)
    if (sum(ok) < 5) NA_real_ else sqrt(mean((y[ok] - f[ok])^2)) }
  dm_rows <- list()
  for (bf in names(bench_map)) {
    if (!file.exists(file.path(DIR_FC, bf))) next
    M <- ld_fc(bf)
    for (h in hz) {
      if (h > ncol(M) || h > ncol(fc_ref)) next
      ok <- complete.cases(yout[, h], fc_ref[, h], M[, h]); if (sum(ok) < 20) next
      d <- tryCatch(forecast::dm.test(ts((yout[ok, h] - fc_ref[ok, h])^2),
                                      ts((yout[ok, h] - M[ok, h])^2),
                                      alternative = "two.sided", h = h),
                    error = function(e) NULL)
      if (is.null(d)) next
      dm_rows[[length(dm_rows) + 1]] <- data.frame(
        benchmark = unname(bench_map[bf]), h = h, DM_p = round(d$p.value, 4),
        ref_wins = isTRUE(rmse_fn(yout[, h], fc_ref[, h]) <
                          rmse_fn(yout[, h], M[, h])))
    }
  }
  dm <- do.call(rbind, dm_rows)
  # green (negative) = 2SRR-AR more accurate; red (positive) = the other model.
  # NOTE: fixes the legend inversion present in the original 04_v2 P_DM figure.
  dm$signed_p  <- ifelse(dm$ref_wins, -dm$DM_p, dm$DM_p)
  dm$benchmark <- factor(dm$benchmark, levels = rev(unname(bench_map)))
  fig4 <- ggplot(dm, aes(factor(h), benchmark,
                         fill = pmin(pmax(signed_p, -0.5), 0.5))) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.2f", DM_p)), size = 3) +
    scale_fill_gradient2(midpoint = 0, low = "#1A7F37", mid = "white",
                         high = "#B0202A", limits = c(-0.5, 0.5),
                         name = "Signed p-value") +
    labs(x = "Horizon (h, months)", y = NULL,
         subtitle = "Reference: 2SRR-AR. Green = 2SRR-AR more accurate; red = the other model. Cell = DM p-value") +
    theme_article + theme(panel.grid = element_blank())
  save_article_fig(fig4, "FIG4_dm_heatmap_2SRR_AR", width = 8.5, height = 6.5)
} else cat("  [skip] forecast pkg or 2SRR_AR.rda not available\n")

# =============================================================================
# FIGURE 5 — Implicit regularization (1/2): how often CV picks the grid ceiling
# =============================================================================
cat("FIGURE 5: lambda saturation...\n")
p5 <- rd(RUN_FINAL, "P5b_lambda_saturation.csv")
p5 <- p5[p5$step == "lambda", ]
p5$case <- factor(p5$case, levels = c("AR", "Factor", "FAVAR"),
                  labels = c("2SRR-AR", "2SRR-Factor", "2SRR-FAVAR"))
fig5 <- ggplot(p5, aes(factor(h), pct_at_top, fill = case)) +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_text(aes(label = sprintf("%.0f%%", pct_at_top)),
            position = position_dodge(0.8), vjust = -0.35, size = 3) +
  scale_fill_manual(values = MODEL_COLORS, name = NULL) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "Horizon (h, months)", y = "% of windows at the grid ceiling",
       subtitle = "In most windows cross-validation picks the largest available penalty, shrinking coefficients toward constancy") +
  theme_article
save_article_fig(fig5, "FIG5_lambda_saturation", width = 9, height = 5.5)

# =============================================================================
# FIGURE 6 — Implicit regularization (2/2): widening the grid hurts accuracy
# =============================================================================
cat("FIGURE 6: original vs extended grid...\n")
g_o <- rd(GRID_ORIG, "P2_tvp_3cases_rmse.csv"); g_o$grid <- "Original grid  exp(-2, 12)"
g_e <- rd(GRID_EXT,  "P2_tvp_3cases_rmse.csv"); g_e$grid <- "Extended grid  exp(-4, 18)"
g <- rbind(g_o[, c("case", "h", "ratio", "grid")],
           g_e[, c("case", "h", "ratio", "grid")])
g$case <- factor(g$case, levels = c("AR", "Factor", "FAVAR"),
                 labels = c("2SRR-AR", "2SRR-Factor", "2SRR-FAVAR"))
g$grid <- factor(g$grid, levels = c("Original grid  exp(-2, 12)",
                                    "Extended grid  exp(-4, 18)"))
fig6 <- ggplot(g, aes(factor(h), ratio, fill = grid)) +
  geom_col(position = position_dodge(0.8), width = 0.72) +
  geom_hline(yintercept = 1, linetype = 2) +
  facet_wrap(~case) +
  scale_fill_manual(values = c("Original grid  exp(-2, 12)" = "#2CA02C",
                               "Extended grid  exp(-4, 18)" = "#B0202A"),
                    name = NULL) +
  labs(x = "Horizon (h, months)", y = "RMSE / RMSE(RW)",
       subtitle = "A wider penalty grid degrades out-of-sample accuracy in 10 of 12 cells") +
  theme_article
save_article_fig(fig6, "FIG6_grid_comparison", width = 10, height = 5)

# =============================================================================
# FIGURE 7 — 2SRR-AR across inflation regimes (RMSE relative to RW)
# =============================================================================
cat("FIGURE 7: subperiods...\n")
p9 <- rd(RUN_FINAL, "P9_subperiods.csv")
p9 <- p9[p9$case == "AR", ]
p9$period <- factor(p9$period,
                    levels = c("Post-GFC/Pre-COVID", "COVID",
                               "High Inflation", "Post-Inflation"))
fig7 <- ggplot(p9, aes(factor(h), ratio_vs_RW, fill = period)) +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_hline(yintercept = 1, linetype = 2) +
  scale_fill_brewer(palette = "Set2", name = NULL) +
  labs(x = "Horizon (h, months)", y = "RMSE / RMSE(RW)",
       subtitle = "2SRR-AR against the random walk (below 1.00 = win); gains peak in the high-inflation and disinflation regimes") +
  theme_article
save_article_fig(fig7, "FIG7_subperiods_2SRR_AR", width = 9.5, height = 5.5)

# =============================================================================
# FIGURE 8 — Mincer-Zarnowitz rationality (joint p-value, selected models)
#   Green when the forecast cannot be rejected as unbiased (p > 0.10).
# =============================================================================
cat("FIGURE 8: Mincer-Zarnowitz heatmap...\n")
mz <- rd(RUN_FINAL, "P14b_MZ_PT_tests.csv")
keep_mz <- c("2SRR_AR", "2SRR_FAVAR", "Random Forest", "LASSO", "ElNET", "AR")
mz <- mz[mz$model %in% c("2SRR_AR", "2SRR_FAVAR", "RF", "LASSO", "ElNET", "AR"),
         c("model", "h", "MZ_p_joint")]
mz$model <- factor(disp(mz$model),
                   levels = rev(c("2SRR-AR", "2SRR-FAVAR", "Random Forest",
                                  "LASSO", "Elastic Net", "AR")))
mz$rational <- ifelse(mz$MZ_p_joint > 0.10, "Rational (p > 0.10)", "Rejected")
fig8 <- ggplot(mz, aes(factor(h), model, fill = rational)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.3f", MZ_p_joint)), size = 3) +
  scale_fill_manual(values = c("Rational (p > 0.10)" = "#2CA02C",
                               "Rejected" = "grey80"), name = NULL) +
  labs(x = "Horizon (h, months)", y = NULL,
       subtitle = "2SRR-AR is the only model that stays unbiased at three consecutive horizons") +
  theme_article + theme(panel.grid = element_blank())
save_article_fig(fig8, "FIG8_mincer_zarnowitz", width = 8.5, height = 4.5)

# =============================================================================
# TABLES (Word-ready .docx via flextable)
#   T1: compact accuracy summary (curated rows, not all 18 models)
#   T2: Mincer-Zarnowitz joint p-values (alternative to FIG 8)
# =============================================================================
if (has_flex) {
  library(flextable)
  set_flextable_defaults(font.family = "Calibri", font.size = 10)

  # T1 — compact accuracy summary
  cat("TABLE 1: compact accuracy summary (.docx)...\n")
  sel_t1 <- c("RW", "LASSO", "ElNET", "RF", "AR", "2SRR_AR", "2SRR_FAVAR")
  t1 <- raw1[match(sel_t1, raw1$model), c("model", "h1", "h3", "h6", "h12")]
  t1$model <- ifelse(t1$model == "RW", "Random walk", disp(t1$model))
  names(t1) <- c("Model", "h = 1", "h = 3", "h = 6", "h = 12")
  ft1 <- flextable(t1) |> colformat_double(digits = 3) |>
    bold(part = "header") |> autofit() |>
    add_header_lines("Table 1 - RMSE relative to the random walk (selected models)")
  save_as_docx(ft1, path = file.path(OUT_TAB, "TABLE1_accuracy.docx"))

  # T2 — Mincer-Zarnowitz
  cat("TABLE 2: Mincer-Zarnowitz (.docx)...\n")
  mzt <- rd(RUN_FINAL, "P14b_MZ_PT_tests.csv")
  mzt <- mzt[mzt$model %in% c("2SRR_AR", "2SRR_FAVAR", "RF", "LASSO", "ElNET", "AR"),
             c("model", "h", "MZ_p_joint")]
  mzt$model <- disp(mzt$model)
  mzt <- tidyr::pivot_wider(mzt, names_from = h, values_from = MZ_p_joint)
  names(mzt) <- c("Model", "h = 1", "h = 3", "h = 6", "h = 12")
  ft2 <- flextable(mzt) |> colformat_double(digits = 3) |>
    bold(part = "header") |> autofit() |>
    add_header_lines("Table 2 - Mincer-Zarnowitz joint p-value (values > 0.10 = unbiased)")
  save_as_docx(ft2, path = file.path(OUT_TAB, "TABLE2_mincer_zarnowitz.docx"))
} else {
  cat("  [skip tables] install 'flextable' for Word-ready .docx output\n")
}

cat("\n== 05_article_figures.R done ==\n")
cat("Figures:", normalizePath(OUT_FIG), "\n")
cat("Tables :", normalizePath(OUT_TAB), "\n")

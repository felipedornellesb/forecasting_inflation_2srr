# ---- 0. Configuration: paths to the result folders --------------------------
# 05 reads only from the two grid runs produced by 04_analysis.R:
RUN_FINAL <- "40_results/run_coulombe_grid_linspace(-2_12_15)"   # main run (default grid)
GRID_ORIG <- RUN_FINAL
GRID_EXT  <- "40_results/run_extended_grid_overexpanded_linspace_-4_18_25"
DIR_BETAS <- "30_output/betas"
DIR_DATA  <- "10_data"

OUT_FIG <- "40_results/article_figures/figures"
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)

# ---- 1. Packages ------------------------------------------------------------
need <- c("ggplot2", "dplyr", "tidyr", "scales")
for (p in need)
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cran.r-project.org")
suppressPackageStartupMessages(invisible(lapply(need, library, character.only = TRUE)))

# ---- 2. House style: palette, theme, helpers --------------------------------
# Fixed colour scheme, consistent across every figure (heatmaps excepted):
#   Realized = black | Ridge benchmark = grey | 2SRR-FAVAR = blue
#   2SRR-Factor = green | 2SRR-AR = red. Remaining models chosen for contrast.
MODEL_COLORS <- c(
  "Realized"           = "black",
  "Ridge"              = "#7F7F7F",   # benchmark: grey
  "2SRR-FAVAR"         = "#1F4E9C",   # blue
  "2SRR-Factor"        = "#2CA02C",   # green
  "2SRR-AR"            = "#D62728",   # red
  "LASSO"              = "#FF7F0E",   # orange
  "Elastic Net"        = "#17BECF",   # cyan
  "Adaptive LASSO"     = "#9467BD",   # purple
  "Adaptive ElNet"     = "#8C564B",   # brown
  "Random Forest"      = "#E377C2",   # pink
  "Bagging"            = "#BCBD22",   # olive
  "Factor"             = "#1A55A3",   # navy
  "Target Factor"      = "#7B4FA3",
  "CSR"                = "#FF1493",
  "AR"                 = "#000080",
  "AR-BIC"             = "#4B0082",
  "Ridge-Step1-AR"     = "#FFA07A",
  "Ridge-Step1-Factor" = "#A0522D",
  "Ridge-Step1-FAVAR"  = "#DEB887")
# One line width for every forecast line (realized included); distinguish series
# of the same colour family with visible dashes (never tiny dotted lines).
LINE_W <- 0.8

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
         "RF" = "Random Forest", "Ridge" = "Ridge")
  out <- m[x]; out[is.na(out)] <- x[is.na(out)]; unname(out)
}

save_article_fig <- function(plt, fname, width = 9, height = 5.5) {
#  ggsave(file.path(OUT_FIG, paste0(fname, ".pdf")), plt, width = width, height = height)
  ggsave(file.path(OUT_FIG, paste0(fname, ".png")), plt, width = width, height = height,
         dpi = 300, bg = "white")
#  cat("  [fig]", fname, ".pdf / .png\n", sep = "")
  cat("  [fig]", fname, ".png\n", sep = "")
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
ord <- c("Random walk", "Ridge", "Factor", "CSR", "Target Factor",
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
# FIGURE 2 — The value of regularized time variation: 2SRR-FAVAR vs two ridges
#   The ordinary high-dimensional ridge (117 series) vs the Step-1 static ridge
#   on the same FAVAR design vs the full time-varying 2SRR-FAVAR. RMSE vs RW.
# =============================================================================
cat("FIGURE 2: 2SRR-FAVAR vs static ridges...\n")
sel <- c("Ridge", "RidgeStep1_FAVAR", "2SRR_FAVAR")
p2 <- p1[FALSE, ]  # not used; rebuild from raw to keep raw codes
raw1 <- rd(RUN_FINAL, "P1_rmsfe_relative_rw_all_h.csv")
d2 <- raw1[raw1$model %in% sel, c("model", "h1", "h3", "h6", "h12")]
d2$model <- factor(disp(d2$model),
                   levels = c("Ridge", "Ridge-Step1-FAVAR", "2SRR-FAVAR"))
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
       subtitle = "2SRR-FAVAR against the ordinary 117-series ridge and its own static (Step-1) counterpart") +
  theme_article
save_article_fig(fig2, "FIG2_2srr_vs_ridges", width = 9, height = 5.5)

# =============================================================================
# FIGURE 3 — Where the parameters move: 2SRR-FAVAR coefficient paths (h = 6).
#   FAVAR carries many coefficients (own lags + PCA factor loadings); the five
#   most time-varying are shown. Per window, the last in-sample beta is used.
# =============================================================================
cat("FIGURE 3: 2SRR-FAVAR coefficient paths (h=6, top-5 by variability)...\n")
load(file.path(DIR_DATA, "data.rda"))            # provides `data` with $date
all_dates <- as.Date(data$date)
fig3 <- NULL
bp <- file.path(DIR_BETAS, "betas_2SRR_FAVAR.rda")
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
  sds  <- apply(mat, 2, sd, na.rm = TRUE)
  topk <- names(sort(sds, decreasing = TRUE))[1:min(5, K)]
  df3  <- data.frame(date = oos_dates, mat[, topk, drop = FALSE], check.names = FALSE)
  long3 <- tidyr::pivot_longer(df3, -date, names_to = "coef", values_to = "value")
  long3$coef <- factor(long3$coef, levels = topk)
  # Static-ridge fallback windows have no TVP betas; bridge the NA gap (a short
  # straight segment) instead of breaking. Disclose in the figure footnote.
  n_fallback <- sum(is.na(long3$value)) / dplyr::n_distinct(long3$coef)
  if (n_fallback > 0)
    cat(sprintf("  note: %d window(s) used ridge fallback (no TVP betas); line bridges them.\n",
                round(n_fallback)))
  long3 <- long3[!is.na(long3$value), ]
  fig3 <- ggplot(long3, aes(date, value, colour = coef)) +
    geom_hline(yintercept = 0, linetype = 3, colour = "grey50") +
    geom_line(linewidth = 0.8, na.rm = TRUE) +
    scale_colour_brewer(palette = "Dark2", name = NULL) +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    labs(x = "Forecast origin", y = "Estimated coefficient",
         subtitle = "Top-5 most time-varying coefficients of 2SRR-FAVAR (own lags + PCA factor loadings)") +
    theme_article
  save_article_fig(fig3, "FIG3_coef_paths_2srrFAVAR_h6", width = 9.5, height = 5.5)
} else cat("  [skip] betas_2SRR_FAVAR.rda not found\n")

# =============================================================================
# FIGURE 4 — Diebold-Mariano heatmap: 2SRR-AR (championed spec) vs each model
#   Computed fresh from the saved forecasts so the REFERENCE is 2SRR-AR (the
#   spec the analysis champions), with the sibling 2SRR specs shown as rows.
#   The DM call uses the squared-error loss differential with a HAC variance,
#   consistent with the Giacomini-White heatmaps above (Figure 9).
# =============================================================================
cat("FIGURE 4: Diebold-Mariano heatmap (reference = 2SRR-FAVAR)...\n")
fig4 <- NULL
DIR_FC <- "30_output/forecasts"
if (requireNamespace("forecast", quietly = TRUE) &&
    file.exists(file.path(DIR_FC, "2SRR_FAVAR.rda"))) {
  ld_fc <- function(f) { e <- new.env(); load(file.path(DIR_FC, f), envir = e)
                         as.matrix(get(ls(e)[1], envir = e)) }
  load(file.path(DIR_FC, "yout.rda"))            # realized cumulative target
  fc_ref <- ld_fc("2SRR_FAVAR.rda")
  bench_map <- c(
    "2SRR_AR.rda" = "2SRR-AR", "2SRR_Factor.rda" = "2SRR-Factor",
    "Ridge_from_2SRR_FAVAR.rda" = "Ridge-Step1-FAVAR", "AR.rda" = "AR",
    "AR_BIC.rda" = "AR-BIC", "RF.rda" = "Random Forest", "LASSO.rda" = "LASSO",
    "ElNET.rda" = "Elastic Net", "AdaLASSO.rda" = "Adaptive LASSO",
    "AdaElNET.rda" = "Adaptive ElNet", "Bagging.rda" = "Bagging",
    "CSR.rda" = "CSR", "Factor.rda" = "Factor", "T.Factor.rda" = "Target Factor",
    "Ridge.rda" = "Ridge")
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
  # green (negative) = 2SRR-FAVAR more accurate; red (positive) = the other model.
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
         subtitle = "Reference: 2SRR-FAVAR. Green = 2SRR-FAVAR more accurate; red = the other model. Cell = DM p-value") +
    theme_article + theme(panel.grid = element_blank())
  save_article_fig(fig4, "FIG4_dm_heatmap_2SRR_FAVAR", width = 8.5, height = 6.5)
} else cat("  [skip] forecast pkg or 2SRR_FAVAR.rda not available\n")

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
# FIGURE 7 — 2SRR-FAVAR across inflation regimes (RMSE relative to RW)
# =============================================================================
cat("FIGURE 7: subperiods...\n")
p9 <- rd(RUN_FINAL, "P9_subperiods.csv")
p9 <- p9[p9$case == "FAVAR", ]
p9$period <- factor(p9$period,
                    levels = c("Post-GFC/Pre-COVID", "COVID",
                               "High Inflation", "Post-Inflation"))
fig7 <- ggplot(p9, aes(factor(h), ratio_vs_RW, fill = period)) +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_hline(yintercept = 1, linetype = 2) +
  scale_fill_brewer(palette = "Set2", name = NULL) +
  labs(x = "Horizon (h, months)", y = "RMSE / RMSE(RW)",
       subtitle = "2SRR-FAVAR against the random walk (below 1.00 = win); it wins outside the pandemic, where the estimated factors falter") +
  theme_article
save_article_fig(fig7, "FIG7_subperiods_2SRR_FAVAR", width = 9.5, height = 5.5)

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
# Giacomini-White (2006) conditional predictive ability test (as provided).
#   x = forecasts of model 1 (reference / benchmark), y = forecasts of model 2
#   (competitor), p = realized OOS values, tau = forecast horizon, method = HAC.
# =============================================================================
suppressPackageStartupMessages({
  if (!requireNamespace("sandwich", quietly = TRUE))
    install.packages("sandwich", repos = "https://cran.r-project.org")
  library(sandwich)
})
gw.test <- function(x, y, p, T, tau,
                    method = c("HAC", "NeweyWest", "Andrews", "LumleyHeagerty"),
                    alternative = c("two.sided", "less", "greater")) {
  if (NCOL(x) > 1) stop("x is not a vector or univariate time series")
  if (tau < 1) stop("Predictive horizon must be a positive integer")
  if (length(x) != length(y)) stop("size of x and y differ")
  method      <- match.arg(method)
  alternative <- match.arg(alternative)
  l1 <- abs(x - p); l2 <- abs(y - p); dif <- l1 - l2
  q  <- length(dif); delta <- mean(dif)
  mod <- lm(dif ~ 0 + rep(1, q))
  if (tau == 1) {
    STATISTIC <- summary(mod)$coefficients[1, 3]
  } else {
    ds <- switch(method,
                 "HAC"       = sqrt(vcovHAC(mod)[1, 1]),
                 "NeweyWest" = sqrt(NeweyWest(mod, tau)[1, 1]),
                 "Andrews"   = sqrt(kernHAC(mod)[1, 1]),
                 sqrt(vcovHAC(mod)[1, 1]))
    STATISTIC <- delta / ds
  }
  PVAL <- switch(alternative,
                 "two.sided" = 2 * pnorm(-abs(STATISTIC)),
                 "less"      = pnorm(STATISTIC),
                 "greater"   = pnorm(STATISTIC, lower.tail = FALSE))
  list(statistic = STATISTIC, p.value = PVAL, method = method,
       alternative = alternative)
}

# =============================================================================
# FIGURE 9 — Giacomini-White heatmaps with the relative RMSE as the focus and
#   the GW p-value (HAC) below as significance. One panel per benchmark:
#   (a) the ordinary RIDGE (the origin of 2SRR), (b) the AR benchmark, and
#   (c) 2SRR-AR. Cell = RMSE(model)/RMSE(benchmark): below 1 (green) the model
#   beats the benchmark; above 1 (red) it does worse. Relative RMSE says who
#   won; the p-value says whether the gain is real.
# =============================================================================
cat("FIGURE 9: GW relative-RMSE heatmaps (Ridge, AR, 2SRR-FAVAR benchmarks)...\n")
DIR_FC <- "30_output/forecasts"
ld_fc  <- function(f) { e <- new.env(); load(file.path(DIR_FC, f), envir = e)
                        as.matrix(get(ls(e)[1], envir = e)) }
if (!exists("yout")) load(file.path(DIR_FC, "yout.rda"))
bigT    <- nrow(yout) + 606L                 # total sample size (tau + n_oos)
rmse_fn <- function(y, f) { ok <- complete.cases(y, f)
  if (sum(ok) < 5) NA_real_ else sqrt(mean((y[ok] - f[ok])^2)) }
hz <- c(1, 3, 6, 12)
gw_all_map <- c(
  "2SRR_AR.rda" = "2SRR-AR", "2SRR_Factor.rda" = "2SRR-Factor",
  "2SRR_FAVAR.rda" = "2SRR-FAVAR", "Ridge_from_2SRR_AR.rda" = "Ridge-Step1-AR",
  "AR.rda" = "AR", "AR_BIC.rda" = "AR-BIC", "RF.rda" = "Random Forest",
  "LASSO.rda" = "LASSO", "ElNET.rda" = "Elastic Net",
  "AdaLASSO.rda" = "Adaptive LASSO", "AdaElNET.rda" = "Adaptive ElNet",
  "Bagging.rda" = "Bagging", "CSR.rda" = "CSR", "Factor.rda" = "Factor",
  "T.Factor.rda" = "Target Factor", "Ridge.rda" = "Ridge")

gw_heatmap <- function(bench_file, bench_label, out_name) {
  if (!file.exists(file.path(DIR_FC, bench_file))) {
    cat("  [skip]", out_name, "(benchmark not found)\n"); return(invisible(NULL)) }
  fc_b <- ld_fc(bench_file)
  rows <- list()
  for (bf in names(gw_all_map)) {
    if (bf == bench_file) next                          # skip the benchmark itself
    if (!file.exists(file.path(DIR_FC, bf))) next
    M <- ld_fc(bf)
    for (h in hz) {
      if (h > ncol(M) || h > ncol(fc_b)) next
      ok <- complete.cases(yout[, h], fc_b[, h], M[, h]); if (sum(ok) < 20) next
      rel <- rmse_fn(yout[, h], M[, h]) / rmse_fn(yout[, h], fc_b[, h])
      pv  <- tryCatch(gw.test(x = fc_b[ok, h], y = M[ok, h], p = yout[ok, h],
                              T = bigT, tau = h, method = "HAC")$p.value,
                      error = function(e) NA_real_)
      rows[[length(rows) + 1]] <- data.frame(
        competitor = unname(gw_all_map[bf]), h = h,
        rel = round(as.numeric(rel), 2), p = round(as.numeric(pv), 3))
    }
  }
  df  <- do.call(rbind, rows)
  lev <- intersect(rev(unname(gw_all_map)), unique(df$competitor))
  df$competitor <- factor(df$competitor, levels = lev)
  df$stars <- ifelse(is.na(df$p), "",
              ifelse(df$p < 0.01, "***",
              ifelse(df$p < 0.05, "**",
              ifelse(df$p < 0.10, "*", ""))))
  df$lab <- sprintf("%.2f\n(%.2f)%s", df$rel, df$p, df$stars)
  g <- ggplot(df, aes(factor(h), competitor, fill = pmin(pmax(rel, 0.5), 1.5))) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = lab), size = 2.7, lineheight = 0.85) +
    scale_fill_gradient2(midpoint = 1, low = "#1A7F37", mid = "white",
                         high = "#B0202A", limits = c(0.5, 1.5),
                         oob = scales::squish,
                         name = sprintf("RMSE(model) / RMSE(%s)", bench_label)) +
    labs(x = "Horizon (h, months)", y = NULL,
         subtitle = sprintf("Below 1 (green) beats %s; GW two-sided HAC p-value below; * p<0.10, ** p<0.05, *** p<0.01", bench_label)) +
    theme_article + theme(panel.grid = element_blank())
  save_article_fig(g, out_name, width = 8.8, height = 7)
}
gw_heatmap("Ridge.rda",    "Ridge",    "FIG9_gw_relrmse_vs_ridge")
gw_heatmap("AR.rda",       "AR",       "FIG9b_gw_relrmse_vs_AR")
gw_heatmap("2SRR_FAVAR.rda", "2SRR-FAVAR", "FIG9c_gw_relrmse_vs_2SRR_FAVAR")

# =============================================================================
# FIGURE 10 — CSFE: cumulative squared forecast-error difference vs a benchmark.
#   D_t = sum_{s<=t} ( e_bench,s^2 - e_2SRR,s^2 ). Rising = the 2SRR spec beats
#   the benchmark; a fall flags where it underperforms (e.g., the pandemic).
#   Produced against TWO benchmarks: the random walk and the ordinary ridge.
#   All three TVP specs are shown (AR red, Factor green, FAVAR blue), h = 1.
# =============================================================================
cat("FIGURE 10: CSFE (vs random walk and vs ridge)...\n")
if (!exists("rw"))   load(file.path(DIR_FC, "rw.rda"))
if (!exists("data")) load(file.path(DIR_DATA, "data.rda"))
n_oos     <- nrow(yout)
oos_dates <- tail(as.Date(data$date), n_oos)
hC     <- 1
fc_tvp <- list("2SRR-AR"     = ld_fc("2SRR_AR.rda"),
               "2SRR-Factor" = ld_fc("2SRR_Factor.rda"),
               "2SRR-FAVAR"  = ld_fc("2SRR_FAVAR.rda"))
make_csfe <- function(bench_mat, bench_label, out_name) {
  e_b <- (yout[, hC] - bench_mat[, hC])^2
  df  <- data.frame(date = oos_dates,
                    ar  = cumsum(e_b - (yout[, hC] - fc_tvp[["2SRR-AR"]][, hC])^2),
                    fac = cumsum(e_b - (yout[, hC] - fc_tvp[["2SRR-Factor"]][, hC])^2),
                    fav = cumsum(e_b - (yout[, hC] - fc_tvp[["2SRR-FAVAR"]][, hC])^2))
  names(df) <- c("date", "2SRR-AR", "2SRR-Factor", "2SRR-FAVAR")
  long <- tidyr::pivot_longer(df, -date, names_to = "Model", values_to = "CSFE")
  long$Model <- factor(long$Model, levels = c("2SRR-AR", "2SRR-Factor", "2SRR-FAVAR"))
  g <- ggplot(long, aes(date, CSFE, colour = Model)) +
    annotate("rect", xmin = as.Date("2020-02-01"), xmax = as.Date("2020-09-01"),
             ymin = -Inf, ymax = Inf, alpha = 0.13, fill = "grey40") +
    annotate("rect", xmin = as.Date("2021-04-01"), xmax = as.Date("2023-01-01"),
             ymin = -Inf, ymax = Inf, alpha = 0.10, fill = "#F4A582") +
    geom_hline(yintercept = 0, linetype = 3, colour = "grey50") +
    geom_line(linewidth = LINE_W, na.rm = TRUE) +
    scale_colour_manual(values = c("2SRR-AR"     = MODEL_COLORS[["2SRR-AR"]],
                                   "2SRR-Factor" = MODEL_COLORS[["2SRR-Factor"]],
                                   "2SRR-FAVAR"  = MODEL_COLORS[["2SRR-FAVAR"]]),
                        name = NULL) +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    labs(x = "Forecast origin",
         y = sprintf("Cumulative squared-error difference vs. %s", bench_label),
         subtitle = sprintf("Rising = beats %s; a fall flags underperformance (grey: COVID 2020; orange: 2021-22 surge), h = 1", bench_label)) +
    theme_article
  save_article_fig(g, out_name, width = 9.5, height = 5.5)
}
make_csfe(rw,                 "the random walk",     "FIG10_csfe_vs_rw")
make_csfe(ld_fc("Ridge.rda"), "the ordinary ridge",  "FIG10b_csfe_vs_ridge")

# =============================================================================
# FIGURE 11 — Forecasts vs realized, per horizon: 2SRR-FAVAR, the Ridge
#   benchmark, and the top-2 Medeiros models (excluding Ridge) ranked by RMSE at
#   that horizon. Equal line width for all; realized = black, 2SRR-FAVAR = blue,
#   Ridge = grey, the Medeiros pair in visible long-dashes. Forecasts are shown
#   as a monthly rate (cumulative forecast / h) so they are comparable.
# =============================================================================
cat("FIGURE 11: forecasts vs realized (2SRR-FAVAR, Ridge, top-2 Medeiros)...\n")
if (!exists("yout")) load(file.path(DIR_FC, "yout.rda"))
if (!exists("data")) load(file.path(DIR_DATA, "data.rda"))
n_oos     <- nrow(yout)
oos_dates <- tail(as.Date(data$date), n_oos)
y_monthly <- tail(data$CPIAUCSL, n_oos)
med_files <- c("LASSO.rda", "ElNET.rda", "AdaLASSO.rda", "AdaElNET.rda", "RF.rda",
               "Bagging.rda", "Factor.rda", "T.Factor.rda", "CSR.rda", "AR.rda",
               "AR_BIC.rda")
med_disp  <- c("LASSO.rda" = "LASSO", "ElNET.rda" = "Elastic Net",
               "AdaLASSO.rda" = "Adaptive LASSO", "AdaElNET.rda" = "Adaptive ElNet",
               "RF.rda" = "Random Forest", "Bagging.rda" = "Bagging",
               "Factor.rda" = "Factor", "T.Factor.rda" = "Target Factor",
               "CSR.rda" = "CSR", "AR.rda" = "AR", "AR_BIC.rda" = "AR-BIC")
fc_fav11 <- ld_fc("2SRR_FAVAR.rda")
fc_rdg11 <- ld_fc("Ridge.rda")
ov_rows <- list()
for (h in c(1, 3, 6, 12)) {
  rmses <- sapply(med_files, function(f) {
    if (!file.exists(file.path(DIR_FC, f))) return(NA_real_)
    M <- ld_fc(f); if (h > ncol(M)) return(NA_real_); rmse_fn(yout[, h], M[, h]) })
  best2 <- names(sort(rmses))[1:2]
  series <- list("Realized"           = y_monthly,
                 "2SRR-FAVAR"          = fc_fav11[, h] / h,
                 "Ridge"  = fc_rdg11[, h] / h)
  for (f in best2) series[[unname(med_disp[f])]] <- ld_fc(f)[, h] / h
  for (nm in names(series))
    ov_rows[[length(ov_rows) + 1]] <- data.frame(
      date = oos_dates, h = h, Series = nm, value = as.numeric(series[[nm]]))
}
ov <- do.call(rbind, ov_rows)
ov$h <- factor(ov$h, levels = c(1, 3, 6, 12), labels = paste0("h = ", c(1, 3, 6, 12)))
key <- c("Realized", "2SRR-FAVAR", "Ridge")
ov$Series <- factor(ov$Series, levels = c(key, setdiff(unique(ov$Series), key)))
ov_cols <- MODEL_COLORS[levels(ov$Series)]; names(ov_cols) <- levels(ov$Series)
ov_cols[is.na(ov_cols)] <- "#444444"
ov_lty  <- setNames(rep("longdash", nlevels(ov$Series)), levels(ov$Series))
ov_lty[key] <- "solid"
fig11 <- ggplot(ov, aes(date, value, colour = Series, linetype = Series)) +
  geom_hline(yintercept = 0, linetype = 3, colour = "grey60") +
  geom_line(linewidth = LINE_W, na.rm = TRUE) +
  facet_wrap(~h, ncol = 2, scales = "free_y") +
  scale_colour_manual(values = ov_cols, name = NULL) +
  scale_linetype_manual(values = ov_lty, name = NULL) +
  scale_x_date(date_breaks = "3 years", date_labels = "%Y") +
  labs(x = NULL, y = "Monthly inflation (%)",
       subtitle = "Realized (black), 2SRR-FAVAR (blue), Ridge benchmark (grey) and the top-2 Medeiros models per horizon (long-dash)") +
  theme_article
save_article_fig(fig11, "FIG11_forecasts_vs_realized", width = 11, height = 7.5)

# Numeric tables are not generated here. The figures above are the article
# outputs; the underlying numbers come from the CSV tables written by
# 04_analysis.R (e.g. P1_rmsfe_relative_rw_all_h.csv, P14b_MZ_PT_tests.csv).

cat("\n== 05_article_figures.R done ==\n")
cat("Figures:", normalizePath(OUT_FIG), "\n")

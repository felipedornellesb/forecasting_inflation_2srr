# =============================================================================
# 05b_article_figures_AR.R
#
# Final AR-benchmark figures + numbers for the Data Analysis (AR-focused article).
# The autoregression (AR) of Medeiros et al. (2021) is the benchmark, and
# 2SRR-AR is the championed specification, so the three TVP specs can be ranked
# and the under-performance of the data-rich 2SRR-FAVAR made explicit.
#
# Produces:
#   FIG4b_dm_heatmap_2SRR_AR  — Diebold-Mariano, reference = 2SRR-AR.
#   FIG10c_csfe_vs_AR         — CSFE of the three TVP specs against AR (h = 1).
# Reuses (from 05_article_figures.R, not redrawn here):
#   FIG9b_gw_relrmse_vs_AR    — Giacomini-White relative-RMSE heatmap vs AR.
#
# Console report (for the prose): RMSE ratios vs AR, GW HAC p-values vs AR,
# DM p-values (2SRR-AR vs each competitor) and final CSFE values.
# =============================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(scales); library(sandwich)
})
OUT_FIG  <- "40_results/article_figures/figures"
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)
DIR_FC   <- "30_output/forecasts"
DIR_DATA <- "10_data"
LINE_W   <- 0.8
MODEL_COLORS <- c("2SRR-AR" = "#D62728", "2SRR-Factor" = "#2CA02C", "2SRR-FAVAR" = "#1F4E9C")
theme_article <- theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        plot.subtitle = element_text(size = 9, colour = "grey30"),
        panel.grid.minor = element_blank())
save_article_fig <- function(plt, fname, width = 9, height = 5.5) {
  ggsave(file.path(OUT_FIG, paste0(fname, ".png")), plt, width = width, height = height,
         dpi = 150, bg = "white")
  ggsave(file.path(OUT_FIG, paste0(fname, ".pdf")), plt, width = width, height = height)
  cat("  [fig]", fname, "\n")
}
ld_fc <- function(f) { e <- new.env(); load(file.path(DIR_FC, f), envir = e)
                       as.matrix(get(ls(e)[1], envir = e)) }
load(file.path(DIR_FC, "yout.rda"))
load(file.path(DIR_FC, "rw.rda"))
load(file.path(DIR_DATA, "data.rda"))
hz      <- c(1, 3, 6, 12)
bigT    <- nrow(yout) + 606L
rmse_fn <- function(y, f) { ok <- complete.cases(y, f)
  if (sum(ok) < 5) NA_real_ else sqrt(mean((y[ok] - f[ok])^2)) }

# Giacomini-White (2006) conditional predictive ability test (professor's HAC variant).
gw.test <- function(x, y, p, T, tau,
                    method = c("HAC", "NeweyWest", "Andrews", "LumleyHeagerty"),
                    alternative = c("two.sided", "less", "greater")) {
  method <- match.arg(method); alternative <- match.arg(alternative)
  l1 <- abs(x - p); l2 <- abs(y - p); dif <- l1 - l2
  q  <- length(dif); delta <- mean(dif); mod <- lm(dif ~ 0 + rep(1, q))
  if (tau == 1) STATISTIC <- summary(mod)$coefficients[1, 3]
  else {
    ds <- switch(method, "HAC" = sqrt(vcovHAC(mod)[1, 1]),
                 "NeweyWest" = sqrt(NeweyWest(mod, tau)[1, 1]),
                 "Andrews" = sqrt(kernHAC(mod)[1, 1]), sqrt(vcovHAC(mod)[1, 1]))
    STATISTIC <- delta / ds
  }
  PVAL <- switch(alternative, "two.sided" = 2 * pnorm(-abs(STATISTIC)),
                 "less" = pnorm(STATISTIC), "greater" = pnorm(STATISTIC, lower.tail = FALSE))
  list(statistic = STATISTIC, p.value = PVAL)
}

AR     <- ld_fc("AR.rda")
fc_tvp <- list("2SRR-AR"     = ld_fc("2SRR_AR.rda"),
               "2SRR-Factor" = ld_fc("2SRR_Factor.rda"),
               "2SRR-FAVAR"  = ld_fc("2SRR_FAVAR.rda"))

# ---------------------------------------------------------------------------- #
# REPORT 1 — RMSE(model)/RMSE(AR) and GW HAC two-sided p-value vs AR
# ---------------------------------------------------------------------------- #
cat("\n=== RMSE relative to AR  and  GW (HAC) p-value vs AR ===\n")
cat(sprintf("%-12s | %s\n", "spec",
            paste(sprintf("h=%-2d rel (gwp)", hz), collapse = "  ")))
for (nm in names(fc_tvp)) {
  M <- fc_tvp[[nm]]; line <- sprintf("%-12s |", nm)
  for (h in hz) {
    ok  <- complete.cases(yout[, h], AR[, h], M[, h])
    rel <- rmse_fn(yout[, h], M[, h]) / rmse_fn(yout[, h], AR[, h])
    gwp <- tryCatch(gw.test(AR[ok, h], M[ok, h], yout[ok, h], bigT, h, "HAC")$p.value,
                    error = function(e) NA_real_)
    line <- paste0(line, sprintf("  %.3f (%.3f)", rel, gwp))
  }
  cat(line, "\n")
}

# ---------------------------------------------------------------------------- #
# FIGURE 4b — Diebold-Mariano heatmap, reference = 2SRR-AR
# ---------------------------------------------------------------------------- #
cat("\nFIGURE 4b: Diebold-Mariano heatmap (reference = 2SRR-AR)...\n")
fc_ref <- fc_tvp[["2SRR-AR"]]
bench_map <- c(
  "2SRR_FAVAR.rda" = "2SRR-FAVAR", "2SRR_Factor.rda" = "2SRR-Factor",
  "Ridge_from_2SRR_AR.rda" = "Ridge-Step1-AR", "AR.rda" = "AR",
  "AR_BIC.rda" = "AR-BIC", "RF.rda" = "Random Forest", "LASSO.rda" = "LASSO",
  "ElNET.rda" = "Elastic Net", "AdaLASSO.rda" = "Adaptive LASSO",
  "AdaElNET.rda" = "Adaptive ElNet", "Bagging.rda" = "Bagging",
  "CSR.rda" = "CSR", "Factor.rda" = "Factor", "T.Factor.rda" = "Target Factor",
  "Ridge.rda" = "Ridge")
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
      ref_wins = isTRUE(rmse_fn(yout[, h], fc_ref[, h]) < rmse_fn(yout[, h], M[, h])))
  }
}
dm <- do.call(rbind, dm_rows)
dm$signed_p  <- ifelse(dm$ref_wins, -dm$DM_p, dm$DM_p)
dm$benchmark <- factor(dm$benchmark, levels = rev(unname(bench_map)))
fig4b <- ggplot(dm, aes(factor(h), benchmark, fill = pmin(pmax(signed_p, -0.5), 0.5))) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", DM_p)), size = 3) +
  scale_fill_gradient2(midpoint = 0, low = "#1A7F37", mid = "white",
                       high = "#B0202A", limits = c(-0.5, 0.5), name = "Signed p-value") +
  labs(x = "Horizon (h, months)", y = NULL,
       subtitle = "Reference: 2SRR-AR. Green = 2SRR-AR more accurate; red = the other model. Cell = DM p-value") +
  theme_article + theme(panel.grid = element_blank())
save_article_fig(fig4b, "FIG4b_dm_heatmap_2SRR_AR", width = 8.5, height = 6.5)

cat("\n=== DM two-sided p-value, 2SRR-AR vs each TVP sibling / key benchmark ===\n")
for (lab in c("AR", "2SRR-Factor", "2SRR-FAVAR", "Random Forest", "LASSO", "Elastic Net")) {
  rr <- dm[dm$benchmark == lab, ]
  if (nrow(rr) == 0) next
  rr <- rr[order(rr$h), ]
  cat(sprintf("  vs %-14s : %s\n", lab,
              paste(sprintf("h=%d p=%.3f%s", rr$h, rr$DM_p,
                            ifelse(rr$ref_wins, "(AR+)", "(rival+)")), collapse = "  ")))
}

# ---------------------------------------------------------------------------- #
# FIGURE 10c — CSFE of the three TVP specs against AR (h = 1)
# ---------------------------------------------------------------------------- #
cat("\nFIGURE 10c: CSFE vs the autoregression (3 TVP specs)...\n")
n_oos     <- nrow(yout)
oos_dates <- tail(as.Date(data$date), n_oos)
hC <- 1
e_b <- (yout[, hC] - AR[, hC])^2
df  <- data.frame(date = oos_dates,
                  a = cumsum(e_b - (yout[, hC] - fc_tvp[["2SRR-AR"]][, hC])^2),
                  f = cumsum(e_b - (yout[, hC] - fc_tvp[["2SRR-Factor"]][, hC])^2),
                  v = cumsum(e_b - (yout[, hC] - fc_tvp[["2SRR-FAVAR"]][, hC])^2))
names(df) <- c("date", "2SRR-AR", "2SRR-Factor", "2SRR-FAVAR")
long <- tidyr::pivot_longer(df, -date, names_to = "Model", values_to = "CSFE")
long$Model <- factor(long$Model, levels = c("2SRR-AR", "2SRR-Factor", "2SRR-FAVAR"))
fig10c <- ggplot(long, aes(date, CSFE, colour = Model)) +
  annotate("rect", xmin = as.Date("2020-02-01"), xmax = as.Date("2020-09-01"),
           ymin = -Inf, ymax = Inf, alpha = 0.13, fill = "grey40") +
  annotate("rect", xmin = as.Date("2021-04-01"), xmax = as.Date("2023-01-01"),
           ymin = -Inf, ymax = Inf, alpha = 0.10, fill = "#F4A582") +
  geom_hline(yintercept = 0, linetype = 3, colour = "grey50") +
  geom_line(linewidth = LINE_W, na.rm = TRUE) +
  scale_colour_manual(values = MODEL_COLORS, name = NULL) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(x = "Forecast origin",
       y = "Cumulative squared-error difference vs. the autoregression",
       subtitle = "Rising = beats AR; a fall flags underperformance (grey: COVID 2020; orange: 2021-22 surge), h = 1") +
  theme_article
save_article_fig(fig10c, "FIG10c_csfe_vs_AR", width = 9.5, height = 5.5)

cat("\n=== Final cumulative CSFE vs AR at end of sample (positive = beats AR) ===\n")
cat(sprintf("  2SRR-AR=%.4f  2SRR-Factor=%.4f  2SRR-FAVAR=%.4f\n",
            tail(df[["2SRR-AR"]], 1), tail(df[["2SRR-Factor"]], 1), tail(df[["2SRR-FAVAR"]], 1)))

# ---------------------------------------------------------------------------- #
# FIGURE 2 (AR) — value of time variation: 2SRR-AR vs the two ridges (RMSE/RW)
# ---------------------------------------------------------------------------- #
cat("\nFIGURE 2 (AR): 2SRR-AR vs static ridges...\n")
RUN_FINAL <- "40_results/run_coulombe_grid_linspace(-2_12_15)"
rd <- function(dir, f) read.csv(file.path(dir, "tables", f), stringsAsFactors = FALSE)
raw1 <- rd(RUN_FINAL, "P1_rmsfe_relative_rw_all_h.csv")
disp2 <- c("Ridge" = "Ridge", "RidgeStep1_AR" = "Ridge-Step1-AR", "2SRR_AR" = "2SRR-AR")
sel  <- c("Ridge", "RidgeStep1_AR", "2SRR_AR")
d2 <- raw1[match(sel, raw1$model), c("model", "h1", "h3", "h6", "h12")]
d2$model <- factor(unname(disp2[d2$model]),
                   levels = c("Ridge", "Ridge-Step1-AR", "2SRR-AR"))
long2 <- tidyr::pivot_longer(d2, -model, names_to = "h", values_to = "ratio")
long2$h <- factor(long2$h, levels = c("h1", "h3", "h6", "h12"), labels = c("1", "3", "6", "12"))
cols2 <- c("Ridge" = "#7F7F7F", "Ridge-Step1-AR" = "#FFA07A", "2SRR-AR" = "#D62728")
fig2 <- ggplot(long2, aes(h, ratio, fill = model)) +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_text(aes(label = sprintf("%.2f", ratio)),
            position = position_dodge(0.8), vjust = -0.35, size = 3) +
  scale_fill_manual(values = cols2, name = NULL) +
  labs(x = "Horizon (h, months)", y = "RMSE / RMSE(RW)",
       subtitle = "2SRR-AR against the ordinary 117-series ridge and its own static (Step-1) counterpart") +
  theme_article
save_article_fig(fig2, "FIG2_2srr_AR_vs_ridges", width = 9, height = 5.5)
cat("  RMSE/RW ratios (Ridge / Ridge-Step1-AR / 2SRR-AR):\n")
print(d2, row.names = FALSE)
cat("done\n")

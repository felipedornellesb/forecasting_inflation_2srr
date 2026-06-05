# ==============================================================================
# 06_robustness.R  --  robustness analyses (advisor's request). Re-evaluates the
# forecasts from 03/03b; no re-forecasting.
#
#   (1) Pre-pandemic sub-sample (target before Jan 2020, 246 of 312 windows):
#       pre-pandemic analogues of Figure 1 (rel-RMSE vs AR + Giacomini-White) and
#       Figure 2 (Diebold-Mariano vs 2SRR-ARDI), with the underlying tables.
#   (2) Grid saturation: share of windows where CV picks the upper bound of
#       Coulombe's default grid (high share = data favour near-constant params).
#
# Inputs : forecasts/{yout,rw,<model>,2SRR_<case>}.rda, betas/betas_2SRR_<case>.rda
# Outputs: 40_results/pre_pandemic/*.csv and article_figures/figures/*prepandemic*
# Run from the project root.
# ==============================================================================
cat("== 06_robustness.R ==\n\n")
suppressPackageStartupMessages({ library(ggplot2); library(scales); library(forecast) })

DIR_FC   <- "30_output/forecasts"
DIR_BETA <- "30_output/betas"
DIR_DATA <- "10_data"
OUT_FIG  <- "40_results/article_figures/figures"
OUT_TAB  <- "40_results/pre_pandemic"
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)

# ---- Helpers ----------------------------------------------------------------
ld <- function(p) { e <- new.env(); load(p, envir = e); get(ls(e)[1], envir = e) }
M  <- function(name) as.matrix(ld(file.path(DIR_FC, paste0(name, ".rda"))))

# Display names: the data-rich case (internal id "FAVAR") is the ARDI
# specification of Goulet Coulombe (2025).
disp <- function(x) {
  map <- c("2SRR_AR" = "2SRR-AR", "2SRR_Factor" = "2SRR-Factor",
           "2SRR_FAVAR" = "2SRR-ARDI", "AR_BIC" = "AR-BIC",
           "T.Factor" = "Target Factor", "ElNET" = "Elastic Net",
           "AdaLASSO" = "Adaptive LASSO", "AdaElNET" = "Adaptive ElNet",
           "RF" = "Random Forest")
  out <- map[x]; out[is.na(out)] <- x[is.na(out)]; unname(out)
}
mod_levels <- rev(c("Ridge", "2SRR-ARDI", "2SRR-Factor", "AR", "AR-BIC",
                    "Bagging", "Factor", "Target Factor", "CSR",
                    "Adaptive ElNet", "Adaptive LASSO", "Elastic Net", "LASSO",
                    "Random Forest", "2SRR-AR"))
theme_article <- theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        plot.subtitle = element_text(size = 9, colour = "grey30"),
        panel.grid.minor = element_blank())
save_png <- function(plt, fname, w, h) {
  ggsave(file.path(OUT_FIG, paste0(fname, ".png")), plt,
         width = w, height = h, dpi = 150, bg = "white")
  cat("  [fig]", fname, "\n")
}

# Giacomini-White (2006) conditional predictive-ability test. Test function
# h_t = (1, d_{t-1})'; statistic ~ chi^2(q = 2) with a Bartlett HAC covariance of
# bandwidth h - 1. Identical to the one used in 04/05.
gw_cond <- function(loss_a, loss_b, h = 1) {
  d <- loss_a - loss_b; d <- d[is.finite(d)]; n <- length(d)
  if (n < 20) return(NA_real_)
  z <- cbind(1, c(NA, d[-n])) * d
  z <- z[complete.cases(z), , drop = FALSE]; T_n <- nrow(z)
  z_bar <- colMeans(z); zc <- z - rep(z_bar, each = T_n)
  Omega <- crossprod(zc) / T_n; bw <- max(0, h - 1)
  if (bw > 0) for (l in 1:bw) {
    w_l <- 1 - l / (bw + 1)
    G <- crossprod(zc[-(1:l), , drop = FALSE],
                   zc[-((T_n - l + 1):T_n), , drop = FALSE]) / T_n
    Omega <- Omega + w_l * (G + t(G))
  }
  stat <- tryCatch(as.numeric(T_n * t(z_bar) %*% solve(Omega) %*% z_bar),
                   error = function(e) NA_real_)
  if (is.na(stat)) NA_real_ else 1 - pchisq(stat, df = ncol(z))
}

# ---- Data and the pre-pandemic window mask ----------------------------------
yout <- M("yout"); rw <- M("rw")
data <- ld(file.path(DIR_DATA, "data.rda")); dates <- as.Date(data$date)
n_obs <- length(dates); n_oos <- nrow(yout); tau <- n_obs - n_oos
target_date <- dates[(tau + 1):n_obs]               # realized target of each window
pre <- which(target_date < as.Date("2020-01-01"))   # pre-pandemic windows
horizons <- c(1, 3, 6, 12)
cat(sprintf("OOS windows: %d total | %d pre-pandemic (target < Jan 2020)\n\n",
            n_oos, length(pre)))

rmse_sub <- function(model_mat, h, idx) {
  ok <- idx[complete.cases(yout[idx, h], model_mat[idx, h])]
  sqrt(mean((yout[ok, h] - model_mat[ok, h])^2))
}
AR_fc <- M("AR")

# ==============================================================================
# (1a) TABLE A + plot data: relative RMSE vs AR and conditional GW p-value
#       (pre-pandemic and, for reference, full sample)
# ==============================================================================
specs <- c("2SRR_FAVAR", "2SRR_AR", "2SRR_Factor", "Ridge", "LASSO", "ElNET", "RF")
tabA <- list()
for (m in specs) {
  Fm <- M(m)
  for (h in horizons) {
    okp <- pre[complete.cases(yout[pre, h], Fm[pre, h], AR_fc[pre, h])]
    okf <- which(complete.cases(yout[, h], Fm[, h], AR_fc[, h]))
    tabA[[length(tabA) + 1]] <- data.frame(
      model = disp(m), h = h,
      relAR_pre  = round(rmse_sub(Fm, h, pre)  / rmse_sub(AR_fc, h, pre), 3),
      relAR_full = round(rmse_sub(Fm, h, 1:n_oos) / rmse_sub(AR_fc, h, 1:n_oos), 3),
      GWp_pre  = round(gw_cond((yout[okp, h] - Fm[okp, h])^2,
                               (yout[okp, h] - AR_fc[okp, h])^2, h), 3),
      GWp_full = round(gw_cond((yout[okf, h] - Fm[okf, h])^2,
                               (yout[okf, h] - AR_fc[okf, h])^2, h), 3))
  }
}
tabA <- do.call(rbind, tabA)
write.csv(tabA, file.path(OUT_TAB, "TabA_prepandemic_relAR_GW.csv"), row.names = FALSE)
cat("  [tbl] TabA_prepandemic_relAR_GW.csv\n")

# ==============================================================================
# (1b) FIGURE 1 (pre-pandemic): relative RMSE vs AR + conditional GW p-value
# ==============================================================================
all_models <- c("Ridge", "LASSO", "ElNET", "AdaLASSO", "AdaElNET", "RF", "Bagging",
                "Factor", "T.Factor", "CSR", "AR", "AR_BIC",
                "2SRR_AR", "2SRR_Factor", "2SRR_FAVAR")
gw_rows <- list()
for (m in setdiff(all_models, "AR")) {
  fp <- file.path(DIR_FC, paste0(m, ".rda")); if (!file.exists(fp)) next
  Fm <- M(m)
  for (h in horizons) {
    ok <- pre[complete.cases(yout[pre, h], AR_fc[pre, h], Fm[pre, h])]
    if (length(ok) < 20) next
    rel <- rmse_sub(Fm, h, pre) / rmse_sub(AR_fc, h, pre)
    pv  <- gw_cond((yout[ok, h] - Fm[ok, h])^2, (yout[ok, h] - AR_fc[ok, h])^2, h)
    gw_rows[[length(gw_rows) + 1]] <- data.frame(model = disp(m), h = h,
                                                 rel = round(rel, 3), p = round(pv, 3))
  }
}
gw_df <- do.call(rbind, gw_rows)
gw_df$stars <- ifelse(is.na(gw_df$p), "",
                ifelse(gw_df$p < 0.01, "***",
                ifelse(gw_df$p < 0.05, "**",
                ifelse(gw_df$p < 0.10, "*", ""))))
gw_df$lab   <- sprintf("%.2f\n(%.2f)%s", gw_df$rel, gw_df$p, gw_df$stars)
gw_df$model <- factor(gw_df$model, levels = intersect(mod_levels, unique(gw_df$model)))
fig1 <- ggplot(gw_df, aes(factor(h), model, fill = pmin(pmax(rel, 0.5), 1.5))) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = lab), size = 2.7, lineheight = 0.85) +
  scale_fill_gradient2(midpoint = 1, low = "#1A7F37", mid = "white", high = "#B0202A",
                       limits = c(0.5, 1.5), oob = scales::squish,
                       name = "RMSE(model) / RMSE(AR)") +
  labs(x = "Horizon (h, months)", y = NULL,
       subtitle = paste0("Pre-pandemic sub-sample (target < Jan 2020, n = ", length(pre),
                         ").\nCell: RMSE relative to AR (green < 1) with Giacomini-White ",
                         "conditional p-value.  * p<.10  ** p<.05  *** p<.01")) +
  theme_article + theme(panel.grid = element_blank())
save_png(fig1, "FIG1_gw_relrmse_vs_AR_prepandemic", 8.8, 7)

# ==============================================================================
# (1c) TABLE B + FIGURE 2 (pre-pandemic): Diebold-Mariano, reference 2SRR-ARDI
# ==============================================================================
ref <- M("2SRR_FAVAR")
bms <- c("AR", "Ridge", "LASSO", "ElNET", "RF", "Factor", "Bagging", "CSR")
dm_rows <- list()
for (m in bms) {
  Fm <- M(m)
  for (h in horizons) {
    ok <- pre[complete.cases(yout[pre, h], ref[pre, h], Fm[pre, h])]
    if (length(ok) < 20) next
    p <- tryCatch(forecast::dm.test(yout[ok, h] - ref[ok, h], yout[ok, h] - Fm[ok, h],
                                    alternative = "two.sided", h = h, power = 2)$p.value,
                  error = function(e) NA_real_)
    dm_rows[[length(dm_rows) + 1]] <- data.frame(
      model = disp(m), h = h, p = round(p, 4),
      ref_wins = isTRUE(rmse_sub(ref, h, pre) < rmse_sub(Fm, h, pre)))
  }
}
dm_df <- do.call(rbind, dm_rows)
write.csv(dm_df[, c("h", "model", "p", "ref_wins")],
          file.path(OUT_TAB, "TabB_prepandemic_DM.csv"), row.names = FALSE)
cat("  [tbl] TabB_prepandemic_DM.csv\n")

dm_df$strength <- ifelse(dm_df$ref_wins, -(1 - dm_df$p), (1 - dm_df$p))
dm_df$model    <- factor(dm_df$model, levels = intersect(mod_levels, unique(dm_df$model)))
fig2 <- ggplot(dm_df, aes(factor(h), model, fill = strength)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", p)), size = 3) +
  scale_fill_gradient2(midpoint = 0, low = "#1A7F37", mid = "white", high = "#B0202A",
                       limits = c(-1, 1), name = "Signed (1 - p)") +
  labs(x = "Horizon (h, months)", y = NULL,
       subtitle = paste0("Pre-pandemic sub-sample (target < Jan 2020, n = ", length(pre),
                         ").\nReference: 2SRR-ARDI (green = more accurate; deeper colour = ",
                         "lower p).  Cell: Diebold-Mariano p-value.")) +
  theme_article + theme(panel.grid = element_blank())
save_png(fig2, "FIG2_dm_heatmap_2SRR_ARDI_prepandemic", 8.5, 6.5)

# ==============================================================================
# (2) CROSS-VALIDATION GRID SATURATION
# ==============================================================================
grid    <- exp(seq(-2, 12, length.out = 15))     # Coulombe default grid
g_max   <- max(grid); g_min <- min(grid)
grid_rows <- list()
cat("\n  Grid saturation (% of windows at the grid bound):\n")
cat(sprintf("    %-12s %5s %9s %9s\n", "case", "h", "at_MAX", "at_min"))
for (cs in c("AR", "Factor", "FAVAR")) {
  bp <- file.path(DIR_BETA, paste0("betas_2SRR_", cs, ".rda")); if (!file.exists(bp)) next
  betas <- ld(bp)
  for (h in horizons) {
    lam <- vapply(betas[[paste0("h", h)]],
                  function(w) if (is.null(w$lambda)) NA_real_ else w$lambda, numeric(1))
    pmax_at <- mean(abs(lam - g_max) < 1e-6, na.rm = TRUE)
    pmin_at <- mean(abs(lam - g_min) < 1e-6, na.rm = TRUE)
    grid_rows[[length(grid_rows) + 1]] <- data.frame(
      case = disp(paste0("2SRR_", cs)), h = h,
      pct_at_max = round(100 * pmax_at, 1), pct_at_min = round(100 * pmin_at, 1))
    cat(sprintf("    %-12s %5d %8.1f%% %8.1f%%\n",
                disp(paste0("2SRR_", cs)), h, 100 * pmax_at, 100 * pmin_at))
  }
}
grid_df <- do.call(rbind, grid_rows)
write.csv(grid_df, file.path(OUT_TAB, "grid_saturation.csv"), row.names = FALSE)
cat(sprintf("\n  Overall share at grid maximum: %.1f%%\n",
            mean(grid_df$pct_at_max)))
cat("  [tbl] grid_saturation.csv\n")

cat("\n== 06_robustness.R done ==\n")

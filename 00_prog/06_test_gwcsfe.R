# Standalone test of the FIG9 (GW relative-RMSE heatmaps) + FIG10 (CSFE) code
# from 05_article_figures.R, run directly against 30_output/forecasts.
#   - GW heatmaps for three benchmarks: Ridge, AR, 2SRR-AR.
#   - CSFE vs the random walk AND vs the ordinary ridge, all three TVP specs.
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(scales); library(sandwich)
})
OUT_FIG  <- "40_results/_test_gwcsfe"
dir.create(OUT_FIG, showWarnings = FALSE, recursive = TRUE)
DIR_FC   <- "30_output/forecasts"
DIR_DATA <- "10_data"
LINE_W   <- 0.8
MODEL_COLORS <- c("2SRR-AR" = "#D62728", "2SRR-Factor" = "#2CA02C", "2SRR-FAVAR" = "#1F4E9C")
theme_article <- theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        plot.subtitle = element_text(size = 9, colour = "grey30"),
        panel.grid.minor = element_blank())
save_article_fig <- function(plt, fname, width = 9, height = 5.5) {
  ggsave(file.path(OUT_FIG, paste0(fname, ".png")), plt, width = width,
         height = height, dpi = 150, bg = "white")
  cat("  [fig]", fname, "\n")
}

# Giacomini-White (2006) conditional predictive ability test (as provided).
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
  list(statistic = STATISTIC, p.value = PVAL, method = method, alternative = alternative)
}

ld_fc <- function(f) { e <- new.env(); load(file.path(DIR_FC, f), envir = e)
                       as.matrix(get(ls(e)[1], envir = e)) }
load(file.path(DIR_FC, "yout.rda"))
bigT    <- nrow(yout) + 606L
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
    cat("  [skip]", out_name, "\n"); return(invisible(NULL)) }
  fc_b <- ld_fc(bench_file); rows <- list()
  for (bf in names(gw_all_map)) {
    if (bf == bench_file) next
    if (!file.exists(file.path(DIR_FC, bf))) next
    M <- ld_fc(bf)
    for (h in hz) {
      if (h > ncol(M) || h > ncol(fc_b)) next
      ok <- complete.cases(yout[, h], fc_b[, h], M[, h]); if (sum(ok) < 20) next
      rel <- rmse_fn(yout[, h], M[, h]) / rmse_fn(yout[, h], fc_b[, h])
      pv  <- tryCatch(gw.test(x = fc_b[ok, h], y = M[ok, h], p = yout[ok, h],
                              T = bigT, tau = h, method = "HAC")$p.value,
                      error = function(e) NA_real_)
      rows[[length(rows) + 1]] <- data.frame(competitor = unname(gw_all_map[bf]),
        h = h, rel = round(as.numeric(rel), 2), p = round(as.numeric(pv), 3))
    }
  }
  df  <- do.call(rbind, rows)
  lev <- intersect(rev(unname(gw_all_map)), unique(df$competitor))
  df$competitor <- factor(df$competitor, levels = lev)
  df$stars <- ifelse(is.na(df$p), "", ifelse(df$p < 0.01, "***",
              ifelse(df$p < 0.05, "**", ifelse(df$p < 0.10, "*", ""))))
  df$lab <- sprintf("%.2f\n(%.2f)%s", df$rel, df$p, df$stars)
  g <- ggplot(df, aes(factor(h), competitor, fill = pmin(pmax(rel, 0.5), 1.5))) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = lab), size = 2.7, lineheight = 0.85) +
    scale_fill_gradient2(midpoint = 1, low = "#1A7F37", mid = "white",
                         high = "#B0202A", limits = c(0.5, 1.5), oob = scales::squish,
                         name = sprintf("RMSE(model) / RMSE(%s)", bench_label)) +
    labs(x = "Horizon (h, months)", y = NULL,
         subtitle = sprintf("Below 1 (green) beats %s; GW two-sided HAC p-value below; * p<0.10, ** p<0.05, *** p<0.01", bench_label)) +
    theme_article + theme(panel.grid = element_blank())
  save_article_fig(g, out_name, width = 8.8, height = 7)
}
cat("FIGURE 9: GW relative-RMSE heatmaps (Ridge, AR, 2SRR-AR)...\n")
gw_heatmap("Ridge.rda",   "Ridge",   "FIG9_gw_relrmse_vs_ridge")
gw_heatmap("AR.rda",      "AR",      "FIG9b_gw_relrmse_vs_AR")
gw_heatmap("2SRR_AR.rda", "2SRR-AR", "FIG9c_gw_relrmse_vs_2SRR_AR")

cat("FIGURE 10: CSFE (vs random walk and vs ridge)...\n")
load(file.path(DIR_FC, "rw.rda")); load(file.path(DIR_DATA, "data.rda"))
n_oos     <- nrow(yout)
oos_dates <- tail(as.Date(data$date), n_oos)
hC     <- 1
fc_tvp <- list("2SRR-AR" = ld_fc("2SRR_AR.rda"), "2SRR-Factor" = ld_fc("2SRR_Factor.rda"),
               "2SRR-FAVAR" = ld_fc("2SRR_FAVAR.rda"))
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
    scale_colour_manual(values = MODEL_COLORS, name = NULL) +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    labs(x = "Forecast origin",
         y = sprintf("Cumulative squared-error difference vs. %s", bench_label),
         subtitle = sprintf("Rising = beats %s; a fall flags underperformance (grey: COVID 2020; orange: 2021-22 surge), h = 1", bench_label)) +
    theme_article
  save_article_fig(g, out_name, width = 9.5, height = 5.5)
}
make_csfe(rw,                 "the random walk",    "FIG10_csfe_vs_rw")
make_csfe(ld_fc("Ridge.rda"), "the ordinary ridge", "FIG10b_csfe_vs_ridge")
cat("done\n")

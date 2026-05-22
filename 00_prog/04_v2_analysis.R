# ==============================================================================
# 04_v2_analysis.R
#
# Supplementary analysis script. Reads existing forecasts/betas and generates
# four additional outputs that the original 04_analysis.R does not produce in
# the form required here:
#
#   (A) P2b: TVP-cases relative performance against the AR benchmark (not RW).
#   (B) P7v2: 2SRR-FAVAR against the 3 best Medeiros models per horizon, with
#       a CONSISTENT model->color palette across all four horizons. Realized
#       is always black.
#   (C) P12v2: MCS heatmap with the legend bug fixed (labels match levels) and
#       a stricter alpha=0.50 run to obtain discrimination.
#   (D) P_DM: full Diebold-Mariano matrix for 2SRR-FAVAR vs every other model
#       at each horizon, exported as CSV table + heatmap (mirrors the GW
#       table structure that already exists).
#
# Outputs land in 40_results/run_v2_<timestamp>/  to keep the original run
# folder intact.
# ==============================================================================

cat("== 04_v2_analysis.R ==\n\n")
# Run this script from the project root directory.
source("00_prog/00_setup.R")

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(scales)
  library(gridExtra); library(reshape2)
})
has_patchwork <- requireNamespace("patchwork", quietly = TRUE)
if (has_patchwork) library(patchwork)
has_mcs <- requireNamespace("MCS", quietly = TRUE)
if (has_mcs) library(MCS)

# Output --------------------------------------------------------------------
TS      <- format(Sys.time(), "%Y%m%d_%H%M")
OUT_DIR <- file.path("40_results", paste0("run_v2_", TS))
FIG_DIR <- file.path(OUT_DIR, "figures")
TAB_DIR <- file.path(OUT_DIR, "tables")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
cat("Output:", OUT_DIR, "\n\n")

save_fig <- function(plt, fname, width = 9, height = 5) {
  if (is.null(plt)) { cat("  [skip]", fname, "\n"); return(invisible(NULL)) }
  pdf_path <- file.path(FIG_DIR, paste0(fname, ".pdf"))
  ggsave(pdf_path, plt, width = width, height = height)
  cat("  [fig]", pdf_path, "\n")
  invisible(plt)
}
save_tbl <- function(df, fname, ...) {
  write.csv(df, file.path(TAB_DIR, paste0(fname, ".csv")), row.names = FALSE)
  cat("  [tbl]", file.path(TAB_DIR, paste0(fname, ".csv")), "\n")
}

# ---------------------------------------------------------------------------
# Global model -> color palette. Used in EVERY plot in this script so that
# the same model always shows up with the same colour.
# ---------------------------------------------------------------------------
MODEL_COLORS <- c(
  "Realized"          = "black",
  "2SRR-FAVAR"        = "#D62728",  # red
  "2SRR-AR"           = "#FF7F0E",  # orange
  "2SRR-Factor"       = "#2CA02C",  # green
  "LASSO"             = "#1F77B4",  # blue
  "ElNET"             = "#17BECF",  # cyan
  "AdaLASSO"          = "#9467BD",  # purple
  "AdaElNET"          = "#8C564B",  # brown
  "RF"                = "#E377C2",  # pink
  "Bagging"           = "#7F7F7F",  # grey
  "Factor"            = "#BCBD22",  # olive
  "T.Factor"          = "#1A55A3",  # navy
  "CSR"               = "#FF1493",  # magenta
  "AR"                = "#000080",  # dark navy
  "AR_BIC"            = "#4B0082",  # indigo
  "Ridge"             = "#8B0000",  # dark red
  "RidgeStep1_AR"     = "#FFA07A",  # light salmon
  "RidgeStep1_Factor" = "#A0522D",  # sienna
  "RidgeStep1_FAVAR"  = "#DEB887"   # burlywood
)
MODEL_LTY <- c(
  "Realized"          = "solid",
  "2SRR-FAVAR"        = "solid",
  "2SRR-AR"           = "dashed",
  "2SRR-Factor"       = "dotted",
  "LASSO"             = "dashed",
  "ElNET"             = "dashed",
  "RF"                = "dashed",
  "AR"                = "dotdash",
  "AR_BIC"            = "dotdash",
  "Ridge"             = "dotted",
  "Bagging"           = "dotted",
  "Factor"            = "dotted",
  "T.Factor"          = "dotted",
  "CSR"               = "dotted",
  "AdaLASSO"          = "dashed",
  "AdaElNET"          = "dashed"
)

# 0. Loading ----------------------------------------------------------------
cat("0. Loading inputs...\n")
load(file.path(DIR_DATA,      "data.rda"))
load(file.path(DIR_FORECASTS, "yout.rda"))
load(file.path(DIR_FORECASTS, "rw.rda"))

horizons   <- c(1, 3, 6, 12)
maxh       <- 12
n_oos      <- nrow(yout)
all_dates  <- as.Date(data$date)
oos_dates  <- tail(all_dates, n_oos)
y_raw_global  <- data$CPIAUCSL
tau_global    <- length(y_raw_global) - n_oos
y_oos_monthly <- y_raw_global[(tau_global + 1):(tau_global + n_oos)]

medeiros_models <- c("Ridge", "LASSO", "ElNET", "AdaLASSO", "AdaElNET",
                      "RF", "Bagging", "Factor", "T.Factor", "CSR",
                      "AR", "AR_BIC")
fc_all <- list()
for (m in medeiros_models) {
  fp <- file.path(DIR_FORECASTS, paste0(m, ".rda"))
  if (file.exists(fp)) {
    e <- new.env(); load(fp, envir = e)
    fc_all[[m]] <- as.matrix(get(ls(e)[1], envir = e))
  }
}
cases_tvp <- c("AR", "Factor", "FAVAR")
fc_2srr <- list(); fc_ridge_step1 <- list()
for (case in cases_tvp) {
  fp <- file.path(DIR_FORECASTS, paste0("2SRR_", case, ".rda"))
  if (file.exists(fp)) {
    e <- new.env(); load(fp, envir = e)
    fc_2srr[[case]] <- as.matrix(get(ls(e)[1], envir = e))
  }
  fp <- file.path(DIR_FORECASTS, paste0("Ridge_from_2SRR_", case, ".rda"))
  if (file.exists(fp)) {
    e <- new.env(); load(fp, envir = e)
    fc_ridge_step1[[case]] <- as.matrix(get(ls(e)[1], envir = e))
  }
}

rmse_fn <- function(y, f) {
  ok <- complete.cases(y, f); if (sum(ok) < 5) return(NA_real_)
  sqrt(mean((y[ok] - f[ok])^2))
}

# ===========================================================================
# (A) P2b: TVP cases vs AR benchmark (not RW)
# ===========================================================================
cat("\n(A) P2b: 3 TVP cases relative to AR benchmark...\n")
rmse_ar <- sapply(1:maxh, function(h)
  if (!is.null(fc_all$AR) && h <= ncol(fc_all$AR))
    rmse_fn(yout[, h], fc_all$AR[, h]) else NA_real_)

tvp_vs_ar <- list()
for (case in names(fc_2srr)) {
  for (h in horizons) {
    if (h > ncol(fc_2srr[[case]])) next
    rmse_s <- rmse_fn(yout[, h], fc_2srr[[case]][, h])
    tvp_vs_ar[[length(tvp_vs_ar) + 1]] <- data.frame(
      case = paste0("2SRR-", case), h = h,
      rmse_2srr = round(rmse_s, 4),
      rmse_AR   = round(rmse_ar[h], 4),
      ratio_AR  = round(rmse_s / rmse_ar[h], 4))
  }
}
tvp_vs_ar_df <- do.call(rbind, tvp_vs_ar)
save_tbl(tvp_vs_ar_df, "P2b_tvp_3cases_vs_AR",
         latex_caption = "RMSE of the 3 TVP cases relative to the AR benchmark")

p <- ggplot(tvp_vs_ar_df, aes(x = factor(h), y = ratio_AR, fill = case)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 1, linetype = 2, color = "black") +
  geom_text(aes(label = sprintf("%.3f", ratio_AR)),
            position = position_dodge(width = 0.9), vjust = -0.3, size = 3) +
  scale_fill_manual(values = c("2SRR-AR"     = MODEL_COLORS["2SRR-AR"],
                                "2SRR-Factor" = MODEL_COLORS["2SRR-Factor"],
                                "2SRR-FAVAR"  = MODEL_COLORS["2SRR-FAVAR"])) +
  labs(x = "Horizon (h)", y = "RMSE(TVP) / RMSE(AR)",
        title = "3 TVP cases vs autoregressive benchmark",
        subtitle = "Below 1.0 = TVP beats AR; above 1.0 = AR beats TVP",
        fill = "Case") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
save_fig(p, "P2b_tvp_3cases_ratio_AR", 9, 5)

# ===========================================================================
# (B) P7v2: 2SRR-FAVAR vs the 3 best Medeiros per horizon.
# Design choices to ensure visual consistency across horizons:
#   1. Compute the UNIVERSE of models that appear in any panel (union of
#      the per-horizon best3). All panels share the SAME factor levels —
#      one model maps to one colour irrespective of the horizon.
#   2. In each panel, columns for models outside that horizon's best3 are
#      filled with NA. The line is skipped (na.rm=TRUE) but the model
#      retains its colour slot in the global legend.
#   3. Realized: solid, thick (1.3). 2SRR-FAVAR: solid, medium (0.95).
#      Other models: dashed, thin (0.7). Visual hierarchy is clear.
#   4. patchwork::wrap_plots(..., guides = "collect") produces a SINGLE
#      legend below the entire figure, eliminating duplication.
# ===========================================================================
cat("\n(B) P7v2: 2SRR-FAVAR vs 3 best Medeiros per horizon...\n")

# Step 1: compute the universe of best-3 models across all horizons
all_best <- list()
for (h in horizons) {
  if (h > ncol(fc_2srr$FAVAR)) next
  rmses <- sapply(fc_all, function(M)
    if (h <= ncol(M)) rmse_fn(yout[, h], M[, h]) else NA)
  rmses <- rmses[!is.na(rmses)]
  if (length(rmses) > 0)
    all_best[[as.character(h)]] <- names(sort(rmses))[1:min(3, length(rmses))]
}
universe_models <- unique(c("Realized", "2SRR-FAVAR", unlist(all_best)))
cat("  Universe of models across panels:",
    paste(universe_models, collapse = ", "), "\n")

# Step 2: per-horizon plot with shared factor levels
plot_p7v2 <- function(h) {
  if (is.null(fc_2srr$FAVAR) || h > ncol(fc_2srr$FAVAR)) return(NULL)
  best3 <- all_best[[as.character(h)]]
  if (is.null(best3)) return(NULL)

  df_p <- data.frame(date = oos_dates, check.names = FALSE)
  df_p[["Realized"]]   <- y_oos_monthly
  df_p[["2SRR-FAVAR"]] <- fc_2srr$FAVAR[, h] / h
  for (m in best3) df_p[[m]] <- fc_all[[m]][, h] / h
  # Fill missing models with NA so every panel has the same column set
  for (m in universe_models)
    if (!(m %in% colnames(df_p))) df_p[[m]] <- NA_real_

  df_long <- pivot_longer(df_p, -date, names_to = "Series", values_to = "v")
  df_long$Series <- factor(df_long$Series, levels = universe_models)
  # Assign a per-series linewidth so Realized stands out from the rest
  df_long$lw <- ifelse(df_long$Series == "Realized", 1.3,
                  ifelse(df_long$Series == "2SRR-FAVAR", 0.95, 0.7))

  ggplot(df_long, aes(date, v, color = Series, linetype = Series,
                       linewidth = lw, group = Series)) +
    geom_line(alpha = 0.92, na.rm = TRUE) +
    geom_hline(yintercept = 0, linetype = 3, alpha = 0.4) +
    scale_color_manual(values = MODEL_COLORS, drop = FALSE,
                        breaks = universe_models) +
    scale_linetype_manual(values = MODEL_LTY,  drop = FALSE,
                           breaks = universe_models) +
    scale_linewidth_identity() +
    labs(title = sprintf("h = %d months", h), x = NULL,
          y = "Monthly inflation (%)",
          color = NULL, linetype = NULL) +
    theme_minimal(base_size = 10)
}
plots_p7v2 <- lapply(setNames(horizons, paste0("h", horizons)), plot_p7v2)

# Save individual panels
for (h in horizons) {
  p <- plots_p7v2[[paste0("h", h)]]
  if (!is.null(p)) save_fig(p +
                              theme(legend.position = "bottom",
                                    legend.text = element_text(size = 8)),
                            sprintf("P7v2_2srr_vs_3best_h%02d", h), 10, 5)
}

# Combined 4h panel with a SINGLE shared legend
if (has_patchwork) {
  panel <- patchwork::wrap_plots(plots_p7v2, ncol = 2, guides = "collect") +
    patchwork::plot_annotation(
      title = "2SRR-FAVAR vs the 3 best Medeiros models — 4 horizons",
      subtitle = paste0("Realized (black) = monthly y_t (identical in all panels). ",
                        "Forecasts h-steps / h = monthly-rate prediction.")
    ) &
    theme(legend.position = "bottom",
          legend.text = element_text(size = 9),
          legend.title = element_blank())
  ggsave(file.path(FIG_DIR, "P7v2_2srr_vs_3best_4h.pdf"),
         panel, width = 14, height = 9)
  cat("  [fig 4h]", file.path(FIG_DIR, "P7v2_2srr_vs_3best_4h.pdf"), "\n")
}

# ===========================================================================
# (C) P12v2: MCS heatmap with legend FIX + alpha = 0.50 for discrimination
# ===========================================================================
cat("\n(C) P12v2: MCS heatmap (legend fix + alpha=0.50)...\n")
if (has_mcs) {
  all_fc_for_mcs <- c(fc_all,
                      setNames(fc_2srr, paste0("2SRR_", names(fc_2srr))),
                      setNames(fc_ridge_step1,
                               paste0("RidgeStep1_", names(fc_ridge_step1))))
  mcs_rows <- list()
  for (h in horizons) {
    models_h <- names(all_fc_for_mcs)[
      sapply(all_fc_for_mcs, function(M) h <= ncol(M))]
    if (length(models_h) < 3) next
    L <- sapply(models_h, function(mn)
      (yout[, h] - all_fc_for_mcs[[mn]][, h])^2)
    cols_all_na <- apply(L, 2, function(x) all(is.na(x)))
    if (any(cols_all_na)) {
      L <- L[, !cols_all_na, drop = FALSE]
      models_h <- models_h[!cols_all_na]
    }
    ok <- complete.cases(L)
    if (sum(ok) < 30) next
    L <- L[ok, , drop = FALSE]
    # Stricter MCS at alpha = 0.50 to obtain discrimination
    mcs_50 <- tryCatch(
      MCS::MCSprocedure(Loss = L, alpha = 0.50, B = 1000,
                         statistic = "Tmax", verbose = FALSE),
      error = function(e) NULL)
    in_50 <- if (!is.null(mcs_50)) rownames(mcs_50@show) else character(0)
    # Also keep the original 90% set
    mcs_10 <- tryCatch(
      MCS::MCSprocedure(Loss = L, alpha = 0.10, B = 1000,
                         statistic = "Tmax", verbose = FALSE),
      error = function(e) NULL)
    in_10 <- if (!is.null(mcs_10)) rownames(mcs_10@show) else character(0)
    for (m in models_h) {
      mcs_rows[[length(mcs_rows) + 1]] <- data.frame(
        model = m, h = h,
        in_MCS_90 = m %in% in_10,
        in_MCS_50 = m %in% in_50)
    }
  }
  mcs_long_v2 <- do.call(rbind, mcs_rows)
  save_tbl(mcs_long_v2, "P12v2_MCS_long_table",
           latex_caption = "MCS membership at alpha=0.10 (90% set) and alpha=0.50 (50% set, stricter)")

  # Heatmap with corrected labels
  mcs_long_v2$Status <- factor(
    ifelse(mcs_long_v2$in_MCS_50, "In 50% MCS (most informative)",
    ifelse(mcs_long_v2$in_MCS_90, "In 90% MCS only", "Out")),
    levels = c("Out", "In 90% MCS only", "In 50% MCS (most informative)"))
  p_mcs <- ggplot(mcs_long_v2, aes(x = factor(h), y = model, fill = Status)) +
    geom_tile(color = "white") +
    scale_fill_manual(values = c(
      "Out" = "grey85",
      "In 90% MCS only" = "#FFD580",
      "In 50% MCS (most informative)" = "#2CA02C"),
      drop = FALSE) +
    labs(title = "Model Confidence Set per horizon",
          subtitle = "Green = survives stricter 50% MCS; orange = only 90% MCS; grey = eliminated",
          x = "Horizon (h)", y = NULL) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom")
  save_fig(p_mcs, "P12v2_MCS_heatmap", 10, 7)
}

# ===========================================================================
# (D) P_DM: full DM matrix for 2SRR-FAVAR vs every benchmark
# ===========================================================================
cat("\n(D) DM matrix for 2SRR-FAVAR vs all benchmarks...\n")
dm_safe <- function(y, f1, f2, h) {
  ok <- complete.cases(y, f1, f2)
  if (sum(ok) < 20) return(list(stat = NA, p = NA))
  tryCatch({
    d <- forecast::dm.test(ts((y[ok] - f1[ok])^2),
                            ts((y[ok] - f2[ok])^2),
                            alternative = "two.sided", h = h)
    list(stat = as.numeric(d$statistic), p = d$p.value)
  }, error = function(e) list(stat = NA, p = NA))
}

benchmarks_dm <- c(intersect(c("Ridge", "LASSO", "ElNET", "AdaLASSO",
                                "AdaElNET", "RF", "Bagging", "Factor",
                                "T.Factor", "CSR", "AR", "AR_BIC"),
                              names(fc_all)),
                    paste0("RidgeStep1_", names(fc_ridge_step1)))

# Loop over BOTH 2SRR references: "AR" is the championed spec; "FAVAR" kept for
# continuity with the earlier run. signed_p sign FIXED so GREEN = the 2SRR
# reference is more accurate (previously inverted relative to the legend).
for (ref_case in c("AR", "FAVAR")) {
  if (is.null(fc_2srr[[ref_case]])) next
  dm_rows <- list()
  for (h in horizons) {
    if (h > ncol(fc_2srr[[ref_case]])) next
    y_h <- yout[, h]
    f_ref <- fc_2srr[[ref_case]][, h]
    rmse_ref <- rmse_fn(y_h, f_ref)
    for (bn in benchmarks_dm) {
      M <- if (startsWith(bn, "RidgeStep1_"))
              fc_ridge_step1[[sub("RidgeStep1_", "", bn)]]
            else fc_all[[bn]]
      if (is.null(M) || h > ncol(M)) next
      f_bn <- M[, h]
      dm <- dm_safe(y_h, f_ref, f_bn, h)
      rmse_bn <- rmse_fn(y_h, f_bn)
      dm_rows[[length(dm_rows) + 1]] <- data.frame(
        h = h, ref = paste0("2SRR-", ref_case), benchmark = bn,
        RMSE_ref = round(rmse_ref, 4),
        RMSE_bn  = round(rmse_bn,  4),
        DM_stat  = round(dm$stat, 3),
        DM_p     = round(dm$p,    4),
        ref_wins = rmse_ref < rmse_bn)
    }
  }
  dm_df <- do.call(rbind, dm_rows)
  save_tbl(dm_df, sprintf("P_DM_2SRR_%s_vs_all", ref_case),
           latex_caption = sprintf("Diebold-Mariano test: 2SRR-%s vs all benchmark models (two-sided)",
                                    ref_case))

  # DM p-value heatmap. GREEN (negative) = 2SRR ref more accurate; RED = benchmark.
  dm_df$signed_p <- ifelse(dm_df$ref_wins, -dm_df$DM_p, dm_df$DM_p)
  p_dm <- ggplot(dm_df, aes(x = factor(h), y = benchmark,
                              fill = pmin(pmax(signed_p, -0.5), 0.5))) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.3f", DM_p)), size = 3) +
    scale_fill_gradient2(midpoint = 0,
                          low = "darkgreen", mid = "white", high = "darkred",
                          limits = c(-0.5, 0.5), name = "Signed p-value") +
    labs(title = sprintf("Diebold-Mariano p-values: 2SRR-%s vs each benchmark", ref_case),
          subtitle = sprintf("Green = 2SRR-%s has lower RMSE; red = benchmark has lower RMSE", ref_case),
          x = "Horizon (h)", y = "Benchmark") +
    theme_minimal(base_size = 10)
  save_fig(p_dm, sprintf("P_DM_2SRR_%s_heatmap", ref_case), 9, 6)

  # Summary of significant DM rejections
  dm_signif <- dm_df %>%
    group_by(h) %>%
    summarise(
      n_total       = n(),
      n_p_below_10  = sum(DM_p < 0.10, na.rm = TRUE),
      n_2srr_wins_signif = sum(DM_p < 0.10 & ref_wins,  na.rm = TRUE),
      n_2srr_loses_signif = sum(DM_p < 0.10 & !ref_wins, na.rm = TRUE),
      .groups = "drop")
  save_tbl(as.data.frame(dm_signif), sprintf("P_DM_%s_summary", ref_case),
           latex_caption = sprintf("Diebold-Mariano significant rejections summary, 2SRR-%s (alpha=0.10)",
                                    ref_case))
}

cat("\n== 04_v2_analysis.R DONE ==\n")
cat(sprintf("All outputs in: %s\n", OUT_DIR))

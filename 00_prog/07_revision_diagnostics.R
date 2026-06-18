# ==============================================================================
# 07_revision_diagnostics.R  --  robustness diagnostics for the major-revision
# points. RE-EVALUATES the saved forecasts/betas only; no re-forecasting.
#
#  (A) MCS power (point 2): re-run the Model Confidence Set with the absolute-error
#      (MAE) loss and on the sample EXCLUDING the 2020-2022 shock, to see whether
#      the test regains power once outliers are tamed.
#  (B) Loss-differential variance inflation by the pandemic (supports A).
#  (C) Time-variation diagnostic (point 1): coefficient drift and 2SRR-ARDI vs its
#      own Step-1 STATIC ridge -- is the time variation active or collapsed?
#  (D) Giacomini-Rossi (2010) Fluctuation Test (point 3): relative predictive
#      ability of 2SRR-ARDI vs AR over time, instead of an ex-post sample split.
#  (E) Mincer-Zarnowitz (point 4): report the slope estimate and re-run pre-2020.
#
# Run from the project root.
# ==============================================================================
cat("== 07_revision_diagnostics.R ==\n\n")
suppressPackageStartupMessages({ library(MCS); library(sandwich); library(ggplot2) })

FC  <- "30_output/forecasts"; BT <- "30_output/betas"
OUT <- "40_results/revision"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
ld  <- function(p){e<-new.env();load(p,envir=e);get(ls(e)[1],envir=e)}
M   <- function(n) as.matrix(ld(file.path(FC, paste0(n,".rda"))))

yout <- M("yout"); rw <- M("rw")
data <- ld("10_data/data.rda"); dts <- as.Date(data$date)
n <- length(dts); nw <- nrow(yout); tau <- n - nw
tgt  <- dts[(tau+1):n]
hz   <- c(1,3,6,12)

# masks over the 312 windows
m_full <- rep(TRUE, nw)
m_pre  <- tgt <  as.Date("2020-01-01")                                   # 246
m_xpan <- tgt <  as.Date("2020-01-01") | tgt >= as.Date("2023-01-01")    # exclude 2020-2022

models <- c("AR","AR_BIC","Ridge","LASSO","ElNET","AdaLASSO","AdaElNET","RF",
            "Bagging","Factor","T.Factor","CSR","2SRR_AR","2SRR_Factor","2SRR_FAVAR")
disp <- function(x){x<-gsub("2SRR_FAVAR","2SRR-ARDI",x);x<-gsub("2SRR_","2SRR-",x);x}
FC_list <- setNames(lapply(models, M), models)

# ============================================================================ #
# (A) MCS under different losses / samples
# ============================================================================ #
cat("############ (A) MCS robustness (squared vs MAE; full vs excl. 2020-22) ####\n")
mcs_survivors <- function(Loss, alpha = 0.25, B = 1000) {
  res <- tryCatch({
    invisible(capture.output(
      out <- MCS::MCSprocedure(Loss = Loss, alpha = alpha, B = B,
                               statistic = "Tmax")))
    out
  }, error = function(e) { message("  MCS error: ", conditionMessage(e)); NULL })
  if (is.null(res)) return(list(surv = NA_character_, k = NA_integer_))
  sm <- tryCatch(slot(res, "show"), error = function(e) NULL)
  if (is.null(sm) || is.null(rownames(sm)))
    return(list(surv = NA_character_, k = NA_integer_))
  list(surv = rownames(sm), k = nrow(sm))
}
build_loss <- function(h, mask, power) {
  L <- sapply(models, function(m) {
    e <- yout[mask, h] - FC_list[[m]][mask, h]
    if (power == 2) e^2 else abs(e)
  })
  colnames(L) <- disp(models)
  L[complete.cases(L), , drop = FALSE]
}
mcs_rows <- list()
for (h in hz) {
  for (cfg in list(list(lab="MSE_full", mask=m_full, pw=2),
                   list(lab="MSE_excl2020-22", mask=m_xpan, pw=2),
                   list(lab="MAE_full", mask=m_full, pw=1))) {
    r <- mcs_survivors(build_loss(h, cfg$mask, cfg$pw))
    surv <- if (length(r$surv) && !all(is.na(r$surv))) paste(r$surv, collapse=", ") else "NA"
    cat(sprintf("  h=%2d | %-16s | survivors: %2s/%d\n", h, cfg$lab,
                ifelse(is.na(r$k), "NA", r$k), length(models)))
    mcs_rows[[length(mcs_rows)+1]] <- data.frame(h=h, config=cfg$lab,
                                                 n_survivors=r$k, survivors=surv)
  }
}
write.csv(do.call(rbind, mcs_rows), file.path(OUT,"A_mcs_robustness.csv"), row.names=FALSE)

# ============================================================================ #
# (B) Variance inflation of the loss differential by the pandemic
#     (ARDI vs the worst model, Bagging; and ARDI vs AR)
# ============================================================================ #
cat("\n############ (B) Loss-differential variance: full vs excl. 2020-22 ########\n")
var_rows <- list()
for (h in hz) for (pair in list(c("Bagging"), c("AR"))) {
  d_full <- (yout[m_full,h]-FC_list[[pair]][m_full,h])^2 - (yout[m_full,h]-FC_list[["2SRR_FAVAR"]][m_full,h])^2
  d_xpan <- (yout[m_xpan,h]-FC_list[[pair]][m_xpan,h])^2 - (yout[m_xpan,h]-FC_list[["2SRR_FAVAR"]][m_xpan,h])^2
  ratio <- var(d_full, na.rm=TRUE)/var(d_xpan, na.rm=TRUE)
  cat(sprintf("  h=%2d | ARDI vs %-7s | var(full)/var(excl)= %.1fx\n", h, pair, ratio))
  var_rows[[length(var_rows)+1]] <- data.frame(h=h, vs=pair,
    var_full=round(var(d_full,na.rm=TRUE),4), var_excl=round(var(d_xpan,na.rm=TRUE),4),
    inflation_x=round(ratio,1))
}
write.csv(do.call(rbind, var_rows), file.path(OUT,"B_variance_inflation.csv"), row.names=FALSE)

# ============================================================================ #
# (C) Time-variation diagnostic: does 2SRR-ARDI actually time-vary?
# ============================================================================ #
cat("\n############ (C) Time variation: coefficient drift + ARDI vs static ridge #\n")
b_fav <- ld(file.path(BT,"betas_2SRR_FAVAR.rda"))
tv_rows <- list()
for (h in hz) {
  bl <- b_fav[[paste0("h",h)]]
  # mean over windows of the within-window temporal SD of each coefficient,
  # scaled by the coefficient level (|mean|): drift-to-level ratio.
  ratios <- sapply(bl, function(w){
    B <- w$betas_tvp; if (is.null(B)) return(NA)
    sds <- apply(B,2,sd); lev <- pmax(colMeans(abs(B)),1e-8); median(sds/lev,na.rm=TRUE)
  })
  # ARDI vs its Step-1 STATIC ridge forecasts
  ardi <- FC_list[["2SRR_FAVAR"]][,h]; stat <- M("Ridge_from_2SRR_FAVAR")[,h]
  ok <- complete.cases(yout[,h],ardi,stat)
  rmse_ardi <- sqrt(mean((yout[ok,h]-ardi[ok])^2)); rmse_stat <- sqrt(mean((yout[ok,h]-stat[ok])^2))
  cat(sprintf("  h=%2d | median drift/level=%.3f | corr(ARDI,staticRidge)=%.3f | RMSE ARDI=%.3f static=%.3f\n",
              h, median(ratios,na.rm=TRUE), cor(ardi[ok],stat[ok]), rmse_ardi, rmse_stat))
  tv_rows[[length(tv_rows)+1]] <- data.frame(h=h,
    drift_to_level=round(median(ratios,na.rm=TRUE),3),
    corr_ardi_static=round(cor(ardi[ok],stat[ok]),3),
    rmse_ardi=round(rmse_ardi,3), rmse_static=round(rmse_stat,3))
}
write.csv(do.call(rbind, tv_rows), file.path(OUT,"C_time_variation.csv"), row.names=FALSE)

# ============================================================================ #
# (D) Giacomini-Rossi (2010) Fluctuation Test: ARDI vs AR over time
# ============================================================================ #
cat("\n############ (D) Giacomini-Rossi Fluctuation Test (ARDI vs AR) ###########\n")
gr_cv <- 3.393  # two-sided 5% critical value for mu = m/P = 0.3 (GR 2010, Table 1)
mu <- 0.3
flu_all <- list()
for (h in c(1,12)) {
  d <- (yout[,h]-FC_list[["AR"]][,h])^2 - (yout[,h]-FC_list[["2SRR_FAVAR"]][,h])^2  # >0: ARDI better
  ok <- which(is.finite(d)); d <- d[ok]; dd <- tgt[ok]; P <- length(d)
  m <- max(20, round(mu*P))
  sig <- sqrt(as.numeric(sandwich::lrvar(d, type="Newey-West") * P))  # HAC long-run sd of the sum-scale
  half <- floor(m/2); Ft <- rep(NA_real_, P)
  for (t in (half+1):(P-half)) {
    s <- sum(d[(t-half):(t-half+m-1)])
    Ft[t] <- s / (sig * sqrt(m))
  }
  flu_all[[as.character(h)]] <- data.frame(date=dd, F=Ft, h=h)
}
flu_df <- do.call(rbind, flu_all)
write.csv(flu_df, file.path(OUT,"D_fluctuation.csv"), row.names=FALSE)
flu_df$h <- factor(paste0("h = ",flu_df$h), levels=c("h = 1","h = 12"))
pfl <- ggplot(flu_df, aes(date, F)) +
  geom_hline(yintercept=c(-gr_cv,gr_cv), linetype=2, colour="#B0202A") +
  geom_hline(yintercept=0, colour="grey70") +
  geom_line(colour="#1F4E9C", linewidth=0.8, na.rm=TRUE) +
  facet_wrap(~h, ncol=1) +
  labs(x="Forecast origin", y="Fluctuation statistic",
       subtitle=paste0("Giacomini-Rossi (2010) test, 2SRR-ARDI vs AR (mu = 0.3).\n",
                       "Above +", gr_cv, " = ARDI significantly better at that date; dashed = 5% bands.")) +
  theme_minimal(base_size=12)+theme(plot.subtitle=element_text(size=9,colour="grey30"))
ggsave(file.path("40_results/article_figures/figures","FIG_fluctuation_ARDI_vs_AR.png"),
       pfl, width=9, height=6, dpi=150, bg="white")
for (h in c(1,12)) { f <- flu_all[[as.character(h)]]$F
  cat(sprintf("  h=%2d | max F=%.2f | %% of dates with F>+CV: %.0f%%\n",
              h, max(f,na.rm=TRUE), 100*mean(f>gr_cv,na.rm=TRUE))) }

# ============================================================================ #
# (E) Mincer-Zarnowitz: slope estimate + pre-2020 sub-sample
# ============================================================================ #
cat("\n############ (E) Mincer-Zarnowitz: slope (bias) + pre-2020 ###############\n")
mz <- function(y,f,h){ok<-complete.cases(y,f); if(sum(ok)<30)return(c(a=NA,b=NA,p=NA))
  y<-y[ok];f<-f[ok]; m<-lm(y~f); cf<-coef(m)
  V<-tryCatch(sandwich::NeweyWest(m,lag=max(0,h-1),prewhite=FALSE,adjust=TRUE),error=function(e)vcov(m))
  d<-cf-c(0,1); w<-tryCatch(as.numeric(t(d)%*%solve(V)%*%d),error=function(e)NA)
  c(a=unname(cf[1]),b=unname(cf[2]),p=if(is.na(w)||w<0)NA else 1-pchisq(w,df=2))}
mz_rows<-list()
for (mdl in c("AR","2SRR_FAVAR","LASSO")) for (h in hz) {
  full<-mz(yout[,h],FC_list[[mdl]][,h],h)
  pre <-mz(yout[m_pre,h],FC_list[[mdl]][m_pre,h],h)
  cat(sprintf("  %-10s h=%2d | full: beta=%5.2f p=%.3f | pre-2020: beta=%5.2f p=%.3f\n",
              disp(mdl),h,full["b"],full["p"],pre["b"],pre["p"]))
  mz_rows[[length(mz_rows)+1]]<-data.frame(model=disp(mdl),h=h,
    beta_full=round(full["b"],3),p_full=round(full["p"],3),
    beta_pre=round(pre["b"],3),p_pre=round(pre["p"],3))
}
write.csv(do.call(rbind, mz_rows), file.path(OUT,"E_mincer_zarnowitz.csv"), row.names=FALSE)

cat("\n== done. CSVs in 40_results/revision/, figure in article_figures/figures/ ==\n")

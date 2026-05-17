# ==============================================================================
# 02_forecast_medeiros.R
#
# Pipeline Medeiros: 12 modelos (Ridge, LASSO, ElNET, AdaLASSO, AdaElNET, RF,
# Bagging, Factor, T.Factor, CSR, AR, AR_BIC) sobre FRED-MD, horizontes 1/3/6/12.
#
# Target: cumulativo trailing Y_h(t) = sum_{j=0}^{h-1} y_{t-j} (id. ao 2SRR).
# Janela: rolling fixo (Medeiros rolling_window.R), 180 janelas POOS.
# ==============================================================================
cat("== 02_forecast_medeiros.R ==\n\n")
source("00_prog/00_setup.R")

load(file.path(DIR_DATA, "data.rda"))
load(file.path(DIR_FORECASTS, "yout.rda"))   # garante n_oos consistente
nwindows <- nrow(yout)                       # 180 — NUNCA hard-codar

variable <- "CPIAUCSL"
horizons <- c(1, 3, 6, 12)
maxh     <- 12

# CRITICAL: remove date column (embed() and model functions require numeric)
dates_col <- data$date
data$date <- NULL
rownames(data) <- as.character(dates_col)
cat("Data prepared:", nrow(data), "obs x", ncol(data),
    "numeric vars | nwindows =", nwindows, "\n")

# ----- Wrappers de modelos -----
# As funcoes runlasso e runar do Medeiros aceitam parametros (alpha, adaptive,
# type) que distinguem variantes. Os wrappers abaixo expoem essas variantes
# como funcoes nomeadas, compativeis com a interface esperada por
# rolling_window(fn = ...).
runridge        <- function(ind, df, variable, horizon)
  runlasso(ind, df, variable, horizon, alpha = 0)
runelnet        <- function(ind, df, variable, horizon)
  runlasso(ind, df, variable, horizon, alpha = 0.5)
runadaptlasso   <- function(ind, df, variable, horizon)
  runlasso(ind, df, variable, horizon, alpha = 1, adaptive = TRUE)
runadaptelnet   <- function(ind, df, variable, horizon)
  runlasso(ind, df, variable, horizon, alpha = 0.5, alpha2 = 0.5, adaptive = TRUE)
runarbic        <- function(ind, df, variable, horizon)
  runar(ind, df, variable, horizon, type = "bic")
# Aliases para os nomes "didaticos" (Medeiros e meio inconsistente nos nomes)
runbag          <- function(...) runbagging(...)
runfactor       <- function(...) runfact(...)

# Patch do runtfact: a versao original do Medeiros (functions.R) so cria a
# variavel `dummy` quando variable == "CPI" ou "PCE" (nomes da base BR do
# paper Medeiros et al. 2021 JBES). Para FRED-MD a variavel se chama
# "CPIAUCSL" e o ramo nao executa, deixando `dummy` indefinido -> tail()
# falha. Aqui criamos uma versao defensiva que SEMPRE define dummy.
runtfact_safe <- function(ind, df, variable, horizon) {
  dfaux <- df[ind, ]
  # Cria dummy 2008-11-01 se data existe E variavel e' inflacao (CPI variants)
  has_2008 <- "2008-11-01" %in% rownames(dfaux)
  is_cpi   <- variable %in% c("CPI", "PCE", "CPIAUCSL", "PCEPI")
  if (has_2008 && is_cpi) {
    dummy <- rep(0, nrow(dfaux))
    dummy[which(rownames(dfaux) == "2008-11-01")] <- 1
  } else {
    dummy <- rep(0, nrow(dfaux))   # default: sem dummy
  }

  index <- which(colnames(dfaux) == variable)
  # Forca matriz numerica (Medeiros original assume base BR ja em matriz; com
  # data.frame da FRED-MD, cbind preserva data.frame -> lm() em tfaux quebra
  # com "tipo invalido (list) para variavel 'w'" porque X[,fixed.controls]
  # vira data.frame em vez de matriz).
  y_emb     <- embed(as.numeric(dfaux[, variable]), 5)
  x_others  <- as.matrix(dfaux[, -index, drop = FALSE])
  n_keep    <- nrow(dfaux) - 4
  mat <- cbind(y_emb,
                tail(dummy,    n_keep),
                tail(x_others, n_keep))
  storage.mode(mat) <- "double"
  pretest <- tfaux(mat, pre.testing = "individual",
                    fixed.controls = 1:4)[-c(1:6)]
  pretest[pretest != 0] <- 1
  aux <- rep(0, ncol(dfaux)); aux[index] <- 1; aux[-index] <- pretest
  # selected = indices em x_others = dfaux[, -index], NAO em df direto.
  # Mapeamos de volta para as colunas reais do df. O Medeiros original
  # esquece esse passo — funciona la' so porque target e' 1a coluna.
  other_cols    <- setdiff(seq_len(ncol(dfaux)), index)
  sel_in_others <- which(pretest == 1)
  if (length(sel_in_others) < 2) {
    sel_in_others <- 1:min(10, length(other_cols))   # fallback: 10 primeiras
  }
  selected_in_df <- other_cols[sel_in_others]
  # dfreduced SEMPRE contem a coluna `variable` (na posicao 1) + as
  # selecionadas. dataprep depois precisa de `variable` por nome.
  dfreduced <- df[, c(index, selected_in_df), drop = FALSE]

  prep_data <- dataprep(ind, dfreduced, variable, horizon,
                          add_dummy = FALSE, factonly = TRUE)
  Xin   <- prep_data$Xin
  yin   <- prep_data$yin
  Xout  <- prep_data$Xout
  dummy_pd <- prep_data$dummy

  bb <- Inf; modelest <- NULL; f.coef <- NULL; coefdum <- 0
  for (i in seq(5, min(20, ncol(Xin)), 5)) {
    m <- tryCatch(lm(yin ~ Xin[, 1:i] + dummy_pd),
                   error = function(e) NULL)
    if (is.null(m)) next
    crit <- BIC(m)
    if (is.finite(crit) && crit < bb) {
      bb <- crit
      modelest <- m
      f.coef <- coef(modelest)
      coefdum <- f.coef[length(f.coef)]
      f.coef  <- f.coef[-length(f.coef)]
    }
  }
  if (is.null(modelest)) return(list(forecast = NA_real_,
                                       outputs = list(coef = NULL)))

  coef <- rep(0, ncol(Xin) + 1)
  coef[1:length(f.coef)] <- f.coef
  coef <- c(coef, coefdum)
  coef[is.na(coef)] <- 0
  forecast <- (cbind(1, Xout, 0) %*% coef)[1]
  list(forecast = forecast, outputs = list(coef = coef))
}
runtargetfactor <- runtfact_safe

# Lista de modelos. fn = nome da funcao (devera existir em escopo global apos
# wrappers acima ou apos source do functions.R do Medeiros).
all_models <- list(
  list(name = "LASSO",    fn = "runlasso"),
  list(name = "Ridge",    fn = "runridge"),
  list(name = "ElNET",    fn = "runelnet"),
  list(name = "AdaLASSO", fn = "runadaptlasso"),
  list(name = "AdaElNET", fn = "runadaptelnet"),
  list(name = "RF",       fn = "runrf"),
  list(name = "Bagging",  fn = "runbag"),
  list(name = "Factor",   fn = "runfactor"),
  list(name = "T.Factor", fn = "runtargetfactor"),
  list(name = "CSR",      fn = "runcsr"),
  list(name = "AR",       fn = "runar"),
  list(name = "AR_BIC",   fn = "runarbic")
)

cat("\nAvailable functions: ")
avail <- sapply(all_models, function(m) exists(m$fn, mode = "function"))
cat(paste(sapply(all_models[avail], "[[", "name"), collapse = ", "), "\n\n")

for (m in all_models) {
  out_path <- file.path(DIR_FORECASTS, paste0(m$name, ".rda"))
  if (file.exists(out_path)) { cat(sprintf("  %-12s exists\n", m$name)); next }
  if (!exists(m$fn, mode = "function")) { next }

  cat(sprintf("  %-12s running...\n", m$name))
  t0 <- Sys.time()
  
  forecasts_mat <- matrix(NA_real_, nwindows, maxh)
  betas_bundle <- list()
  
  for (h in horizons) {
    cat(sprintf("    h=%2d...", h))
    
    # Direct forecasting: create cumulative target for horizon h
    data_h <- data
    y_h <- as.numeric(stats::filter(data[[variable]], rep(1, h), sides = 1))
    if (h > 1) y_h[1:(h-1)] <- y_h[h]
    data_h[[variable]] <- y_h
    
    tryCatch({
      result <- rolling_window(get(m$fn), data_h, nwindows, h, variable)
      forecasts_mat[, h] <- result$forecast
      betas_bundle[[paste0("h", h)]] <- result$outputs
      cat(" done.\n")
    }, error = function(e) cat(sprintf(" FAILED: %s\n", e$message)))
  }
  
  tryCatch({
    forecasts <- forecasts_mat
    colnames(forecasts) <- paste0("h", 1:maxh)
    save(forecasts, file = out_path)
    
    # Save the betas and lambdas bundle
    if (length(betas_bundle) > 0) {
      save(betas_bundle, file = file.path(DIR_BETAS, paste0("betas_", m$name, ".rda")))
    }
    
    cat(sprintf("  %-12s %.1f min\n", m$name, difftime(Sys.time(), t0, units = "mins")))
  }, error = function(e) cat(sprintf(" FAILED to save: %s\n", e$message)))
}

cat("== done ==\n")

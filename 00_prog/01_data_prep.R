# ==============================================================================
# 01_data_prep.R
#
# Gera matrizes de realizados (yout) e benchmark Random Walk (rw) para POOS.
#
# IMPORTANTE: alinhamento com o rolling_window() do Medeiros
#   - Medeiros usa janelas indexadas por j=1..nwindows, cada uma cobrindo
#     obs j..(j + window_size - 1) onde window_size = n - nwindows.
#   - A ultima obs IN-SAMPLE da janela j e portanto: T_j = j + window_size - 1
#                                                       = (n - nwindows) + j - 1
#                                                       = tau + j - 1
#   - O Medeiros prediz Y_h(T_j + h) = sum(y[(T_j+1):(T_j+h)]).
#   - Logo: yout[j, h] = sum(y[(tau+j):(tau+j+h-1)])  com t_end = tau + j - 1.
#
# Versao anterior usava t_end <- tau + i (um passo a mais), o que causava
# desalinhamento de 1 mes entre as previsoes Medeiros e yout, inflando RMSE.
# ==============================================================================
cat("== 01_data_prep.R ==\n\n")
source("00_prog/00_setup.R")

load(file.path(DIR_DATA, "data.rda"))
cat("Data:", nrow(data), "obs x", ncol(data), "vars\n")

variable <- "CPIAUCSL"
nwindows <- 180
maxh     <- 12
y   <- data[[variable]]
n   <- length(y)
tau <- n - nwindows
cat(sprintf("Janelas POOS: %d | tau=%d | window_size=%d\n",
            nwindows, tau, tau))

yout <- matrix(NA_real_, nwindows, maxh)
rw   <- matrix(NA_real_,  nwindows, maxh)
for (h in 1:maxh) for (i in 1:nwindows) {
  t_end <- tau + i - 1                                # última obs in-sample
  if ((t_end + h) <= n) yout[i, h] <- y[t_end + h]    # taxa em t+h  (π_{t+h})
  rw[i, h] <- y[t_end]                                # RW = última taxa observada
}
colnames(yout) <- colnames(rw) <- paste0("h", 1:maxh)

save(yout, file = file.path(DIR_FORECASTS, "yout.rda"))
save(rw,   file = file.path(DIR_FORECASTS, "rw.rda"))
cat("Saved yout e rw:", nwindows, "x", maxh, "\n")
cat("yout[1, 1] = sum(y[", tau + 1, ":", tau + 1,
    "]) =", yout[1, 1], "\n", sep = "")
cat("yout[", nwindows, ", 1] = sum(y[", n, ":", n,
    "]) =", yout[nwindows, 1], "\n", sep = "")
cat("== done ==\n")

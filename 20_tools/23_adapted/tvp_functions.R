# ==============================================================================
# tvp_functions.R
#
# Implementacao do Algoritmo 1 de Coulombe (IJF 2025, sec. 2.4) — Two-Step
# Ridge Regression (2SRR) — em formulacao dual rapida.
#
# DERIVACAO MATEMATICA (compativel bit-a-bit com o paper):
#
#   Modelo:        y_t = X_t' beta_t + eps_t,  beta_t = beta_{t-1} + u_t
#   Reparametriza: beta = C theta, com C = (I_K kron C_T), C_T triangular
#                  inferior de 1's. Para a variavel k: beta_k = C_T theta_k.
#   Design:        Z = W C, com W = [diag(X_1) | ... | diag(X_K)] (T x KT).
#                  Cada bloco de Z e Z_k = diag(X_k) C_T (T x T).
#
#   Ridge primal (Eq. 4):  min ||y - Z theta||^2 + lambda ||theta||^2
#       theta_hat = (Z'Z + lambda I)^-1 Z'y
#       beta_hat  = C theta_hat
#
#   Ridge dual (Eq. 9, mais rapido qd. KT > T):
#       alpha_hat = (ZZ' + lambda I_T)^-1 y          (NOTE: y, NAO Z'y)
#       beta_hat  = C Z' alpha_hat
#                  para a variavel k:
#                  beta_k = C_T * Z_k' * alpha = C_T * C_T' * (X_k o alpha)
#                  i.e. (beta_k)_t = sum_s min(t,s) X_{s,k} alpha_s.
#
#   Algorithm 1 (sec. 2.4, pag. 987):
#       Step 1: aprox. homogenea (Omega_u=sigma^2_u I, Omega_eps=sigma^2_eps I).
#               Obter beta_hat_1 via Eq. (9), com lambda escolhido por CV.
#       Step 2: ajustar volatilidade (GARCH ou rolling-var) nos residuos do
#               Step 1, obter sigma^2_eps,t. Normalizar para media 1.
#       Step 3: estimar sigma^2_u,k = (1/T) sum_t u_hat^2_{t,k}. Normalizar.
#       Step 4: estimacao final reweighted (Eq. 11), re-rodando CV. Equivale
#               a rodar o ridge dual em variaveis transformadas:
#                  X_tilde[t,k] = (sigma_{u,k}/sigma_{eps,t}) X[t,k]
#                  y_tilde[t]   = y[t]/sigma_{eps,t}
#               apos resolver, desfaz a transformacao: beta_k = sigma_{u,k} * beta_tilde_k
#
# VERIFICACAO: rodar .self_test() ao fim do arquivo verifica equivalencia
# numerica entre dual e primal em exemplo sintetico (tolerancia 1e-8).
# ==============================================================================


# ============================================================================ #
#  KERNEL DUAL E HELPERS BASICOS
# ============================================================================ #

# C_T: matriz T x T triangular inferior (com diagonal) de 1's.
make_C0 <- function(T_obs) {
  C0 <- matrix(0, T_obs, T_obs)
  C0[lower.tri(C0, diag = TRUE)] <- 1
  C0
}

# Constroi o kernel ZZ' (T x T). Operacao matematica:
#   ZZ' = sum_k Z_k Z_k', com Z_k = diag(X_k) C_T.
#   (Z_k Z_k')[i,j] = X_{i,k} X_{j,k} * min(i,j)
make_ZZt <- function(X) {
  T_obs <- nrow(X); K <- ncol(X)
  C0 <- make_C0(T_obs)
  ZZt <- matrix(0, T_obs, T_obs)
  for (k in 1:K) {
    Zk  <- X[, k] * C0                   # T x T: linha i = X[i,k] * C0[i,]
    ZZt <- ZZt + tcrossprod(Zk)
  }
  list(ZZt = ZZt, C0 = C0, K = K, T_obs = T_obs)
}

# Resolucao do sistema dual: (ZZ' + lambda I) alpha = y.
dual_solve <- function(ZZt, y, lam, eps = 1e-8) {
  T_obs <- nrow(ZZt)
  M <- ZZt + (lam + eps) * diag(T_obs)
  tryCatch(solve(M, y),
           error = function(e) solve(M + 1e-4 * diag(T_obs), y))
}

# Eigendecomp de ZZ' (PSD, simetrica). Usada para resolver multiplos lambdas
# em O(T^2) por lambda em vez de O(T^3).
eig_decomp <- function(ZZt) {
  ed <- eigen(ZZt, symmetric = TRUE)
  list(U = ed$vectors, d = pmax(ed$values, 0))
}

# Usando a eigendecomp: alpha = (ZZ' + lambda I)^-1 y = U diag(1/(d+lambda)) U' y.
solve_from_eig <- function(eig, y, lam, eps = 1e-8) {
  Uty <- as.numeric(crossprod(eig$U, y))
  as.numeric(eig$U %*% (Uty / (eig$d + lam + eps)))
}

# Recupera betas TVP a partir do alpha dual.
# Formula (derivada acima): beta_k = C_T C_T' * (X_k o alpha).
# Implementacao em O(T) por variavel via duas operacoes cumsum/reverse:
#   v       = X[,k] o alpha           (T x 1)
#   (C_T' v) = rev(cumsum(rev(v)))     (T x 1, soma "para tras")
#   beta_k  = cumsum(C_T' v)           (T x 1)
recover_beta <- function(X, alpha) {
  T_obs <- nrow(X); K <- ncol(X)
  beta <- matrix(NA_real_, T_obs, K)
  for (k in 1:K) {
    v <- X[, k] * alpha
    beta[, k] <- cumsum(rev(cumsum(rev(v))))
  }
  beta
}

# Versao rapida: somente beta_T. Identidade analitica derivada:
#   beta_T,k = (C_T C_T' v)_T = sum_s s * v[s] = sum_s s * X[s,k] alpha[s]
recover_beta_last <- function(X, alpha) {
  as.numeric(crossprod(X, alpha * seq_len(nrow(X))))
}


# ============================================================================ #
#  CROSS-VALIDATION DUAL COM EIGENDECOMP
# ============================================================================ #

# CV em blocos para series temporais. Cada fold faz UMA eigendecomp e itera
# sobre lambdas em O(T^2). Predicao OOS usa beta_T_train aplicado as observacoes
# de teste, replicando a logica de "ajustar com o passado, prever o futuro".
cv_ridge_dual <- function(X, y, lambdas, kfold = 5, block_size = NULL) {
  T_obs <- nrow(X)
  if (is.null(block_size)) block_size <- max(6, round(T_obs / kfold))
  folds <- rep(1:kfold, each = block_size, length.out = T_obs)

  cv_sse <- numeric(length(lambdas))
  n_pred <- 0L

  for (f in 1:kfold) {
    tr <- which(folds != f); te <- which(folds == f)
    if (length(tr) < 30 || length(te) < 2) next
    X_tr  <- X[tr, , drop = FALSE]; y_tr <- y[tr]
    zz    <- make_ZZt(X_tr)
    eig   <- eig_decomp(zz$ZZt)
    for (li in seq_along(lambdas)) {
      alpha <- solve_from_eig(eig, y_tr, lambdas[li])
      beta_T <- recover_beta_last(X_tr, alpha)
      pred_te <- as.numeric(X[te, , drop = FALSE] %*% beta_T)
      cv_sse[li] <- cv_sse[li] + sum((y[te] - pred_te)^2)
    }
    n_pred <- n_pred + length(te)
  }
  if (n_pred == 0L) return(lambdas[ceiling(length(lambdas) / 2)])
  lambdas[which.min(cv_sse / n_pred)]
}


# ============================================================================ #
#  STEP 1: 1SRR (RIDGE TVP HOMOGENEO)
# ============================================================================ #

tvp_1srr_standalone <- function(X, y, kfold = 5,
                                lambdas = exp(seq(-4, 20, length.out = 25))) {
  lam   <- cv_ridge_dual(X, y, lambdas, kfold)
  zz    <- make_ZZt(X)
  alpha <- dual_solve(zz$ZZt, y, lam)
  beta  <- recover_beta(X, alpha)
  yhat  <- rowSums(X * beta)
  list(beta = beta, resid = y - yhat, lambda = lam, alpha = alpha, yhat = yhat)
}


# ============================================================================ #
#  STEP 2 e 3: ESTIMACAO DE sigma^2_t (residuos) E sigma^2_u,k (parametros)
# ============================================================================ #

# sigma^2_t via GARCH(1,1). Coulombe recomenda GARCH no Step 2. Se falhar
# ou demorar demais, caimos para variancia rolling 12m. Normalizado a media 1.
# Rolling-variance fallback (deterministic, crash-safe) for Step 2.
.sigma2_rolling <- function(resid, window = 12) {
  T_obs <- length(resid); s <- rep(NA_real_, T_obs)
  for (t in seq_len(T_obs)) {
    w <- max(1, t - (window - 1)):t
    s[t] <- if (length(w) >= 2) var(resid[w], na.rm = TRUE) else NA_real_
  }
  bad <- !is.finite(s) | s <= 0
  if (all(bad)) return(rep(1, T_obs))
  if (any(bad)) s[bad] <- s[which(!bad)[1]]
  pmax(s, 1e-8) / mean(pmax(s, 1e-8))
}

# Persistent, crash-isolated GARCH back-end --------------------------------- #
# rugarch's "hybrid" solver cascades to gosolnp (random restarts) and, on
# Windows, INTERMITTENTLY crashes R at the C level (segmentation fault). The
# crash is not a catchable R error, so in a long-lived process (a parallel
# worker, or a serial loop over 180 windows) it eventually kills the session
# and leaves a whole 2SRR case unsaved (this is what removed 2SRR-AR).
# To keep the GARCH(1,1) step (faithful to Coulombe) AND guarantee the run
# completes, the rugarch fit is executed inside a persistent callr subprocess:
# if that subprocess crashes, callr reports it as an ordinary, catchable error
# in the parent, we drop the dead session, and the affected window falls back
# to a rolling variance. The parent process never dies.
.garch_state <- new.env(parent = emptyenv())

.garch_session <- function() {
  rs <- .garch_state$session
  if (!is.null(rs) && inherits(rs, "r_session") &&
      tryCatch(rs$is_alive(), error = function(e) FALSE)) return(rs)
  rs <- tryCatch(callr::r_session$new(), error = function(e) NULL)
  if (!is.null(rs))
    tryCatch(rs$run(function() { suppressMessages(requireNamespace("rugarch")); TRUE }),
             error = function(e) NULL)
  .garch_state$session <- rs
  rs
}

# Close the GARCH subprocess (call at the end of a run; optional — it also dies
# automatically when the R session exits).
garch_session_close <- function() {
  rs <- .garch_state$session
  if (!is.null(rs)) try(rs$kill(), silent = TRUE)
  .garch_state$session <- NULL
  invisible(NULL)
}

# Self-contained GARCH(1,1) fit (runs in the subprocess; uses only its args).
.garch_fit_fun <- function(resid, solver) {
  spec <- rugarch::ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
    mean.model     = list(armaOrder = c(0, 0), include.mean = TRUE))
  fit <- rugarch::ugarchfit(spec, resid, solver = solver,
                            solver.control = list(trace = 0))
  as.numeric(rugarch::sigma(fit))^2
}

# Time-varying residual variance for Step 2 of the 2SRR. solver defaults to
# "hybrid" (the rugarch default, matching the Factor/FAVAR runs); isolation is
# on by default. Set options(garch_isolate = FALSE) to run rugarch in-process.
estimate_sigma2_standalone <- function(resid, window = 12,
                                        solver  = getOption("garch_solver",  "hybrid"),
                                        isolate = getOption("garch_isolate", TRUE), ...) {
  T_obs <- length(resid)
  if (!requireNamespace("rugarch", quietly = TRUE)) return(.sigma2_rolling(resid, window))

  out <- NULL
  if (isTRUE(isolate) && requireNamespace("callr", quietly = TRUE)) {
    rs <- .garch_session()
    if (!is.null(rs))
      out <- tryCatch(
        rs$run(.garch_fit_fun, list(resid = as.numeric(resid), solver = solver)),
        error = function(e) { try(rs$kill(), silent = TRUE)
                              .garch_state$session <- NULL; NULL })
  } else {
    out <- tryCatch(.garch_fit_fun(resid, solver), error = function(e) NULL)
  }
  if (is.null(out) || length(out) != T_obs || any(!is.finite(out)) || all(out <= 0))
    return(.sigma2_rolling(resid, window))
  pmax(out, 1e-8) / mean(pmax(out, 1e-8))
}

# sigma^2_u,k = (1/T) sum_t (delta beta_t,k)^2. Normalizado a media 1.
estimate_omega_standalone <- function(beta) {
  omega <- colMeans(diff(beta)^2, na.rm = TRUE)
  pmax(omega, 1e-12) / mean(pmax(omega, 1e-12))
}


# ============================================================================ #
#  ENGINE: 2SRR STANDALONE (Algorithm 1 completo)
# ============================================================================ #

tvp_2srr_standalone <- function(X, y, kfold = 5,
                                lambdas = exp(seq(-4, 20, length.out = 25)),
                                eps = 1e-8) {
  # Step 1
  s1 <- tvp_1srr_standalone(X, y, kfold, lambdas)

  # Step 2, 3
  sigma2 <- estimate_sigma2_standalone(s1$resid)
  omega  <- estimate_omega_standalone(s1$beta)

  # Step 4: variaveis transformadas
  isig <- 1 / sqrt(pmax(sigma2, eps))      # 1/sigma_{eps,t}
  som  <- sqrt(pmax(omega, eps))            # sigma_{u,k}
  X_t  <- sweep(sweep(X, 2, som, "*"), 1, isig, "*")
  y_t  <- y * isig

  lam2  <- cv_ridge_dual(X_t, y_t, lambdas, kfold)
  zz2   <- make_ZZt(X_t)
  a2    <- dual_solve(zz2$ZZt, y_t, lam2)
  bt    <- recover_beta(X_t, a2)
  beta  <- sweep(bt, 2, som, "*")            # desfaz transformacao A_theta
  yhat  <- rowSums(X * beta)

  list(beta         = beta,
       resid        = y - yhat,
       lambda       = lam2,
       omega        = omega,
       sigma2       = sigma2,
       lambda_step1 = s1$lambda,
       beta_step1   = s1$beta,
       yhat         = yhat)
}


# ============================================================================ #
#  ENGINE PRINCIPAL: COULOMBE-FAST
#
#  Implementacao RAPIDA do TVPRR_cosso do Coulombe (type=2). Reproduz a EXATA
#  matematica do paper (sec. 2.4 + Algorithm 1 + dualGRR), mas evita inflar
#  Zprime para T x KT. Usa fato: o kernel ZZ' de Coulombe tem forma fechada
#  somavel em T x T:
#
#      Kmat_du(lambda1, lambda2) = (1/lambda1) * sum_k sw[k] X_k X_k' o M_inn
#                                 + (1/lambda2) X_full X_full'                (*)
#
#  onde M_inn[i,j] = min(i,j) - 1, X_full = [1 | X_std], sw[k] = sigma^2_u,k.
#  As betas se reconstroem por
#      uhat_innov_k[s] = (sw[k]/lambda1) * sum_{t>s} X_full[t,k] alpha[t]
#      uhat_int_k      = (1/lambda2)     * sum_t X_full[t,k] alpha[t]
#      beta_t,k        = uhat_int_k + sum_{s<t} uhat_innov_k[s]
#  (com correcao OLS-prior identica a dualGRR do Coulombe quando olsprior=1).
#
#  CV sobre lambda1 (lambda2 fixo) usa decomposicao generalizada via Cholesky:
#      Q = (1/lambda2) X_full X_full' + diag(eweights)
#      M(lambda1) = (1/lambda1) K_innov + Q
#      L = chol(Q), H = L^-1 K_innov L^-T, H = U D U'
#      M(lambda1)^-1 y = L^-T U ((1/lambda1) D + I)^-1 U' L^-1 y
#
#  Custo por fold: O(T^3) decomp + O(T^2) por lambda. Total ~ 3 horas em 180
#  janelas x 3 casos x 4 horizontes.
# ============================================================================ #

# Constroi K_innov = sum_k sweights[k] (X_full[,k] X_full[,k]') o M_inn
# onde M_inn[i,j] = min(i,j) - 1.
build_K_innov <- function(X_full, sweights, M_inn = NULL) {
  T_obs <- nrow(X_full); Kp <- ncol(X_full)
  if (is.null(M_inn))
    M_inn <- outer(1:T_obs, 1:T_obs, function(i, j) pmin(i, j) - 1)
  K_innov <- matrix(0, T_obs, T_obs)
  for (k in 1:Kp)
    K_innov <- K_innov + sweights[k] * (X_full[, k] %o% X_full[, k]) * M_inn
  K_innov
}

# Reconstroi betas (T x Kp) a partir de alpha (T-vec), via formulas acima.
# Se olsprior=1, soma beta_ols no bloco constante (igual ao dualGRR original).
recover_beta_coulombe <- function(X_full, alpha, lambda1, lambda2, sweights,
                                    beta_ols = NULL) {
  T_obs <- nrow(X_full); Kp <- ncol(X_full)
  beta  <- matrix(NA_real_, T_obs, Kp)
  for (k in 1:Kp) {
    v        <- X_full[, k] * alpha
    rev_csum <- rev(cumsum(rev(v)))        # rev_csum[s] = sum_{t>=s} v[t]
    excl     <- c(rev_csum[-1], 0)         # excl[s] = sum_{t > s} v[t]
    uhat_inn <- (sweights[k] / lambda1) * excl[1:(T_obs - 1)]
    uhat_int <- (1 / lambda2) * sum(v)
    if (!is.null(beta_ols)) uhat_int <- uhat_int + beta_ols[k]
    beta[, k] <- uhat_int + c(0, cumsum(uhat_inn))
  }
  beta
}

# Fit unico (lambda1 fixo), equivalente a dualGRR(type=2) do Coulombe.
dualGRR_fast <- function(X_full, y, lambda1, lambda2 = 0.1,
                          sweights = NULL, eweights = NULL,
                          olsprior = 1, M_inn = NULL) {
  T_obs <- nrow(X_full); Kp <- ncol(X_full)
  if (is.null(sweights) || length(sweights) == 1) sweights <- rep(1, Kp)
  if (is.null(eweights) || length(eweights) == 1) eweights <- rep(1, T_obs)

  K_innov <- build_K_innov(X_full, sweights, M_inn)
  K_int   <- tcrossprod(X_full)
  Kmat_du <- (1 / lambda1) * K_innov + (1 / lambda2) * K_int
  Lambda_T <- diag(as.numeric(eweights))

  beta_ols <- NULL
  if (olsprior == 1) {
    # OLS no bloco constante (X_full), com penalty 1e-7 para estabilidade
    beta_ols <- as.numeric(solve(crossprod(X_full) + 1e-7 * diag(Kp),
                                  crossprod(X_full, y)))
    y_resid  <- y - X_full %*% beta_ols
    alpha    <- as.numeric(solve(Kmat_du + Lambda_T, y_resid))
  } else {
    alpha    <- as.numeric(solve(Kmat_du + Lambda_T, y))
  }

  beta <- recover_beta_coulombe(X_full, alpha, lambda1, lambda2, sweights, beta_ols)
  yhat <- rowSums(X_full * beta)
  list(beta = beta, alpha = alpha, yhat = yhat,
       lambdas = c(lambda1, lambda2), sweights = sweights)
}

# CV rapida sobre lambda1, com decomposicao generalizada via Cholesky.
# Reproduz a CV de cvgs.bhk2015 (Frisch-Waugh) com mesma logica de blocos.
cv_dualGRR_fast <- function(X_full, y, lambda1_vec, lambda2 = 0.1,
                              sweights = NULL, eweights = NULL,
                              olsprior = 1, kfold = 5) {
  T_obs <- nrow(X_full); Kp <- ncol(X_full)
  if (is.null(sweights) || length(sweights) == 1) sweights <- rep(1, Kp)
  if (is.null(eweights) || length(eweights) == 1) eweights <- rep(1, T_obs)

  block_size <- max(6, round(T_obs / kfold))
  folds <- rep(1:kfold, each = block_size, length.out = T_obs)

  cv_sse <- numeric(length(lambda1_vec))
  n_pred <- 0L

  for (f in 1:kfold) {
    tr <- which(folds != f); te <- which(folds == f)
    if (length(tr) < 30 || length(te) < 2) next

    X_tr <- X_full[tr, , drop = FALSE]
    y_tr <- y[tr]; T_tr <- length(tr)

    K_innov_tr <- build_K_innov(X_tr, sweights)
    K_int_tr   <- tcrossprod(X_tr)
    Q_tr       <- (1 / lambda2) * K_int_tr + diag(eweights[tr])

    # Cholesky generalizada: Q = L L'; H = L^-1 K_innov L^-T; H = U D U'
    Lc <- tryCatch(chol(Q_tr + 1e-8 * diag(T_tr)),
                    error = function(e) { warning("chol failed"); NULL })
    if (is.null(Lc)) next
    L     <- t(Lc)                                       # lower
    Linv  <- backsolve(Lc, diag(T_tr), transpose = TRUE) # L^-1 = L_low^-1
    H     <- Linv %*% K_innov_tr %*% t(Linv)
    eig   <- eigen((H + t(H)) / 2, symmetric = TRUE)
    U     <- eig$vectors;  d <- pmax(eig$values, 0)

    # OLS prior no bloco constante (para os folds de treino)
    beta_ols_tr <- NULL
    if (olsprior == 1) {
      beta_ols_tr <- as.numeric(solve(crossprod(X_tr) + 1e-7 * diag(Kp),
                                       crossprod(X_tr, y_tr)))
      y_eff <- y_tr - X_tr %*% beta_ols_tr
    } else y_eff <- y_tr
    Liny <- as.numeric(Linv %*% y_eff)
    UtLy <- as.numeric(crossprod(U, Liny))

    # Para cada lambda1: alpha = L^-T U ((1/lambda1) D + I)^-1 U' L^-1 y
    for (li in seq_along(lambda1_vec)) {
      lam1   <- lambda1_vec[li]
      coefs  <- UtLy / ((1 / lam1) * d + 1)
      alpha  <- as.numeric(t(Linv) %*% (U %*% coefs))
      beta_tr <- recover_beta_coulombe(X_tr, alpha, lam1, lambda2, sweights,
                                        beta_ols_tr)
      # Predicao no holdout: usa beta no instante T_tr (mais recente do treino)
      beta_T <- beta_tr[T_tr, ]
      pred_te <- as.numeric(X_full[te, , drop = FALSE] %*% beta_T)
      cv_sse[li] <- cv_sse[li] + sum((y[te] - pred_te)^2)
    }
    n_pred <- n_pred + length(te)
  }
  if (n_pred == 0L) return(lambda1_vec[ceiling(length(lambda1_vec) / 2)])
  lambda1_vec[which.min(cv_sse / n_pred)]
}

# Engine principal: Algorithm 1 completo, FIEL ao TVPRR_cosso(type=2).
tvp_2srr_coulombe_fast <- function(X, y, kfold = 5,
                                     lambdavec = exp(seq(-2, 12, length.out = 15)),
                                     lambda2 = 0.1,
                                     homo.param = 0.75,
                                     sv.param   = 0.75,
                                     olsprior   = 1,
                                     eps        = 1e-8) {
  T_obs <- nrow(X); K <- ncol(X)

  # 1. Padronizacao (igual TVPRR_cosso: dividir cada coluna pelo proprio sd,
  #    e y por sd(y); ao final, re-escalar betas).
  sdy           <- sd(y);            if (sdy < eps) sdy <- 1
  X_sd          <- apply(X, 2, sd);  X_sd[X_sd < eps] <- 1
  scalingfactor <- sdy / X_sd
  X_std         <- sweep(X, 2, X_sd, "/")
  y_std         <- y / sdy

  # 2. Adiciona intercepto, monta X_full = [1 | X_std], T x (K+1)
  X_full <- cbind(intercept = 1, X_std)
  Kp     <- ncol(X_full)

  # Step 1: homogeneous CV (sweights=1, eweights=1)
  lam1_s1   <- cv_dualGRR_fast(X_full, y_std, lambdavec, lambda2,
                                 sweights = rep(1, Kp),
                                 eweights = rep(1, T_obs),
                                 olsprior = olsprior, kfold = kfold)
  fit_s1    <- dualGRR_fast(X_full, y_std, lam1_s1, lambda2,
                              sweights = rep(1, Kp),
                              eweights = rep(1, T_obs), olsprior = olsprior)
  betas_s1  <- fit_s1$beta
  resid_s1  <- y_std - fit_s1$yhat

  # Step 2: sigma^2_eps,t via GARCH (com queda para rolling-var se falhar)
  sigma2_eps <- estimate_sigma2_standalone(resid_s1)
  EW <- sigma2_eps ^ sv.param;  EW <- EW / mean(EW)

  # Step 3: sigma^2_u,k via media dos quadrados das diferencas de beta
  umat       <- diff(betas_s1)               # (T-1) x Kp
  sigma2_u   <- colMeans(umat^2) ^ homo.param
  sigma2_u   <- sigma2_u / mean(sigma2_u)
  sigma2_u   <- pmax(sigma2_u, eps)          # estabilidade

  # Step 4: CV ponderada + fit final
  lam1_s4    <- cv_dualGRR_fast(X_full, y_std, lambdavec, lambda2,
                                  sweights = sigma2_u, eweights = EW,
                                  olsprior = olsprior, kfold = kfold)
  fit_s4     <- dualGRR_fast(X_full, y_std, lam1_s4, lambda2,
                               sweights = sigma2_u, eweights = EW,
                               olsprior = olsprior)

  # 3. Re-escala betas: beta_intercept * sdy, beta_j * scalingfactor[j]
  beta_unscaled <- matrix(NA_real_, T_obs, Kp)
  beta_unscaled[, 1] <- fit_s4$beta[, 1] * sdy
  for (j in 1:K) beta_unscaled[, j + 1] <- fit_s4$beta[, j + 1] * scalingfactor[j]
  yhat_unscaled <- rowSums(cbind(1, X) * beta_unscaled)

  list(beta         = beta_unscaled,
       resid        = y - yhat_unscaled,
       lambda       = lam1_s4,
       lambda_step1 = lam1_s1,
       omega        = sigma2_u,
       sigma2       = sigma2_eps,
       beta_step1   = betas_s1,
       yhat         = yhat_unscaled,
       lambda2      = lambda2)
}


# ============================================================================ #
#  ENGINE LEGADO: COULOMBE ORIGINAL (wrapper para TVPRR_cosso, sanity check)
# ============================================================================ #

tvp_2srr_coulombe <- function(X, y, kfold = 5,
                              lambdavec = exp(pracma::linspace(-2, 12, n = 15)),
                              lambda2 = 0.1) {
  if (!exists("TVPRR_cosso", mode = "function"))
    stop("TVPRR_cosso nao carregado. Rode 00_setup.R primeiro.")
  # IMPORTANTE: sweigths em TVPRR_cosso/dualGRR e' POR-PARAMETRO (length nf = K+1),
  # NAO por-observacao. Default scalar 1 e' expandido internamente para length nf.
  # Passar rep(1, nrow(X)) (length T) causava "indice fora de limites" no loop
  # do dualGRR que faz Kmat_half[begin:end,] = sweigths[m] * (...).
  res <- TVPRR_cosso(
    X = X, y = y, type = 2,
    lambdavec = lambdavec, lambda2 = lambda2,
    kfold = kfold,
    sweigths = 1,                    # scalar -> expandido para length nf=K+1
    silent = 1, tol = 1e-3, maxit = 5)

  # TVPRR_cosso (type=2) devolve betas em res$grrats$betas_grr (1 x (K+1) x T)
  # (K+1 porque TVPRR_cosso adiciona um intercepto internamente; primeira coluna
  #  de betas e' o TVP intercept)
  betas <- if (!is.null(res$grrats$betas_grr))
             res$grrats$betas_grr
           else if (!is.null(res$grr$betas_grr))
             res$grr$betas_grr
           else NULL
  if (is.null(betas)) return(NULL)

  # Reformatar para T x (K+1), consistente com coulombe_fast
  Kp <- dim(betas)[2]; T_obs <- dim(betas)[3]
  beta_TxKp <- matrix(NA_real_, T_obs, Kp)
  for (k in 1:Kp) beta_TxKp[, k] <- betas[1, k, ]

  # yhat = sum_k X_full[, k] * beta_t,k, com X_full = [1 | X]
  X_full <- cbind(1, X)
  yhat <- rowSums(X_full * beta_TxKp)
  list(beta   = beta_TxKp,
       resid  = y - yhat,
       lambda = res$lambda1,
       yhat   = yhat)
}


# ============================================================================ #
#  OUTLIER FILTER (Coulombe, sec. 4 do paper)
# ============================================================================ #

outlier_filter <- function(pred, y_is, pred_fallback) {
  y_mean <- mean(y_is, na.rm = TRUE)
  y_max  <- max(y_is,  na.rm = TRUE)
  y_min  <- min(y_is,  na.rm = TRUE)
  bound  <- 2 * max(abs(y_max - y_mean), abs(y_min - y_mean))
  if (!is.finite(pred) || abs(pred - y_mean) > bound) return(pred_fallback)
  pred
}


# ============================================================================ #
#  CONSTRUCAO DO DESIGN MATRIX PARA OS 3 CASOS TVP
# ============================================================================ #

build_lags <- function(mat, start_lag, n_lags) {
  if (n_lags <= 0) return(matrix(0, nrow = nrow(mat), ncol = 0))
  mat <- as.matrix(mat)
  T_obs <- nrow(mat); p <- ncol(mat)
  out <- matrix(NA_real_, T_obs, n_lags * p)
  cn  <- character(n_lags * p)
  base_names <- if (is.null(colnames(mat))) paste0("v", 1:p) else colnames(mat)
  for (i in 1:n_lags) {
    lag_i <- start_lag + i - 1
    out[(lag_i + 1):T_obs, ((i - 1) * p + 1):(i * p)] <-
      mat[1:(T_obs - lag_i), , drop = FALSE]
    cn[((i - 1) * p + 1):(i * p)] <- paste0(base_names, "_L", lag_i)
  }
  colnames(out) <- cn
  out
}

# Casos TVP (Hudson Request):
#   "AR"    : univar=T, nofact=T   -> apenas lags de Y_h (mais intercept)
#   "Factor": factonly=T           -> apenas lags dos fatores PCA
#   "FAVAR" : padrao Coulombe      -> lags de Y_h + lags dos fatores
# Target em todos: a TAXA de inflacao h passos a frente (pi_{t+h}), igual ao
# 01_data_prep.R e ao 02_forecast_medeiros.R, garantindo comparabilidade direta.
build_design_tvp <- function(y_is, X_is_raw, h, case = c("AR", "Factor", "FAVAR"),
                              ly = 2, lf = 2, nf = 4, include_intercept = TRUE) {
  case  <- match.arg(case)
  T_obs <- length(y_is)

  # Target: TAXA de inflacao h passos a frente (direct multi-step). Com os
  # preditores defasados em h (start_lag = h, abaixo), a regressao in-sample e
  #   y[t] ~ y[t-h], y[t-h-1], ..., F[t-h], ...  e a linha OOS preve y[T_obs+h].
  # Isso casa com yout[i,h] = pi_{t+h} (01_data_prep.R) e com as previsoes
  # diretas do rolling_window do Medeiros. (Antes: soma movel = acumulado, o que
  # punha o 2SRR numa escala diferente do alvo e inflava o RMSE em h>1.)
  y_h <- as.numeric(y_is)

  # Fatores PCA (se aplicavel)
  factors_full <- NULL; nf_eff <- 0
  if (case != "AR" && !is.null(X_is_raw) && ncol(X_is_raw) > 0) {
    Xc <- X_is_raw; Xc[is.na(Xc)] <- 0
    pc <- tryCatch(prcomp(Xc, center = FALSE, scale. = FALSE),
                   error = function(e) NULL)
    if (!is.null(pc)) {
      nf_eff <- min(nf, ncol(pc$x))
      factors_full <- pc$x[, 1:nf_eff, drop = FALSE]
      colnames(factors_full) <- paste0("F", 1:nf_eff)
    }
  }

  ly_eff <- if (case == "Factor") 0 else ly
  lags_y <- build_lags(matrix(y_h, ncol = 1, dimnames = list(NULL, "Yh")),
                       start_lag = h, n_lags = ly_eff)
  lags_f <- if (!is.null(factors_full))
              build_lags(factors_full, start_lag = h, n_lags = lf)
            else matrix(0, nrow = T_obs, ncol = 0)

  X_full <- cbind(lags_y, lags_f)
  if (ncol(X_full) == 0) stop("Design vazio para case=", case)
  if (include_intercept) X_full <- cbind(intercept = 1, X_full)

  # Trim NAs (proveniente dos lags)
  ok   <- complete.cases(X_full) & !is.na(y_h)
  X_in <- X_full[ok, , drop = FALSE]
  y_in <- y_h[ok]

  # Regressor OOS (instante T+h): aplica modelo treinado a (Y_h(T), Y_h(T-1), F(T), ...)
  x_out <- numeric(ncol(X_full))
  names(x_out) <- colnames(X_full)
  pos <- 1
  if (include_intercept) { x_out[pos] <- 1; pos <- pos + 1 }
  if (ly_eff > 0) {
    x_out[pos:(pos + ly_eff - 1)] <- y_h[T_obs:(T_obs - ly_eff + 1)]
    pos <- pos + ly_eff
  }
  if (!is.null(factors_full) && lf > 0) {
    for (i in 1:lf) {
      block <- (pos + (i - 1) * nf_eff):(pos + i * nf_eff - 1)
      x_out[block] <- factors_full[T_obs - i + 1, ]
    }
  }

  list(X_in = X_in, y_in = y_in, x_out = x_out,
       var_names = colnames(X_full),
       nf_eff = nf_eff, case = case, h = h, T_in = nrow(X_in))
}


# ============================================================================ #
#  WRAPPER: AJUSTA 2SRR (UM CASO, UM HORIZONTE, UMA JANELA)
# ============================================================================ #

fit_2srr_window <- function(y_is, X_is_raw, h, case = "FAVAR",
                             ly = 2, lf = 2, nf = 4, kfold = 5,
                             lambdas = exp(seq(-4, 20, length.out = 20)),
                             engine = c("coulombe_fast", "standalone", "coulombe"),
                             lambda2 = 0.1,
                             coulombe_lambdavec = NULL) {
  engine <- match.arg(engine)

  design <- tryCatch(
    build_design_tvp(y_is, X_is_raw, h, case, ly, lf, nf,
                     include_intercept = FALSE),  # intercept entra dentro da engine
    error = function(e) NULL)
  if (is.null(design) || design$T_in < 30 || ncol(design$X_in) < 2) {
    return(list(forecast = NA_real_, ridge_forecast = NA_real_,
                betas_tvp = NULL, lambda = NA_real_, omega = NULL,
                sigma2 = NULL, var_names = NULL, engine = engine,
                case = case, h = h, status = "design_failed"))
  }

  X_in  <- design$X_in;  y_in <- design$y_in;  x_out <- design$x_out
  # design foi construido SEM intercepto. Cada engine adiciona o que precisa:
  # - coulombe_fast: intercepto entra como 1a coluna do X_full, com lambda2
  # - standalone:    intercepto entra como 1a coluna, com mesmo lambda dos demais
  # - coulombe (original TVPRR_cosso): adiciona internamente

  lvec <- if (is.null(coulombe_lambdavec))
            exp(pracma::linspace(-2, 12, n = 15))
          else coulombe_lambdavec

  fit <- switch(engine,
    "coulombe_fast" = tryCatch(
      tvp_2srr_coulombe_fast(X_in, y_in, kfold = kfold,
                              lambdavec = lvec, lambda2 = lambda2),
      error = function(e) { message("  coulombe_fast falhou: ", e$message); NULL }),
    "coulombe" = tryCatch(
      tvp_2srr_coulombe(X_in, y_in, kfold = kfold, lambdavec = lvec,
                         lambda2 = lambda2),
      error = function(e) { message("  Coulombe orig falhou: ", e$message); NULL }),
    "standalone" = {
      X_std_int <- cbind(intercept = 1, X_in)
      tryCatch(
        tvp_2srr_standalone(X_std_int, y_in, kfold = kfold, lambdas = lambdas),
        error = function(e) { message("  Standalone falhou: ", e$message); NULL })
    }
  )

  # Ridge baseline (Step 1 do Coulombe / baseline Medeiros-Ridge equivalente).
  # glmnet adiciona seu proprio intercepto, entao passamos X_in sem intercepto.
  ridge_fc <- tryCatch({
    cv_r <- glmnet::cv.glmnet(X_in, y_in, alpha = 0,
                              nfolds = min(kfold, floor(nrow(X_in) / 3)))
    as.numeric(predict(cv_r, newx = matrix(x_out, nrow = 1), s = "lambda.min"))
  }, error = function(e) mean(y_in, na.rm = TRUE))

  if (is.null(fit)) {
    return(list(forecast = ridge_fc, ridge_forecast = ridge_fc,
                betas_tvp = NULL, lambda = NA_real_, omega = NULL,
                sigma2 = NULL, var_names = design$var_names, engine = engine,
                case = case, h = h, status = "fit_failed_fallback_ridge"))
  }

  # fit$beta tem dim T x (K+1): primeira coluna = intercepto, outras = regressores.
  # x_out tem dim K (sem intercepto), entao prefixamos 1 para o forecast.
  beta_T   <- fit$beta[nrow(fit$beta), ]
  x_out_full <- c(1, x_out)
  fcast    <- sum(x_out_full * beta_T)
  fcast_f  <- outlier_filter(fcast, y_in, ridge_fc)

  # Nomes das variaveis incluindo intercepto (para diagnostico posterior)
  var_names_full <- c("intercept", design$var_names)

  list(
    forecast       = fcast_f,
    forecast_raw   = fcast,
    ridge_forecast = ridge_fc,
    betas_tvp      = fit$beta,
    lambda         = fit$lambda,
    lambda_step1   = if (!is.null(fit$lambda_step1)) fit$lambda_step1 else NA,
    lambda2        = if (!is.null(fit$lambda2)) fit$lambda2 else NA,
    omega          = fit$omega,
    sigma2         = fit$sigma2,
    var_names      = var_names_full,
    engine         = engine,
    case           = case, h = h,
    n_obs          = design$T_in, n_vars = ncol(fit$beta),
    status         = "ok"
  )
}


# ============================================================================ #
#  AUTO-TESTE: equivalencia dual <-> primal direto (Eq. 5 do paper)
# ============================================================================ #

# Verifica que para um exemplo sintetico (T=30, K=3), os betas recuperados
# via dual (make_ZZt + dual_solve + recover_beta) batem com o primal direto
# theta = (Z'Z + lam I)^-1 Z'y, beta = C theta. Se a diferenca max > 1e-6,
# emite warning. Garante a correcao matematica da implementacao standalone.
tvp_self_test <- function(seed = 7, T_obs = 30, K = 3, lam = 1, verbose = TRUE) {
  set.seed(seed)
  X <- matrix(rnorm(T_obs * K), T_obs, K)
  y <- rnorm(T_obs)

  # Primal direto: theta_hat = (Z'Z + lam I)^-1 Z'y, beta = (I_K kron C) theta
  C   <- make_C0(T_obs)
  W   <- do.call(cbind, lapply(1:K, function(k) diag(X[, k])))           # T x KT
  IkC <- kronecker(diag(K), C)                                            # KT x KT
  Z   <- W %*% IkC                                                        # T x KT
  theta_p <- as.numeric(solve(crossprod(Z) + lam * diag(K * T_obs),
                              crossprod(Z, y)))
  beta_p  <- matrix(NA_real_, T_obs, K)
  for (k in 1:K) {
    th_k <- theta_p[((k - 1) * T_obs + 1):(k * T_obs)]
    beta_p[, k] <- as.numeric(C %*% th_k)
  }

  # Dual via standalone
  zz    <- make_ZZt(X)
  alpha <- dual_solve(zz$ZZt, y, lam)
  beta_d <- recover_beta(X, alpha)

  diff <- max(abs(beta_p - beta_d))
  if (verbose) {
    cat(sprintf("[tvp_self_test] T=%d K=%d lam=%.2f | max|beta_primal - beta_dual| = %.3e\n",
                T_obs, K, lam, diff))
  }
  if (diff > 1e-6)
    warning(sprintf("Equivalencia dual<->primal FALHOU: diff=%.3e", diff))
  invisible(diff)
}

# Roda o auto-teste automaticamente na carga (silencioso). Se descomentado,
# o sourceamento aborta com warning visivel caso a matematica esteja errada.
# Recomendo deixar comentado para nao bagular o run em producao.
# tvp_self_test()


# Auto-teste 2: equivalencia entre dualGRR_fast (engine coulombe_fast) e a
# implementacao ORIGINAL dualGRR de Coulombe, com mesmos lambda1/lambda2 e
# pesos. Se dualGRR existir (sourceado pelo 00_setup) e a matematica
# do fast estiver correta, betas devem coincidir ate precisao numerica.
tvp_coulombe_fast_vs_original_test <- function(seed = 11, T_obs = 25, K = 4,
                                                lambda1 = 5, lambda2 = 0.1,
                                                tol = 1e-4, verbose = TRUE) {
  if (!exists("dualGRR", mode = "function") ||
      !exists("Zfun", mode = "function")) {
    if (verbose) cat("[coulombe_fast_test] dualGRR/Zfun nao carregados — skip\n")
    return(invisible(NA))
  }
  set.seed(seed)
  X <- matrix(rnorm(T_obs * K), T_obs, K)
  y <- as.numeric(rnorm(T_obs))

  # ORIGINAL: Zfun + dualGRR
  Zp <- Zfun(X)                                   # T x (K+1)*T
  res_orig <- dualGRR(Zprime = Zp, y = y, dimX = K + 1,
                      lambda1 = lambda1, lambda2 = lambda2,
                      olsprior = 1, sweigths = rep(1, K + 1), eweigths = 1,
                      calcul_beta = 1)
  # betas_grr orig: array 1 x (K+1) x T -> reformata para T x (K+1)
  beta_orig <- t(res_orig$betas_grr[1, , ])

  # FAST: dualGRR_fast com mesmos parametros
  X_full <- cbind(intercept = 1, X)
  res_fast <- dualGRR_fast(X_full, y, lambda1 = lambda1, lambda2 = lambda2,
                            sweights = rep(1, K + 1), eweights = rep(1, T_obs),
                            olsprior = 1)
  beta_fast <- res_fast$beta

  diff <- max(abs(beta_orig - beta_fast))
  if (verbose) {
    cat(sprintf("[coulombe_fast_test] T=%d K=%d lam1=%.2f lam2=%.2f | max|orig - fast| = %.3e %s\n",
                T_obs, K, lambda1, lambda2, diff,
                ifelse(diff < tol, "[OK]", "[FALHA]")))
  }
  invisible(diff)
}

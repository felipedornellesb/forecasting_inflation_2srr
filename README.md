# Inflation Forecasting with Two-Step Ridge Regression

**Two-Step Ridge Regression** (2SRR; Goulet Coulombe, *International Journal of
Forecasting*, 2025) applied to U.S. CPI inflation forecasting (CPIAUCSL,
FRED-MD), benchmarked against twelve econometric and machine-learning models —
**Ridge, LASSO, Elastic Net, Adaptive LASSO, Adaptive Elastic Net, Random
Forest, Bagging, Factor, Target Factor, Complete Subset Regression, AR and
AR-BIC** — from the Medeiros et al. (*JBES*, 2021) framework.

Reproducible R pipeline: from the FRED-MD panel to forecasts, formal
predictive-accuracy tests, and publication-ready figures and tables.

---

## Pipeline overview

```
01_data_prep.R         -> yout.rda (realized cumulative target) + rw.rda (random-walk benchmark)
02_forecast_medeiros.R -> 12 econometric / ML benchmarks (Medeiros et al., 2021)
03_forecast_2srr.R     -> 2SRR in 3 TVP specifications (AR, Factor, FAVAR)
04_analysis.R          -> RMSE, Diebold-Mariano, Giacomini-White, Model Confidence Set,
                          Mincer-Zarnowitz; beta trajectories; figures + CSV tables
05_figures_tables.R    -> publication-ready figures (PNG/PDF) + Word (.docx) tables
```

Run the scripts in order. `04_analysis.R` is robust to missing inputs — any
benchmark model whose forecast file does not exist is silently skipped.

---

## How to run

From the **project root**, in R:

```r
source("00_prog/00_setup.R")             # install packages; download Coulombe + Medeiros functions
source("00_prog/01_data_prep.R")         # build yout.rda and rw.rda (180 windows x 12 horizons)
source("00_prog/02_forecast_medeiros.R") # 12 benchmark models (rolling window)
source("00_prog/03_forecast_2srr.R")     # 2SRR via coulombe_fast (parallel over windows)
source("00_prog/04_analysis.R")          # tests, figures and CSV tables
source("00_prog/05_figures_tables.R")    # final figures (PNG/PDF) and Word tables
```

Prerequisite: place `data.rda` in `10_data/`. Target variable: `CPIAUCSL`.

---

## Changing the penalty grid

The cross-validation penalty grid for 2SRR is defined in the **Hyperparameters
block of `00_prog/03_forecast_2srr.R`**. To test other grids, edit:

```r
n_lambda    <- 15
lambda_vec  <- exp(pracma::linspace(-2, 12, n = n_lambda))   # active grid
```

A wider grid is provided just above this line, commented out, for experiments.
The default `exp(linspace(-2, 12))` (15 points) is intentional — see the grid
decision below.

---

## Design decisions

### 2SRR — main engine: `coulombe_fast`

The original `TVPRR_cosso` materializes the `T × KT` kernel explicitly. For
monthly FRED-MD with `T ≈ 600+`, that explodes in memory and time. The
`coulombe_fast` engine implements **the same mathematics** in closed form on the
`T × T` kernel:

$$K_{du}(\lambda_1, \lambda_2) = \frac{1}{\lambda_1} \sum_k \omega_k\, (X_k X_k') \odot M_{inn} + \frac{1}{\lambda_2} X_{full} X_{full}'$$

where `M_inn[i,j] = min(i,j) - 1`. Cross-validation via a generalized Cholesky
decomposition sweeps `n_lambda = 15` in O(T²) per λ — about 10–50× faster, with
no mathematical loss. Two self-tests at the start of `03` (standalone vs. direct
primal; `coulombe_fast` vs. the original `dualGRR`) abort the run if they fail.

### Penalty-grid decision

The default grid is `lambda_vec = exp(linspace(-2, 12, 15))`. Cross-validation
on this grid often saturates at the maximum value for monthly data, yet
empirical testing shows that an **expanded** grid (allowing larger λ) *degrades*
out-of-sample accuracy. The restricted grid acts as an implicit regularization
that prevents cross-validation from eliminating the time-varying structure,
preserving the residual parameter variation that lets 2SRR outperform a
constant-parameter ridge.

### Three TVP specifications

| Case | Regressors |
|------|------------|
| 2SRR-AR | `ly` lags of cumulative inflation only |
| 2SRR-Factor | `lf` lags of PCA factors only |
| 2SRR-FAVAR | inflation lags + PCA factor lags |

### h-step target: trailing cumulative

$$Y_h(t) = \sum_{j=0}^{h-1} y_{t-j}$$

Identical to Medeiros et al. (2021), which guarantees cell-by-cell
comparability between 2SRR and the twelve benchmarks.

### Out-of-sample design: 180 windows

OOS from **Jul/2010 to Jun/2025**, horizons `h ∈ {1, 3, 6, 12}`. The benchmarks
use a rolling window fixed at 606 observations (Medeiros `rolling_window.R`);
2SRR uses an expanding window starting at the same 606-observation training set.

### Parallelization

`parallel::makeCluster(N_CORES, type = "PSOCK")` in `03` parallelizes the 12
jobs (3 cases × 4 horizons) per window. Default: `detectCores() - 1`.

---

## Repository structure

```
00_prog/                Scripts (numbered, run in order)
  00_setup.R               Installs packages, downloads external functions, defines factor() PCA
  01_data_prep.R           yout, rw (180 x 12)
  02_forecast_medeiros.R   12 benchmark models
  03_forecast_2srr.R       3 TVP cases x 4 horizons x 180 windows
  04_analysis.R            Tests, figures, CSV tables
  05_figures_tables.R      Publication figures (PNG/PDF) and Word tables

10_data/                data.rda (FRED-MD, 786 obs x 117 vars + date) — not tracked

20_tools/
  21_coulombe/            Estimation functions from Goulet Coulombe's repository (auto-downloaded)
  22_medeiros/            functions.R + rolling_window.R (auto-downloaded)
  23_adapted/             tvp_functions.R (main engine + self-tests)

30_output/                forecasts/ and betas/ (.rda) — not tracked

40_results/               Figures (PDF/PNG) and CSV tables of the most recent run
```

---

## Estimation functions (`20_tools/21_coulombe/`)

| File | Key function | Role |
|------|--------------|------|
| `dualGRRmdA_v190215.R` | `dualGRR()` | Dual ridge solver |
| `TVPRR_v181111.R` | `TVPRR()` | TVP ridge core |
| `TVPRRcosso_v181120.R` | `TVPRR_cosso()` | Algorithm orchestrator (Steps 1–4) |
| `CVGSBHK_v181127.R` | `cvgs.bhk2015()` | Cross-validation (Bergmeir, Hyndman e Koo, 2018) |
| `zfun_v190304.R` | `Zfun()`, `make_reg_matrix()` | Z matrix construction |
| `fastZrot_v181125.R` | `fastZrot()` | Fast Z rotation for the factor version |
| `EM_sw.R` | `EM_sw()` | Stock-Watson EM imputation (FRED-MD NAs) |
| `ICp2.R` | `ICp2()` | Information criterion for factors |

The `factor()` PCA function is redefined in `00_setup.R` **before** sourcing the
external code, because `EM_sw()` calls `factor(X, n_fac = n)`, which clashes with
`base::factor()`.

## Adapted functions (`20_tools/23_adapted/tvp_functions.R`)

| Function | Role |
|----------|------|
| `tvp_2srr_coulombe_fast()` | Main engine: full algorithm in closed form |
| `cv_dualGRR_fast()` | CV via generalized Cholesky |
| `fit_2srr_window()` | Unified wrapper: selects the engine per window |
| `build_design_tvp()` | Builds the design matrices for the 3 TVP cases |
| `tvp_self_test()`, `tvp_coulombe_fast_vs_original_test()` | Numerical self-tests |

---

## Required packages

`glmnet`, `pracma`, `randomForest`, `forecast`, `lmtest`, `sandwich`,
`ggplot2`, `reshape2`, `HDeconometrics` (GitHub), `rugarch`, `fGarch`, `dplyr`,
`tidyr`, `patchwork`, `scales`, `RColorBrewer`, `gridExtra`, `MCS`, `flextable`.
All auto-install via `00_setup.R` and the script headers.

---

## References

- **Goulet Coulombe, P.** (2025). Time-Varying Parameters as Ridge Regressions.
  *International Journal of Forecasting*, 41(3), 982–1002.
  https://doi.org/10.1016/j.ijforecast.2024.08.006
  (Replication code: github.com/hugocout/Replication-codes-for-Time-Varying-Parameters-as-Ridge-Regressions)
- **Medeiros, M. C., Vasconcelos, G. F. R., Veiga, Á., & Zilberman, E.** (2021).
  Forecasting Inflation in a Data-Rich Environment: The Benefits of Machine
  Learning Methods. *Journal of Business & Economic Statistics*, 39(1), 98–119.
  (Code: github.com/gabrielrvsc/ForecastingInflation)
- **McCracken, M. W., & Ng, S.** (2016). FRED-MD: A Monthly Database for
  Macroeconomic Research. *Journal of Business & Economic Statistics*, 34(4),
  574–589.
- **Hansen, P. R., Lunde, A., & Nason, J. M.** (2011). The Model Confidence Set.
  *Econometrica*, 79(2), 453–497.
- **Giacomini, R., & White, H.** (2006). Tests of Conditional Predictive
  Ability. *Econometrica*, 74(6), 1545–1578.
- **De Mol, C., Giannone, D., & Reichlin, L.** (2008). Forecasting Using a Large
  Number of Predictors. *Journal of Econometrics*, 146(2), 318–328.
- **Bergmeir, C., Hyndman, R. J., & Koo, B.** (2018). A Note on the Validity of
  Cross-Validation for Evaluating Autoregressive Time Series Prediction.
  *Computational Statistics & Data Analysis*, 120, 70–83.

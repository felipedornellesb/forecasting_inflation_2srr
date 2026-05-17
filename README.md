# Inflation Forecasting with 2SRR

**Two-Step Ridge Regression** (Coulombe, IJF 2025) applied to US CPI inflation
forecasting (CPIAUCSL, FRED-MD), benchmarked against **Ridge, LASSO, Random
Forest, Bagging, Factor, T.Factor, CSR, and adaptive variants** from the
Medeiros et al. (JBES 2021) framework.

Undergraduate thesis in Economics, UFRGS — Felipe Dornelles, 2026.

---

## Pipeline overview

```
01_data_prep.R       -> yout.rda (realized) + rw.rda (Random Walk benchmark)
02_forecast_medeiros -> 12 ML / econometric models from Medeiros et al. (2021)
03_forecast_2srr.R   -> 2SRR in 3 TVP specifications (AR, Factor, FAVAR)
04_analysis.R        -> consolidated analysis (RMSFE, DM, GW, MCS, plots,
                        validation)
```

Run sequentially, in order. Script `04` is robust to missing inputs — any
Medeiros model that does not exist is silently ignored.

---

## How to run

```r
setwd("~/tcc/forecasting_inflation_2srr")
source("00_prog/00_setup.R")             # install packages, download Coulombe + Medeiros functions
source("00_prog/01_data_prep.R")         # generates yout.rda and rw.rda (180 windows x 12 horizons)
source("00_prog/02_forecast_medeiros.R") # 12 baseline models (~3-4h sequential)
source("00_prog/03_forecast_2srr.R")     # 2SRR via coulombe_fast (~2.5h with 7 cores parallel)
source("00_prog/04_analysis.R")          # final analysis + plots + LaTeX tables
```

Prerequisite: `data.rda` placed in `10_data/` (from the Medeiros repository
or from the course). Target variable: `CPIAUCSL`.

---

## Design decisions

### 2SRR — main engine: `coulombe_fast`

Coulombe's original `TVPRR_cosso` materializes the `T × KT` kernel explicitly
in `Zprime`. For monthly FRED-MD with T ≈ 600+, this explodes in memory and
time (~50-100h sequential).

The `coulombe_fast` engine implements **the same mathematics** but in closed
form on the `T × T` kernel:

$$K_{du}(\lambda_1, \lambda_2) = \frac{1}{\lambda_1} \sum_k \omega_k\, (X_k X_k') \odot M_{inn} + \frac{1}{\lambda_2} X_{full} X_{full}'$$

where `M_inn[i,j] = min(i,j) - 1`. CV via generalized Cholesky decomposition
sweeps `n_lambda = 15` in O(T²) per λ. **Result: ~10-50x faster with no
mathematical loss.**

#### Validations at the start of `03`

1. **Validation 1 — standalone vs direct primal (Eq. 5 of the paper)**: 3
   synthetic seeds, max diff < 1e-6.
2. **Validation 2 — `coulombe_fast` vs original `dualGRR`**: 2 synthetic
   configurations, max diff < 1e-10 (essentially numerical equality in
   float64).

If either validation fails, the script aborts before spending hours
computing.

### 3 TVP specifications compared

| Case | univar | factonly | nofact | Regressors |
|------|--------|----------|--------|------------|
| TVP-AR | T | — | T | Only `ly` lags of cumulative inflation Y_h |
| TVP-Factor | F | T | F | Only `lf` lags of PCA factors |
| TVP-FAVAR | F | F | F | Y_h lags + PCA factor lags |

### h-step target: trailing cumulative

$$Y_h(t) = \sum_{j=0}^{h-1} y_{t-j}$$

Identical to Medeiros et al. (2021) -> guarantees cell-by-cell comparability
between 2SRR and the 12 baseline models.

### Rolling window: 180 POOS windows

OOS from **Jul/2010 to Jun/2025**. In-sample window with `T = 606 to 785`
(grows by one observation at each step, the Medeiros `rolling_window.R`
convention).

### Parallelization

`parallel::makeCluster(N_CORES, type = "PSOCK")` in `03` parallelizes 12 jobs
(3 cases × 4 horizons) per window. Default: `detectCores() - 1`.

---

## Structure

```
00_prog/                Scripts (numbered, run in order)
  00_setup.R               Installs packages, downloads Coulombe+Medeiros,
                           defines factor() PCA
  01_data_prep.R           yout, rw (180 x 12)
  02_forecast_medeiros.R   12 ML models (Ridge, LASSO, RF, ...)
  03_forecast_2srr.R       3 TVP cases x 4 horizons x 180 windows
  04_analysis.R            Complete analysis: tables, figures, tests

10_data/                data.rda (FRED-MD, 786 obs x 117 vars + date)

20_tools/
  21_coulombe/            10 functions from Coulombe's repo (auto-downloaded)
  22_medeiros/            functions.R + rolling_window.R (auto-downloaded)
  23_adapted/             tvp_functions.R (main engine + self-tests)

30_output/
  forecasts/              <Model>.rda (180x12), yout.rda, rw.rda
  betas/                  betas_<Model>.rda (bundle per horizon)
  checkpoints/            ckpt_2SRR_w030.rda, ckpt_2SRR_w060.rda, ...

40_results/run_final_<timestamp>/
  figures/                Individual PDFs + combined 4h panels + plotly HTMLs
  tables/                 CSVs + .tex (LaTeX)
  final_narrative.txt     Automatic synthesis of the 12 sections
```

---

## Coulombe functions (`20_tools/21_coulombe/`)

| File | Key function | Role |
|------|--------------|------|
| `dualGRRmdA_v190215.R` | `dualGRR()` | Dual ridge solver (Eq. 9, 11 of the paper) |
| `TVPRR_v181111.R` | `TVPRR()` | TVP ridge core |
| `TVPRRcosso_v181120.R` | `TVPRR_cosso()` | Algorithm 1 orchestrator (Steps 1-4) |
| `CVGSBHK_v181127.R` | `cvgs.bhk2015()` | Cross-validation (Bergmeir, Hyndman, Koo 2018) |
| `CVKFMV_v190214.R` | CV helpers | Additional variants |
| `zfun_v190304.R` | `Zfun()`, `make_reg_matrix()` | Z matrix construction |
| `fastZrot_v181125.R` | `fastZrot()` | Fast Z rotation for the factor version |
| `EM_sw.R` | `EM_sw()` | Stock-Watson EM imputation (FRED-MD NAs) |
| `ICp2.R` | `ICp2()` | Information criterion for factors |
| `Xgenerators_v190127.R` | (generators) | Simulation support |

The `factor()` PCA function is redefined in `00_setup.R` **before** sourcing
the Coulombe code — because `EM_sw()` calls `factor(X, n_fac=n)` which
clashes with `base::factor()` (categorical factor levels).

---

## Adapted functions (`20_tools/23_adapted/tvp_functions.R`)

| Function | Role |
|----------|------|
| `tvp_2srr_coulombe_fast()` | **Main engine**: complete Algorithm 1 in O(T³) closed form |
| `dualGRR_fast()` | Fast equivalent of `dualGRR` (direct T × T kernel) |
| `cv_dualGRR_fast()` | CV via generalized Cholesky (sweeps n_lambda in O(T²)) |
| `tvp_2srr_standalone()` | Simplified version (single λ, no OLS prior) — alternative |
| `tvp_2srr_coulombe()` | Wrapper around the original `TVPRR_cosso` — used only for the sanity check |
| `fit_2srr_window()` | Unified wrapper: picks the engine via parameter |
| `build_design_tvp()` | Builds X_in, y_in, x_out for the 3 TVP cases |
| `tvp_self_test()` | Self-test: standalone vs direct primal (Eq. 5) |
| `tvp_coulombe_fast_vs_original_test()` | Self-test: `coulombe_fast` vs original `dualGRR` |
| `clark_west()` | Clark-West test (not used — non-nested models) |

---

## Analytical outputs of `04_analysis.R`

`04_analysis.R` covers **the advisor's 10 requests** + extras:

| Part | Analysis | Output |
|------|----------|--------|
| 0b | yout audit (cumulative consistency) | P0b_yout_*, audit table |
| 1 | RMSFE relative to Random Walk | P1_rmsfe_*, LaTeX |
| 2 | RMSE of the 3 TVP cases (AR vs Factor vs FAVAR) | P2_tvp_3cases_* |
| 3 | TVP beta trajectories per horizon | P3_betas_*_h<H>.pdf, 4h panel |
| 4 | Correlation between betas of the 3 TVP cases | P4_betas_cross_AR_FAVAR |
| 5 | Step 1 vs Step 4 lambdas (recalibration) | P5_lambdas_*, P5b_saturation |
| 6 | 2SRR vs Ridge (DM, CSSED, rolling RMSE) | P6_2srr_vs_ridge.csv + figs |
| 7 | 2SRR vs Medeiros (2 best + worst) | P7_2srr_vs_med_h<H>, 4h |
| 8 | Parsimony (HHI, near-zero, σ²_u) | P8_parsimony |
| 9 | Sub-periods (Pre-GFC, GFC, COVID, ...) | P9_subperiods, heatmap |
| 10 | Interactive plotly (HTML) | P10_heatmap_FAVAR_h<H>.html |
| 11 | Sanity check original Coulombe vs fast | P11_sanity_coulombe |
| 11b | 2SRR-FAVAR vs other TVPs (monthly scale) | P11b_2srr_vs_tvps_* |
| 11c | 2SRR-FAVAR vs classical Ridge (monthly scale) | P11c_2srr_vs_ridgemed_* |
| 11d | TVP betas vs constant Ridge | P11d_betas_tvp_vs_ridge_* |
| 11e | Evolution of TVP betas (3 cases × 4h) | P11e_betas_evolution_<case>_4h |
| 11f | Combined 4h panels (CSSED, rolling, λ) | P6_*_4h, P5_*_4h |
| 12 | **MCS** — Model Confidence Set (Hansen-Lunde-Nason 2011) | P12_MCS_* |
| 13 | **GW** — Giacomini-White (2006) | P13_GW_test, heatmap |
| 13b | PhD-level econometric validation (12 auditable items) | P13b_econometric_validation |
| 14 | Consolidated final narrative | final_narrative.txt |

All figures are exported as individual PDF **and** as a 2×2 combined panel
(`_4h.pdf`) across the 4 horizons.

---

## Required packages

`glmnet`, `pracma`, `randomForest`, `forecast`, `lmtest`, `sandwich`,
`ggplot2`, `reshape2`, `xtable`, `HDeconometrics` (from GitHub), `rugarch`,
`fGarch`, `dplyr`, `tidyr`, `patchwork`, `scales`, `RColorBrewer`, `plotly`,
`htmlwidgets`, `knitr`, `gridExtra`, `MCS`.

All auto-install via `00_setup.R` and the header of `04_analysis.R`.

---

## References

- **Coulombe, P. G.** (2025). Time-Varying Parameters as Ridge Regressions.
  *International Journal of Forecasting* 41(3), 982–1002.
  https://doi.org/10.1016/j.ijforecast.2024.08.006
  - Replication code:
    github.com/hugocout/Replication-codes-for-Time-Varying-Parameters-as-Ridge-Regressions

- **Medeiros, M. C., Vasconcelos, G. F. R., Veiga, Á., & Zilberman, E.**
  (2021). Forecasting Inflation in a Data-Rich Environment: The Benefits of
  Machine Learning Methods. *Journal of Business & Economic Statistics*,
  39(1), 98–119. https://doi.org/10.1080/07350015.2019.1637745
  - Code: github.com/gabrielrvsc/ForecastingInflation

- **Hansen, P. R., Lunde, A., & Nason, J. M.** (2011). The Model Confidence
  Set. *Econometrica*, 79(2), 453–497.

- **Giacomini, R. & White, H.** (2006). Tests of Conditional Predictive
  Ability. *Econometrica*, 74(6), 1545–1578.

- **McCracken, M. W., & Ng, S.** (2016). FRED-MD: A Monthly Database for
  Macroeconomic Research. *Journal of Business & Economic Statistics*, 34(4),
  574–589. (Underlying dataset.)

---

**Felipe Dornelles** — Universidade Federal do Rio Grande do Sul (UFRGS) —
Economics Undergraduate Thesis 2026
Advisor: Prof. Hudson Chaves Costa

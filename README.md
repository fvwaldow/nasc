
<!-- README.md is generated from README.Rmd. Please edit that file -->

# nasc - Network-Aware Synthetic Control R Package

<!-- badges: start -->

<!-- badges: end -->

## Introduction

`nasc` is an R package that implements the **Network-Aware Synthetic
Control (NASC)** estimator introduced by von Waldow (2026) that
addresses the limitations of conventional synthetic control methods when
SUTVA is violated and interference between units prevails.

The estimator follows a **modular Bayesian (cut-posterior) via multiple
imputation**: a spatial autocorrelation parameter *ρ* is first estimated
from a spatial panel model, and then passed into a synthetic-control
likelihood that estimates donor weights, the optional NASC penalty, and
the bias-correction factor. The cut prevents feedback from the
synthetic-control fit back into the spatial model.

## Installation

The current version of `nasc` can be installed from
[GitHub](https://github.com/) with:

``` r
pak::pkg_install("fvwaldow/nasc")
```

As the package compiles Stan models,
[**Rtools**](https://cran.r-project.org/bin/windows/Rtools/) and
[**RStan**](https://mc-stan.org/rstan/) are needed. The `nasc` package
was built under R 4.4.2.

## Quick Start

``` r
library(nasc)

# View basic package information
packageDescription("nasc")

# List all functions in the package
ls("package:nasc")
```

A minimal workflow

``` r
library(nasc)

# `panel`        long-format panel data with columns: id, year, y, treated
# `W`            row-standardized spatial-weights matrix (units x units)
# `panel_covs`   long-format covariates with columns: id, year, x1, x2, ...

mod <- nascSynth$new(
  data            = panel,
  time            = year,
  id              = id,
  treated         = treated,
  outcome         = y,
  W               = W,
  spatial_model   = "SAR",
  bias_correction = TRUE,
  nasc_penalty    = TRUE,
  ci_width        = 0.90,

  # Module 2
  covariates            = panel_covs,
  predictors.op         = "mean",
  time.predictors.prior = 1995:1999,            # pre-treatment window
  special.predictors    = list(
    list("y", c(1992, 1995, 1998), "mean"),     # lagged-outcome matching rows
    list("x1", 1999, "mean")
  ),
  # Importance weights V
  predictor_weights = c(0.8, 0.1),

  # Module 1
  rho.covariates  = c("x1", "x2"),              # estimate rho from these covariates only
  rho.time.window = 1996:1999                   # pre-treatment window
)

mod$fit(cores = 12, n_samples = "auto")

mod$summary(ci_width = 0.90)     # ATT, per-period TE, weights, diagnostics
mod$indirectEffect()             # posterior draws/summaries of the spillover effect
mod$updateWidth(0.80)            # re-summarise at a new CI width without re-running MCMC

mod$syntheticPlot()              # observed vs. synthetic outcome path
mod$effectPlot()                 # per-period TE with CrI
mod$attPlot()                    # ATT density 
mod$tauPlot()                    # per-period TE densities 
mod$posteriorPlot()              # rho, beta, lambda, sigma, bias_correction, ...
mod$weightDraws()                # donor weight densities
mod$weightCorr()                 # donor weight correlation

contaminationGraph(mod)          # network graph colored by contamination value
effectGraph(mod)                 # network graph colored by indirect effect
```

## Usage

### Core functions

| Function / method | Purpose |
|----|----|
| `nascSynth$new(...)` | Constructs an estimator. Selects a Stan model based on `spatial_model`, `bias_correction`, and `nasc_penalty`. Predictors used in Module 2 covariate matching are specified in **Synth**-compatible style (`covariates`, `predictors.op`, `special.predictors`, `time.predictors.prior`), and the resulting predictor rows can be reweighted in the matching loss via `predictor_weights` (the diagonal V matrix). The Module 1 spatial regression that estimates *ρ* can be restricted, independently of the Module 2 predictors, to a subset of covariates via `rho.covariates` and to a pre-treatment window via `rho.time.window`. |
| `mod$fit(n_samples, n_samples_cap, n_samples_min, cores, worker_iter, worker_warmup, ...)` | Runs the (optional) Module 1 spatial model and the Module 2 synthetic-control model via MCMC |
| `mod$summary(ci_width, print)` | Posterior summary: ATT, per-period τₜ, donor weights, model parameters, and MCMC diagnostics. When the model uses a network, also reports per-period and per-donor indirect (spillover) effects and the average indirect effect. Also dispatched via `summary(mod)` |
| `mod$indirectEffect()` | Posterior draws and summaries of the indirect (spillover) effect at the per-period and per-donor level. Returns `NULL` for non-network configurations |
| `mod$updateWidth(ci_width)` | Recomputes credible intervals at a different width without re-running MCMC |
| `mod$plotData` | Tibble of observed and counterfactual outcomes with credible-interval bounds |
| `mod$interventionTime` | The first treatment period |

#### `nascSynth$new()` arguments

| Argument | Default | Description |
|----|----|----|
| `data`, `time`, `id`, `treated`, `outcome` | — | Long-format panel and column accessors. |
| `ci_width` | `0.75` | Width of the credible interval. |
| `covariates` | `NULL` | Long-format covariate panel. Required for `spatial_model = "SAR"` or `"SDM"`. |
| `W` | `NULL` | Row-standardized spatial-weights matrix. Required when `bias_correction = TRUE` or `nasc_penalty = TRUE`. |
| `spatial_model` | `"none"` | One of `"none"`, `"SAR"`, `"SDM"`, `"exogenous"`. |
| `rho` | `NULL` | Optional scalar in (-1, 1). Required when `spatial_model = "exogenous"`; supplying it for any other choice skips Module 1. |
| `bias_correction` | model-aware | `TRUE` for SAR/SDM/exogenous, `FALSE` for `"none"`. |
| `nasc_penalty` | model-aware | Same default rule as `bias_correction`. |
| `predictor_weights` | `NULL` | Optional non-negative numeric vector with one entry per predictor row, equivalent to the diagonal of Abadie-Diamond-Hainmueller’s V matrix. Accepts either a named vector (matched by predictor row label; every predictor row must be named) or an unnamed vector of the right length (matched positionally; regular covariates first, then `special.predictors` in user order). `NULL` gives equal weight to every predictor row; setting an entry to 0 drops that row from the matching loss. |
| `predictors.op` | `"mean"` | Character string naming the operator applied to every regular predictor (every non-`id`, non-`time` column in `covariates`) over the pre-treatment period. Any name resolvable by `match.fun()` works (`"mean"`, `"median"`, etc.). Mirrors `predictors.op` in the **Synth** package. |
| `special.predictors` | `NULL` | Synth-style list of length-3 specifications `list(<predictor>, <time periods>, <operator>)`, one per extra matching row. The predictor may be any column in `data` or `covariates`, so lagged outcomes (e.g. `list("y", 1993, "mean")`) are supported without restating them as covariates. Mirrors `special.predictors` in the **Synth** package. |
| `time.predictors.prior` | `NULL` | Optional numeric vector of pre-treatment periods over which regular predictors are aggregated. Default `NULL` uses all pre-intervention periods. Does not affect `special.predictors`. Mirrors `time.predictors.prior` in the **Synth** package. |
| `rho.covariates` | `NULL` | Optional character vector naming which columns of `covariates` enter the Module 1 *ρ* regression. `NULL` uses every covariate column. Only consulted when Module 1 actually runs (`spatial_model = "SAR"` or `"SDM"` with no exogenous `rho`); the Module 2 predictor set is left unchanged. |
| `rho.time.window` | `NULL` | Optional numeric vector of pre-treatment periods to which the Module 1 *ρ* regression is restricted. `NULL` uses all pre-treatment periods. Periods that are post-treatment or absent from the panel are dropped with a warning. Does not affect the Module 2 outcome match. |

#### `fit()` arguments

| Argument | Default | Description |
|----|----|----|
| `n_samples` | `100` | Number of *ρ* draws propagated from Module 1 to Module 2. May be a positive integer or the string `"auto"` (see below). Ignored when *ρ* is exogenous or unused. |
| `n_samples_cap` | `500` | Upper bound on the auto-selected `n_samples`. Only used when `n_samples = "auto"`. |
| `n_samples_min` | `30` | Lower bound on the auto-selected `n_samples`. Only used when `n_samples = "auto"`. |
| `cores` | `detectCores() - 1` | Number of CPU cores for parallel execution. |
| `worker_iter` | `2000` | Iterations per worker chain in the multi-*ρ* parallel loop. |
| `worker_warmup` | `1000` | Warmup per worker chain. |
| `...` |  | Additional arguments forwarded to `rstan::sampling()`. |

### Plot functions

| Method / function | Output |
|----|----|
| `mod$syntheticPlot()` | Observed vs. synthetic outcome path with credible band |
| `mod$effectPlot()` | Direct treatment effect τₜ over time with credible band and intervention line |
| `mod$attPlot(indirect)` | Posterior density of the ATT. When the model uses a network and `indirect = TRUE` (the default in that case), a second panel adds the average indirect (spillover) effect |
| `mod$tauPlot(indirect)` | Stacked posterior densities of period-by-period τₜ. When the model uses a network and `indirect = TRUE` (default), a second column shows the per-period donor-average indirect effect |
| `mod$posteriorPlot()` | Posterior densities of the estimated scalar parameters: *ρ*, *θ*ₖ (SDM only), *β*ₖ (SAR only, time-invariant covariates skipped), `sigma_step1`, `lambda` (when `nasc_penalty = TRUE`), `sigma_step2`, and `bias_correction` (when `bias_correction = TRUE`). Each parameter only appears if the configuration actually estimated it. |
| `mod$weightDraws()` | Ridge-style densities of donor weights across posterior draws |
| `mod$weightCorr()` | Heatmap of correlations between donor weights across draws |
| `contaminationGraph(model, ...)` | Spatial-weights graph with donor nodes coloured by posterior-mean contamination \|sⱼ\| (requires `igraph`) |
| `effectGraph(model, ...)` | Spatial-weights graph with donor nodes coloured by posterior-mean indirect effect δ̄ⱼ and the treated node coloured by the ATT, on a shared diverging palette (requires `igraph`) |
| `nascPlot(models, show_ci)` | Overlays the synthetic and direct-effect plots of a named list of fitted `nascSynth` objects |
| `nascWeight(models, ...)` | Multi-model ridgeline of posterior donor-weight densities, one row per donor, one density per model, sharing the colour palette of `nascPlot()` |

## Citation

If `nasc` is implemented in academic work, please cite

> von Waldow, F. (2026). *Network-Aware Synthetic Control - Bias
> Correction and Regularization under Interference.* Master Thesis.

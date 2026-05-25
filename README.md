
<!-- README.md is generated from README.Rmd. Please edit that file -->

# nasc - Network-Aware Synthetic Control R Package

<!-- badges: start -->

<!-- badges: end -->

## Introduction

`nasc` is an R package that implements the **Network-Aware Synthetic
Control (NASC)** estimator. It is based on the methodology of von Waldow
(2026) and addresses the limitations of conventional synthetic control
methods when SUTVA is violated and interference between units prevails.

The estimator follows a **modular Bayesian (cut-posterior) via multiple
imputation**: a spatial autocorrelation parameter *ρ* is first estimated
from a SAR or SDM panel model, and then propagated as fixed data into a
synthetic-control likelihood that estimates donor weights, the optional
NASC penalty, and the bias-correction factor. The cut prevents feedback
from the synthetic-control fit back into the spatial model.

## Installation

The development version of `nasc` can be installed from
[GitHub](https://github.com/) with:

``` r
pak::pkg_install("fvwaldow/nasc")
```

As the package compiles Stan models at install time, a working C++
toolchain and a recent version of
[**rstan**](https://mc-stan.org/rstan/) is needed. The contamination and
effect network plots additionally depend on **igraph**, which can be
installed with `install.packages("igraph")`.

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
  covariates      = panel_covs,
  W               = W,
  spatial_model   = "SAR",
  bias_correction = TRUE,
  nasc_penalty    = TRUE
)

mod$fit(cores = 12)              # default n_samples = 100
# or, sized from the Step-1 ESS of rho:
# mod$fit(cores = 12, n_samples = "auto")

mod$summary()                    # ATT, per-period effects, donor weights, diagnostics
mod$syntheticPlot()              # observed vs. synthetic outcome path
mod$effectPlot()                 # tau_t with credible interval
```

## Usage

### Core functions

| Function / method | Purpose |
|----|----|
| `nascSynth$new(...)` | Constructs an estimator. Selects a Stan model based on `spatial_model`, `bias_correction`, and `nasc_penalty`. Predictors used in Step-2 covariate matching are specified in **Synth**-compatible style (`covariates`, `predictors.op`, `special.predictors`, `time.predictors.prior`), and the resulting predictor rows can be reweighted in the matching loss via `predictor_weights` (the diagonal V matrix). |
| `mod$fit(n_samples, n_samples_cap, n_samples_min, cores, worker_iter, worker_warmup, ...)` | Runs the (optional) Step 1 spatial model and the Step 2 synthetic-control model via MCMC |
| `mod$summary(ci_width, print)` | Posterior summary: ATT, per-period τₜ, donor weights, model parameters, and MCMC diagnostics. When the model uses a network, also reports per-period and per-donor indirect (spillover) effects and the average indirect effect. Also dispatched via `summary(mod)` |
| `mod$indirectEffect()` | Posterior draws and summaries of the indirect (spillover) effect at the per-period and per-donor level. Returns `NULL` for non-network configurations |
| `mod$updateWidth(ci_width)` | Recomputes credible intervals at a different width without re-running MCMC |
| `mod$plotData` | Tibble of observed and counterfactual outcomes with credible-interval bounds |
| `mod$interventionTime` | The first treatment period |

#### `nascSynth$new()` arguments at a glance

| Argument | Default | Description |
|----|----|----|
| `data`, `time`, `id`, `treated`, `outcome` | — | Long-format panel and column accessors. |
| `ci_width` | `0.75` | Width of the credible interval. |
| `covariates` | `NULL` | Long-format covariate panel. Required for `spatial_model = "SAR"` or `"SDM"`. |
| `W` | `NULL` | Row-standardized spatial-weights matrix. Required when `bias_correction = TRUE` or `nasc_penalty = TRUE`. |
| `spatial_model` | `"none"` | One of `"none"`, `"SAR"`, `"SDM"`, `"exogenous"`. |
| `rho` | `NULL` | Optional scalar in (-1, 1). Required when `spatial_model = "exogenous"`; supplying it for any other choice skips Step 1. |
| `bias_correction` | model-aware | `TRUE` for SAR/SDM/exogenous, `FALSE` for `"none"`. |
| `nasc_penalty` | model-aware | Same default rule as `bias_correction`. |
| `predictor_weights` | `NULL` | Optional non-negative numeric vector with one entry per predictor row, equivalent to the diagonal of Abadie-Diamond-Hainmueller’s V matrix. Accepts either a named vector (matched by predictor row label) or an unnamed vector of the right length (matched positionally; regular covariates first, then `special.predictors` in user order). `NULL` gives equal weight to every predictor row; setting an entry to 0 drops that row from the matching loss. |
| `predictors.op` | `"mean"` | Character string naming the operator applied to every regular predictor (every non-`id`, non-`time` column in `covariates`) over the pre-treatment period. Any name resolvable by `match.fun()` works (`"mean"`, `"median"`, etc.). Mirrors `predictors.op` in the **Synth** package. |
| `special.predictors` | `NULL` | Synth-style list of length-3 specifications `list(<predictor>, <time periods>, <operator>)`, one per extra matching row. The predictor may be any column in `data` or `covariates`, so lagged outcomes (e.g. `list("y", 1993, "mean")`) are supported without restating them as covariates. Mirrors `special.predictors` in the **Synth** package. |
| `time.predictors.prior` | `NULL` | Optional numeric vector of pre-treatment periods over which regular predictors are aggregated. Default `NULL` uses all pre-intervention periods. Does not affect `special.predictors`. Mirrors `time.predictors.prior` in the **Synth** package. |

#### `fit()` arguments

| Argument | Default | Description |
|----|----|----|
| `n_samples` | `100` | Number of *ρ* draws propagated from Step 1 to Step 2. May be a positive integer or the string `"auto"` (see below). Ignored when *ρ* is exogenous or unused. |
| `n_samples_cap` | `500` | Upper bound on the auto-selected `n_samples`. Only used when `n_samples = "auto"`. |
| `n_samples_min` | `30` | Lower bound on the auto-selected `n_samples`. Only used when `n_samples = "auto"`. |
| `cores` | `detectCores() - 1` | Number of CPU cores for parallel execution. |
| `worker_iter` | `2000` | Iterations per worker chain in the multi-*ρ* parallel loop. |
| `worker_warmup` | `1000` | Warmup per worker chain. |
| `...` |  | Additional arguments forwarded to `rstan::sampling()`. |

When `n_samples = "auto"`, the package reads the effective sample size
of *ρ* from the Step-1 fit and sets the propagated subsample to
`min(round(n_eff_rho), n_samples_cap)`, floored at `n_samples_min`. This
avoids oversampling a low-information posterior of *ρ*.

### Visualization functions

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
| `contaminationScatter(model, ...)` | Per-donor scatter of posterior-mean weight against posterior-mean contamination |
| `nascPlot(models, show_ci)` | Overlays the synthetic and direct-effect plots of a named list of fitted `nascSynth` objects |
| `nascWeight(models, ...)` | Multi-model ridgeline of posterior donor-weight densities, one row per donor, one density per model, sharing the colour palette of `nascPlot()` |

## Citation

If `nasc` is implemented in academic work, please cite

> von Waldow, F. (2026). *Network-Aware Synthetic Control - Bias
> Correction and Regularization under Interference.* Master Thesis.

## License

MIT License. See `LICENSE` for details.

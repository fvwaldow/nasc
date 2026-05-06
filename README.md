
<!-- README.md is generated from README.Rmd. Please edit that file -->

# nasc - Network-Aware Synthetic Control R Package

<!-- badges: start -->

<!-- badges: end -->

## Introduction

`nasc` is an R package that implements the **Network-Aware Synthetic
Control (NASC)** estimator. It is based on the methodology of von Waldow
(2026) and addresses the limitations of conventional synthetic control
methods when SUTVA is violated and interference between units prevails.

The estimator follows a **modular Bayesian (cut-posterior) workflow**: a
spatial autocorrelation parameter *ρ* is first estimated from a SAR or
SDM panel model, and then propagated as fixed data into a
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

A minimal end-to-end workflow looks like this:

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

The package exposes a single user-facing [**R6**](https://r6.r-lib.org/)
class, `nascSynth`, several stand-alone diagnostic plot functions, and
comparison helpers for overlaying several fitted models.

| Function / method | Purpose |
|----|----|
| `nascSynth$new(...)` | Constructs an estimator. Selects a Stan model based on `spatial_model`, `bias_correction`, and `nasc_penalty`. Optionally accepts `predictor_weights` to set a diagonal V matrix on the covariate-matching loss in Step 2. |
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
| `predictor_weights` | `NULL` | Optional non-negative numeric vector with one entry per covariate, equivalent to the diagonal of Abadie-Diamond-Hainmueller’s V matrix. `NULL` gives equal weight to every covariate; setting an entry to 0 drops that covariate from the matching loss. |

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

All plotting methods produce base R graphics.

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

## Methodological notes

The estimator decomposes into two steps that the package dispatches
automatically.

**Step 1 — spatial model.** When `spatial_model` is `"SAR"` or `"SDM"`
and no exogenous `rho` is supplied, a Bayesian SAR or Spatial Durbin
Model is fit on the pre-treatment panel using the user-supplied
covariates. Posterior draws of *ρ* are stratified across chains and
propagated to Step 2 as fixed data (`n_samples` controls the size of the
propagated subsample, or `"auto"` selects it from the Step-1 ESS). When
`spatial_model = "exogenous"`, Step 1 is skipped and the user-supplied
scalar `rho` is used directly. When neither the bias correction nor the
NASC penalty is active, no spatial parameter is needed and Step 1 is
skipped entirely.

**Step 2 — NASC weight estimation.** Donor weights *w* are estimated on
the pre-treatment outcome panel (and, when applicable, on covariate
moments weighted by `predictor_weights`). The post-treatment
counterfactual is reconstructed as a convex combination of the donor
pool; if `bias_correction = TRUE`, it is rescaled by the bias-correction
factor 1 / (1 − ⟨w, s⟩), so that contamination of the donor pool by
spillovers from the treated unit is accounted for. If
`nasc_penalty = TRUE`, the log-likelihood includes a penalty term −λ⟨w,
\|s\|⟩, which discourages weight on units with strong network exposure
to the treated unit. When several *ρ* draws are propagated from Step 1,
Step 2 is run in parallel across them with `furrr` (one Stan fit per
*ρ*, single-chain `worker_iter` iterations each) and the resulting
posteriors are pooled to approximate the cut posterior.

**Cut posterior.** *ρ* enters Step 2 as data, not as a parameter.
Conditional on each *ρ* draw, Step 2’s likelihood and posterior are
sealed off from Step 1, so the synthetic-control fit cannot inform the
spatial model. The final approximation to the cut posterior is the
equally-weighted Monte Carlo mixture over the propagated *ρ* draws.
Per-worker MCMC diagnostics (split-Rhat, n_eff, divergent transitions,
max-treedepth saturation) are tracked and summarized by `mod$summary()`.

**Indirect effects.** When the model uses a network, the per-donor
spillover at post-period *t* is δᵢₜᴺᴬˢᶜ = sᵢ · τ₁ₜᴺᴬˢᶜ (Proposition
6.2), with the contamination vector s = ρ (Iⱼ − ρWⱼ)⁻¹ wⱼ₁ pulled from
each posterior draw. `mod$indirectEffect()` returns the full posterior
tensor, while `summary()`, `attPlot()` and `tauPlot()` surface the
per-period donor-average and the scalar average indirect effect (the
spillover analog of the ATT). For a network-graph view, see
`effectGraph()`.

**Spatial weights.** The matrix `W` must be row-standardized and known a
priori. The code checks this explicitly and aligns row/column names with
the unit identifier so that donor and treated rows can be picked
unambiguously.

**Identification scope.** The current implementation supports a single
treated unit. Multiple treated units, staggered adoption, and
non-stationary network structures are subject to further research.

## Citation

If `nasc` is implemented in academic work, please cite

> von Waldow, F. (2026). *Network-Aware Synthetic Control - Bias
> Correction and Regularization under Interference.* Master Thesis.

## License

MIT License. See `LICENSE` for details.

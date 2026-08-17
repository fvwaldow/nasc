\name{nascSynth}
\alias{nascSynth}
\title{Network-Aware Synthetic Control}
\description{
An R6 class implementing a Bayesian Synthetic Control estimator with optional
Network-Aware penalty and correction. Module 1 estimates the spatial
autocorrelation parameter \eqn{\rho}. By multiple imputations, draws of
\eqn{\rho} are passed into Module 2, a SC likelihood over the simplex of donor
weights. Both network components are driven by the contamination vector
\eqn{s = \rho (I_J - \rho W_J)^{-1} w_{J1}}
\itemize{
\item \code{bias_correction}: factor that corrects for bias induced by
contaminated donors \eqn{1/(1 - \gamma's)}.
\item \code{nasc_penalty}: contamination regularization through penalty
\eqn{\lambda \gamma'|s|} that enters the log-likelihood. \eqn{\lambda} is
pre-specified (\code{lambda}) or selected by hold-out CV in \code{$fit()}.
}
}
\usage{
nascSynth$new(data, time, id, treated, outcome, ci_width = 0.75,
              covariates = NULL, W = NULL, spatial_model = "none", rho = NULL,
              bias_correction = NULL, nasc_penalty = NULL,
              predictor_weights = NULL, predictors.op = "mean",
              special.predictors = NULL, time.predictors.prior = NULL,
              rho.covariates = NULL, rho.time.window = NULL,
              lambda = NULL, lambda_cv_grid = seq(0, 100, by = 1),
              lambda_train_frac = 0.7)
}
\arguments{
\item{data}{Long-format panel with unit, time, outcome and a binary
treatment column (1 = treated unit in treated periods).}

\item{time, id, treated, outcome}{Unquoted column names.}

\item{ci_width}{Credible-interval width.}

\item{covariates}{Long-format covariate panel. Required for
\code{spatial_model = "SAR"} or \code{"SDM"}.}

\item{W}{Row-standardized weights matrix, dimnames matching \code{id}.
Required when \code{bias_correction} or \code{nasc_penalty} is TRUE.}

\item{spatial_model}{\code{"none"}, \code{"SAR"}, \code{"SDM"} or
\code{"exogenous"}.}

\item{rho}{Scalar in (-1, 1). Skips Module 1 whatever \code{spatial_model}
is; required for \code{"exogenous"}.}

\item{bias_correction, nasc_penalty}{Logical. Default TRUE for
SAR/SDM/exogenous, FALSE for \code{"none"}.}

\item{predictor_weights}{Non-negative vector, one entry per predictor row
(the diagonal of the V matrix). Named (by row label) or unnamed
(positional: regular covariates, then \code{special.predictors}).
NULL weights every row 1; 0 drops a row from the matching loss.}

\item{predictors.op}{Operator aggregating each regular predictor over the
pre-period; any name resolvable by \code{\link[base]{match.fun}}. NAs
dropped.}

\item{special.predictors}{List of \code{list(<predictor>, <periods>, <op>)}
triples, one extra matching row each, e.g. \code{list("Y", 1993, "mean")}.
The predictor may be any column of \code{data} or \code{covariates}.}

\item{time.predictors.prior}{Pre-periods for regular predictors. NULL uses
all; no effect on \code{special.predictors}.}

\item{rho.covariates}{Covariate columns entering Module 1. NULL uses all;
Module 2 predictors unchanged.}

\item{rho.time.window}{Pre-periods for Module 1. NULL uses all;
post-treatment or absent periods are dropped with a warning.}

\item{lambda}{Fixed penalty strength. NULL selects it by CV in
\code{$fit()}.}

\item{lambda_cv_grid}{Candidate \eqn{\lambda} values. Default
\code{seq(0, 100, 1)}.}

\item{lambda_train_frac}{Share of the pre-period used for CV training; the
adjacent block before treatment is the validation window. Default 0.7.}
}
\section{Active bindings}{
\describe{
\item{\code{plotData}}{Tibble of observed and counterfactual outcomes with CrI bounds.}
\item{\code{interventionTime}}{First treatment period.}
}
}
\section{Methods}{
\describe{
\item{\code{$fit(n_samples = 100, n_samples_cap = 500L, n_samples_min = 30L,
cores = parallel::detectCores() - 1, worker_iter = 2000L,
worker_warmup = 1000L, ...)}}{Runs Module 1 (if specified) and Module 2 via
MCMC. \code{n_samples} is the number of \eqn{\rho} draws cut into Module 2:
a positive integer, or \code{"auto"}, which sets it to the STAN \code{n_eff}
of \eqn{\rho} clamped to \code{[n_samples_min, n_samples_cap]} and to the
available draws. \code{...} goes to
\code{\link[rstan:stanmodel-method-sampling]{rstan::sampling}}.}
\item{\code{$summary(ci_width = NULL, print = TRUE)}}{ATT, per-period
direct and indirect \eqn{\tau_t}, donor weights, parameters and MCMC
diagnostics. Prints, and invisibly returns a \code{summary.nascSynth} list.
Also reached via \code{summary(obj)}.}
\item{\code{$indirectEffect()}}{Draws and summaries of the spillover
\eqn{\delta_{j,t} = s_j \tau_{1t}}. Returns a list with
\code{donor_names}, \code{time_post}, \code{delta_arr}
\eqn{[draws \times T_{post} \times J]}, \code{delta_total},
\code{avg_per_donor} and \code{avg_total}; NULL without a network.}
\item{\code{$lambdaCV()}}{List with the deployed \code{lambda} and the CV
\code{table} (\code{lambda}, \code{rmse}, \code{se}, \code{contamination},
\code{picked}). \code{table} is NULL when \code{lambda} was user-fixed.}
\item{\code{$updateWidth(ci_width = 0.75)}}{Recomputes credible intervals and
plot data without re-running MCMC.}
\item{\code{$syntheticPlot()}}{Observed vs. synthetic path with credible band.}
\item{\code{$effectPlot(indirect = NULL, show_avg = TRUE, max_legend = 12L,
indirect_pre = c("zero", "drop", "placebo"))}}{Direct effect
\eqn{\tau_{1t}} over time with credible band. With a network and
\code{indirect = TRUE} (the default there) a second panel draws the
per-donor spillover, its donor average if \code{show_avg}, and a legend
for up to \code{max_legend} donors. \code{indirect_pre} sets the
pre-treatment display: zero, dropped, or placebo path.}
\item{\code{$attPlot(indirect = NULL)}}{Posterior density of the ATT; with a
network and \code{indirect = TRUE} (default) a second panel adds the
average indirect effect.}
\item{\code{$tauPlot(indirect = NULL)}}{Period-by-period densities of
\eqn{\tau_{1t}}; with a network and \code{indirect = TRUE} (default) a
second column adds the donor-average \eqn{\bar\delta_t}.}
\item{\code{$posteriorPlot()}}{Densities of the estimated scalars, panel per
parameter: \code{rho}, \code{theta[k]} (SDM), \code{beta[k]} (SAR;
time-invariant covariates are unidentified under the within transform and
are skipped), \code{sigma_step1}, \code{sigma_step2},
\code{bias_correction}, \code{lambda_tilde}. Constant quantities are
skipped, so \code{lambda_tilde} normally does not appear and
\code{bias_correction} only when it is on.}
\item{\code{$weightDraws(overlap = 0.5, scale = 1.4, fill_alpha = 0.85)}}{Ridgeline
of posterior donor weights. \code{overlap} in [0, 1) is the vertical
overlap between ridges, \code{scale} their height, \code{fill_alpha} the
fill transparency.}
\item{\code{$weightCorr()}}{Heatmap of pairwise correlations between donor
weights across draws; the diagonal is NA.}
\item{\code{$clone(deep = FALSE)}}{Clone the object.}
}
}
\seealso{
\code{\link{nascPlot}}, \code{\link{nascWeight}},
\code{\link{contaminationGraph}}, \code{\link{effectGraph}}.
}

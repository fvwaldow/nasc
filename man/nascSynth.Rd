% Handwritten Manual
\name{nascSynth}
\alias{nascSynth}
\title{Bayesian Network-Aware Synthetic Control}
\description{
An R6 class implementing Bayesian Synthetic Control estimators with optional
network/spatial awareness.

Three orthogonal options govern how the spatial structure enters the
estimator:

\itemize{
\item \code{spatial_model}: \code{"none"}, \code{"SAR"}, \code{"SDM"}, or
\code{"exogenous"}. The first three behave as before; the last skips
Step 1 and uses a user-supplied scalar \code{rho}.
\item \code{rho}: optional numeric in \code{(-1, 1)}. When provided,
Step 1 is skipped regardless of \code{spatial_model} and the supplied
value is used directly.
\item \code{bias_correction} (logical): whether the post-treatment effect
is rescaled by \eqn{1/(1 - w's)}.
\item \code{nasc_penalty} (logical): whether the NASC penalty
\eqn{-\lambda \langle w, |s| \rangle} enters the likelihood.
}

Model dispatch (single source of truth):
\itemize{
\item \code{nasc_penalty = TRUE}  -> stan_2_NASC.stan (penalty always on
in this file, bias_correction toggleable).
\item \code{nasc_penalty = FALSE} -> model1.stan (no penalty; bias
correction toggleable through the same generated-quantities mechanism).
}

Step 1 (SAR or SDM rho estimation) runs only when (i) the chosen Step 2
configuration uses a non-trivial rho (penalty on, or bias correction on)
AND (ii) no exogenous rho was supplied.
}
\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\describe{
\item{\code{plotData}}{Tibble with observed and counterfactual outcomes.}

\item{\code{interventionTime}}{The first treatment period.}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-nascSynth-new}{\code{nascSynth$new()}}
\item \href{#method-nascSynth-fit}{\code{nascSynth$fit()}}
\item \href{#method-nascSynth-updateWidth}{\code{nascSynth$updateWidth()}}
\item \href{#method-nascSynth-summary}{\code{nascSynth$summary()}}
\item \href{#method-nascSynth-indirectEffect}{\code{nascSynth$indirectEffect()}}
\item \href{#method-nascSynth-syntheticPlot}{\code{nascSynth$syntheticPlot()}}
\item \href{#method-nascSynth-effectPlot}{\code{nascSynth$effectPlot()}}
\item \href{#method-nascSynth-posteriorPlot}{\code{nascSynth$posteriorPlot()}}
\item \href{#method-nascSynth-attPlot}{\code{nascSynth$attPlot()}}
\item \href{#method-nascSynth-tauPlot}{\code{nascSynth$tauPlot()}}
\item \href{#method-nascSynth-weightDraws}{\code{nascSynth$weightDraws()}}
\item \href{#method-nascSynth-weightCorr}{\code{nascSynth$weightCorr()}}
\item \href{#method-nascSynth-clone}{\code{nascSynth$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-nascSynth-new"></a>}}
\if{latex}{\out{\hypertarget{method-nascSynth-new}{}}}
\subsection{Method \code{new()}}{
Create a new nascSynth object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{nascSynth$new(
  data,
  time,
  id,
  treated,
  outcome,
  ci_width = 0.75,
  covariates = NULL,
  W = NULL,
  spatial_model = "none",
  rho = NULL,
  bias_correction = NULL,
  nasc_penalty = NULL,
  predictor_weights = NULL,
  predictors.op = "mean",
  special.predictors = NULL,
  time.predictors.prior = NULL
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{data}}{A data frame in long format containing the panel data.
Must include columns for the unit identifier, time, outcome, and
a binary treatment indicator (1 = treated, 0 = untreated).}

\item{\code{time}}{Unquoted column name for the time variable.}

\item{\code{id}}{Unquoted column name for the unit identifier.}

\item{\code{treated}}{Unquoted column name for the binary treatment indicator
(1 in treated periods for the treated unit, 0 otherwise).}

\item{\code{outcome}}{Unquoted column name for the outcome variable.}

\item{\code{ci_width}}{Numeric in \code{(0, 1)}. Width of the credible interval.
Default is \code{0.75}.}

\item{\code{covariates}}{Optional data frame in long format with columns for
\code{id}, \code{time}, and one or more predictor variables. Required
when \code{spatial_model} is \code{"SAR"} or \code{"SDM"}.}

\item{\code{W}}{Optional row-standardized spatial weights matrix. Row and column
names should match the unit identifiers in \code{data}. Required when
\code{bias_correction = TRUE} or \code{nasc_penalty = TRUE}.}

\item{\code{spatial_model}}{One of \code{"none"}, \code{"SAR"}, \code{"SDM"},
or \code{"exogenous"}.}

\item{\code{rho}}{Optional scalar in \code{(-1, 1)}. When supplied, Step 1 is
skipped and this value is used in Step 2 regardless of
\code{spatial_model}. Required when \code{spatial_model = "exogenous"}.}

\item{\code{bias_correction}}{Logical. If \code{TRUE}, the post-treatment
counterfactual is rescaled by \eqn{1/(1 - w's)}. Default is \code{TRUE}
for SAR/SDM/exogenous and \code{FALSE} for
\code{spatial_model = "none"}.}

\item{\code{nasc_penalty}}{Logical. If \code{TRUE}, the NASC penalty
\eqn{-\lambda \langle w, |s| \rangle} enters the likelihood. Same
model-aware default as \code{bias_correction}.}

\item{\code{predictor_weights}}{Optional non-negative numeric vector giving the
importance weight for each \emph{predictor row} in the Step-2 SCM
matching loss (the diagonal of Abadie-Diamond-Hainmueller's V matrix;
same role as the \code{vs} argument in the bsynth package). With the
default \code{predictors.op = "mean"} and no \code{special.predictors},
predictor rows are exactly the columns of \code{covariates} and the
length must match the number of covariate columns. When
\code{special.predictors} adds extra rows, the vector must cover all
predictor rows that result from the combined regular and special
specification. Two input forms are accepted: a \emph{named} numeric
vector (matched by row label, e.g. \code{c(x1 = 1, x2 = 0.5,
"Y_mean_1993" = 2)}) or an \emph{unnamed} numeric vector of the right
length (matched positionally, regular predictors first in the order
given in \code{covariates}, then special predictors in the order
given). \code{NULL} (default) gives equal weight (= 1) to every
predictor row. Setting an entry to 0 removes that row from the
matching loss without dropping it from the data. Outcome rows always
carry implicit weight 1; \code{predictor_weights} only re-scales the
predictor-matching contribution. Has no effect when neither
\code{covariates} nor \code{special.predictors} is supplied.}

\item{\code{predictors.op}}{Character string naming the operator
applied to every regular predictor (i.e. every non-id, non-time
column in \code{covariates}) over the pre-treatment period. Any
function name resolvable by \code{\link[base]{match.fun}} is
accepted; \code{"mean"} (default) and \code{"median"} are the
common choices. The chosen operator is applied unit-by-unit and
\code{NA}s are dropped before aggregation. Mirrors
\code{predictors.op} in the \pkg{Synth} package's
\code{dataprep()}.}

\item{\code{special.predictors}}{Optional list giving Synth-style
\dQuote{special predictors}: each predictor specified here
contributes one extra matching row, computed by applying a chosen
operator over a chosen subset of pre-treatment periods. Each entry
must be a list of length 3:
\enumerate{
\item the predictor, given either as a column name (character)
present in \code{data} or in \code{covariates}, or as a column
number referencing \code{data};
\item a numeric vector of time periods to aggregate over
(e.g. \code{1990:2000} or \code{1980});
\item a character string naming the operator (e.g. \code{"mean"},
\code{"median"}).
}
For example,
\preformatted{
special.predictors = list(
  list("x1", 1990:2000, "mean"),
  list("x2", 1980:1983, "median"),
  list("x3", 1980,      "mean")
)
}
specifies that \code{x1} is matched on its pre-period mean over
1990-2000, \code{x2} on its 1980-1983 median, and \code{x3} on its
1980 value alone. Lagged outcomes (\code{list("Y", 1993, "mean")})
are a typical use even when the outcome variable is not part of
\code{covariates}. Each entry produces a row labelled
\code{<var>_<op>_<period>} in the matching matrix; these labels are
also used to match \code{predictor_weights}. Mirrors
\code{special.predictors} in the \pkg{Synth} package's
\code{dataprep()}. Default \code{NULL} (no special predictors).}

\item{\code{time.predictors.prior}}{Optional numeric vector of
pre-treatment time periods over which regular predictors are
aggregated by \code{predictors.op}. Default \code{NULL} aggregates
over all pre-intervention periods (the previous behaviour). Has no
effect on \code{special.predictors}, whose time windows are set
per-entry. Mirrors \code{time.predictors.prior} in the \pkg{Synth}
package's \code{dataprep()}.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-nascSynth-fit"></a>}}
\if{latex}{\out{\hypertarget{method-nascSynth-fit}{}}}
\subsection{Method \code{fit()}}{
Fit the Stan model via MCMC.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{nascSynth$fit(
  n_samples = 100,
  n_samples_cap = 500L,
  n_samples_min = 30L,
  cores = parallel::detectCores() - 1,
  worker_iter = 2000L,
  worker_warmup = 1000L,
  ...
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n_samples}}{Number of rho draws to propagate to Step 2 when Step 1
runs. Either a positive integer (default 100), or the string
\code{"auto"}. When \code{"auto"}, the number is set to
\code{min(round(n_eff_rho), n_samples_cap)} after Step 1 completes,
floored at \code{n_samples_min}, so that the cut-posterior
approximation does not oversample a low-information rho posterior.
Ignored when an exogenous rho was supplied or when no rho is in use.}

\item{\code{n_samples_cap}}{Upper bound on the auto-selected \code{n_samples}.
Default 500. Only used when \code{n_samples = "auto"}.}

\item{\code{n_samples_min}}{Lower bound on the auto-selected \code{n_samples}.
Default 30. Only used when \code{n_samples = "auto"}.}

\item{\code{cores}}{Number of CPU cores for parallel execution.}

\item{\code{worker_iter}}{Iterations per worker chain in the multi-rho parallel
loop. Default 2000. Increase if Stan reports low Bulk/Tail ESS.}

\item{\code{worker_warmup}}{Warmup per worker chain in the multi-rho parallel
loop. Default 1000.}

\item{\code{...}}{Additional arguments forwarded to \code{\link[rstan:stanmodel-method-sampling]{rstan::sampling()}}.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-nascSynth-updateWidth"></a>}}
\if{latex}{\out{\hypertarget{method-nascSynth-updateWidth}{}}}
\subsection{Method \code{updateWidth()}}{
Update the credible interval width and recompute plot data.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{nascSynth$updateWidth(ci_width = 0.75)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{ci_width}}{Numeric in \code{(0, 1)}. New credible-interval width.
Default \code{0.75}.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-nascSynth-summary"></a>}}
\if{latex}{\out{\hypertarget{method-nascSynth-summary}{}}}
\subsection{Method \code{summary()}}{
Posterior summary of the fit: ATT, per-period treatment effects,
donor weights, model parameters and MCMC diagnostics. When the model
uses a network, the summary additionally reports the per-period
donor-average indirect (spillover) effect, the per-donor average
spillover, and the scalar average indirect effect (the spillover
analog of the ATT). Prints a formatted summary and invisibly returns
a \code{summary.nascSynth} list. Also dispatched from
\code{summary(obj)} via the S3 method \code{summary.nascSynth}.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{nascSynth$summary(ci_width = NULL, print = TRUE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{ci_width}}{Optional override for the credible interval width.}

\item{\code{print}}{Logical; if \code{TRUE} (default), print the formatted
summary to the console. Set to \code{FALSE} to retrieve the
structured list silently for programmatic use.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-nascSynth-indirectEffect"></a>}}
\if{latex}{\out{\hypertarget{method-nascSynth-indirectEffect}{}}}
\subsection{Method \code{indirectEffect()}}{
Posterior draws and summaries of the indirect (spillover) treatment
effect. By Proposition 6.2, the per-donor spillover at post-period
\eqn{t} is \eqn{\delta_{i,t}^{NASC} = s_i \cdot \tau_{1t}^{NASC}},
where \eqn{s = \rho (I_J - \rho W_J)^{-1} w_{J1}}.

Returns \code{NULL} when the model is configured without a network
(\code{uses_rho = FALSE}), since spillover is then identically zero
by construction.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{nascSynth$indirectEffect()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
Either \code{NULL} (no rho) or a named list with components:
\itemize{
\item \code{donor_names} -- character vector of donor IDs.
\item \code{time_post}  -- post-period time vector.
\item \code{delta_arr}  -- \code{[n_draws x T_post x J]} array
of per-draw, per-period, per-donor spillover.
\item \code{delta_total} -- \code{[n_draws x T_post]} matrix of
per-period total spillover (\eqn{\sum_i \delta_{i,t}}).
\item \code{avg_per_donor} -- \code{[n_draws x J]} matrix of
per-donor average spillover across post-periods.
\item \code{avg_total} -- \code{[n_draws]} vector of average
total spillover (the spillover analog of ATT).
}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-nascSynth-syntheticPlot"></a>}}
\if{latex}{\out{\hypertarget{method-nascSynth-syntheticPlot}{}}}
\subsection{Method \code{syntheticPlot()}}{
Base R plot of the observed and synthetic-control outcome trajectories,
with a shaded credible-interval band around the synthetic series and a
dotted vertical line at the intervention time.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{nascSynth$syntheticPlot()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-nascSynth-effectPlot"></a>}}
\if{latex}{\out{\hypertarget{method-nascSynth-effectPlot}{}}}
\subsection{Method \code{effectPlot()}}{
Plot the estimated direct treatment effect (\eqn{\tau_{1t}}) over
time as a line with a shaded credible-interval band, with a horizontal
zero reference and a dotted vertical line at the intervention time.
The indirect-effect overlay has been removed; for the per-period or
per-donor spillover see \code{$tauPlot()}, \code{$attPlot()}, or
\code{$indirectEffect()}.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{nascSynth$effectPlot(indirect = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{indirect}}{Ignored. Retained in the signature so that
existing scripts calling \code{effectPlot(indirect = ...)} do not
break; no indirect line will be drawn.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-nascSynth-posteriorPlot"></a>}}
\if{latex}{\out{\hypertarget{method-nascSynth-posteriorPlot}{}}}
\subsection{Method \code{posteriorPlot()}}{
Posterior density plots of the main scalar Bayesian parameters that
the model actually estimated. Depending on the configuration, panels
may include:
\itemize{
\item \code{rho} -- Step-1 spatial autocorrelation (SAR/SDM).
\item \code{theta[k]} -- Step-1 SDM spillover coefficients (one panel
per covariate, when \code{spatial_model = "SDM"}). Labelled by
covariate name when available.
\item \code{beta[k]} -- Step-1 SAR within-model covariate coefficients
(one panel per covariate). Labelled by covariate name when
available. Time-invariant covariates are not identified under the
within transform and are skipped.
\item \code{sigma_step1} -- Step-1 residual SD (\code{sigma_sar} or
\code{sigma_sdm}, surfaced under a unified label).
\item \code{lambda} -- NASC penalty strength (when
\code{nasc_penalty = TRUE}).
\item \code{sigma_step2} -- Step-2 residual SD on the standardized
outcome scale.
\item \code{bias_correction} -- Posterior of the multiplicative bias
factor \eqn{1/(1 - \langle w, s \rangle)} (when
\code{bias_correction = TRUE}). Skipped when bias correction is
off, since draws are constant at 1.
}
Each parameter is shown in its own panel as a color-filled density
curve. For more than four panels the layout switches to a 2-column
grid so densities stay readable.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{nascSynth$posteriorPlot()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-nascSynth-attPlot"></a>}}
\if{latex}{\out{\hypertarget{method-nascSynth-attPlot}{}}}
\subsection{Method \code{attPlot()}}{
Posterior density plot of the average treatment effect on the treated
(ATT), pooled across all post-treatment periods. When the model uses
a network, a second panel is added showing the density of the average
indirect (spillover) effect, defined as the per-draw average across
both donors and post-periods of \eqn{\delta_{i,t}^{NASC}} (the
spillover analog of the ATT, Proposition 6.2). This matches the
"Average indirect Treatment Effect" reported in \code{summary()}. Each
density is drawn as a color-filled curve with a dotted zero reference
line.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{nascSynth$attPlot(indirect = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{indirect}}{Logical or \code{NULL}. If \code{TRUE} (default
when the model uses a network), append a second panel for the
average indirect effect. Set to \code{FALSE} to suppress the
indirect panel even when a network is in use. \code{NULL} (default)
uses the model-aware default.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-nascSynth-tauPlot"></a>}}
\if{latex}{\out{\hypertarget{method-nascSynth-tauPlot}{}}}
\subsection{Method \code{tauPlot()}}{
Posterior density plots of the period-by-period treatment effects
for all post-treatment periods. When the model uses a network (i.e.
\code{uses_rho = TRUE}), two columns are drawn side by side: the
left column is the direct effect \eqn{\tau_{1t}} and the right
column is the donor-average indirect effect
\eqn{\bar\delta_t = (1/J)\sum_i \delta_{i,t}^{NASC}}, with one row
per post-period. When no network is in use, only the direct-effect
column is drawn (single-column layout). Panels within a column share
a common x-axis; the two columns use independent x-axes since direct
and indirect effects are typically on very different scales.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{nascSynth$tauPlot(indirect = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{indirect}}{Logical or \code{NULL}. If \code{TRUE} (default
when the model uses a network), also draw the per-period
donor-average indirect-effect column. Set to \code{FALSE} to
suppress the indirect column even when a network is in use.
\code{NULL} (default) uses the model-aware default.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-nascSynth-weightDraws"></a>}}
\if{latex}{\out{\hypertarget{method-nascSynth-weightDraws}{}}}
\subsection{Method \code{weightDraws()}}{
Ridgeline plot of the posterior weight density per donor. Densities
are stacked vertically, each labelled by donor name, and adjacent
ridges overlap by a configurable fraction so that each individual
density can be drawn taller without compressing the layout. Fills
match the steelblue palette of \code{$posteriorPlot()}.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{nascSynth$weightDraws(overlap = 0.5, scale = 1.4, fill_alpha = 0.85)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{overlap}}{Numeric in \code{[0, 1)}. Fraction of vertical
overlap between adjacent ridges. \code{0} reproduces the old
non-overlapping layout; \code{0.5} (default) makes each ridge
reach halfway into the row above; values close to \code{1}
approach a fully-stacked look. Default \code{0.5}.}

\item{\code{scale}}{Numeric > 0. Multiplicative height of every ridge
relative to the baseline step. Combined with \code{overlap} this
controls how prominent each density is. Default \code{1.4}.}

\item{\code{fill_alpha}}{Numeric in \code{(0, 1]}. Multiplier on the
default fill transparency (\eqn{\approx 0.4} alpha at
\code{fill_alpha = 1}). Lower values make overlapping ridges more
transparent. Default \code{0.85}, matching the look of
\code{$posteriorPlot()}.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-nascSynth-weightCorr"></a>}}
\if{latex}{\out{\hypertarget{method-nascSynth-weightCorr}{}}}
\subsection{Method \code{weightCorr()}}{
Heatmap of pairwise Pearson correlations between donor weights across
posterior draws. The diagonal is set to \code{NA} so the off-diagonal
structure is not visually dominated by the trivial unit correlations.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{nascSynth$weightCorr()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-nascSynth-clone"></a>}}
\if{latex}{\out{\hypertarget{method-nascSynth-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{nascSynth$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}

\name{nascWeight}
\alias{nascWeight}
\title{Posterior Donor-Weight Distributions Across Models Plot}
\usage{
nascWeight(
  models,
  overlap = 0.5,
  scale = 1.4,
  fill_alpha = 0.45,
  max_donors = NULL
)
}
\arguments{
\item{models}{Named list of fitted \code{nascSynth} objects; names must be
unique and are used as legend labels.}

\item{overlap}{Numeric in \code{[0, 1)}. Vertical overlap between adjacent
donor rows. Default \code{0.5}.}

\item{scale}{Numeric > 0. Ridge height relative to the row step. Default
\code{1.4}.}

\item{fill_alpha}{Numeric in \code{(0, 1]}. Fill transparency; lower values
help when many models overlap. Default \code{0.45}.}

\item{max_donors}{Integer or \code{NULL}. Keep only the \code{max_donors}
donors at the top of the plot (the largest IDs). Default \code{NULL} (all).}
}
\value{
Invisibly, a tibble with one row per donor-model pair giving the
posterior mean weight, its SD, and the share of draws above \code{1e-3},
ordered as plotted.
}
\description{
Ridgeline plot of posterior donor-weight densities: one row per donor, one
translucent density per model. Multi-model counterpart of
\code{$weightDraws()}, sharing the colour palette of \code{\link{nascPlot}}.
}
\details{
Donors are ordered by identifier with the largest ID on top, sorted
numerically where possible and lexicographically otherwise. Across differing
donor pools the union is shown; a model missing a donor leaves that row empty.
}
\seealso{
\code{\link{nascPlot}} for the time-series counterpart.
}

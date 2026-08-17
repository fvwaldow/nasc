\name{nascPlot}
\alias{nascPlot}
\title{Multiple NASC Model Plot}
\usage{
nascPlot(
  models,
  show_ci = FALSE,
  indirect = TRUE,
  show_avg = TRUE,
  indirect_pre = c("zero", "drop", "placebo")
)
}
\arguments{
\item{models}{Named list of fitted \code{nascSynth} objects; the names are the
legend labels.}

\item{show_ci}{Adds shaded credible-interval bands. Default
\code{FALSE}, since overlapping bands are hard to read.}

\item{indirect}{Adds a third panel with the indirect (spillover)
effects. Models without a network are skipped with a message. Default
\code{TRUE}.}

\item{show_avg}{Draws each model's donor average as a bold line in the
spillover panel. Default \code{TRUE}.}

\item{indirect_pre}{Pre-treatment display in the spillover panel:
\code{"zero"} holds the pre-period at 0 (default, since
\eqn{\delta = s\tau} and \eqn{\tau = 0} pre-treatment), \code{"placebo"}
shows the rescaled pre-treatment residual, \code{"drop"} starts at the
intervention.}
}
\value{
Invisibly \code{NULL}; called for its plots.
}
\description{
Overlays the outcome paths, direct effects and (optionally) spillovers of
several fitted \code{nascSynth} objects, e.g. conventional SC,
bias-corrected SC and NASC.
}
\details{
Panels are drawn in sequence on the active device. Set
\code{par(mfrow = ...)} beforehand to arrange them. Model colours come from
\code{hcl.colors(palette = "Dark 3")}, shared with \code{\link{nascWeight}}.
}
\seealso{
\code{\link{nascWeight}} for the donor-weight counterpart.
}

\name{effectGraph}
\alias{effectGraph}
\title{Network graph coloured by direct and indirect treatment effects}
\usage{
effectGraph(
  model,
  signed = TRUE,
  edge_threshold = 0.001,
  layout = NULL,
  vertex_size = 12,
  label_cex = 0.8,
  directed = FALSE,
  show_values = FALSE,
  digits = 2,
  seed = 1L
)
}
\arguments{
\item{model}{A fitted \code{nascSynth} object using a network; without one
the spillovers are zero by construction.}

\item{signed}{\code{TRUE} (default) colours signed effects on a
diverging palette centred at zero; \code{FALSE} uses absolute effects on a
sequential palette.}

\item{edge_threshold}{Edges with absolute weight below this are dropped.
Default \code{1e-3}.}

\item{layout}{An \code{igraph} layout function or a two-column coordinate
matrix. Default \code{igraph::layout_with_fr}.}

\item{vertex_size}{Node size; the treated unit is drawn 1.2x larger. Default
\code{12}.}

\item{label_cex}{Label size. Default \code{0.8}.}

\item{directed}{Treat \eqn{W} as directed and draw arrows. Default
\code{FALSE}.}

\item{show_values}{Append the numeric effect to each node label
(\eqn{\bar\delta_j} for donors, ATT for the treated unit). Default
\code{FALSE}.}

\item{digits}{Significant digits used when \code{show_values = TRUE}. Default
\code{2}.}

\item{seed}{Integer or \code{NULL}. Seed for stochastic layouts, applied
locally with the user's RNG state restored on exit. Ignored for a
precomputed layout. Default \code{1L}.}
}
\value{
Invisibly, a tibble with one row per node giving \code{effect_mean}
(\eqn{\bar\delta_j} for donors, ATT for the treated unit),
\code{abs_effect_mean} and \code{is_treated}.
}
\description{
Draws the spatial-weights network of a fitted \code{nascSynth} object with
each donor coloured by its posterior mean spillover \eqn{\bar\delta_j} over
the post-treatment periods and the treated node by the ATT. Since
\eqn{\delta_{j,t} = s_j \tau_{1t}}. Both effects share one colour scale, so the treated node is typically
strong and peripheral donors pale, and the contrast itself shows the direct
share of the total effect. Requires \pkg{igraph}.
}
\seealso{
\code{\link{contaminationGraph}} for the exposure-side analog.
}

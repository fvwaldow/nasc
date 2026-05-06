% Handwritten Manual
\name{contaminationGraph}
\alias{contaminationGraph}
\title{Network graph coloured by donor contamination}
\usage{
contaminationGraph(
  model,
  signed = FALSE,
  edge_threshold = 0.001,
  layout = NULL,
  vertex_size = 12,
  label_cex = 0.8,
  directed = FALSE,
  seed = 1L
)
}
\arguments{
\item{model}{A fitted \code{nascSynth} object.}

\item{signed}{Logical. If \code{TRUE}, nodes are coloured by signed
\eqn{s_j} on a diverging palette; if \code{FALSE} (default), by
\eqn{|s_j|} on a sequential palette.}

\item{edge_threshold}{Numeric. Edges with absolute weight below this
value are not drawn. Default \code{1e-3}.}

\item{layout}{An \code{igraph} layout function or a two-column matrix
of node coordinates. Default \code{igraph::layout_with_fr}.}

\item{vertex_size}{Numeric scalar for node size. Default \code{12}.}

\item{label_cex}{Numeric scalar for label size. Default \code{0.8}.}

\item{directed}{Logical. If \code{TRUE}, treats \eqn{W} as a directed
graph and draws arrows. Default \code{FALSE}.}

\item{seed}{Integer or \code{NULL}. Seed used to make stochastic
layouts (e.g. \code{layout_with_fr}) reproducible across calls.
The seed is applied locally and the user's RNG state is restored
on exit. Set to \code{NULL} for the original stochastic behaviour.
Ignored when \code{layout} is a precomputed coordinate matrix.
Default \code{1L}.}
}
\value{
Invisibly returns a tibble with one row per donor giving the
posterior mean of \eqn{s_j} and \eqn{|s_j|}.
}
\description{
Plots the spatial-weights network underlying a fitted \code{nascSynth}
object as a graph, with each donor node coloured by its posterior mean
contamination \eqn{|s_j|}. The treated unit is drawn as a distinct
marker and is not coloured on the contamination scale (its row of
\eqn{W} feeds donors but it has no weight of its own).

Edges follow the row-standardised weights matrix \eqn{W}. By default the
graph is treated as undirected and edges weaker than \code{edge_threshold}
are dropped to keep the plot readable.
}
\seealso{
\code{\link{effectGraph}} for the effect-side analog
(direct effect at the treated unit, indirect effect at every donor),
and \code{\link{contaminationScatter}} for a per-donor weight-vs-
contamination scatter plot.
}

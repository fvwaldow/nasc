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

\item{signed}{Colour nodes by signed \eqn{s_j} on a diverging
palette; \code{FALSE} (default) uses \eqn{|s_j|} on a sequential palette.}

\item{edge_threshold}{Edges with absolute weight below this are dropped.
Default \code{1e-3}.}

\item{layout}{An \code{igraph} layout function or a two-column coordinate
matrix. Default \code{igraph::layout_with_fr}.}

\item{vertex_size}{Node size. Default \code{12}.}

\item{label_cex}{Label size. Default \code{0.8}.}

\item{directed}{Treat \eqn{W} as directed and draw arrows. Default
\code{FALSE}.}

\item{seed}{Integer or \code{NULL}. Seed for stochastic layouts, applied
locally with the user's RNG state restored on exit. Ignored for a
precomputed layout. Default \code{1L}.}
}
\value{
Invisibly, a tibble with one row per donor giving the posterior means
of \eqn{s_j} and \eqn{|s_j|}.
}
\description{
Draws the spatial-weights network of a fitted \code{nascSynth} object, each
donor coloured by its posterior mean contamination \eqn{|s_j|}. The treated
unit gets a distinct marker and no contamination colour: its row of \eqn{W}
feeds the donors but it carries no weight itself. Requires \pkg{igraph}.
}
\seealso{
\code{\link{effectGraph}} for the effect-side analog.
}

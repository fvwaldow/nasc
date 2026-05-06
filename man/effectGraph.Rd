% Handwritten Manual
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
\item{model}{A fitted \code{nascSynth} object. The model must use a
network (\code{uses_rho = TRUE}); otherwise spillover is identically
zero by construction and only the treated unit would carry a
non-trivial colour.}

\item{signed}{Logical. If \code{TRUE} (default), nodes are coloured
by signed effects on a diverging blue-white-red palette centered at
zero -- the natural choice for treatment effects, where the sign
distinguishes positive from negative effects. If \code{FALSE}, by
absolute effects on a sequential palette.}

\item{edge_threshold}{Numeric. Edges with absolute weight below this
value are not drawn. Default \code{1e-3}.}

\item{layout}{An \code{igraph} layout function or a two-column matrix
of node coordinates. Default \code{igraph::layout_with_fr}.}

\item{vertex_size}{Numeric scalar for node size. The treated unit is
drawn 1.2x this size for visual identification. Default \code{12}.}

\item{label_cex}{Numeric scalar for label size. Default \code{0.8}.}

\item{directed}{Logical. If \code{TRUE}, treats \eqn{W} as a directed
graph and draws arrows. Default \code{FALSE}.}

\item{show_values}{Logical. If \code{TRUE}, append the numeric
effect to each node label (\eqn{\bar\delta_j} for donors, ATT for
the treated unit). Useful for small networks; can clutter larger
ones. Default \code{FALSE}.}

\item{digits}{Integer. Number of significant digits used when
\code{show_values = TRUE}. Default \code{2}.}

\item{seed}{Integer or \code{NULL}. Seed used to make stochastic
layouts (e.g. \code{layout_with_fr}) reproducible across calls.
The seed is applied locally and the user's RNG state is restored
on exit. Set to \code{NULL} for the original stochastic behaviour.
Ignored when \code{layout} is a precomputed coordinate matrix.
Default \code{1L}.}
}
\value{
Invisibly returns a tibble with one row per node giving the
posterior mean effect (\code{effect_mean} -- \eqn{\bar\delta_j}
for donors, ATT for the treated unit), its absolute value
(\code{abs_effect_mean}), and a logical flag \code{is_treated}.
}
\description{
Plots the spatial-weights network underlying a fitted \code{nascSynth}
object as a graph, with each donor node coloured by its posterior
mean indirect (spillover) effect \eqn{\bar\delta_j} averaged over
the post-treatment periods, and the treated node coloured by the
posterior-mean ATT (the direct effect averaged over post-treatment
periods). By Proposition 6.2,
\eqn{\delta_{j,t}^{NASC} = s_j \cdot \tau_{1t}^{NASC}}, so this plot
is the effect-side analog of \code{\link{contaminationGraph}}: where
\code{contaminationGraph} shows how \emph{exposed} each donor is to
the treated unit through the network, this plot shows the
\emph{realized} treatment effect at every node -- direct at the
treated unit, indirect at every donor.

Direct and indirect effects share a single colour scale so the
relative magnitudes are directly comparable on the figure. Because
the ATT is typically larger than any single donor's spillover, the
treated node usually appears strongly coloured while donors near
the periphery of the network appear pale -- the contrast itself
communicates the share of the total effect that is direct vs
indirect. The treated unit is also drawn as a slightly larger
square so it remains visually identifiable.

Edges follow the row-standardised weights matrix \eqn{W}. By default
the graph is treated as undirected and edges weaker than
\code{edge_threshold} are dropped to keep the plot readable.
}
\seealso{
\code{\link{contaminationGraph}} for the exposure-side
analog (donors coloured by \eqn{|s_j|}).
}

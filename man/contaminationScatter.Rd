% Handwritten Manual
\name{contaminationScatter}
\alias{contaminationScatter}
\title{Mean contamination vs mean weight scatter plot}
\usage{
contaminationScatter(model, signed = FALSE, label = TRUE, top_n = 10L)
}
\arguments{
\item{model}{A fitted \code{nascSynth} object.}

\item{signed}{Logical. If \code{TRUE}, plot signed \eqn{\bar s_j} on the
x-axis with a vertical reference line at zero; if \code{FALSE}
(default), plot \eqn{\overline{|s_j|}}.}

\item{label}{Logical. If \code{TRUE} (default), every marker is labelled
with its donor name. Set to \code{FALSE} to suppress all labels (useful
for very dense donor pools where labels would be unreadable).}

\item{top_n}{Integer. Retained for backward compatibility with earlier
versions of this function that highlighted the top contributors.
Currently unused -- all labels are drawn at uniform size and weight.
Default \code{10}.}
}
\value{
Invisibly returns a tibble with columns \code{donor},
\code{w_mean}, \code{s_mean}, \code{abs_s_mean}, sorted by descending
contribution \eqn{\bar w_j |\bar s_j|}.
}
\description{
Per-donor scatter plot of posterior-mean synthetic-control weight
(\eqn{\bar w_j}) against posterior-mean contamination
(\eqn{\bar s_j} or \eqn{\overline{|s_j|}}). Donors that combine large
weight with large contamination contribute the most to the NASC penalty
\eqn{\langle w, |s|\rangle} and to the bias-correction term \eqn{w's}
and are the diagnostic targets of this plot.
}
\details{
To reduce label overlap without adding a dependency on \pkg{ggrepel},
points in the upper half of the y-range are labelled below the marker
and points in the lower half above. For large donor pools where labels
still collide, set \code{label = FALSE}.
}
\seealso{
\code{\link{contaminationGraph}} for the network-graph view of
the same posterior-mean contamination values.
}

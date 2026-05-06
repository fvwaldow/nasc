% Handwritten Manual
\docType{package}
\name{nasc-package}
\alias{nasc}
\alias{nasc-package}
\title{nasc: Network-Aware Synthetic Control Estimator}
\description{
Implements the Network-Aware Synthetic Control (NASC) estimator for causal inference in panel data where the SUTVA is violated due to network or spatial interference. The package utilizes a modular Bayesian framework to estimate treatment effects and network dependencies while accounting for spillover magnitude and direction. Besides bias correction it features a novel penalized-weight regularization that discounts donor units based on their contamination risk. Supports Stan-based MCMC sampling for robust inference, covariate adjustment, and automated construction of network-aware counterfactuals using row-standardized spatial weight matrices.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/fvwaldow/nasc}
  \item Report bugs at \url{https://github.com/fvwaldow/nasc/issues}
}

}
\author{
\strong{Maintainer}: Frederik von Waldow \email{f.vonwaldow@uni-bonn.de} (\href{https://orcid.org/0009-0009-4269-2088}{ORCID})

}
\keyword{internal}

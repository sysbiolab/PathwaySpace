% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pspaceDecay.R
\name{weibullDecay}
\alias{weibullDecay}
\title{Weibull decay function}
\usage{
weibullDecay(x, signal = 1, decay = 0.001, shape = 1.05)
}
\arguments{
\item{x}{A numeric vector of distances (in [0,1]).}

\item{signal}{A single numeric value representing a signal.}

\item{decay}{A decay factor (in [0,1]). This term indicates how much the
\code{signal} decreases as a function of distance in pathway space. 
For example, at a specific distance defined by the \code{pdist} parameter 
(see \code{\link{circularProjection}}), the signal intensity will
be the initial signal multiplied by \code{decay}.}

\item{shape}{A parameter (>=1) of a Weibull function. When \code{shape=1}
the Weibull decay follows an exponential decay. When \code{shape>1}
the function is first convex, then concave with an inflection point.}
}
\value{
A numeric vector; if missing 'x', it will return decay function.
}
\description{
The \code{weibullDecay} function is used by
PathwaySpace's methods for signal convolution and projection.
}
\examples{
x <- seq(0, 2, 0.01)
y <- weibullDecay(x, signal = 1, decay = 0.5, shape = 4)
plot(x, y)

}
\seealso{
\code{\link{expDecay}}, \code{\link{linearDecay}}
}
\author{
Mauro Castro.
}

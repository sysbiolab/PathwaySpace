% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pspaceMethods.R
\docType{methods}
\name{vertexWeight,PathwaySpace-method}
\alias{vertexWeight,PathwaySpace-method}
\alias{vertexWeight}
\alias{vertexWeight<-}
\alias{vertexWeight<-,PathwaySpace-method}
\title{Accessor functions for fetching slots from a PathwaySpace object.}
\usage{
\S4method{vertexWeight}{PathwaySpace}(pts)

\S4method{vertexWeight}{PathwaySpace}(pts) <- value
}
\arguments{
\item{pts}{A \linkS4class{PathwaySpace} class object.}

\item{value}{A numeric vector with values representing vertex weights.
This vector should be aligned to the "vertex" slot.}
}
\value{
A numeric vector.
}
\description{
Get or set 'weights' for a \linkS4class{PathwaySpace}
class object.
}
\examples{
data('gtoy1', package = 'RGraphSpace')
pts <- buildPathwaySpace(gtoy1, nrc = 100)
vertexWeight(pts)

}

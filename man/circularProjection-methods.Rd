% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pspaceMethods.R
\docType{methods}
\name{circularProjection,PathwaySpace-method}
\alias{circularProjection,PathwaySpace-method}
\alias{circularProjection}
\title{Creating 2D-landscape images from graph objects.}
\usage{
\S4method{circularProjection}{PathwaySpace}(
  pts,
  kns = 8,
  pdist = 0.15,
  rescale = TRUE,
  verbose = TRUE,
  .decay_fun = weibullDecay,
  knn = deprecated(),
  ...
)
}
\arguments{
\item{pts}{A \linkS4class{PathwaySpace} class object.}

\item{kns}{A single positive integer determining the k-nearest signal 
sources used in the signal convolution operation.}

\item{pdist}{A term (in \code{[0,1]}) determining a distance unit for the
signal convolution related to the image space. This distance will affect the
extent over which the convolution operation projects the signal between
source- and destination points. The signal projection tends to zero at
`pdist`, and when `pdist = 1` the projection will extend to its maximum 
distance, defined by the diameter of the inscribed circle within the 
pathway space.}

\item{rescale}{A single logical value indicating whether to rescale 
the signal. If the signal \code{>=0}, then it will be rescaled to 
\code{[0,1]}; if the signal \code{<=0}, then it will be rescaled to 
\code{[-1,0]}; and if the signal in \code{(-Inf,+Inf)}, then it will be 
rescaled to \code{[-1,1]}.}

\item{verbose}{A single logical value specifying to display detailed 
messages (when \code{verbose=TRUE}) or not (when \code{verbose=FALSE}).}

\item{.decay_fun}{A signal decay function. Available: 'Weibull',
'exponential', and 'linear' functions (see \code{\link{weibullDecay}}).}

\item{knn}{Deprecated from PathwaySpace 1.0; use `kns` instead.}

\item{...}{Additional arguments passed to the decay function.}
}
\value{
A preprocessed \linkS4class{PathwaySpace} class object.
}
\description{
\code{circularProjection} implements a convolution
algorithm to project signal across a 2D-coordinate system.
}
\examples{
# Load a demo igraph
data('gtoy1', package = 'RGraphSpace')

# Create a new PathwaySpace object
pts <- buildPathwaySpace(gtoy1, nrc = 100)
# note: adjust 'nrc' to increase image resolution

# Create a 2D-landscape image
pts <- circularProjection(pts)

}
\seealso{
\code{\link{buildPathwaySpace}}
}
\author{
Victor Apolonio, Vinicius Chagas, Mauro Castro, 
and TCGA Network.
}

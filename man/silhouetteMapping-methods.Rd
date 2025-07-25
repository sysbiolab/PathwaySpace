% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pspaceMethods.R
\docType{methods}
\name{silhouetteMapping,PathwaySpace-method}
\alias{silhouetteMapping,PathwaySpace-method}
\alias{silhouetteMapping}
\title{Decorating PathwaySpace Images with Graph Silhouettes}
\usage{
\S4method{silhouetteMapping}{PathwaySpace}(
  ps,
  baseline = 0.01,
  pdist = 0.05,
  fill.cavity = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{ps}{A \linkS4class{PathwaySpace} class object.}

\item{baseline}{A fraction (in \code{[0,1]}) of the signal scale of a 
PathwaySpace image. This term only affects the image baseline projection, 
which represents a silhouette of the graph's layout outlined in the 
resulting image. When \code{baseline = 0} (i.e. lower level of the signal 
scale), the baseline will extend over the entire image space, so no 
silhouette will be visible.}

\item{pdist}{A term (in \code{[0,1]}) determining a distance unit for the
signal decay function. This distance will affect the extent over which 
the convolution operation projects the image baseline.}

\item{fill.cavity}{A single logical value specifying to fill cavities 
in the silhouette mask (when \code{verbose=TRUE}) or not 
(when \code{verbose=FALSE}).}

\item{verbose}{A single logical value specifying to display detailed 
messages (when \code{verbose=TRUE}) or not (when \code{verbose=FALSE}).}
}
\value{
A preprocessed \linkS4class{PathwaySpace} class object.
}
\description{
\code{silhouetteMapping} constructs an image baseline used
to outline the graph layout in a PathwaySpace image.
}
\examples{
# Load a demo igraph
data('gtoy1', package = 'RGraphSpace')

# Create a new PathwaySpace object
ps <- buildPathwaySpace(gtoy1, nrc = 100)
# note: adjust 'nrc' to increase image resolution

# Set '1s' as vertex signal
vertexSignal(ps) <- 1

# Map graph silhouette
ps <- silhouetteMapping(ps, pdist = 0.1)

}
\seealso{
\code{\link{circularProjection}}
}
\author{
Mauro Castro and TCGA Network.
}

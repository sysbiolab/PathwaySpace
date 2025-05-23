% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pspacePlots.R
\docType{methods}
\name{plotPathwaySpace,PathwaySpace-method}
\alias{plotPathwaySpace,PathwaySpace-method}
\alias{plotPathwaySpace}
\title{Plotting 2D-landscape images for the PathwaySpace package.}
\usage{
\S4method{plotPathwaySpace}{PathwaySpace}(
  pts,
  colors = pspace.cols(),
  trim.colors = c(3, 2, 1, 2, 3),
  bg.color = "grey85",
  theme.name = c("th0", "th1", "th2", "th3"),
  title = "PathwaySpace",
  font.size = 1,
  font.color = "white",
  xlab = "Pathway coordinates 1",
  ylab = "Pathway coordinates 2",
  zlab = "Density",
  zlim = NULL,
  slices = 25,
  add.grid = TRUE,
  grid.color = "white",
  add.contour = TRUE,
  contour.color = "white",
  label.summits = TRUE,
  marks = FALSE,
  mark.size = 3,
  mark.color = "white",
  mark.padding = 0.5,
  mark.line.width = 0.5,
  use.dotmark = FALSE
)
}
\arguments{
\item{pts}{A \linkS4class{PathwaySpace} class object.}

\item{colors}{A vector of colors. Each color is a specific tone used to 
create a customized color palette that matches the signal type projected 
in the image space, such as binary or numeric values in (-Inf, +Inf). The 
palette will be created with 5-color tones, interpolated according to the
provided sequence in the vector of colors. The proportion of each color 
tone can be adjusted by the 'trim.colors' argument. The primary objective 
of this color palette construction is to enhance the visibility of summits 
and valleys in the image space. If you wish to bypass the automatic palette 
generation and use the 'colors' input as-is, simply set 'trim.colors' 
to NULL.}

\item{trim.colors}{An vector with 5 positive integer numbers. This argument
can be used to adjust the proportion of each color tone in the palette.}

\item{bg.color}{A single color for background.}

\item{theme.name}{Name of a custom PathwaySpace theme. These themes 
(from 'th0' to 'th3') consist mainly of preconfigured ggplot settings, 
which the user can subsequently fine-tune within the resulting 
ggplot object.}

\item{title}{A string for the title.}

\item{font.size}{A single numeric value passed to ggplot themes.}

\item{font.color}{A single color passed to ggplot themes.}

\item{xlab}{The title for the 'x' axis of a 2D-image space.}

\item{ylab}{The title for the 'y' axis of a 2D-image space.}

\item{zlab}{The title for the 'z' axis of the image signal.}

\item{zlim}{The 'z' limits of the plot (a numeric vector with two numbers).
If NULL, limits are determined from the range of the input values.}

\item{slices}{A single positive integer value used to split 
the image signal into equally-spaced intervals.}

\item{add.grid}{A logical value indicating whether to add gridlines to 
the image space. However, gridlines will only appear when the image 
is decorated with graph silhouettes (see \code{\link{silhouetteMapping}}).}

\item{grid.color}{A color passed to \code{\link[ggplot2]{geom_point}}.}

\item{add.contour}{A logical value indicating whether to add contour 
lines to 'summits' (see \code{\link{summitMapping}}).}

\item{contour.color}{A color passed to \code{\link[ggplot2]{geom_tile}}.}

\item{label.summits}{A logical value indicating whether to label summits,
(when summits are available).}

\item{marks}{A logical value indicating whether to add 'marks' to vertex 
positions. Alternatively, this could be a vector listing vertex names.}

\item{mark.size}{A font size argument passed to 
\code{\link[ggplot2]{geom_text}}.}

\item{mark.color}{A color passed to \code{\link[ggrepel]{geom_text_repel}}.}

\item{mark.padding}{A box padding argument passed to 
\code{\link[ggrepel]{geom_text_repel}}.}

\item{mark.line.width}{A line width argument passed to 
\code{\link[ggrepel]{geom_text_repel}}.}

\item{use.dotmark}{A logical value indicating whether "marks" should be 
represented as dots.}
}
\value{
A ggplot-class object.
}
\description{
\code{plotPathwaySpace} is a wrapper function to 
create dedicated ggplot graphics for PathwaySpace-class objects.
}
\examples{
# Load a demo igraph
data('gtoy1', package = 'RGraphSpace')

# Create a new PathwaySpace object
pts <- buildPathwaySpace(gtoy1, nrc = 100)
# note: adjust 'nrc' to increase image resolution

# Create a 2D-landscape image
pts <- circularProjection(pts)

# Map graph silhouette
pts <- silhouetteMapping(pts)

# Plot a 2D-landscape image
plotPathwaySpace(pts)

}
\seealso{
\code{\link{circularProjection}}
}
\author{
Vinicius Chagas, Victor Apolonio, Mauro Castro, 
and TCGA Network.
}

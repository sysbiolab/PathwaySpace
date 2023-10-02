## ----Load packages - quick start, eval=TRUE, message=FALSE--------------------
#--- Load required packages for this section
library(PathwaySpace)
library(igraph)

## ----Making a toy igraph - 1, eval=TRUE, message=FALSE------------------------
# Make a 'toy' undirected igraph
gtoy1 <- make_star(5, mode="undirected")

# Assign xy coordinates to each vertex
V(gtoy1)$x <- c(0, 1.5, -4, -4, -9)
V(gtoy1)$y <- c(0, 0,  4, -4,  0)

# Assign a name to each vertex (here, from n1 to n5)
V(gtoy1)$name <- paste0("n", 1:5)

## ----node size, eval=TRUE, message=FALSE, echo=FALSE--------------------------
V(gtoy1)$size <- 20

## ----Making a toy igraph - 2, eval=TRUE, message=FALSE, out.width="30%"-------
# Check the graph layout
plot.igraph(gtoy1)

## ----PathwaySpace constructor, eval=TRUE, message=FALSE-----------------------
# Run the PathwaySpace constructor
pspace1 <- buildPathwaySpace(gtoy1, mar = 0.2)

## ----Setting vertex signal - 1, eval=TRUE, message=FALSE----------------------
# Check the number of vertices in the PathwaySpace object
length(pspace1)

# Check vertex names
names(pspace1)

# Check signal (initialized with '0')
vertexSignal(pspace1)

## ----Setting vertex signal - 2, eval=TRUE, message=FALSE----------------------
# Set new signal to all vertices
vertexSignal(pspace1) <- c(1, 3, 2, 3, 2)

# Set a new signal to the 1st vertex
vertexSignal(pspace1)[1] <- 2

# Set a new signal to vertex "n1"
vertexSignal(pspace1)["n1"] <- 4

# Check updated signal values
vertexSignal(pspace1)

## ----Circular projection - 1, eval=FALSE, message=FALSE, out.width="50%"------
#  # Run network signal projection
#  pspace1 <- circularProjection(pspace1, knn = 1, pdist = 0.4)
#  
#  # Plot a PathwaySpace image
#  plotPathwaySpace(pspace1, marks = TRUE)

## ----Circular projection - 2, eval=FALSE, message=FALSE, out.width="50%"------
#  # Re-run the network signal projection with 'knn = 2'
#  pspace1 <- circularProjection(pspace1, knn = 2, pdist = 0.4)
#  
#  # Plot the PathwaySpace image
#  plotPathwaySpace(pspace1, marks = c("n3","n4"), theme = "th2")

## ----Circular projection - 3, eval=FALSE, message=FALSE, out.width="50%"------
#  # Re-run the network signal projection, passing 'shape' to the decay function
#  pspace1 <- circularProjection(pspace1, knn = 2, pdist = 0.2, shape = 2)
#  
#  # Plot the PathwaySpace image
#  plotPathwaySpace(pspace1, marks = "n1", theme = "th2")

## ----Polar projection - 1, eval=TRUE, message=FALSE, out.width="25%"----------
# Load a pre-processed directed igraph object
data("gtoy2")

# Check the graph layout
plot.igraph(gtoy2)

## ----Polar projection - 2, eval=TRUE, message=FALSE---------------------------
# Build a PathwaySpace for the 'gtoy2' igraph
pspace2 <- buildPathwaySpace(gtoy2, mar = 0.2)

# Set '1s' as vertex signal
vertexSignal(pspace2) <- 1

## ----Polar projection - 3, eval=FALSE, message=FALSE, out.width="50%"---------
#  # Run the network signal projection using polar coordinates
#  pspace2 <- polarProjection(pspace2, knn = 2, theta = 45, shape = 2)
#  
#  # Plot the PathwaySpace image
#  plotPathwaySpace(pspace2, theme = "th2", marks = TRUE)

## ----Polar projection - 4, eval=FALSE, message=FALSE, out.width="50%"---------
#  # Re-run the network signal projection using 'directional = TRUE'
#  pspace2 <- polarProjection(pspace2, knn = 2, theta = 45, shape = 2,
#    directional = TRUE)
#  
#  # Plot the PathwaySpace image
#  plotPathwaySpace(pspace2, theme = "th2", marks = c("n1","n3","n4","n5"))

## ----Signal types, eval=FALSE, message=FALSE, out.width="50%"-----------------
#  # Set a negative signal to vertices "n3" and "n4"
#  vertexSignal(pspace1)[c("n3","n4")] <- c(-2, -4)
#  
#  # Check updated signal vector
#  vertexSignal(pspace1)
#  #  4  3 -2 -4  2
#  
#  # Re-run the network signal projection
#  pspace1 <- circularProjection(pspace1, knn = 2, shape = 2)
#  
#  # Plot the PathwaySpace image
#  plotPathwaySpace(pspace1, bg.color = "white", font.color = "grey20",
#    marks = TRUE, mark.color = "magenta")

## ----label='Session information', eval=TRUE, echo=FALSE-----------------------
sessionInfo()


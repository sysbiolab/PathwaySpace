## ----Load packages - quick start, eval=TRUE, message=FALSE--------------------
#--- Load required packages for this section
library(PathwaySpace)
library(igraph)
library(ggplot2)

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

## ----Making a toy igraph - 2, eval=FALSE, message=FALSE, out.width="100%"-----
#  # Check the graph layout
#  plot.igraph(gtoy1)

## ----PathwaySpace constructor, eval=TRUE, message=FALSE-----------------------
# Run the PathwaySpace constructor
pspace1 <- buildPathwaySpace(gtoy1, mar = 0.2)

## ----Setting vertex signal - 1, eval=TRUE, message=FALSE, results='hide'------
# Check the number of vertices in the PathwaySpace object
length(pspace1)
## [1] 5

# Check vertex names
names(pspace1)
## [1] "n1" "n2" "n3" "n4" "n5"

# Check signal (initialized with '0')
vertexSignal(pspace1)
## n1 n2 n3 n4 n5 
##  0  0  0  0  0

## ----Setting vertex signal - 2, eval=TRUE, message=FALSE, results='hide'------
# Set new signal to all vertices
vertexSignal(pspace1) <- c(1, 3, 2, 3, 2)

# Set a new signal to the 1st vertex
vertexSignal(pspace1)[1] <- 2

# Set a new signal to vertex "n1"
vertexSignal(pspace1)["n1"] <- 4

# Check updated signal values
vertexSignal(pspace1)
## n1 n2 n3 n4 n5 
##  4  3  2  3  2

## ----Circular projection - 1, eval=FALSE, message=FALSE, out.width="70%"------
#  # Run network signal projection
#  pspace1 <- circularProjection(pspace1, knn = 1, pdist = 0.4)
#  
#  # Plot a PathwaySpace image
#  plotPathwaySpace(pspace1, marks = TRUE)

## ----Circular projection - 2, eval=FALSE, message=FALSE, out.width="70%"------
#  # Re-run the network signal projection with 'knn = 2'
#  pspace1 <- circularProjection(pspace1, knn = 2, pdist = 0.4)
#  
#  # Plot the PathwaySpace image
#  plotPathwaySpace(pspace1, marks = c("n3","n4"), theme = "th2")

## ----Circular projection - 3, eval=FALSE, message=FALSE, out.width="70%"------
#  # Re-run the network signal projection, passing 'shape' to the decay function
#  pspace1 <- circularProjection(pspace1, knn = 2, pdist = 0.2, shape = 2)
#  
#  # Plot the PathwaySpace image
#  plotPathwaySpace(pspace1, marks = "n1", theme = "th2")

## ----Polar projection - 1, eval=TRUE, message=FALSE, out.width="100%"---------
# Load a pre-processed directed igraph object
data("gtoy2", package = "PathwaySpace")

## ----Polar projection - 2, eval=FALSE, message=FALSE, out.width="100%"--------
#  # Check the graph layout
#  plot.igraph(gtoy2)

## ----Polar projection - 3, eval=TRUE, message=FALSE---------------------------
# Build a PathwaySpace for the 'gtoy2' igraph
pspace2 <- buildPathwaySpace(gtoy2, mar = 0.2)

# Set '1s' as vertex signal
vertexSignal(pspace2) <- 1

## ----Polar projection - 4, eval=FALSE, message=FALSE, out.width="70%"---------
#  # Run the network signal projection using polar coordinates
#  pspace2 <- polarProjection(pspace2, knn = 2, theta = 45, shape = 2)
#  
#  # Plot the PathwaySpace image
#  plotPathwaySpace(pspace2, theme = "th2", marks = TRUE)

## ----Polar projection - 6, eval=FALSE, message=FALSE, out.width="70%"---------
#  # Re-run the network signal projection using 'directional = TRUE'
#  pspace2 <- polarProjection(pspace2, knn = 2, theta = 45, shape = 2,
#    directional = TRUE)
#  
#  # Plot the PathwaySpace image
#  plotPathwaySpace(pspace2, theme = "th2", marks = c("n1","n3","n4","n5"))

## ----Signal types, eval=FALSE, message=FALSE, out.width="70%"-----------------
#  # Set a negative signal to vertices "n3" and "n4"
#  vertexSignal(pspace1)[c("n3","n4")] <- c(-2, -4)
#  
#  # Check updated signal vector
#  vertexSignal(pspace1)
#  # n1 n2 n3 n4 n5
#  #  4  3 -2 -4  2
#  
#  # Re-run the network signal projection
#  pspace1 <- circularProjection(pspace1, knn = 2, shape = 2)
#  
#  # Plot the PathwaySpace image
#  plotPathwaySpace(pspace1, bg.color = "white", font.color = "grey20",
#    marks = TRUE, mark.color = "magenta", theme = "th2")

## ----Load packages - case study, eval=FALSE, message=FALSE--------------------
#  #--- Load required packages for this section
#  library(PathwaySpace)
#  library(RGraphSpace)
#  library(igraph)
#  library(ggplot2)

## ----PathwaySpace decoration - 1, eval=TRUE, message=FALSE, results='hide'----
# Load a large igraph object
data("PCv12_pruned_igraph", package = "PathwaySpace")

# Check number of vertices
length(PCv12_pruned_igraph)
# [1] 12990

# Check vertex names
head(V(PCv12_pruned_igraph)$name)
# [1] "A1BG" "AKT1" "CRISP3" "GRB2" "PIK3CA" "PIK3R1"

# Get top-connected nodes for visualization
top10hubs <- igraph::degree(PCv12_pruned_igraph)
top10hubs <- names(sort(top10hubs, decreasing = TRUE)[1:10])
head(top10hubs)
# [1] "GNB1" "TRIM28" "RPS27A" "CTNNB1" "TP53" "ACTB"

## ----PathwaySpace decoration - 2, eval=FALSE, message=FALSE-------------------
#  ## Visualize the graph layout labeled with 'top10hubs' nodes
#  plotGraphSpace(PCv12_pruned_igraph, marks = top10hubs,
#    mark.color = "blue", theme = "th3")

## ----PathwaySpace decoration - 3, eval=FALSE, message=FALSE-------------------
#  # Load a list with Hallmark gene sets
#  data("Hallmarks_v2023_1_Hs_symbols", package = "PathwaySpace")
#  
#  # There are 50 gene sets in "hallmarks"
#  length(hallmarks)
#  # [1] 50
#  
#  # We will use the 'HALLMARK_P53_PATHWAY' (n=200 genes) for demonstration
#  length(hallmarks$HALLMARK_P53_PATHWAY)
#  # [1] 200

## ----PathwaySpace decoration - 4, eval=FALSE, message=FALSE-------------------
#  # Run the PathwaySpace constructor
#  pspace_PCv12 <- buildPathwaySpace(g=PCv12_pruned_igraph, nrc=500)
#  # Note: 'nrc' sets the number of rows and columns of the
#  # image space, which will affect the image resolution (in pixels)

## ----PathwaySpace decoration - 5, eval=FALSE, message=FALSE-------------------
#  # Intersect Hallmark genes with the PathwaySpace
#  hallmarks <- lapply(hallmarks, intersect, y = names(pspace_PCv12) )
#  
#  # After intersection, the 'HALLMARK_P53_PATHWAY' dropped to n=173 genes
#  length(hallmarks$HALLMARK_P53_PATHWAY)
#  # [1] 173
#  
#  # Set a binary signal (1s) to 'HALLMARK_P53_PATHWAY' genes
#  vertexSignal(pspace_PCv12) <- 0
#  vertexSignal(pspace_PCv12)[ hallmarks$HALLMARK_P53_PATHWAY ] <- 1

## ----PathwaySpace decoration - 6, eval=FALSE, message=FALSE-------------------
#  # Run network signal projection
#  pspace_PCv12 <- circularProjection(pspace_PCv12)
#  plotPathwaySpace(pspace_PCv12, title="HALLMARK_P53_PATHWAY",
#    marks = top10hubs, mark.size = 2, theme = "th3")

## ----PathwaySpace decoration - 7, eval=FALSE, message=FALSE-------------------
#  # Add silhouettes
#  pspace_PCv12 <- silhouetteMapping(pspace_PCv12)
#  plotPathwaySpace(pspace_PCv12, title="HALLMARK_P53_PATHWAY",
#    marks = top10hubs, mark.size = 2, theme = "th3")

## ----PathwaySpace decoration - 9, eval=FALSE, message=FALSE-------------------
#  # Mapping summits
#  pspace_PCv12 <- summitMapping(pspace_PCv12, minsize = 50)
#  plotPathwaySpace(pspace_PCv12, title="HALLMARK_P53_PATHWAY", theme = "th3")

## ----PathwaySpace decoration - 10, eval=FALSE, message=FALSE------------------
#  # Extracting summits from a PathwaySpace
#  summits <- getPathwaySpace(pspace_PCv12, "summits")
#  class(summits)
#  # [1] "list"

## ----label='Session information', eval=TRUE, echo=FALSE-----------------------
sessionInfo()


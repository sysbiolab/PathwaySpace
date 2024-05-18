#--- Load required packages for this section
library(PathwaySpace)
library(igraph, include.only = c("make_star", "V"))
library(ggplot2, include.only = c("plot.igraph"))

test_that("Execution test of circularProjection", {

#--- Build a toy graph
  # Make a 'toy' undirected igraph
  toy_graph <- make_star(5, mode="undirected")
  # Assign xy coordinates to each vertex
  V(toy_graph)$x <- c(0, 1.5, -4, -4, -9)
  V(toy_graph)$y <- c(0, 0,  4, -4,  0)
  # Assign a name to each vertex (here, from n1 to n5)
  V(toy_graph)$name <- paste0("n", 1:5)
#--- Build a pathwayspace
  pspace_toy <- buildPathwaySpace(toy_graph, mar = 0.2, verbose = FALSE)
  pspace_projection <- circularProjection(pspace_toy, knn = 1, pdist = 0.4, verbose = FALSE)

  expect_equal(class(pspace_projection@gxyz), c("matrix", "array"))
})

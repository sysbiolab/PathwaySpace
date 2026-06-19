# Accessors for Fetching Slots from a PathwaySpace Object

`getPathwaySpace` retrives information from individual slots available
in a PathwaySpace object.

## Usage

``` r
# S4 method for class 'PathwaySpace'
getPathwaySpace(ps, what = "status")
```

## Arguments

- ps:

  A preprocessed
  [PathwaySpace](https://github.com/sysbiolab/PathwaySpace/reference/PathwaySpace-class.md)
  class object

- what:

  A character value specifying which information should be retrieved
  from the slots. Options: "nodes", "edges", "graph", "image", "pars",
  "misc", "signal","projection", "status", "silhouette", "summits",
  "summit_mask", "summit_contour"

## Value

Content from slots in the
[PathwaySpace](https://github.com/sysbiolab/PathwaySpace/reference/PathwaySpace-class.md)
object.

## Examples

``` r
library(PathwaySpace)

# Load a demo igraph
data('gtoy1', package = 'RGraphSpace')

# Create a new PathwaySpace object
ps <- buildPathwaySpace(gtoy1, nrc = 100)
#> Validating arguments...
#> Validating the 'igraph' object...
#> Normalizing node coordinates to graph space...
#> Creating a 'PathwaySpace' object...
# note: adjust 'nrc' to increase image resolution

# Get the 'status' slot in ps
status <- getPathwaySpace(ps, what = 'status')
```

# Accessor Functions for PathwaySpace Objects

Get or set edge and vertex attributes in
[PathwaySpace](https://github.com/sysbiolab/PathwaySpace/reference/PathwaySpace-class.md)
class object.

## Usage

``` r
# S4 method for class 'PathwaySpace'
gs_vertex_attr(x, name, ...) <- value

# S4 method for class 'PathwaySpace'
gs_edge_attr(x, name, ...) <- value
```

## Arguments

- x:

  A
  [PathwaySpace](https://github.com/sysbiolab/PathwaySpace/reference/PathwaySpace-class.md)
  class object.

- name:

  Name of the attribute.

- ...:

  Additional arguments passed to igraph methods.

- value:

  The new value of the attribute.

## Value

Updated
[PathwaySpace](https://github.com/sysbiolab/PathwaySpace/reference/PathwaySpace-class.md)
object.

## Examples

``` r
library(PathwaySpace)

# Load a demo igraph
data('gtoy1', package = 'RGraphSpace')
ps <- buildPathwaySpace(gtoy1, nrc = 100)
#> Validating arguments...
#> Validating the 'igraph' object...
#> Normalizing node coordinates to graph space...
#> Creating a 'PathwaySpace' object...

# Get vertex count
gs_vcount(ps)
#> [1] 5

# Get edge count
gs_ecount(ps)
#> [1] 4

# Access a specific vertex attribute
gs_vertex_attr(ps, "signal")
#> n1 n2 n3 n4 n5 
#>  0  0  0  0  0 

# Replace an entire vertex attribute
gs_vertex_attr(ps, "signal") <- 1

# Modify a single value within a vertex attribute
gs_vertex_attr(ps, "signal")["n1"] <- 1

# Access a specific edge attribute
gs_edge_attr(ps, "weight")
#> [1] 1 1 1 1

# Replace an entire edge attribute
gs_edge_attr(ps, "weight") <- 1
```

# Constructor of PathwaySpace-class Objects

`buildPathwaySpace` is a constructor of PathwaySpace-class objects.

## Usage

``` r
buildPathwaySpace(gs, nrc = 500, verbose = TRUE)
```

## Arguments

- gs:

  A
  [`GraphSpace`](https://sysbiolab.github.io/RGraphSpace/reference/GraphSpace-methods.html)
  object. Alternatively, an
  [`igraph`](https://r.igraph.org/reference/aaa-igraph-package.html)
  object with node coordinates assigned to `x` and `y` vertex
  attributes, and node labels assigned to `name` vertex attribute.

- nrc:

  A single positive integer indicating the number of rows and columns
  (in pixels) for a square image matrix. This argument will affect the
  resulting image size and resolution.

- verbose:

  A logical value specifying to display detailed messages (when
  `verbose=TRUE`) or not (when `verbose=FALSE`).

## Value

A pre-processed
[PathwaySpace](https://github.com/sysbiolab/PathwaySpace/reference/PathwaySpace-class.md)
class object.

## See also

[`circularProjection`](https://github.com/sysbiolab/PathwaySpace/reference/circularProjection-methods.md),
[`polarProjection`](https://github.com/sysbiolab/PathwaySpace/reference/polarProjection-methods.md)

## Author

Sysbiolab Team

## Examples

``` r
library(PathwaySpace)

# Load a demo igraph
data('gtoy1', package = 'RGraphSpace')

# Check graph validity
gs <- GraphSpace(gtoy1)
#> Validating the 'igraph' object...
#> Creating a 'GraphSpace' object...

gs <- normalizeGraphSpace(gs)
#> Normalizing node coordinates to graph space...

# Create a new PathwaySpace object
ps <- buildPathwaySpace(gs, nrc = 100)
#> Validating arguments...
#> Creating a 'PathwaySpace' object...
# note: adjust 'nrc' to increase image resolution
```

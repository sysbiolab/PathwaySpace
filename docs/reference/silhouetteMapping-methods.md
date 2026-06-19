# Decorating PathwaySpace Images with Graph Silhouettes

`silhouetteMapping` constructs an image baseline used to outline the
graph layout in a PathwaySpace image.

## Usage

``` r
# S4 method for class 'PathwaySpace'
silhouetteMapping(
  ps,
  pdist = 0.05,
  baseline = 0.01,
  fill.cavity = TRUE,
  verbose = TRUE
)
```

## Arguments

- ps:

  A
  [PathwaySpace](https://github.com/sysbiolab/PathwaySpace/reference/PathwaySpace-class.md)
  class object.

- pdist:

  A term (in `[0,1]`) determining a distance unit for the silhouette
  projection.

- baseline:

  A fraction (in `[0,1]`) of the silhouette projection, representing the
  level over which a silhouette will outline the graph layout. When
  `baseline = 0` (i.e. lower level of the projection), the silhouette
  will extend over the entire image space, so no outline will be
  visible.

- fill.cavity:

  A logical value specifying to fill cavities in the silhouette mask
  (when `fill.cavity=TRUE`) or not (when `fill.cavity=FALSE`).

- verbose:

  A logical value specifying to display detailed messages (when
  `verbose=TRUE`) or not (when `verbose=FALSE`).

## Value

A preprocessed
[PathwaySpace](https://github.com/sysbiolab/PathwaySpace/reference/PathwaySpace-class.md)
class object.

## See also

[`circularProjection`](https://github.com/sysbiolab/PathwaySpace/reference/circularProjection-methods.md)

## Author

Sysbiolab Team

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

# Set '1s' as vertex signal
vertexSignal(ps) <- 1

# Map graph silhouette
ps <- silhouetteMapping(ps, pdist = 0.1)
#> Validating arguments...
#> Mapping graph silhouette...
#> Mapping 'x' and 'y' coordinates...
#> Silhouette: 8.12% of the landscape area!
```

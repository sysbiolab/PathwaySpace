# Mapping Summits on PathwaySpace Images

The `summitMapping` method implements a segmentation strategy to
identify summits on a 2D-landscape image (see
[`summitWatershed`](https://github.com/sysbiolab/PathwaySpace/reference/summitWatershed.md)).

## Usage

``` r
# S4 method for class 'PathwaySpace'
summitMapping(
  ps,
  maxset = 30,
  minsize = 30,
  threshold = 0.5,
  segm.fun = summitWatershed,
  ...
)
```

## Arguments

- ps:

  A
  [PathwaySpace](https://github.com/sysbiolab/PathwaySpace/reference/PathwaySpace-class.md)
  class object.

- maxset:

  A single positive integer indicating the maximum number of summits to
  be returned by the segmentation function.

- minsize:

  A single positive integer indicating the minimum size of the summits.

- threshold:

  A threshold provided as a fraction (in `[0,1]`) of the max signal
  intensity.

- segm.fun:

  A segmentation function used to detect summits (see
  [`summitWatershed`](https://github.com/sysbiolab/PathwaySpace/reference/summitWatershed.md)).

- ...:

  Additional arguments passed to the segmentation function.

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

# Load a large igraph
data("PCv12_pruned_igraph", package = "PathwaySpace")

# Continue this example from the PathwaySpace vignette,
# in the 'PathwaySpace decoration' section
```

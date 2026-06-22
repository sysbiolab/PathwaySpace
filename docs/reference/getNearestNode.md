# getNearestNode

Retrieves the nearest neighbor for each node from a
[PathwaySpace](https://github.com/sysbiolab/PathwaySpace/reference/PathwaySpace-class.md)
object using Euclidean distance.

## Usage

``` r
getNearestNode(ps)
```

## Arguments

- ps:

  Either a
  [PathwaySpace](https://github.com/sysbiolab/PathwaySpace/reference/PathwaySpace-class.md)
  or
  [GraphSpace](https://sysbiolab.github.io/RGraphSpace/reference/GraphSpace-methods.html)
  object.

## Value

A `data.frame` with columns `from`, `to`, and `dist`, listing each
node's nearest neighbor and the Euclidean distance between them.

## See also

[`nn2`](https://jefferislab.github.io/RANN/reference/nn2.html)

## Examples

``` r
# See examples in the PathwaySpace's tutorials:
# https://sysbiolab.github.io/PathwaySpace/
```

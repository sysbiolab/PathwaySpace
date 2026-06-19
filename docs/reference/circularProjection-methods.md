# Circular Projection of Graph-Associated Signals

`circularProjection()` implements a convolution algorithm to project
vertex-associated signals onto a 2D image space using a circular decay
function.

## Usage

``` r
# S4 method for class 'PathwaySpace'
circularProjection(
  ps,
  feature = activeFeature(ps),
  decay.fun = weibullDecay(),
  aggregate.fun = signalAggregation(),
  k = gs_vcount(ps),
  rescale = TRUE,
  verbose = TRUE,
  pdist = deprecated()
)
```

## Arguments

- ps:

  A
  [PathwaySpace](https://github.com/sysbiolab/PathwaySpace/reference/PathwaySpace-class.md)
  class object.

- feature:

  A single string specifying the feature to project as a signal. Must
  match either a feature name (see `gs_features(ps)`) or a node
  attribute (see `gs_names(ps)`). If a node attribute, make sure it is
  of numeric type. If the signal does not come from internal features,
  assign it directly using the
  [`vertexSignal`](https://github.com/sysbiolab/PathwaySpace/reference/vertexSignal-accessors.md)
  accessor.

- decay.fun:

  A signal decay function. Available options include 'Weibull',
  'exponential', and 'linear' (see
  [`weibullDecay`](https://github.com/sysbiolab/PathwaySpace/reference/weibullDecay.md)).
  Users may also define a custom decay model with at least two
  arguments, e.g., `function(x, signal) { ... }`, which should returns a
  vector of projected signals of the same length as `x`. Additional
  arguments may include any variable available as a graph vertex
  attribute.

- aggregate.fun:

  A function used to aggregate the projected signals. It must be
  provided as a unary function, e.g., `function(x) { ... }`, which
  should aggregate a vector of signals to a scalar value. Available
  options include 'mean', 'wmean', 'log.wmean', and 'exp.wmean' (See
  [`signalAggregation`](https://github.com/sysbiolab/PathwaySpace/reference/signalAggregation.md)).

- k:

  A single positive integer specifying the maximum number of vertices
  whose signals contribute to the projection. Defaults to
  `gs_vcount(ps)`, i.e. all vertices are considered. Specifically, at
  each point in space, the *k*-top decayed signals are retained prior to
  aggregation. Reducing *k* focuses the projection on the strongest
  local signals, filtering out weaker contributions.

- rescale:

  A logical value indicating whether to rescale the signal. If the
  signal `>=0`, then it will be rescaled to `[0, 1]`; if the signal
  `<=0`, then it will be rescaled to `[-1, 0]`; and if the signal in
  `(-Inf, +Inf)`, then it will be rescaled to `[-1, 1]`.

- verbose:

  A logical value specifying to display detailed messages (when
  `verbose=TRUE`) or not (when `verbose=FALSE`).

- pdist:

  Deprecated as of PathwaySpace 1.0.2; this parameter is now passed
  internally through `decay.fun`.

## Value

A preprocessed
[PathwaySpace](https://github.com/sysbiolab/PathwaySpace/reference/PathwaySpace-class.md)
class object.

## See also

[`buildPathwaySpace`](https://github.com/sysbiolab/PathwaySpace/reference/buildPathwaySpace.md)

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

# Create a 2D-landscape image
ps <- circularProjection(ps)
#> Validating arguments...
#> Using circular projection...
#> Mapping 'x' and 'y' coordinates...
#> Running signal convolution...
```

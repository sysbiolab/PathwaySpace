# Polar Projection of Graph-Associated Signals

`polarProjection()` implements a convolution algorithm to project
vertex-associated signals onto a 2D image space along graph edges, using
a polar decay function.

## Usage

``` r
# S4 method for class 'PathwaySpace'
polarProjection(
  ps,
  feature = activeFeature(ps),
  decay.fun = weibullDecay(pdist = 1),
  aggregate.fun = signalAggregation(),
  polar.fun = polarDecay(),
  k = gs_vcount(ps),
  beta = 10,
  directional = FALSE,
  edge.norm = TRUE,
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

- polar.fun:

  A polar decay function (see
  [`polarDecay`](https://github.com/sysbiolab/PathwaySpace/reference/polarDecay.md)).

- k:

  A single positive integer specifying the maximum number of vertices
  whose signals contribute to the projection. Defaults to
  `gs_vcount(ps)`, i.e. all vertices are considered. Specifically, at
  each point in space, the *k*-top decayed signals are retained prior to
  aggregation. Reducing *k* focuses the projection on the strongest
  local signals, filtering out weaker contributions.

- beta:

  An exponent (in `>=0)`) used in the polar projection functions (see
  [`polarDecay`](https://github.com/sysbiolab/PathwaySpace/reference/polarDecay.md)).
  It controls the shape of the polar projection by modulating the
  angular span. For example, \\beta = 0\\ yields a circular projection,
  \\beta = 1\\ produces a cardioid-like shape, and `beta > 1`
  progressively narrows the projection along a reference edge axis.

- directional:

  If directional edges are available, this argument can be used to
  orientate the signal projection on directed graphs.

- edge.norm:

  Scale distances based on edge lengths (when `edge.norm=TRUE`) or based
  on full coordinate space (when `edge.norm=FALSE`).

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
# Load a demo igraph
data('gtoy2', package = 'RGraphSpace')

# Create a new PathwaySpace object
ps <- buildPathwaySpace(gtoy2, nrc = 100)
#> Validating arguments...
#> Validating the 'igraph' object...
#> Normalizing node coordinates to graph space...
#> Creating a 'PathwaySpace' object...
# note: adjust 'nrc' to increase image resolution

# Set '1s' as vertex signal
vertexSignal(ps) <- 1

# Set edge weight
# gs_edge_attr(ps, "weight") <- c(-1, 1, 1, 1, 1, 1)

# Create a 2D-landscape image
ps <- polarProjection(ps)
#> Validating arguments...
#> Using polar projection on undirected graph...
#> Mapping 'x' and 'y' coordinates...
#> Computing linear and angular distances...
#> Running signal convolution...
```

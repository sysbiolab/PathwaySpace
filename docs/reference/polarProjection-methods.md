# Polar Projection of Graph-Associated Signals

`polarProjection` implements a convolution algorithm to project signals
across a 2D-coordinate system.

## Usage

``` r
# S4 method for class 'PathwaySpace'
polarProjection(
  ps,
  k = 2,
  beta = 10,
  decay.fun = weibullDecay(pdist = 1),
  aggregate.fun = signalAggregation(),
  polar.fun = polarDecay(),
  directional = FALSE,
  edge.norm = TRUE,
  rescale = TRUE,
  verbose = TRUE,
  theta = deprecated(),
  pdist = deprecated()
)
```

## Arguments

- ps:

  A
  [PathwaySpace](https://github.com/sysbiolab/PathwaySpace/reference/PathwaySpace-class.md)
  class object.

- k:

  A single positive integer determining the k-top signals for the
  convolution operation.

- beta:

  An exponent (in `[0, +Inf)`) used in the polar projection functions
  (see
  [`polarDecay`](https://github.com/sysbiolab/PathwaySpace/reference/polarDecay.md)).
  It controls the shape of the polar projection by modulating the
  angular span. For example, \\beta = 0\\ yields a circular projection,
  \\beta = 1\\ produces a cardioid-like shape, and `beta > 1`
  progressively narrows the projection along a reference edge axis.

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

- theta:

  Deprecated as of PathwaySpace 1.0.2; use 'beta' instead.

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
ps <- polarProjection(ps, pdist=1)
#> Validating arguments...
#> Using polar projection on undirected graph...
#> Mapping 'x' and 'y' coordinates...
#> Computing linear and angular distances...
#> Running signal convolution...
```

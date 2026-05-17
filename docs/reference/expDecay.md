# Constructor of exponential decay functions

The `expDecay()` constructor either creates a decay function or returns
a `ggplot` object for visualizing the decay model. It is a utility
function used internally by
[`circularProjection`](https://github.com/sysbiolab/PathwaySpace/reference/circularProjection-methods.md)
and
[`polarProjection`](https://github.com/sysbiolab/PathwaySpace/reference/polarProjection-methods.md).

## Usage

``` r
expDecay(decay = 0.001, pdist = 0.15, plot = FALSE, demo.signal = 1)
```

## Arguments

- decay:

  A decay factor (in `[0,1]`). This term indicates how much a `signal`
  decreases as a function of distance in pathway space. For example, at
  a specific distance defined by the `pdist` parameter, the signal
  intensity will be the initial signal multiplied by `decay`.

- pdist:

  A distance normalization term (in (0, 1\]) at which the signal reaches
  `signal * decay`. This parameter is used to anchor the decay to a
  meaningful distance (see `details`). Also, when `pdist = 1`, it will
  represent the diameter of the inscribed circle within the coordinate
  space of a `PathwaySpace` object.

- plot:

  A logical value indicating whether to return a `ggplot` object.

- demo.signal:

  A numeric value in `[-Inf, Inf]`, only passed when `plot = TRUE` to
  visualize the decay curve with a specific signal intensity. The value
  is ignored by the function constructor, as the decay function itself
  is returned without using an initial signal.

## Value

Returns either a function of the form `function(x, signal) { ... }` or,
if `plot = TRUE`, a `ggplot` object illustrating the decay model.

## Details

The `expDecay()` constructor creates an exponential decay model. It
describes how a signal decreases as a function of distance, controlled
by a decay rate parameter.

The decay function is defined as:

\$\$y = signal \times decay^{\left(\frac{x}{pdist}\right)}\$\$

where \\signal\\ represents the initial intensity, \\decay\\ controls
the rate of attenuation, and \\x\\ is a vector of normalized distances.
The \\pdist\\ parameter anchors the model such that:

- \\y = signal\\ when \\x = 0\\

- \\y = signal \times decay\\ when \\x = pdist\\

## See also

[`linearDecay`](https://github.com/sysbiolab/PathwaySpace/reference/linearDecay.md),
[`weibullDecay`](https://github.com/sysbiolab/PathwaySpace/reference/weibullDecay.md)

## Author

Sysbiolab Team

## Examples

``` r
# Return a decay function
decay_fun <- expDecay(decay = 0.25, pdist = 0.5)

# Plot decay model parameters
# expDecay(decay = 0.25, pdist = 0.5, plot = TRUE)
```

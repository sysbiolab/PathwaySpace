# Constructor of linear decay functions

The `linearDecay()` constructor either creates a decay function or
returns a `ggplot` object for visualizing the decay model. It is a
utility function used internally by
[`circularProjection`](https://github.com/sysbiolab/PathwaySpace/reference/circularProjection-methods.md)
and
[`polarProjection`](https://github.com/sysbiolab/PathwaySpace/reference/polarProjection-methods.md).

## Usage

``` r
linearDecay(decay = 0.001, pdist = 0.15, plot = FALSE, demo.signal = 1)
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

The `linearDecay()` constructor creates a simple linear decay model. It
describes how a signal decreases proportionally with distance.

The decay function is defined as: \$\$y = signal \times \left(1 - (1 -
decay) \times \frac{x}{pdist}\right)\$\$

where \\signal\\ represents the initial intensity, \\decay\\ defines the
relative signal level at \\pdist\\, and \\x\\ is a vector of normalized
distances. The signal decreases uniformly from its initial value to
\\pdist\\, which is a reference distance that anchors the model such
that:

- \\y = signal\\ when \\x = 0\\

- \\y = signal \times decay\\ when \\x = pdist\\

This makes the linear form consistent with the exponential and Weibull
decay functions, both of which also reach \\signal \times decay\\ at the
reference distance.

## See also

[`expDecay`](https://github.com/sysbiolab/PathwaySpace/reference/expDecay.md),
[`weibullDecay`](https://github.com/sysbiolab/PathwaySpace/reference/weibullDecay.md)

## Author

Sysbiolab Team

## Examples

``` r
# Return a decay function
decay_fun <- linearDecay(decay = 0.5, pdist = 0.25)

# Plot decay model parameters
# linearDecay(decay = 0.5, pdist = 0.25, plot = TRUE)
```

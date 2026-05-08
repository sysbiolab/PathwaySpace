# Polar transformation functions

Creates polar transformation functions for
[`polarProjection`](https://github.com/sysbiolab/PathwaySpace/reference/polarProjection-methods.md)
internal calls. These functions are used to adjusts signal decay
according to point-to-edge angular distances, with options to attenuate
angular shapes.

## Usage

``` r
polarDecay(
  method = c("power", "gaussian", "logistic"),
  s = 0.5,
  k = 10,
  m = 0.5
)
```

## Arguments

- method:

  String indicating the transformation to apply. Must be one of:
  "power", "gaussian", or "logistic".

- s:

  Single numeric value in `[0, 1]`. Controls the spread around the `x`
  mean of the Gaussian function.

- k:

  Single numeric value `>=1`. Controls the steepness of the logistic
  function.

- m:

  Single numeric value in `[0, 1]`. Specifies the midpoint of the
  logistic function.

## Value

Returns a function of the form: `function(x, beta) { ... }`, that
applies the specified shape-based transformation.

## Details

The polar transformation controls how much the projected signal decays
as a function of the angular distance between a point in pathway space
and a reference edge axis. The function returned by `polarDecay()`
expects two arguments, with the following signature:
`function(x, beta) { ... }`.

**Power:** \$\$x^{\beta}\$\$ where \\x\\ is a vector of normalized
angular distances (in `[0, 1]`) and \\beta\\ is a non-negative exponent
that controls the rate of signal decay. Increasing \\beta\\ results in a
steeper decay rate, modulating the angular span of the projection.

**Gaussian:**
\$\$\exp\left(-\frac{(1-x)^2}{2\sigma^2}\right)^{\beta}\$\$ where
\\sigma\\ controls the spread around the mean, creating fuzzier effect
on projections.

**Logistic:** \$\$(1 / (1 + \exp(k (x - m))))^{\beta}\$\$ where \\k\\ is
the steepness and \\m\\ is the function's midpoint, making more gradual
transitions.

These transformations are intended to be plugged into the higher-level
[`polarProjection`](https://github.com/sysbiolab/PathwaySpace/reference/polarProjection-methods.md)
function, allowing user control over the polar projection profiles.

## See also

[`polarProjection`](https://github.com/sysbiolab/PathwaySpace/reference/polarProjection-methods.md)

## Author

Sysbiolab Team

## Examples

``` r
polar.fun <- polarDecay("power")
```

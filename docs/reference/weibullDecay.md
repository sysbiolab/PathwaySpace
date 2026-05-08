# Constructor of Weibull decay functions

The \`weibullDecay()\` constructor either creates a decay function or
returns a \`ggplot\` object for visualizing the decay model. It is a
utility function used internally by
[`circularProjection`](https://github.com/sysbiolab/PathwaySpace/reference/circularProjection-methods.md)
and
[`polarProjection`](https://github.com/sysbiolab/PathwaySpace/reference/polarProjection-methods.md).

## Usage

``` r
weibullDecay(
  decay = 0.001,
  pdist = 0.15,
  shape = 1.05,
  plot = FALSE,
  demo.signal = 1
)
```

## Arguments

- decay:

  A decay factor (in \[0,1\]). This term indicates how much a `signal`
  decreases as a function of distance in pathway space. For example, at
  a specific distance defined by the `pdist` parameter, the signal
  intensity will be the initial signal multiplied by `decay`.

- pdist:

  A distance normalization term (in (0, 1\]) at which the signal reaches
  \`signal \* decay\`. This parameter is used to anchor the decay to a
  meaningful distance (see \`details\`). Also, when `pdist = 1`, it will
  represent the diameter of the inscribed circle within the coordinate
  space of a \`PathwaySpace\` object.

- shape:

  A parameter (\>=1) of a Weibull function. When `shape=1` the Weibull
  decay follows an exponential decay. When `shape>1` the function is
  first convex, then concave with an inflection point.

- plot:

  A logical value indicating whether to return a \`ggplot\` object.

- demo.signal:

  A numeric value in \`\[-Inf, Inf\]\`, only passed when `plot = TRUE`
  to visualize the decay curve with a specific signal intensity. The
  value is ignored by the function constructor, as the decay function
  itself is returned without using an initial signal.

## Value

Returns either a function of the form `function(x, signal) { ... }` or,
if `plot = TRUE`, a \`ggplot\` object illustrating the decay model.

## Details

The \`weibullDecay()\` constructor creates a decay model based on the
Weibull distribution. It describes how a signal decreases as a function
of distance, controlled by both a decay rate and a shape parameter.

The decay function is defined as:

\$\$y = signal \times decay^{\left(\frac{x}{pdist}\right)^{shape}}\$\$

where \\signal\\ represents the initial intensity, \\decay\\ controls
the rate of attenuation, \\x\\ is a vector of normalized distances, and
\\shape\\ adjusts the curvature of the decay. When \\shape = 1\\, the
function follows an exponential decay. For \\shape \> 1\\, the curve
transitions from convex to concave, exhibiting an inflection point. The
\\pdist\\ parameter anchors the model such that:

- \\y = signal\\ when \\x = 0\\

- \\y = signal \times decay\\ when \\x = pdist\\

## See also

[`linearDecay`](https://github.com/sysbiolab/PathwaySpace/reference/linearDecay.md),
[`expDecay`](https://github.com/sysbiolab/PathwaySpace/reference/expDecay.md)

## Author

Sysbiolab Team

## Examples

``` r
# Return a decay function
decay_fun <- weibullDecay(decay = 0.5, pdist = 0.4, shape = 2)

# Plot decay model parameters
# weibullDecay(decay = 0.5, pdist = 0.4, shape = 2, plot = TRUE)
```

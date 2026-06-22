# Signal aggregation functions

Signal aggregation functions for
[`circularProjection`](https://github.com/sysbiolab/PathwaySpace/reference/circularProjection-methods.md)
and
[`polarProjection`](https://github.com/sysbiolab/PathwaySpace/reference/polarProjection-methods.md)
internal calls. The aggregation should be symmetric with respect to
signal polarity, ensuring that opposite signals produce corresponding
outputs.

## Usage

``` r
signalAggregation(method = c("mean", "wmean", "log.wmean", "exp.wmean"))
```

## Arguments

- method:

  A character string specifying the method for signal aggregation,
  returning either a customized
  [`mean`](https://rdrr.io/r/base/mean.html) or
  [`weighted.mean`](https://rdrr.io/r/stats/weighted.mean.html)
  function.

## Value

Returns a function of the form: `function(x) { ... }`

## Details

At each point in pathway space, multiple vertices may each contribute a
decayed signal value; `method` controls how these contributions are
combined into a single value:

- **mean**: a plain, unweighted average. Because the same number of
  potential contributors is assumed throughout the image, points reached
  by few vertices can show a diluted value compared to points reached by
  many, even when the underlying signals are equally strong.

- **wmean**: each contribution is weighted by its own magnitude, so the
  strongest nearby signal dominates the result regardless of how many
  (or how weak) the other contributions are.

- **log.wmean**: like `wmean`, but the weighting is compressed, giving
  moderate signals comparatively more influence relative to the single
  strongest one.

- **exp.wmean**: like `wmean`, but the weighting is sharpened, so the
  strongest nearby signal dominates the result even more than under
  `wmean`.

Unlike `mean`, the weighted variants are not affected by the dilution
described above, since contributions with zero weight do not affect
their result. `mean` is a reasonable default for simple or binary
signals, such as those used in introductory examples; for continuous,
more nuanced analyses, `wmean` (or one of its variants) is generally
preferable. For other aggregation rules, see *fuzzy logic* functions in
the online tutorials: https://sysbiolab.github.io/PathwaySpace/

## See also

[`circularProjection`](https://github.com/sysbiolab/PathwaySpace/reference/circularProjection-methods.md),
[`polarProjection`](https://github.com/sysbiolab/PathwaySpace/reference/polarProjection-methods.md),
[`weighted.mean`](https://rdrr.io/r/stats/weighted.mean.html)

## Author

Sysbiolab Team

## Examples

``` r
aggregate.fun <- signalAggregation()
```

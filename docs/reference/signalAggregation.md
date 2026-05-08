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

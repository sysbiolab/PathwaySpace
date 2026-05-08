# Create interpolated color palettes for PathwaySpace images

Creates mixed color palettes by interpolating and offsetting hues,
useful for generating transitions between hues.

## Usage

``` r
pspace.pals(
  colors = c("#303f9d", "#578edb", "#63b946", "#f3930c", "#a60d0d"),
  trim.colors = c(3, 2, 1, 2, 3),
  offset = 0.5,
  n = 25
)
```

## Arguments

- colors:

  A vector of five base colors used to construct the custom diverging
  palette. These colors are interpolated according to the
  \`trim.colors\` values.

- trim.colors:

  A vector of five positive integers that control the relative weight of
  each hue in the five-color diverging palette.

- offset:

  Adjusts brightness by shifting hues toward the center, either brighter
  (\`offset \> 0\`) or darker (\` offset \< 0\`).

- n:

  The number of colors to generate in the output palette.

## Value

A vector with hexadecimal color codes.

## See also

[`plotPathwaySpace`](https://github.com/sysbiolab/PathwaySpace/reference/plotPathwaySpace-methods.md)

## Examples

``` r
pspace.pals()
#>  [1] "#253799" "#3141A1" "#3E4CAA" "#4A57B3" "#5661BC" "#626CC5" "#6D77CE"
#>  [8] "#7982D8" "#5681D1" "#568BD3" "#719CDD" "#8DADE8" "#4FA72B" "#E6AF89"
#> [15] "#E39E4D" "#DF8C10" "#DE7A25" "#E66464" "#DB5656" "#D04747" "#C63939"
#> [22] "#BB2A2A" "#B11919" "#A60404" "#990101"
```

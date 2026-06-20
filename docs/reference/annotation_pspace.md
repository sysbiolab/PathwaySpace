# Annotation Functions for PathwaySpace Plots

`annotation_pspace_signal()` annotates a `ggplot`-based `PathwaySpace`
plot with the projected signal layer, rendered as a raster heatmap
bounded to the normalized unit space `[0, 1]`.

## Usage

``` r
annotation_pspace_signal(
  ps,
  title = activeFeature(ps),
  colors = pspace.cols(),
  bg.color = "grey95",
  si.color = "grey85",
  si.alpha = 1,
  zlab = "Density",
  zlim = NULL,
  slices = 25,
  interpolate = FALSE,
  ...
)
```

## Arguments

- ps:

  A
  [PathwaySpace](https://github.com/sysbiolab/PathwaySpace/reference/PathwaySpace-class.md)
  object containing a valid signal projection. Run
  [`circularProjection`](https://github.com/sysbiolab/PathwaySpace/reference/circularProjection-methods.md)
  before calling this function.

- title:

  A string used as the plot title. Defaults to the active feature name
  returned by
  [`activeFeature`](https://github.com/sysbiolab/PathwaySpace/reference/vertexSignal-accessors.md).

- colors:

  A character vector of colors used to build the signal palette.
  Defaults to
  [`pspace.cols`](https://github.com/sysbiolab/PathwaySpace/reference/pspace.cols.md).

- bg.color:

  A string specifying the background color, used for zero-signal or
  masked regions. Defaults to `"grey95"`.

- si.color:

  A single color for silhouette. (see
  [`silhouetteMapping`](https://github.com/sysbiolab/PathwaySpace/reference/silhouetteMapping-methods.md)).

- si.alpha:

  A transparency level in `[0, 1]`, used to adjust the opacity of the
  silhouette. This parameter is useful for improving the perception of a
  background image, when one is available.

- zlab:

  The title for the 'z' axis of the image signal.

- zlim:

  The 'z' limits of the plot (a numeric vector with two numbers). If
  NULL, limits are determined from the range of the input values.

- slices:

  An integer specifying the number of discrete color levels used to
  quantize the signal. For `"negpos"` scale types, this is rounded up to
  the nearest even number.

- interpolate:

  A logical value indicating whether to apply linear interpolation with
  [`geom_raster`](https://ggplot2.tidyverse.org/reference/geom_tile.html).

- ...:

  Additional parameters passed to the title annotation; see
  [`geom_text`](https://ggplot2.tidyverse.org/reference/geom_text.html).

## Value

A list of `ggplot2` layer objects that can be added to a
[`ggplot()`](https://ggplot2.tidyverse.org/reference/ggplot.html) call
with `+`.

## See also

[`circularProjection`](https://github.com/sysbiolab/PathwaySpace/reference/circularProjection-methods.md)

## Examples

``` r
data("gtoy1", package = "RGraphSpace")
ps <- buildPathwaySpace(gtoy1, nrc = 100)
#> Validating arguments...
#> Validating the 'igraph' object...
#> Normalizing node coordinates to graph space...
#> Creating a 'PathwaySpace' object...
vertexSignal(ps) <- 1
ps <- circularProjection(ps)
#> Validating arguments...
#> Using circular projection...
#> Mapping 'x' and 'y' coordinates...
#> Running signal convolution...

if (FALSE) { # \dontrun{
ggplot(ps) +
  annotation_pspace_signal(ps, si.alpha = 0.8) +
  theme_gspace_coords(is_norm = TRUE)
} # }
```

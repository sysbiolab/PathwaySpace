# PathwaySpace: Spatial Projection of Network Signals along Geodesic Paths

For a given graph containing vertices, edges, and a signal associated
with the vertices, the PathwaySpace package performs a convolution
operation, which involves a weighted combination of neighboring vertices
and their associated signals. The package then uses a decay function to
propagate these signals, creating geodesic paths on a 2D-image space.

## Details

|             |                                         |
|-------------|-----------------------------------------|
| Package:    | PathwaySpace                            |
| Type:       | Software                                |
| License:    | Artistic-2.0                            |
| Maintainer: | Mauro Castro <mauro.a.castro@gmail.com> |

## Index

|  |  |
|----|----|
| [PathwaySpace-class](https://github.com/sysbiolab/PathwaySpace/reference/PathwaySpace-class.md): | An S4 class for signal propagation on pathway spaces. |
| [buildPathwaySpace](https://github.com/sysbiolab/PathwaySpace/reference/buildPathwaySpace.md): | Constructor of PathwaySpace-class objects. |
| [circularProjection](https://github.com/sysbiolab/PathwaySpace/reference/circularProjection-methods.md): | Creating 2D-landscape images from graph objects. |
| [polarProjection](https://github.com/sysbiolab/PathwaySpace/reference/polarProjection-methods.md): | Creating 2D-landscape images from graph objects. |
| [silhouetteMapping](https://github.com/sysbiolab/PathwaySpace/reference/silhouetteMapping-methods.md): | Mapping graph silhouettes on PathwaySpace images. |
| [summitMapping](https://github.com/sysbiolab/PathwaySpace/reference/summitMapping-methods.md): | Mapping summits on a 2D-landscape image. |
| [getPathwaySpace](https://github.com/sysbiolab/PathwaySpace/reference/getPathwaySpace-methods.md): | Accessory method for fetching slots from a PathwaySpace object. |
| [plotPathwaySpace](https://github.com/sysbiolab/PathwaySpace/reference/plotPathwaySpace-methods.md): | Plotting 2D-landscape images for the PathwaySpace package. |

Further information is available in the vignettes by typing
`vignette('PathwaySpace')`. Documented topics are also available in HTML
by typing [`help.start()`](https://rdrr.io/r/utils/help.start.html) and
selecting the PathwaySpace package from the menu.

## References

The Cancer Genome Atlas Analysis Network (2023). PathwaySpace: Spatial
propagation of network signals along geodesic paths. R package version
0.99.

## See also

Useful links:

- <https://sysbiolab.github.io/PathwaySpace/>

- <https://github.com/sysbiolab/PathwaySpace>

- Report bugs at <https://github.com/sysbiolab/PathwaySpace/issues>

## Author

The Cancer Genome Atlas Analysis Network

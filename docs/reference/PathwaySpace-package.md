# PathwaySpace: Spatial Projection of Network Signals along Geodesic Paths

For a given graph containing vertices, edges, and a signal associated
with the vertices, the 'PathwaySpace' package performs a convolution
operation, which involves a weighted combination of neighboring vertices
and their associated signals. The package uses a decay function to
project these signals, creating geodesic paths on a 2D-image space.
'PathwaySpace' has various applications, such as visualizing network
data in a graphical format that highlights the relationships and signal
strengths between vertices. By combining graph theory, signal
processing, and visualization, 'PathwaySpace' provides a way of
representing graph data on a continuous projection space. Based on
methods introduced in Tercan et al. (2025)
[doi:10.1016/j.xpro.2025.103681](https://doi.org/10.1016/j.xpro.2025.103681)
and Ellrott et al. (2025)
[doi:10.1016/j.ccell.2024.12.002](https://doi.org/10.1016/j.ccell.2024.12.002)
.

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
[`vignette('PathwaySpace')`](https://github.com/sysbiolab/PathwaySpace/articles/PathwaySpace.md).
Documented topics are also available in HTML by typing
[`help.start()`](https://rdrr.io/r/utils/help.start.html) and selecting
the PathwaySpace package from the menu.

## References

Sysbiolab Team (2026). *RGraphSpace: A lightweight interface between
'igraph' and 'ggplot2' graphics*. R package version 1.3.1 (Doi:
10.32614/CRAN.package.RGraphSpace),
<https://CRAN.R-project.org/package=RGraphSpace>.

## See also

Useful links:

- <https://sysbiolab.github.io/PathwaySpace/>

- <https://github.com/sysbiolab/PathwaySpace>

- Report bugs at <https://github.com/sysbiolab/PathwaySpace/issues>

## Author

**Maintainer**: Mauro Castro <mauro.a.castro@gmail.com>
([ORCID](https://orcid.org/0000-0003-4942-8131))

Authors:

- Sysbiolab Team

Other contributors:

- Victor Apolonio \[contributor\]

- Jonathan Back \[contributor\]

- Lana Querne \[contributor\]

- Vinicius Chagas \[contributor\]

- Bahar Tercan \[contributor\]

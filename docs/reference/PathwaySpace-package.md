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

## Details

For a hands-on introduction, see the vignette:
[`vignette("PathwaySpace")`](https://github.com/sysbiolab/PathwaySpace/articles/PathwaySpace.md).

The full set of documented topics can also be browsed in HTML by running
[`help.start()`](https://rdrr.io/r/utils/help.start.html) and selecting
the PathwaySpace package from the package list.

## References

Tercan B, Apolonio VH, Chagas VS, Wong CK, Lee JA, Yau C, Benz CC,
Stuart JM, Karlberg BJ, Ellrott K, Grewal JK, Jones SJ, Network TCGAA,
Zenklusen JC, Robertson AG, Laird PW, Cherniack AD, Castro MA (2025).
"Protocol for assessing distances in pathway space for classifier
feature sets from machine learning methods." *STAR Protocols*, *6*(2),
103681. doi:10.1016/j.xpro.2025.103681
<https://doi.org/10.1016/j.xpro.2025.103681>.

Ellrott K, Wong CK, Yau C, Castro MA, Lee JA, Karlberg BJ, Grewal JK,
Lagani V, Tercan B, al. e (2025). "Classification of non-TCGA cancer
samples to TCGA molecular subtypes using compact feature sets." *Cancer
Cell*, *43*(1), 1. doi:10.1016/j.ccell.2024.12.002
<https://doi.org/10.1016/j.ccell.2024.12.002>.

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

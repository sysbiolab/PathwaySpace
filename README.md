### *PathwaySpace*: Spatial projection of network signals along geodesic paths
  <!-- badges: start -->
  [![](https://www.r-pkg.org/badges/version/PathwaySpace)](https://cran.r-project.org/package=PathwaySpace)
  [![](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
  [![](https://cranlogs.r-pkg.org/badges/PathwaySpace)](https://cran.r-project.org/package=PathwaySpace)
  [![](https://img.shields.io/badge/license-Artistic--2.0-blue.svg)](https://cran.r-project.org/web/licenses/Artistic-2.0)
  [![](https://img.shields.io/badge/doi-10.32614/CRAN.package.PathwaySpace-blue.svg)](https://doi.org/10.32614/CRAN.package.PathwaySpace)
  <!-- badges: end -->
*PathwaySpace* is an R package that creates landscape images from graphs containing vertices (nodes), edges (lines), and a signal associated with the vertices. The package processes the signal using a convolution algorithm that considers the graph's topology to project the signal on a 2D space. 

*PathwaySpace* could have various applications, such as visualizing network data in a graphical format that highlights the relationships and signal strengths between vertices. 

### Installation in R (>=4.4)

##### Install dependencies to build the package's vignettes

```r
install.packages("knitr")
install.packages("rmarkdown")
```

##### Install the PathwaySpace package

```r
install.packages("remotes")
remotes::install_github("sysbiolab/RGraphSpace", build_vignettes=TRUE)
remotes::install_github("sysbiolab/PathwaySpace", build_vignettes=TRUE)
```

### Examples

Follow the *PathwaySpace* vignette and try to make some *brain plots*!

```r
library(PathwaySpace)
vignette("PathwaySpace")
```

### Licenses

The *PathwaySpace* package is distributed under [Artistic-2.0](https://www.r-project.org/Licenses/Artistic-2.0)

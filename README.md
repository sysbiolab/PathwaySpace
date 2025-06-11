### *PathwaySpace*: Spatial projection of network signals along geodesic paths
  <!-- badges: start -->
  [![](https://www.r-pkg.org/badges/version/PathwaySpace)](https://cran.r-project.org/package=PathwaySpace)
  [![](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
  [![](https://cranlogs.r-pkg.org/badges/PathwaySpace)](https://cranlogs.r-pkg.org/badges/PathwaySpace)
  [![](https://img.shields.io/badge/license-Artistic--2.0-blue.svg)](https://cran.r-project.org/web/licenses/Artistic-2.0)
  [![](https://img.shields.io/badge/doi-10.32614/CRAN.package.PathwaySpace-blue.svg)](https://doi.org/10.32614/CRAN.package.PathwaySpace)
  <!-- badges: end -->
*PathwaySpace* is an R package that creates landscape images from graphs 
containing vertices (nodes), edges (lines), and a signal associated with 
the vertices. The package processes the signal using a convolution algorithm 
that considers the graph's topology to project the signal on a 2D space.
**Figure 1** illustrates the convolution operation problem addressed by 
the *PathwaySpace* package. For detailed documentation and usage examples, 
see the package's vignettes and workflows.

*PathwaySpace* could have various applications, such as visualizing network 
data in a graphical format that highlights the relationships and signal 
strengths between vertices.

![Alt text](vignettes/figures/fig1.png?raw=true)

**Figure 1.** Signal processing addressed by the *PathwaySpace* package. 
**A**) Graph overlaid on a 2D coordinate system. Each projection cone represents 
the signal associated with a graph vertex (referred to as *vertex-signal positions*), 
while question marks indicate positions with no signal information (referred to 
as *null-signal positions*). *Inset*: Graph layout of the toy example used in 
the *quick start* section of the package's vignette. **B**) Illustration of signal 
projection from two neighboring vertices, simplified to one dimension. 
*Right*: Signal profiles from aggregation and decay functions.

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

### Citation

If you use *PathwaySpace*, please cite:

* Tercan & Apolonio *et al.* Protocol for assessing distances in pathway space for classifier feature sets from machine learning methods. *STAR Protocols*, 2025. https://doi.org/10.1016/j.xpro.2025.103681

* Ellrott *et al.* Classification of non-TCGA cancer samples to TCGA molecular subtypes using compact feature sets. *Cancer Cell*, 2025. https://doi.org/10.1016/j.ccell.2024.12.002

#### Supporting Material for Tercan *et al.* (2025)

Download and uncompress *Tercan_et_al_20250112.zip*, then follow the instructions in the *pspace_perturbation.R* script. This R script has been developed to reproduce the results presented in *Figure S1* of Tercan *et al.* (2025).

* [Tercan_et_al_20250112.zip](https://github.com/sysbiolab/PathwaySpace/blob/main/Tercan_et_al_20250112.zip)

### Licenses

The *PathwaySpace* package is distributed under [Artistic-2.0](https://www.r-project.org/Licenses/Artistic-2.0)

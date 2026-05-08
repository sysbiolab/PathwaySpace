# Installation Instructions

To install *PathwaySpace*, R version [4.5](https://www.r-project.org/)
or greater is required.

## Release version

``` r

# Release version from CRAN
install.packages("PathwaySpace")
```

## Development version

``` r

# Dependencies to build the vignettes
install.packages("knitr")
install.packages("rmarkdown")
install.packages("remotes")

# Package source
remotes::install_github("sysbiolab/RGraphSpace", build_vignettes=TRUE)
remotes::install_github("sysbiolab/PathwaySpace", build_vignettes=TRUE)

# Dependencies used in the spatial transcriptomics tutorial
remotes::install_github("satijalab/Seurat")
remotes::install_github("satijalab/seurat-data")
```

## Citation

- Tercan & Apolonio *et al.* Protocol for assessing distances in pathway
  space for classifier feature sets from machine learning methods. *STAR
  Protocols*, 2025. <https://doi.org/10.1016/j.xpro.2025.103681>

- Ellrott *et al.* Classification of non-TCGA cancer samples to TCGA
  molecular subtypes using compact feature sets. *Cancer Cell*, 2025.
  <https://doi.org/10.1016/j.ccell.2024.12.002>

## Licenses

The *PathwaySpace* package is distributed under
[Artistic-2.0](https://www.r-project.org/Licenses/Artistic-2.0)

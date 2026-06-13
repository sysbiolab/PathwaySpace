# Modeling signal decay functions

**Package**: PathwaySpace 1.3.1  

## Overview

In this tutorial, we demonstrate how to model signal decay functions
targeted at specific nodes. Using a simple lattice graph, we describe
the steps to create and apply different decay functions for
*PathwaySpace* projections. Nodes in this example represent spots,
providing an intuitive starting point for understanding how these
functions can be used to capture behaviors in larger or more complex
graphs, such as from [**spatial
transcriptomics**](https://github.com/sysbiolab/PathwaySpace/articles/spatial-transcriptomics.md)
data.

## Required packages

``` r

# Check required packages for this vignette
if (!require("remotes", quietly = TRUE)){
  install.packages("remotes")
}
if (!require("RGraphSpace", quietly = TRUE)){
  remotes::install_github("sysbiolab/RGraphSpace")
}
if (!require("PathwaySpace", quietly = TRUE)){
  remotes::install_github("sysbiolab/PathwaySpace")
}
```

``` r

# Check versions
if (packageVersion("RGraphSpace") < "1.3.1"){
  message("Need to update 'RGraphSpace' for this vignette")
  remotes::install_github("sysbiolab/RGraphSpace")
}
if (packageVersion("PathwaySpace") < "1.3.1"){
  message("Need to update 'PathwaySpace' for this vignette")
  remotes::install_github("sysbiolab/PathwaySpace")
}
```

``` r

# Load packages
library("igraph")
library("ggplot2")
library("RGraphSpace")
library("PathwaySpace")
library("patchwork")
```

## Setting basic input data

``` r

# Create a lattice graph with igraph
g <- make_lattice(c(9, 9), directed = FALSE)
V(g)$name <- paste0("n",1:vcount(g))
gs <- GraphSpace(g, layout = layout_on_grid(g))
gs <- normalizeGraphSpace(gs)
plotGraphSpace(gs, add.labels = TRUE)
```

![](modeling-signal-decay_files/figure-html/Setting%20basic%20input%20data%20-%201-1.png)

``` r

# Build a PathwaySpace object
ps <- buildPathwaySpace(gs)
```

## Creating a decay function

In the context of *PathwaySpace*, a decay function describes how a
signal decreases as a function of distance, modeling the gradual loss of
intensity. These functions are used to attenuate signals over
graph-based domains, so contributions from distant vertices are weighted
less than those nearby. *PathwaySpace* provides three built-in decay
function constructors,
[`weibullDecay()`](https://github.com/sysbiolab/PathwaySpace/reference/weibullDecay.md),
[`expDecay()`](https://github.com/sysbiolab/PathwaySpace/reference/expDecay.md),
and
[`linearDecay()`](https://github.com/sysbiolab/PathwaySpace/reference/linearDecay.md),
all sharing a consistent argument structure for practical use and
comparison:

``` math
y = signal \times \text{decay}^{\left(\frac{x}{\text{pdist}}\right)^{\text{shape}}} \quad \text{(Weibull)}
```

``` math
y = signal \times \text{decay}^{\left(\frac{x}{\text{pdist}}\right)} \quad \text{(Exponential)}
```

``` math
y = signal \times \left(1 - (1 - \text{decay}) \times \frac{x}{\text{pdist}}\right) \quad \text{(Linear*)}
```

``` math
\text{*output clipped to prevent flipping sign}
```

where $`\textbf{signal}`$ represents the initial intensity,
$`\textbf{decay}`$ controls the rate of attenuation, $`\textbf{x}`$ is a
vector of normalized distances, $`\textbf{shape}`$ adjusts the curvature
of the decay, $`\textbf{pdist}`$ is a normalization term, and
$`\textbf{y}`$ is the resulting signal decay values. Next, the main
arguments are illustrated using the
[`weibullDecay()`](https://github.com/sysbiolab/PathwaySpace/reference/weibullDecay.md)
constructor, which returns a decay function with customized parameters.

``` r

weibullDecay(decay = 0.25, shape = 2, pdist = 0.75)
#> function (x, signal) 
#> {
#>     y <- signal * 0.25^((x/0.75)^2)
#>     return(y)
#> }
#> <environment: 0x5fef69c4b638>
#> attr(,"name")
#> [1] "weibullDecay"
```

… and to visualize how different parameters affect the signal
attenuation, rerun the
[`weibullDecay()`](https://github.com/sysbiolab/PathwaySpace/reference/weibullDecay.md)
constructor with `plot = TRUE`.

``` r

# Run Weibull constructor with decay = 0.5, shape = 1, and pdist = 0.25
p1 <- weibullDecay(decay = 0.4, shape = 1, pdist = 0.2, plot = TRUE)

# Run Weibull constructor with decay = 0.5, shape = 2, and pdist = 0.50
p2 <- weibullDecay(decay = 0.4, shape = 2, pdist = 0.4, plot = TRUE)

# Run Weibull constructor with decay = 0.25, shape = 3, and pdist = 0.75
p3 <- weibullDecay(decay = 0.2, shape = 4, pdist = 0.8, plot = TRUE)

p1 + p2 + p3
```

![](modeling-signal-decay_files/figure-html/Creating%20a%20decay%20function%20-%204.2-1.png)

Note that the normalization term $`\textbf{pdist}`$ anchors the decay to
a reference distance where the initial signal $`S_0`$ decreases to
$`S_0 * \textbf{decay}`$, controlling the extent over which the signal
is projected along the x-axis. The $`\textbf{shape}`$ parameter, on the
other hand, controls the curvature of the decay function: When
$`\textbf{shape} = 1`$, the function follows an exponential decay; and
for $`\textbf{shape} > 1`$, the curve transitions from convex to
concave, increasingly sigmoidal.

Next, we demonstrate the
[`expDecay()`](https://github.com/sysbiolab/PathwaySpace/reference/expDecay.md)
and
[`linearDecay()`](https://github.com/sysbiolab/PathwaySpace/reference/linearDecay.md)
constructors:

``` r

# Run the exp expDecay constructor with decay = 0.4 and pdist = 0.2
p1 <- expDecay(decay = 0.4, pdist = 0.2, plot = TRUE)

# Run the linearDecay constructor with decay = 0.5,and pdist = 0.25
p2 <- linearDecay(decay = 0.4, pdist = 0.3, plot = TRUE)

p1 + p2
```

![](modeling-signal-decay_files/figure-html/Creating%20a%20decay%20function%20-%204.3-1.png)

The exponential model follows an asymptotic decay, while the linear
model decreases proportionally with distance and is clipped at zero to
prevent flipping sign. Together, these decay models produce distinct
geometric profiles when projecting signals onto the 2D coordinate space,
ranging from cone-like projections in the linear model (see the
conceptual representation in [**Figure
1A**](https://github.com/sysbiolab/PathwaySpace/articles/overview.md))
to radially symmetric peaks in the exponential model, with the Weibull
model transitioning between them to form dome-like surfaces depending on
its $`\textbf{shape}`$ parameter.

## Assigning a decay model

Returning to our lattice example, next we assign a linear decay model to
the *PathwaySpace* object, setting $`\textbf{pdist}`$ as the average
center-to-center distance between vertices. For more details, refer to
the documentation of the
[`vertexDecay()`](https://github.com/sysbiolab/PathwaySpace/reference/vertexSignal-accessors.md)
accessor.

``` r

# Get distance to the nearest vertex
near_df <- getNearestNode(ps)
pdist <- mean(near_df$dist)
# 'pdist' set as the average center-to-center distance between vertices
pdist
#> [1] 0.1
```

``` r

# Setting a linear decay model for all vertices
vertexDecay(ps) <- linearDecay(pdist = pdist)
```

## Running *PathwaySpace*

Now we project a random binary signal using the
[`circularProjection()`](https://github.com/sysbiolab/PathwaySpace/reference/circularProjection-methods.md)
function.

``` r

# Add a random binary signal
set.seed(10)
vertexSignal(ps) <- sample(c(0,1), gs_vcount(ps), replace = TRUE)
```

``` r

# Running and plotting projections
ps <- circularProjection(ps, k = 1)
plotPathwaySpace(ps, marks = "n34")
```

![](modeling-signal-decay_files/figure-html/Running%20PathwaySpace%20-%202-1.png)

By assigning a new decay model to specific vertices, we can explore how
signals propagate across the network under different rules. To
illustrate, we assign a Weibull decay function to vertex `n34`.

``` r

# Changing decay model of vertex 'n34'
vertexDecay(ps)[["n34"]] <- weibullDecay(shape = 2, pdist = pdist*10)
```

``` r

# Running and plotting projections
ps <- circularProjection(ps, k = 1)
plotPathwaySpace(ps, marks = "n34")
```

![](modeling-signal-decay_files/figure-html/Running%20PathwaySpace%20-%204-1.png)

## Session information

    #> R version 4.6.0 (2026-04-24)
    #> Platform: x86_64-pc-linux-gnu
    #> Running under: Ubuntu 24.04.4 LTS
    #> 
    #> Matrix products: default
    #> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    #> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    #> 
    #> locale:
    #>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    #>  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    #>  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    #>  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    #>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    #> [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    #> 
    #> time zone: America/Sao_Paulo
    #> tzcode source: system (glibc)
    #> 
    #> attached base packages:
    #> [1] stats     graphics  grDevices utils     datasets  methods   base     
    #> 
    #> other attached packages:
    #> [1] patchwork_1.3.2    igraph_2.3.2       PathwaySpace_1.3.1 RGraphSpace_1.4.0 
    #> [5] ggplot2_4.0.3      remotes_2.5.0     
    #> 
    #> loaded via a namespace (and not attached):
    #>  [1] sass_0.4.10        generics_0.1.4     tidyr_1.3.2        lattice_0.22-9    
    #>  [5] digest_0.6.39      magrittr_2.0.5     evaluate_1.0.5     grid_4.6.0        
    #>  [9] RColorBrewer_1.1-3 fastmap_1.2.0      jsonlite_2.0.0     Matrix_1.7-5      
    #> [13] ggrepel_0.9.8      ggnewscale_0.5.2   ggrastr_1.0.2      purrr_1.2.2       
    #> [17] scales_1.4.0       textshaping_1.0.5  jquerylib_0.1.4    cli_3.6.6         
    #> [21] rlang_1.2.0        tidygraph_1.3.1    withr_3.0.2        RANN_2.6.2        
    #> [25] cachem_1.1.0       yaml_2.3.12        otel_0.2.0         ggbeeswarm_0.7.3  
    #> [29] tools_4.6.0        dplyr_1.2.1        colorspace_2.1-2   vctrs_0.7.3       
    #> [33] R6_2.6.1           lifecycle_1.0.5    fs_2.1.0           htmlwidgets_1.6.4 
    #> [37] vipor_0.4.7        ragg_1.5.2         pkgconfig_2.0.3    beeswarm_0.4.0    
    #> [41] desc_1.4.3         pkgdown_2.2.0      pillar_1.11.1      bslib_0.11.0      
    #> [45] gtable_0.3.6       Rcpp_1.1.1-1.1     glue_1.8.1         systemfonts_1.3.2 
    #> [49] xfun_0.58          tibble_3.3.1       tidyselect_1.2.1   rstudioapi_0.18.0 
    #> [53] knitr_1.51         farver_2.1.2       htmltools_0.5.9    rmarkdown_2.31    
    #> [57] compiler_4.6.0     S7_0.2.2

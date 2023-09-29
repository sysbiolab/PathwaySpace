### PathwaySpace: Spatial projection of network signals along geodesic paths

The *PathwaySpace* is an R package that creates landscape images from graphs containing vertices (nodes), edges (lines), and a signal associated with the vertices. The package processes the signal using a convolution algorithm that considers the graph's topology to project the signal on a 2D space. 

The *PathwaySpace* could have various applications, such as visualizing network data in a graphical format that highlights the relationships and signal strengths between vertices. 


#### Installation

```r
install.packages("remotes")
```

```r
library(remotes)
install_github("sysbiolab/PathwaySpace")
```

#### Example use

Try the examples available in the *PathwaySpace* vignettes.

```r
library(PathwaySpace)
vignette("PathwaySpace")
```

#### Licenses

The *PathwaySpace* package is distributed under [Artistic-2.0](https://www.r-project.org/Licenses/Artistic-2.0)

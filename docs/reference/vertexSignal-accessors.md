# Accessor Functions for PathwaySpace Objects

Get or set vertex signals, decay functions, and the active feature in a
[PathwaySpace](https://github.com/sysbiolab/PathwaySpace/reference/PathwaySpace-class.md)
object.

`vertexSignal()` gets or sets the numeric signal assigned to each
vertex, used as input for spatial projection.

`vertexDecay()` gets or sets the decay function assigned to each vertex,
controlling how the signal attenuates with distance.

`activeFeature()` gets or sets the active feature name, which
automatically extracts the corresponding signal from the `fdata` slot or
node attributes and assigns it to `vertexSignal()`.

## Usage

``` r
# S4 method for class 'PathwaySpace'
vertexSignal(x)

# S4 method for class 'PathwaySpace'
vertexSignal(x) <- value

# S4 method for class 'PathwaySpace'
vertexDecay(x)

# S4 method for class 'PathwaySpace'
vertexDecay(x) <- value

# S4 method for class 'PathwaySpace'
activeFeature(x)

# S4 method for class 'PathwaySpace'
activeFeature(x) <- value
```

## Arguments

- x:

  A
  [PathwaySpace](https://github.com/sysbiolab/PathwaySpace/reference/PathwaySpace-class.md)
  class object.

- value:

  The new value to assign:

  - For `vertexSignal()`: a numeric vector or scalar.

  - For `vertexDecay()`: a decay function or list of decay functions
    (see
    [`linearDecay`](https://github.com/sysbiolab/PathwaySpace/reference/linearDecay.md),
    [`weibullDecay`](https://github.com/sysbiolab/PathwaySpace/reference/weibullDecay.md)).

  - For `activeFeature()`: a single string matching a feature name (see
    [`gs_features`](https://sysbiolab.github.io/RGraphSpace/reference/GraphSpace-accessors.html))
    or a node attribute (see
    [`gs_names`](https://sysbiolab.github.io/RGraphSpace/reference/GraphSpace-accessors.html)).

## Value

The updated
[PathwaySpace](https://github.com/sysbiolab/PathwaySpace/reference/PathwaySpace-class.md)
object.

## Examples

``` r
library(PathwaySpace)

# Load a demo igraph
data('gtoy1', package = 'RGraphSpace')
ps <- buildPathwaySpace(gtoy1, nrc = 100)
#> Validating arguments...
#> Validating the 'igraph' object...
#> Normalizing node coordinates to graph space...
#> Creating a 'PathwaySpace' object...

# Check vertex names
names(ps)
#> [1] "n1" "n2" "n3" "n4" "n5"

##--------------------------------------
## 'vertexSignal' accessor

# Access signal values from all vertices
vertexSignal(ps)
#> n1 n2 n3 n4 n5 
#>  0  0  0  0  0 

# Modify signal value of a specific vertex
vertexSignal(ps)[1] <- 1

# Modify signal value of specific vertices
vertexSignal(ps)[c("n2","n3")] <- 1

# Set '1s' to all vertices
vertexSignal(ps) <- 1

##--------------------------------------
## 'activeFeature' accessor

# Assign a signal feature matrix
signal_mtx <- matrix(
  rep(rnorm(gs_vcount(ps)), 2),
  ncol = 2,
  dimnames = list(names(ps), c("feature1", "feature2"))
)
gs_fdata(ps) <- signal_mtx

# Set the active feature — automatically updates vertexSignal()
activeFeature(ps) <- "feature1"
#> Setting active feature 'feature1' from feature matrix...

##--------------------------------------
## 'vertexDecay' accessor

# Access decay function of a specific vertex
vertexDecay(ps)[["n3"]]
#> function (x, signal) 
#> {
#>     y <- signal * 0.001^((x/0.15)^1.05)
#>     return(y)
#> }
#> <environment: 0x609f31ddcc10>
#> attr(,"name")
#> [1] "weibullDecay"

# Modify decay function of a specific vertex
vertexDecay(ps)[["n3"]] <- linearDecay()

# Modify decay functions of two vertices
vertexDecay(ps)[c("n1","n3")] <- list( weibullDecay() )

# Modify decay functions of all vertices
vertexDecay(ps) <- weibullDecay(shape = 2)
```

# Accessor Functions for PathwaySpace Objects

Get or set 'signal' and 'decay' functions in a
[PathwaySpace](https://github.com/sysbiolab/PathwaySpace/reference/PathwaySpace-class.md)
class object.

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
```

## Arguments

- x:

  A
  [PathwaySpace](https://github.com/sysbiolab/PathwaySpace/reference/PathwaySpace-class.md)
  class object.

- value:

  The new value of the attribute.

## Value

Updated
[PathwaySpace](https://github.com/sysbiolab/PathwaySpace/reference/PathwaySpace-class.md)
object.

## Examples

``` r
data('gtoy1', package = 'RGraphSpace')
ps <- buildPathwaySpace(gtoy1, nrc = 100)
#> Validating arguments...
#> Validating the 'igraph' object...
#> Normalizing node coordinates to graph space...
#> Creating a 'PathwaySpace' object...

# Check vertex names
names(ps)
#> [1] "n1" "n2" "n3" "n4" "n5"

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

#----

# Access decay function of a specific vertex
vertexDecay(ps)[["n3"]]
#> function (x, signal) 
#> {
#>     y <- signal * 0.001^((x/0.15)^1.05)
#>     return(y)
#> }
#> <environment: 0x566836c26760>
#> attr(,"name")
#> [1] "weibullDecay"

# Modify decay function of a specific vertex
vertexDecay(ps)[["n3"]] <- linearDecay()

# Modify decay functions of two vertices
vertexDecay(ps)[c("n1","n3")] <- list( weibullDecay() )

# Modify decay functions of all vertices
vertexDecay(ps) <- weibullDecay(shape = 2)
```

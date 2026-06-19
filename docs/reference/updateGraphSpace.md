# Update a PathwaySpace object

Updates outdated `PathwaySpace` objects serialized from previous package
versions, adding any missing slots with default values.

## Usage

``` r
# S4 method for class 'PathwaySpace'
updateGraphSpace(x, verbose = FALSE)
```

## Arguments

- x:

  A `PathwaySpace` object.

- verbose:

  Logical; if `TRUE`, reports which slots were added.

## Value

An updated `PathwaySpace` object.

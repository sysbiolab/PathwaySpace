# PathwaySpace: An S4 class for signal propagation on image spaces

PathwaySpace: An S4 class for signal propagation on image spaces

## Value

An S4 class object.

## Slots

- `nodes`:

  A data frame with xy-vertex coordinates.

- `edges`:

  A data frame with edges.

- `graph`:

  An igraph object.

- `image`:

  A raster background image matrix.

- `pars`:

  A list inherited GraphSpace parameters.

- `misc`:

  A list with intermediate objects for downstream methods.

- `projections`:

  A list with processed objects for downstream methods.

- `pars_ps`:

  A list with PathwaySpace parameters.

- `status`:

  A vector containing the processing status of the PathwaySpace object.

## Constructor

see
[`buildPathwaySpace`](https://github.com/sysbiolab/PathwaySpace/reference/buildPathwaySpace.md)
constructor.

## Author

Sysbiolab Team, Mauro Castro (<mauro.castro@ufpr.br>)

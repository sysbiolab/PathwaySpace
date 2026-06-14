# PathwaySpace: An S4 class for signal projection on image spaces

`PathwaySpace` extends the
[GraphSpace](https://sysbiolab.github.io/RGraphSpace/reference/GraphSpace-class.html)
class with signal projection slots. It stores projected signal matrices,
projection parameters, and workflow status, and is the main object used
by the PathwaySpace package.

## Value

A `PathwaySpace` object.

## Slots

- `nodes`:

  Inherited from
  [GraphSpace](https://sysbiolab.github.io/RGraphSpace/reference/GraphSpace-class.html).

- `edges`:

  Inherited from
  [GraphSpace](https://sysbiolab.github.io/RGraphSpace/reference/GraphSpace-class.html).

- `graph`:

  Inherited from
  [GraphSpace](https://sysbiolab.github.io/RGraphSpace/reference/GraphSpace-class.html).

- `image`:

  Inherited from
  [GraphSpace](https://sysbiolab.github.io/RGraphSpace/reference/GraphSpace-class.html).

- `fdata`:

  Inherited from
  [GraphSpace](https://sysbiolab.github.io/RGraphSpace/reference/GraphSpace-class.html).

- `pars`:

  Inherited from
  [GraphSpace](https://sysbiolab.github.io/RGraphSpace/reference/GraphSpace-class.html).

- `misc`:

  Inherited from
  [GraphSpace](https://sysbiolab.github.io/RGraphSpace/reference/GraphSpace-class.html).

- `uuid`:

  Inherited from
  [GraphSpace](https://sysbiolab.github.io/RGraphSpace/reference/GraphSpace-class.html).

- `projection`:

  A
  [SpaceProjection](https://github.com/sysbiolab/PathwaySpace/reference/SpaceProjection-class.md)
  object storing the intermediate and final matrices produced by a
  projection method.

- `pars_ps`:

  A list with PathwaySpace parameters.

- `status`:

  A vector containing the processing status of the PathwaySpace object.

## Constructor

See
[`buildPathwaySpace`](https://github.com/sysbiolab/PathwaySpace/reference/buildPathwaySpace.md).

## See also

[GraphSpace](https://sysbiolab.github.io/RGraphSpace/reference/GraphSpace-class.html),
[SpaceProjection](https://github.com/sysbiolab/PathwaySpace/reference/SpaceProjection-class.md),
[`buildPathwaySpace`](https://github.com/sysbiolab/PathwaySpace/reference/buildPathwaySpace.md),
[`circularProjection`](https://github.com/sysbiolab/PathwaySpace/reference/circularProjection-methods.md)

## Author

Sysbiolab Team, Mauro Castro (<mauro.castro@ufpr.br>)

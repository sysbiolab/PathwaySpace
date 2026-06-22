# Projection Data Container for PathwaySpace Objects

`SpaceProjection` is an S4 class that stores the intermediate and final
matrices produced during signal projection in a
[PathwaySpace](https://github.com/sysbiolab/PathwaySpace/reference/PathwaySpace-class.md)
object. It is created internally by
[`circularProjection`](https://github.com/sysbiolab/PathwaySpace/reference/circularProjection-methods.md)
or
[`polarProjection`](https://github.com/sysbiolab/PathwaySpace/reference/polarProjection-methods.md)
and is not intended to be constructed directly by the user.

## Value

A `SpaceProjection` object.

## Slots

- `coordinates`:

  A numeric matrix with one row per graph node and four columns (`X`,
  `Y`, `Xint`, `Yint`), storing continuous and integer grid coordinates
  for each node.

- `floor`:

  A numeric matrix of dimensions `nrc x nrc` representing the projection
  floor, used to mask regions outside the graph silhouette.

- `signal`:

  A numeric matrix of dimensions `nrc x nrc` storing the smoothed signal
  values before final scaling.

- `result`:

  A numeric matrix of dimensions `nrc x nrc` containing the final
  projected signal, scaled to `[0, 1]` (or `[-1, 1]` for `"negpos"`
  scale types).

## See also

[PathwaySpace](https://github.com/sysbiolab/PathwaySpace/reference/PathwaySpace-class.md)

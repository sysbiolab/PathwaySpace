# Package index

## PathwaySpace

- [`PathwaySpace-class`](https://github.com/sysbiolab/PathwaySpace/reference/PathwaySpace-class.md)
  : PathwaySpace: An S4 class for signal projection on image spaces
- [`buildPathwaySpace()`](https://github.com/sysbiolab/PathwaySpace/reference/buildPathwaySpace.md)
  : Constructor of PathwaySpace-class Objects
- [`getPathwaySpace(`*`<PathwaySpace>`*`)`](https://github.com/sysbiolab/PathwaySpace/reference/getPathwaySpace-methods.md)
  : Accessors for Fetching Slots from a PathwaySpace Object
- [`updatePathwaySpace(`*`<PathwaySpace>`*`)`](https://github.com/sysbiolab/PathwaySpace/reference/updatePathwaySpace.md)
  : Update a PathwaySpace object

## Accessors

- [`` `gs_vertex_attr<-`( ``*`<PathwaySpace>`*`)`](https://github.com/sysbiolab/PathwaySpace/reference/PathwaySpace-accessors.md)
  [`` `gs_edge_attr<-`( ``*`<PathwaySpace>`*`)`](https://github.com/sysbiolab/PathwaySpace/reference/PathwaySpace-accessors.md)
  : Accessor Functions for PathwaySpace Objects
- [`vertexSignal(`*`<PathwaySpace>`*`)`](https://github.com/sysbiolab/PathwaySpace/reference/vertexSignal-accessors.md)
  [`` `vertexSignal<-`( ``*`<PathwaySpace>`*`)`](https://github.com/sysbiolab/PathwaySpace/reference/vertexSignal-accessors.md)
  [`vertexDecay(`*`<PathwaySpace>`*`)`](https://github.com/sysbiolab/PathwaySpace/reference/vertexSignal-accessors.md)
  [`` `vertexDecay<-`( ``*`<PathwaySpace>`*`)`](https://github.com/sysbiolab/PathwaySpace/reference/vertexSignal-accessors.md)
  [`activeFeature(`*`<PathwaySpace>`*`)`](https://github.com/sysbiolab/PathwaySpace/reference/vertexSignal-accessors.md)
  [`` `activeFeature<-`( ``*`<PathwaySpace>`*`)`](https://github.com/sysbiolab/PathwaySpace/reference/vertexSignal-accessors.md)
  : Accessor Functions for PathwaySpace Objects

## Projections

Spatial projection methods and the SpaceProjection class

- [`SpaceProjection-class`](https://github.com/sysbiolab/PathwaySpace/reference/SpaceProjection-class.md)
  : Projection Data Container for PathwaySpace Objects
- [`circularProjection(`*`<PathwaySpace>`*`)`](https://github.com/sysbiolab/PathwaySpace/reference/circularProjection-methods.md)
  : Circular Projection of Graph-Associated Signals
- [`polarProjection(`*`<PathwaySpace>`*`)`](https://github.com/sysbiolab/PathwaySpace/reference/polarProjection-methods.md)
  : Polar Projection of Graph-Associated Signals

## Decorations

- [`silhouetteMapping(`*`<PathwaySpace>`*`)`](https://github.com/sysbiolab/PathwaySpace/reference/silhouetteMapping-methods.md)
  : Decorating PathwaySpace Images with Graph Silhouettes
- [`summitMapping(`*`<PathwaySpace>`*`)`](https://github.com/sysbiolab/PathwaySpace/reference/summitMapping-methods.md)
  : Mapping Summits on PathwaySpace Images
- [`summitWatershed()`](https://github.com/sysbiolab/PathwaySpace/reference/summitWatershed.md)
  : Variation of the watershed algorithm for summit detection

## Decay functions

Signal decay models used in spatial projection

- [`linearDecay()`](https://github.com/sysbiolab/PathwaySpace/reference/linearDecay.md)
  : Constructor of linear decay functions
- [`expDecay()`](https://github.com/sysbiolab/PathwaySpace/reference/expDecay.md)
  : Constructor of exponential decay functions
- [`polarDecay()`](https://github.com/sysbiolab/PathwaySpace/reference/polarDecay.md)
  : Polar transformation functions
- [`weibullDecay()`](https://github.com/sysbiolab/PathwaySpace/reference/weibullDecay.md)
  : Constructor of Weibull decay functions

## Aggregation functions

- [`signalAggregation()`](https://github.com/sysbiolab/PathwaySpace/reference/signalAggregation.md)
  : Signal aggregation functions

## Path distances

- [`pathDistances()`](https://github.com/sysbiolab/PathwaySpace/reference/pathDistances.md)
  : Calculate a pathway space distance between two vectors
- [`plotPathDistances()`](https://github.com/sysbiolab/PathwaySpace/reference/plotPathDistances.md)
  : Accessory function to plot pathway space distances
- [`getNearestNode()`](https://github.com/sysbiolab/PathwaySpace/reference/getNearestNode.md)
  : getNearestNode

## Plotting

- [`plotPathwaySpace(`*`<PathwaySpace>`*`)`](https://github.com/sysbiolab/PathwaySpace/reference/plotPathwaySpace-methods.md)
  : Plotting 2D-landscape images for the PathwaySpace package
- [`annotation_pspace_signal()`](https://github.com/sysbiolab/PathwaySpace/reference/annotation_pspace_signal.md)
  : Annotation Functions for PathwaySpace Plots

## Utilities

- [`pspace.cols()`](https://github.com/sysbiolab/PathwaySpace/reference/pspace.cols.md)
  : A simple vector of colors for PathwaySpace images
- [`pspace.pals()`](https://github.com/sysbiolab/PathwaySpace/reference/pspace.pals.md)
  : Create interpolated color palettes for PathwaySpace images
- [`gimage`](https://github.com/sysbiolab/PathwaySpace/reference/gimage.md)
  : An image matrix

## Data

- [`CGC_20211118`](https://github.com/sysbiolab/PathwaySpace/reference/CGC_20211118.md)
  : COSMIC-CGC genes mapped to PathwaySpace images
- [`Hallmarks_v2023_1_Hs_symbols`](https://github.com/sysbiolab/PathwaySpace/reference/Hallmarks_v2023_1_Hs_symbols.md)
  : A list with Hallmark gene sets (v2023.1)
- [`PCv12_pruned_igraph`](https://github.com/sysbiolab/PathwaySpace/reference/PCv12_pruned_igraph.md)
  : A pruned and laid out igraph object from Pathway Commons V12

## Package

- [`PathwaySpace`](https://github.com/sysbiolab/PathwaySpace/reference/PathwaySpace-package.md)
  [`PathwaySpace-package`](https://github.com/sysbiolab/PathwaySpace/reference/PathwaySpace-package.md)
  : PathwaySpace: Spatial Projection of Network Signals along Geodesic
  Paths

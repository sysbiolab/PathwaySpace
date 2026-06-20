
#' @title  Constructor of PathwaySpace-class Objects
#' 
#' @description \code{buildPathwaySpace} is a constructor of
#' PathwaySpace-class objects.
#' 
#' @param gs A \code{\link[RGraphSpace]{GraphSpace}} object. Alternatively, 
#' an \code{\link[igraph]{igraph}} object with node coordinates assigned 
#' to \code{x} and \code{y} vertex attributes, and node labels assigned 
#' to \code{name} vertex attribute.
#' @param nrc A single positive integer indicating the number of rows and 
#' columns (in pixels) for a square image matrix. This argument will 
#' affect the resulting image size and resolution.
#' @param verbose A logical value specifying to display detailed 
#' messages (when \code{verbose=TRUE}) or not (when \code{verbose=FALSE}).
#' @return A pre-processed \linkS4class{PathwaySpace} class object.
#' @author Sysbiolab Team
#' @seealso \code{\link{circularProjection}}, \code{\link{polarProjection}}
#' @examples
#' library(PathwaySpace)
#' 
#' # Load a demo igraph
#' data('gtoy1', package = 'RGraphSpace')
#' 
#' # Check graph validity
#' gs <- GraphSpace(gtoy1)
#' 
#' gs <- normalizeGraphSpace(gs)
#' 
#' # Create a new PathwaySpace object
#' ps <- buildPathwaySpace(gs, nrc = 100)
#' # note: adjust 'nrc' to increase image resolution
#' 
#' @importFrom igraph degree vcount ecount which_mutual is_igraph
#' @importFrom igraph as_edgelist as_adjacency_matrix
#' @importFrom igraph simplify V E 'V<-' 'E<-' is_directed
#' @importFrom stats quantile sd
#' @importFrom scales rescale
#' @importFrom RGraphSpace GraphSpace getGraphSpace normalizeGraphSpace
#' @importFrom RANN nn2
#' @importFrom rlang abort warn
#' @aliases buildPathwaySpace
#' @export
#' 
buildPathwaySpace <- function(gs, nrc = 500, verbose = TRUE) {
  .validate.ps.args("singleInteger", "nrc", nrc)
  .validate.ps.args("singleLogical", "verbose", verbose)
  if(verbose) rlang::inform("Validating arguments...")
  #--- validate argument types
  if (nrc < 2) {
    rlang::abort("'nrc' must be >=2")
  }
  #--- validate the graph object
  if(is_igraph(gs)){
    if(verbose) rlang::inform("Validating the 'igraph' object...")
    gs <- GraphSpace(gs, verbose=FALSE)
    gs <- normalizeGraphSpace(gs)
  }
  .validate.gspace(gs)
  if(!gs@pars$is.normalized){
    gs <- normalizeGraphSpace(gs)
  }
    
  #--- build PathwaySpace-class
  ps <- .buildPathwaySpace(gs, nrc, verbose)
  ps <- .updateStatus(ps, "Preprocess")
  return(ps)
}

#-------------------------------------------------------------------------------
#' @title Circular Projection of Graph-Associated Signals
#'
#' @description
#' \code{circularProjection()} implements a convolution algorithm to project
#' vertex-associated signals onto a 2D image space using a circular decay
#' function.
#' 
#' @param ps A \linkS4class{PathwaySpace} class object.
#' @param feature A single string specifying the feature to project as a
#' signal. Must match either a feature name (see \code{gs_features(ps)}) or 
#' a node attribute (see \code{gs_names(ps)}). If a node attribute, make sure
#' it is of numeric type. If the signal does not come from internal features, 
#' assign it directly using the \code{\link{vertexSignal}} accessor.
#' @param decay.fun A signal decay function. Available options include 
#' 'Weibull', 'exponential', and 'linear' (see \code{\link{weibullDecay}}).
#' Users may also define a custom decay model with at least two arguments, 
#' e.g., \code{function(x, signal) \{ ... \}}, which should returns a vector of  
#' projected signals of the same length as \code{x}. Additional arguments may  
#' include any variable available as a graph vertex attribute.
#' @param aggregate.fun A function used to aggregate the projected signals. 
#' It must be provided as a unary function, e.g., \code{function(x) { ... }}, 
#' which should aggregate a vector of signals to a scalar value. 
#' Available options include 'mean', 'wmean', 'log.wmean', and 'exp.wmean' 
#' (See \code{\link{signalAggregation}}).
#' @param k A single positive integer specifying the maximum number of 
#' vertices whose signals contribute to the projection. Defaults to 
#' \code{gs_vcount(ps)}, i.e. all vertices are considered. Specifically, 
#' at each point in space, the \emph{k}-top decayed signals are retained 
#' prior to aggregation. Reducing \emph{k} focuses the projection on the 
#' strongest local signals, filtering out weaker contributions.
#' @param rescale A logical value indicating whether to rescale 
#' the signal. If the signal \code{>=0}, then it will be rescaled to 
#' \code{[0, 1]}; if the signal \code{<=0}, then it will be rescaled to 
#' \code{[-1, 0]}; and if the signal in \code{(-Inf, +Inf)}, then it will be 
#' rescaled to \code{[-1, 1]}.
#' @param verbose A logical value specifying to display detailed 
#' messages (when \code{verbose=TRUE}) or not (when \code{verbose=FALSE}).
#' @param pdist Deprecated as of PathwaySpace 1.0.2; this parameter is now 
#' passed internally through \code{decay.fun}.
#' @return A preprocessed \linkS4class{PathwaySpace} class object.
#' @author Sysbiolab Team
#' @seealso \code{\link{buildPathwaySpace}}
#' @examples
#' library(PathwaySpace)
#' 
#' # Load a demo igraph
#' data('gtoy1', package = 'RGraphSpace')
#'
#' # Create a new PathwaySpace object
#' ps <- buildPathwaySpace(gtoy1, nrc = 100)
#' # note: adjust 'nrc' to increase image resolution
#' 
#' # Set '1s' as vertex signal
#' vertexSignal(ps) <- 1
#' 
#' # Create a 2D-landscape image
#' ps <- circularProjection(ps)
#' 
#' @import methods
#' @importFrom lifecycle deprecated deprecate_soft is_present deprecate_stop
#' @importFrom RGraphSpace gs_vcount
#' @docType methods
#' @rdname circularProjection-methods
#' @aliases circularProjection
#' @export
#'
setMethod("circularProjection", "PathwaySpace", function(ps, 
  feature = activeFeature(ps), 
  decay.fun = weibullDecay(), 
  aggregate.fun = signalAggregation(), 
  k = gs_vcount(ps), 
  rescale = TRUE, verbose = TRUE, 
  pdist = deprecated()) {
  
  ps <- updateGraphSpace(ps)
  
  ### deprecate
  if (lifecycle::is_present(pdist)) {
    deprecate_soft("1.0.2", "circularProjection(pdist)", 
      "circularProjection(decay.fun)")
  }
  #--- validate the pipeline status
  if (!.checkStatus(ps, "Preprocess")) {
    rlang::abort("The 'ps' object has not been preprocessed.")
  }
  #--- validate argument types
  .validate.ps.args("singleInteger", "k", k)
  .validate.ps.args("function", "aggregate.fun", aggregate.fun)
  .validate.ps.args("function", "decay.fun", decay.fun)
  .validate.ps.args("singleLogical", "rescale", rescale)
  .validate.ps.args("singleLogical", "verbose", verbose)
  if(verbose) rlang::inform("Validating arguments...")
  #--- validate argument values
  if (k < 1) {
    rlang::abort("'k' must be >=1")
  }
  n <- gs_vcount(ps)
  if (k > n) k <- n
  #--- validate functions
  if(!missing(decay.fun)){
    gs_vertex_attr(ps, "decayFunction") <- decay.fun
  }
  .validate_aggregate_fun(aggregate.fun)
  
  #--- add feature signal
  if(!is.null(feature)){
    .validate.ps.args("singleString", "feature", feature)
    if (!identical(feature, activeFeature(ps))) {
      activeFeature(ps) <- feature
    }
  }
  
  #--- pack args
  pars <- list(k = k, rescale = rescale, 
    aggregate.fun = aggregate.fun, 
    feature = feature, 
    projection = "Circular")
  for (nm in names(pars)) {
    ps@pars_ps[[nm]] <- pars[[nm]]
  }
  #--- run ps pipeline
  ps <- .circularProjection(ps, verbose)
  ps <- .updateStatus(ps, "CircularProjection")
  if (.checkStatus(ps, "PolarProjection")) {
    if(verbose) rlang::inform("-- polar projection replaced by circular.")
    ps <- .updateStatus(ps, "PolarProjection", FALSE)
  }
  return(ps)
})

#-------------------------------------------------------------------------------
#' @title Polar Projection of Graph-Associated Signals
#'
#' @description
#' \code{polarProjection()} implements a convolution algorithm to project
#' vertex-associated signals onto a 2D image space along graph edges, using
#' a polar decay function.
#'
#' @param ps A \linkS4class{PathwaySpace} class object.
#' @param feature A single string specifying the feature to project as a
#' signal. Must match either a feature name (see \code{gs_features(ps)}) or 
#' a node attribute (see \code{gs_names(ps)}). If a node attribute, make sure
#' it is of numeric type. If the signal does not come from internal features, 
#' assign it directly using the \code{\link{vertexSignal}} accessor.
#' @param decay.fun A signal decay function. Available options include 
#' 'Weibull', 'exponential', and 'linear' (see \code{\link{weibullDecay}}).
#' Users may also define a custom decay model with at least two arguments, 
#' e.g., \code{function(x, signal) \{ ... \}}, which should returns a vector of  
#' projected signals of the same length as \code{x}. Additional arguments may  
#' include any variable available as a graph vertex attribute.
#' @param aggregate.fun A function used to aggregate the projected signals. 
#' It must be provided as a unary function, e.g., \code{function(x) { ... }}, 
#' which should aggregate a vector of signals to a scalar value. 
#' Available options include 'mean', 'wmean', 'log.wmean', and 'exp.wmean' 
#' (See \code{\link{signalAggregation}}).
#' @param polar.fun A polar decay function (see \code{\link{polarDecay}}).
#' @param k A single positive integer specifying the maximum number of 
#' vertices whose signals contribute to the projection. Defaults to 
#' \code{gs_vcount(ps)}, i.e. all vertices are considered. Specifically, 
#' at each point in space, the \emph{k}-top decayed signals are retained 
#' prior to aggregation. Reducing \emph{k} focuses the projection on the 
#' strongest local signals, filtering out weaker contributions.
#' @param beta An exponent (in \code{>=0)}) used in the polar projection 
#' functions (see \code{\link{polarDecay}}). It controls the shape of the 
#' polar projection by modulating the angular span. For example, 
#' \eqn{beta = 0} yields a circular projection, \eqn{beta = 1} produces 
#' a cardioid-like shape, and \code{beta > 1} progressively narrows 
#' the projection along a reference edge axis.
#' @param directional If directional edges are available, this argument can 
#' be used to orientate the signal projection on directed graphs.
#' @param edge.norm Scale distances based on edge lengths 
#' (when \code{edge.norm=TRUE}) or based on full coordinate space 
#' (when \code{edge.norm=FALSE}).
#' @param rescale A logical value indicating whether to rescale 
#' the signal. If the signal \code{>=0}, then it will be rescaled to 
#' \code{[0, 1]}; if the signal \code{<=0}, then it will be rescaled to 
#' \code{[-1, 0]}; and if the signal in \code{(-Inf, +Inf)}, then it will be 
#' rescaled to \code{[-1, 1]}.
#' @param verbose A logical value specifying to display detailed 
#' messages (when \code{verbose=TRUE}) or not (when \code{verbose=FALSE}).
#' @param pdist Deprecated as of PathwaySpace 1.0.2; this parameter is now 
#' passed internally through \code{decay.fun}.
#' @return A preprocessed \linkS4class{PathwaySpace} class object.
#' @author Sysbiolab Team
#' @seealso \code{\link{buildPathwaySpace}}
#' @examples
#' library(PathwaySpace)
#' 
#' # Load a demo igraph
#' data('gtoy2', package = 'RGraphSpace')
#' 
#' # Create a new PathwaySpace object
#' ps <- buildPathwaySpace(gtoy2, nrc = 100)
#' # note: adjust 'nrc' to increase image resolution
#' 
#' # Set '1s' as vertex signal
#' vertexSignal(ps) <- 1
#' 
#' # Set edge weight
#' # gs_edge_attr(ps, "weight") <- c(-1, 1, 1, 1, 1, 1)
#' 
#' # Create a 2D-landscape image
#' ps <- polarProjection(ps)
#' 
#' @import methods
#' @docType methods
#' @rdname polarProjection-methods
#' @aliases polarProjection
#' @export
#'
setMethod("polarProjection", "PathwaySpace", function(ps, 
  feature = activeFeature(ps), 
  decay.fun = weibullDecay(pdist = 1),
  aggregate.fun = signalAggregation(), 
  polar.fun = polarDecay(), 
  k = gs_vcount(ps), 
  beta = 10,
  directional = FALSE,
  edge.norm = TRUE,
  rescale = TRUE, 
  verbose = TRUE, 
  pdist = deprecated()) {
  
  ps <- updateGraphSpace(ps)
  
  #--- validate the pipeline status
  if (!.checkStatus(ps, "Preprocess")) {
    rlang::abort("The 'ps' object has not been preprocessed.")
  }
  ### deprecated
  if (lifecycle::is_present(pdist)) {
    deprecate_soft("1.0.2", "polarProjection(pdist)", 
      "polarProjection(decay.fun)")
  }
  ###
  .validate.ps.args("singleInteger", "k", k)
  .validate.ps.args("singleNumber", "beta", beta)
  .validate.ps.args("function", "decay.fun", decay.fun)
  .validate.ps.args("function", "aggregate.fun", aggregate.fun)
  .validate.ps.args("function", "polar.fun", polar.fun)
  .validate.ps.args("singleLogical", "directional", directional)
  .validate.ps.args("singleLogical", "edge.norm", edge.norm)
  .validate.ps.args("singleLogical", "rescale", rescale)
  .validate.ps.args("singleLogical", "verbose", verbose)
  if(verbose) rlang::inform("Validating arguments...")
  if (k < 1) {
    rlang::abort("'k' must be >=1")
  }
  n <- gs_vcount(ps)
  if (k > n) k <- n
  if (beta < 0) {
    rlang::abort("'beta' must be an exponent in [0,+Inf)")
  }
  
  #--- validate functions
  if(!missing(decay.fun)){
    gs_vertex_attr(ps, "decayFunction") <- decay.fun
  }
  .validate_aggregate_fun(aggregate.fun)
  .validate_polar_fun(polar.fun)
  
  #--- add feature signal
  if(!is.null(feature)){
    .validate.ps.args("singleString", "feature", feature)
    if (!identical(feature, activeFeature(ps))) {
      activeFeature(ps) <- feature
    }
  }
  
  #--- pack args
  pars <- list(k = k, beta = beta,
    edge.norm = edge.norm, rescale = rescale, directional = directional,
    polar.fun = polar.fun, aggregate.fun = aggregate.fun, 
    projection = "Polar")
  for (nm in names(pars)) {
    ps@pars_ps[[nm]] <- pars[[nm]]
  }
  ps <- .polarProjection(ps, verbose)
  ps <- .updateStatus(ps, "PolarProjection")
  if (.checkStatus(ps, "CircularProjection")) {
    if(verbose) rlang::inform("-- circular projection replaced by polar.")
    ps <- .updateStatus(ps, "CircularProjection", FALSE)
  }
  return(ps)
})

#' @title Decorating PathwaySpace Images with Graph Silhouettes
#'
#' @description \code{silhouetteMapping} constructs an image baseline used
#' to outline the graph layout in a PathwaySpace image.
#'
#' @param ps A \linkS4class{PathwaySpace} class object.
#' @param pdist A term (in \code{(0,1]}) determining a distance unit for the
#' silhouette projection.
#' @param baseline A fraction (in \code{[0,1]}) of the silhouette projection,
#' representing the level over which a silhouette will outline the graph layout.
#' When \code{baseline = 0} (i.e. lower level of the projection), the 
#' silhouette will extend over the entire image space, so no outline will 
#' be visible.
#' @param fill.cavity A logical value specifying to fill cavities 
#' in the silhouette mask (when \code{fill.cavity=TRUE}) or not 
#' (when \code{fill.cavity=FALSE}).
#' @param verbose A logical value specifying to display detailed 
#' messages (when \code{verbose=TRUE}) or not (when \code{verbose=FALSE}).
#' @return A preprocessed \linkS4class{PathwaySpace} class object.
#' @author Sysbiolab Team
#' @seealso \code{\link{circularProjection}}
#' @examples
#' library(PathwaySpace)
#' 
#' # Load a demo igraph
#' data('gtoy1', package = 'RGraphSpace')
#'
#' # Create a new PathwaySpace object
#' ps <- buildPathwaySpace(gtoy1, nrc = 100)
#' # note: adjust 'nrc' to increase image resolution
#' 
#' # Set '1s' as vertex signal
#' vertexSignal(ps) <- 1
#'
#' # Map graph silhouette
#' ps <- silhouetteMapping(ps, pdist = 0.1)
#'
#' @import methods
#' @docType methods
#' @rdname silhouetteMapping-methods
#' @aliases silhouetteMapping
#' @export
#'
setMethod("silhouetteMapping", "PathwaySpace", function(ps,
  pdist = 0.05, baseline = 0.01, fill.cavity = TRUE, verbose = TRUE) {
  
  ps <- updateGraphSpace(ps)

  #--- validate the pipeline status
  if (!.checkStatus(ps, "Preprocess")) {
    rlang::abort("The 'ps' object has not been preprocessed.")
  }
  
  #--- validate argument types
  .validate.ps.args("singleNumber", "pdist", pdist)
  .validate.ps.args("singleNumber", "baseline", baseline)
  .validate.ps.args("singleLogical", "fill.cavity", fill.cavity)
  .validate.ps.args("singleLogical", "verbose", verbose)
  if(verbose) rlang::inform("Validating arguments...")
  
  #--- validate argument values
  if (baseline < 0 || baseline > 1) {
    rlang::abort("'baseline' must be in [0,1]")
  }
  if (pdist <= 0 || pdist > 1) {
    rlang::abort("'pdist' must be in (0,1]")
  }
  
  #--- pack args (for default projection)
  k <- min(8, gs_vcount(ps))
  pars <- list(baseline = baseline, pdist = pdist, k = k, 
    fill.cavity = fill.cavity, 
    decay.fun = weibullDecay(pdist=1))
  for (nm in names(pars)) {
    ps@pars_ps$silh[[nm]] <- pars[[nm]]
  }
  #--- run ps pipeline
  if(verbose) rlang::inform("Mapping graph silhouette...")
  ps <- .silhouetteCircular(ps, verbose)
  ps <- .updateStatus(ps, "Silhouette")
  return(ps)
})

#' @title Mapping Summits on PathwaySpace Images
#'
#' @description The \code{summitMapping} method implements a segmentation
#' strategy to identify summits on a 2D-landscape image 
#' (see \code{\link{summitWatershed}}).
#'
#' @param ps A \linkS4class{PathwaySpace} class object.
#' @param maxset A single positive integer indicating the maximum number 
#' of summits to be returned by the segmentation function.
#' @param minsize A single positive integer indicating the minimum size 
#' of the summits.
#' @param threshold A threshold provided as a fraction (in \code{[0,1]}) of the
#' max signal intensity.
#' @param segm.fun A segmentation function used to detect summits
#' (see \code{\link{summitWatershed}}).
#' @param ... Additional arguments passed to the segmentation function.
#' @return A preprocessed \linkS4class{PathwaySpace} class object.
#' @author Sysbiolab Team
#' @seealso \code{\link{circularProjection}}
#' @examples
#' library(PathwaySpace)
#' 
#' # Load a large igraph
#' data("PCv12_pruned_igraph", package = "PathwaySpace")
#' 
#' # Continue this example from the PathwaySpace vignette,
#' # in the 'PathwaySpace decoration' section
#'
#' @import methods
#' @docType methods
#' @rdname summitMapping-methods
#' @aliases summitMapping
#' @export
#'
setMethod("summitMapping", "PathwaySpace", function(ps, maxset = 30, 
  minsize = 30, threshold = 0.5, segm.fun = summitWatershed, ...) {
  
  ps <- updateGraphSpace(ps)

  #--- validate the pipeline status
  if (!.checkStatus(ps, "Projection")) {
    rlang::abort(c(
      "The 'ps' object has not been evaluated by a 'projection' method.",
      "i" = "Run a projection method on 'ps' before calling this function."
    ))
  }
  
  #--- validate argument types
  .validate.ps.args("singleInteger", "maxset", maxset)
  .validate.ps.args("singleInteger", "minsize", minsize)
  .validate.ps.args("singleNumber", "threshold", threshold)
  .validate.ps.args("function", "segm.fun", segm.fun)
  
  #--- validate argument values
  if (maxset < 1) {
    rlang::abort("'maxset' must be >=1")
  }
  if (minsize < 1) {
    rlang::abort("'minsize' must be >=1")
  }
  if (threshold < 0 || threshold > 1) {
    rlang::abort("'threshold' must be in [0,1]")
  }
  
  #--- pack args
  pars <- list(maxset = maxset, minsize = minsize,
    summit_threshold = threshold, segm.fun = segm.fun,
    segm_arg = list(...=...))
  for (nm in names(pars)) {
    ps@pars_ps$summit[[nm]] <- pars[[nm]]
  }
  #--- run ps pipeline
  ps <- .summitMapping(ps, ...=...)
  ps <- .updateStatus(ps, "Summits")
  return(ps)
})

#' @title Accessors for Fetching Slots from a PathwaySpace Object
#'
#' @description \code{getPathwaySpace} retrives information from
#' individual slots available in a PathwaySpace object.
#'
#' @param ps A preprocessed \linkS4class{PathwaySpace} class object
#' @param what A character value specifying which information should 
#' be retrieved from the slots.
#' Options: "nodes", "edges", "graph", "image", "pars", "misc", 
#' "signal","projection", "status", "silhouette", "summits", 
#' "summit_mask", "summit_contour"
#' @return Content from slots in the \linkS4class{PathwaySpace} object.
#' @examples
#' library(PathwaySpace)
#' 
#' # Load a demo igraph
#' data('gtoy1', package = 'RGraphSpace')
#'
#' # Create a new PathwaySpace object
#' ps <- buildPathwaySpace(gtoy1, nrc = 100)
#' # note: adjust 'nrc' to increase image resolution
#'
#' # Get the 'status' slot in ps
#' status <- getPathwaySpace(ps, what = 'status')
#'
#' @import methods
#' @docType methods
#' @rdname getPathwaySpace-methods
#' @aliases getPathwaySpace
#' @export
setMethod("getPathwaySpace", "PathwaySpace", function(ps, what = "status") {
  opts <- c("nodes", "edges", "graph", "image", "pars", "misc", 
    "projection", "status", "signal", "silhouette", "summits", 
    "summit_mask", "summit_contour")
  if (!what %in% opts) {
    opts <- paste0(opts, collapse = ", ")
    rlang::abort(sprintf("'what' must be one of: %s", opts))
  }
  
  .check_updated_ps(ps)
  
  if (what == "nodes") {
    obj <- ps@nodes
  } else if (what == "edges") {
    obj <- ps@edges
  } else if (what == "graph") {
    obj <- ps@graph
  } else if (what == "image") {
    obj <- ps@image
  } else if (what == "pars") {
    obj <- ps@pars_ps
  } else if (what == "projection") {
    obj <- ps@projection
  } else if (what == "misc") {
    obj <- ps@misc
  } else if (what == "status") {
    obj <- ps@status
  } else if (what == "signal") {
    obj <- gs_vertex_attr(ps, "signal")
  } else if (what == "silhouette") {
    obj <- ps@projection@floor
  } else if (what == "summits") {
    obj <- ps@misc$summits$lset
  } else if (what == "summit_mask") {
    obj <- ps@misc$summits$mset
  } else if (what == "summit_contour") {
    obj <- ps@misc$summits$cset
  }
  return(obj)
})

################################################################################
### Accessors
################################################################################

#-------------------------------------------------------------------------------
#' @title Accessor Functions for PathwaySpace Objects
#'
#' @description 
#' Get or set vertex signals, decay functions, and the active feature in a
#' \linkS4class{PathwaySpace} object.
#'
#' \code{vertexSignal()} gets or sets the numeric signal assigned to each
#' vertex, used as input for spatial projection.
#' 
#' \code{vertexDecay()} gets or sets the decay function assigned to each
#' vertex, controlling how the signal attenuates with distance.
#' 
#' \code{activeFeature()} gets or sets the active feature name, which
#' automatically extracts the corresponding signal from the \code{fdata} slot
#' or node attributes and assigns it to \code{vertexSignal()}.
#' 
#' @param x A \linkS4class{PathwaySpace} class object.
#' @param value The new value to assign:
#'   \itemize{
#'     \item For \code{vertexSignal()}: a numeric vector or scalar.
#'     \item For \code{vertexDecay()}: a decay function or list of decay
#'       functions (see \code{\link{linearDecay}}, \code{\link{weibullDecay}}).
#'     \item For \code{activeFeature()}: a single string matching a feature
#'       name (see \code{\link[RGraphSpace]{gs_features}}) or a node attribute
#'       (see \code{\link[RGraphSpace]{gs_names}}).
#'   }
#' @return The updated \linkS4class{PathwaySpace} object.
#' 
#' @examples
#' library(PathwaySpace)
#' 
#' # Load a demo igraph
#' data('gtoy1', package = 'RGraphSpace')
#' ps <- buildPathwaySpace(gtoy1, nrc = 100)
#' 
#' # Check vertex names
#' names(ps)
#' 
#' ##--------------------------------------
#' ## 'vertexSignal' accessor
#' 
#' # Access signal values from all vertices
#' vertexSignal(ps)
#' 
#' # Modify signal value of a specific vertex
#' vertexSignal(ps)[1] <- 1
#' 
#' # Modify signal value of specific vertices
#' vertexSignal(ps)[c("n2","n3")] <- 1
#' 
#' # Set '1s' to all vertices
#' vertexSignal(ps) <- 1
#' 
#' ##--------------------------------------
#' ## 'activeFeature' accessor
#' 
#' # Assign a signal feature matrix
#' signal_mtx <- matrix(
#'   rep(rnorm(gs_vcount(ps)), 2),
#'   ncol = 2,
#'   dimnames = list(names(ps), c("feature1", "feature2"))
#' )
#' gs_fdata(ps) <- signal_mtx
#' 
#' # Set the active feature — automatically updates vertexSignal()
#' activeFeature(ps) <- "feature1"
#' 
#' ##--------------------------------------
#' ## 'vertexDecay' accessor
#' 
#' # Access decay function of a specific vertex
#' vertexDecay(ps)[["n3"]]
#' 
#' # Modify decay function of a specific vertex
#' vertexDecay(ps)[["n3"]] <- linearDecay()
#' 
#' # Modify decay functions of two vertices
#' vertexDecay(ps)[c("n1","n3")] <- list( weibullDecay() )
#' 
#' # Modify decay functions of all vertices
#' vertexDecay(ps) <- weibullDecay(shape = 2)
#' 
#' @import methods
#' @docType methods
#' @rdname vertexSignal-accessors
#' @aliases vertexSignal
#' @aliases vertexSignal<-
#' @aliases vertexDecay
#' @aliases vertexDecay<-
#' @aliases activeFeature
#' @aliases activeFeature<-
#' @export
setMethod("vertexSignal", "PathwaySpace", function(x){
  gs_vertex_attr(x, "signal")
})

#' @rdname vertexSignal-accessors
#' @export
setMethod("vertexSignal<-", "PathwaySpace",
  function(x, value) {
    
    .check_updated_ps(x)
    
    if (!is.numeric(value) || !is.vector(value)){
      rlang::abort("'signal' must be a numeric vector or scalar.")
    }
    gs_vertex_attr(x, "signal") <- value
    return(x)
  }
)

#' @rdname vertexSignal-accessors
#' @export
setMethod("vertexDecay", "PathwaySpace", function(x){
  gs_vertex_attr(x, "decayFunction")
})

#' @rdname vertexSignal-accessors
#' @export
setMethod("vertexDecay<-", "PathwaySpace",
  function(x, value) {
    
    .check_updated_ps(x)
    
    gs_vertex_attr(x, "decayFunction") <- value
    
    return(x)
    
  }
)

#' @rdname vertexSignal-accessors
#' @export
setMethod("activeFeature", "PathwaySpace", function(x) {
  
  .check_updated_ps(x)
  
  feat <- x@pars_ps$active.feature
  
  if (is.null(feat) || length(feat) == 0) {
    return(NULL)
  }
  
  x@pars_ps$active.feature
  
})

#' @importFrom RGraphSpace gs_features gs_names gs_fdata gs_nodes
#' @rdname vertexSignal-accessors
#' @export
setReplaceMethod("activeFeature", "PathwaySpace", function(x, value) {
  
  .check_updated_ps(x)
  
  .validate.ps.args("singleString", "value", value)
  
  b1 <- value %in% gs_features(x)
  b2 <- value %in% gs_names(x)

  if(!b1 && !b2){
    rlang::abort(c(
      "x" = sprintf("Feature '%s' not found.", value),
      "i" = paste("Use `gs_features()` to list available features",
        "or `gs_names()` for node attribute names.")
    ))
  }
  if (b1 && b2) {
    rlang::warn(c(
      sprintf("Feature '%s' found in both feature matrix and node attributes.", value),
      "i" = "Using the feature matrix."
    ))
  }
  if(b1){
    rlang::inform(sprintf("Setting active feature '%s' from feature matrix...", value))
    vertexSignal(x) <- gs_fdata(x)[, value]
  } else {
    rlang::inform(sprintf("Setting active feature '%s' from node attributes...", value))
    vertexSignal(x) <- gs_nodes(x)[, value]
  }
  x@pars_ps$active.feature <- value
  x
})

#-------------------------------------------------------------------------------
#' @title Accessor Functions for PathwaySpace Objects
#'
#' @description Get or set edge and vertex attributes in
#' \linkS4class{PathwaySpace} class object.
#'
#' @param x A \linkS4class{PathwaySpace} class object.
#' @param name Name of the attribute.
#' @param value The new value of the attribute.
#' @param ... Additional arguments passed to igraph methods.
#' @return Updated \linkS4class{PathwaySpace} object.
#' @examples
#' library(PathwaySpace)
#' 
#' # Load a demo igraph
#' data('gtoy1', package = 'RGraphSpace')
#' ps <- buildPathwaySpace(gtoy1, nrc = 100)
#' 
#' # Get vertex count
#' gs_vcount(ps)
#' 
#' # Get edge count
#' gs_ecount(ps)
#' 
#' # Access a specific vertex attribute
#' gs_vertex_attr(ps, "signal")
#' 
#' # Replace an entire vertex attribute
#' gs_vertex_attr(ps, "signal") <- 1
#' 
#' # Modify a single value within a vertex attribute
#' gs_vertex_attr(ps, "signal")["n1"] <- 1
#'
#' # Access a specific edge attribute
#' gs_edge_attr(ps, "weight")
#' 
#' # Replace an entire edge attribute
#' gs_edge_attr(ps, "weight") <- 1
#' 
#' @rdname PathwaySpace-accessors
#' @importFrom RGraphSpace gs_vertex_attr<- gs_edge_attr<-
#' @importFrom RGraphSpace gs_vertex_attr gs_edge_attr
#' @importFrom RGraphSpace gs_ecount gs_vcount
#' @aliases gs_vertex_attr<-
#' @export
setReplaceMethod(
  "gs_vertex_attr","PathwaySpace", function(x, name, ..., value) {
    
    .check_updated_ps(x)
    
    # Call the GraphSpace method first
    x <- callNextMethod()
    
    x <- .validate_ps_containers(x, changed = name)
    
    return(x)
  }
)

#' @rdname PathwaySpace-accessors
#' @aliases gs_edge_attr<-
#' @export
setReplaceMethod(
  "gs_edge_attr","PathwaySpace", function(x, name, ..., value) {
    
    .check_updated_ps(x)
    
    # Call the GraphSpace method first
    x <- callNextMethod()
    
    x <- .validate_ps_containers(x, changed = name)
    
    return(x)
  }
)

#-------------------------------------------------------------------------------
.validate_ps_containers <- function(ps, changed = NULL) {
  if (is.null(changed) || changed == "signal") {
    ps <- .validate_signal(ps)
  }
  if (is.null(changed) || changed == "weight") {
    ps <- .validate_weights(ps)
  }
  if (is.null(changed) || changed == "decayFunction") {
    ps <- .validate_decayFunction(ps)
  }
  return(ps)
}

#-------------------------------------------------------------------------------
.validate_signal <- function(ps) {
  signal <- gs_vertex_attr(ps, "signal")
  if(is.null(signal)){
    rlang::abort("'signal' vertex attribute must be available.")
  }
  if (!is.numeric(signal)){
    rlang::abort("vertex 'signal' variable must be numeric.")
  }
  ps@nodes$signal <- .revise_signal(signal)
  return(ps)
}
.revise_signal <- function(x){
  x[is.nan(x)] <- NA
  x[x == Inf] <- NA
  x[x == -Inf] <- NA
  if (all(is.na(x))) x[] <- 0
  return(x)
}

#-------------------------------------------------------------------------------
.validate_weights <- function(ps) {
  if(gs_ecount(ps)>0){
    weight <- gs_edge_attr(ps, "weight")
    if(is.null(weight)){
      rlang::abort("'weight' edge attribute must be available.")
    }
    if (!is.numeric(weight)){
      rlang::abort("edge 'weight' variable must be numeric.")
    }
    ps@edges$weight <- .revise_weights(weight)
  }
  return(ps)
}
.revise_weights <- function(wt){
  if (all(is.na(wt))) wt[] <- 1
  s <- sd(wt, na.rm = TRUE)
  if (!is.na(s) && s != 0) {
    wt <- wt/max(abs(wt), na.rm = TRUE)
    wt[is.na(wt)] <- 0
  } else {
    wt[] <- 1
  }
  return(wt)
}

#-------------------------------------------------------------------------------
.validate_decayFunction <- function(ps){
  att <- names(gs_vertex_attr(ps))
  if(! "decayFunction" %in% att){
    rlang::abort("Missing a vertex 'decayFunction' attribute.")
  }
  decayFunction <- gs_vertex_attr(ps, "decayFunction")
  # check function, args, and vertex attributes
  lg <- unlist(lapply(decayFunction, is.function))
  if(!all(lg)){
    rlang::abort(c(
      "Vertex 'decayFunction' attribute is invalid.",
      "i" = "Each vertex 'decay function' must be a function.",
      "i" = "e.g. function(x, signal) { ... }"
    ))
  }
  not_used <- lapply(decayFunction, .check_decay_args, nodes=ps@nodes)
  if(.all_equal_fun(decayFunction)){
    dfun <- attributes(decayFunction[[1]])$name
    dfun <- ifelse(.is_singleString(dfun), dfun, "customized")
    ps@pars_ps$decay$fun <- dfun
    ps@pars_ps$decay$info <- "global-defined-decay"
  } else {
    ps@pars_ps$decay$fun <- "customized"
    ps@pars_ps$decay$info <- "local-defined-decay"
  }
  ps@pars_ps$decay$is_default_args <- .is_default_args(decayFunction)
  return(ps)
}
.check_decay_args <- function(decay_fun, nodes, args = c("x","signal")){
  fargs <- formalArgs(args(decay_fun))
  missing_args <- setdiff(args, fargs)
  if (length(missing_args) > 0) {
    rlang::abort(c(
      "Invalid 'decay function'.",
      "i" = "Expected 'x' and 'signal' arguments not found.",
      "i" = "e.g. function(x, signal) { ... }"
    ))
  }
  
  decay_args <- c(args, setdiff(fargs, args))
  
  if ("..." %in% decay_args) {
    rlang::abort(c(
      "Invalid 'decay function'.",
      "i" = "'...' is not supported; every argument must be named and match a vertex attribute.",
      "i" = "e.g. function(x, signal, weight) { ... }"
    ))
  }
  
  if(!all(decay_args %in% colnames(nodes))){
    extra_args <- decay_args[!decay_args %in% colnames(nodes)]
    extra_args <- paste0(sQuote(extra_args, q=FALSE), collapse = ", ")
    rlang::abort(c(
      "Each 'decay function' argument must correspond to a vertex attribute.",
      "i" = paste0("The following argument(s) do not match any vertex attribute: ", 
        extra_args)
    ))
  }
  
  TRUE
  
}

.is_default_args <- function(decayFunction, args = c("x","signal")){
  fargs <- lapply(decayFunction, formalArgs)
  n_args <- unlist(lapply(fargs, length))
  if(all(n_args==length(args))){
    fargs <- as.character(unlist(fargs))
    is_default <- all( fargs %in% args )
  } else {
    is_default <- FALSE
  }
  return(is_default)
}
.all_equal_fun <- function(lst){
  all(unlist(lapply(lst[-1], function(x) identical(x, lst[[1]]))))
}

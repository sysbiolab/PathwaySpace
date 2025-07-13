
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
#' @param verbose A single logical value specifying to display detailed 
#' messages (when \code{verbose=TRUE}) or not (when \code{verbose=FALSE}).
#' @param g Deprecated from PathwaySpace 1.0.1; use 'gs' instead.
#' @param mar Deprecated from PathwaySpace 1.0.1; use 'mar' in  
#' \code{\link[RGraphSpace]{GraphSpace}} instead.
#' @return A pre-processed \linkS4class{PathwaySpace} class object.
#' @author Victor Apolonio, Vinicius Chagas, Mauro Castro,
#' and TCGA Network.
#' @seealso \code{\link[igraph]{undirected_graph}}
#' @examples
#' # Load a demo igraph
#' data('gtoy1', package = 'RGraphSpace')
#' 
#' # Check graph validity
#' gs <- GraphSpace(gtoy1, mar = 0.1)
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
#' @importFrom RGraphSpace GraphSpace getGraphSpace
#' @importFrom RANN nn2
#' @aliases buildPathwaySpace
#' @export
#' 
buildPathwaySpace <- function(gs, nrc = 500, verbose = TRUE, 
  g = deprecated(), mar = deprecated()) {
  if(verbose) message("Validating arguments...")
  #--- validate argument types
  .validate.args("singleInteger", "nrc", nrc)
  .validate.args("singleLogical", "verbose", verbose)
  ### deprecate
  if (lifecycle::is_present(g)) {
    deprecate_soft("1.0.1", "buildPathwaySpace(g)", 
      "buildPathwaySpace(gs)")
    gs <- g
  }
  if (lifecycle::is_present(mar)) {
    deprecate_soft("1.0.1", "buildPathwaySpace(mar)", 
      "GraphSpace(mar)")
  }
  ###
  #--- validate argument values
  if (nrc < 2) {
    stop("'nrc' should be >=2", call. = FALSE)
  }
  #--- validate the graph object
  if(is_igraph(gs)){
    if(verbose) message("Validating the 'igraph' object...")
    gs <- GraphSpace(gs, verbose=FALSE)
  }
  .validate.gspace(gs)
  
  #--- build PathwaySpace-class
  ps <- .buildPathwaySpace(gs, nrc, verbose)
  ps <- .updateStatus(ps, "Preprocess")
  return(ps)
}

#' @title Spatial Projection of Graph-Associated Signals
#'
#' @description \code{circularProjection} implements a convolution
#' algorithm to project signals onto a 2D-coordinate system.
#'
#' @param ps A \linkS4class{PathwaySpace} class object.
#' @param k A single positive integer determining the k-top signals in the 
#' signal convolution operation.
#' @param pdist A term (in \code{[0, 1]}) determining a distance unit for the
#' signal decay function. When `pdist = 1` it will represent the diameter of 
#' the inscribed circle within the pathway space. This distance will affect the
#' extent over which the convolution operation projects the signal between
#' source- and destination points.
#' @param rescale A single logical value indicating whether to rescale 
#' the signal. If the signal \code{>=0}, then it will be rescaled to 
#' \code{[0, 1]}; if the signal \code{<=0}, then it will be rescaled to 
#' \code{[-1, 0]}; and if the signal in \code{(-Inf, +Inf)}, then it will be 
#' rescaled to \code{[-1, 1]}.
#' @param verbose A single logical value specifying to display detailed 
#' messages (when \code{verbose=TRUE}) or not (when \code{verbose=FALSE}).
#' @param decay.fun A signal decay function. Available options include 
#' 'Weibull', 'exponential', and 'linear' (see \code{\link{signalDecay}}).
#' Users may also define a custom decay model with at least two arguments, 
#' e.g., \code{function(x, signal) \{ ... \}}, which should returns a vector of  
#' projected signals of the same length as \code{x}. Additional arguments may  
#' include any variable available as graph vertex attribute.
#' @param aggregate.fun A function used to aggregate the projected signals. 
#' It must be provided as a unary function, e.g., \code{function(x) { ... }}, 
#' which should aggregate a vector of signals to a single scalar value. 
#' Available options include 'mean', 'wmean', 'log.wmean', and 'exp.wmean' 
#' (See \code{\link{signalAggregation}}).
#' @param kns Deprecated from PathwaySpace 1.0.1; use 'k' instead.
#' @return A preprocessed \linkS4class{PathwaySpace} class object.
#' @author Victor Apolonio, Vinicius Chagas, Mauro Castro, 
#' and TCGA Network.
#' @seealso \code{\link{buildPathwaySpace}}
#' @examples
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
#' @importFrom lifecycle deprecated deprecate_soft is_present
#' @docType methods
#' @rdname circularProjection-methods
#' @aliases circularProjection
#' @export
#'
setMethod("circularProjection", "PathwaySpace", function(ps, k = 8,
  pdist = 0.15, rescale = TRUE, verbose = TRUE, 
  decay.fun = signalDecay(), aggregate.fun = signalAggregation(),
  kns = deprecated()) {
  #--- validate the pipeline status
  if (!.checkStatus(ps, "Preprocess")) {
    stop("NOTE: the 'ps' object needs preprocessing!", call. = FALSE)
  }
  ### deprecate
  if (lifecycle::is_present(kns)) {
    deprecate_soft("1.0.1", "circularProjection(kns)", 
      "circularProjection(k)")
    k <- kns
  }
  ###
  if(verbose) message("Validating arguments...")
  #--- validate argument types
  .validate.args("singleInteger", "k", k)
  .validate.args("singleNumber", "pdist", pdist)
  .validate.args("singleLogical", "rescale", rescale)
  .validate.args("singleLogical", "verbose", verbose)
  .validate.args("function", "aggregate.fun", aggregate.fun)
  .validate.args("function", "decay.fun", decay.fun)
  #--- validate argument values
  if (k < 1) {
    stop("'k' should be >=1", call. = FALSE)
  }
  n <- gs_vcount(ps)
  if (k > n) k <- n
  if (pdist < 0 || pdist > 1) {
    stop("'pdist' should be in [0,1]", call. = FALSE)
  }
  #--- validate functions
  if(!missing(decay.fun)){
    gs_vertex_attr(ps, "decayFunction") <- decay.fun
  }
  .validate_aggregate_fun(aggregate.fun)
  
  #--- pack args
  pars <- list(k = k, pdist = pdist, rescale = rescale, 
    aggregate.fun = aggregate.fun, projection="Circular")
  for (nm in names(pars)) {
    ps@pars$ps[[nm]] <- pars[[nm]]
  }
  #--- run ps pipeline
  ps <- .circularProjection(ps, verbose)
  ps <- .updateStatus(ps, "CircularProjection")
  if (.checkStatus(ps, "PolarProjection")) {
    if(verbose) message("-- polar projection replaced by circular.")
    ps <- .updateStatus(ps, "PolarProjection", FALSE)
  }
  return(ps)
})

#' @title Spatial Projection of Graph-Associated Signals
#'
#' @description \code{polarProjection} implements a convolution algorithm
#' to project signals across a 2D-coordinate system.
#'
#' @param ps A \linkS4class{PathwaySpace} class object.
#' @param k A single positive integer determining the k-top signals in the 
#' signal convolution operation.
#' @param pdist A term (in \code{[0, 1]}) determining a distance unit for the
#' signal decay function. When `pdist = 1` it will represent the diameter of 
#' the inscribed circle within the pathway space. This distance will affect the
#' extent over which the convolution operation projects the signal between
#' source- and destination points.
#' @param theta Angle of projection (degrees in \code{(0,360]}).
#' @param directional If directional edges are available, this argument can 
#' be used to orientate the signal projection on directed graphs.
#' @param rescale A single logical value indicating whether to rescale 
#' the signal. If the signal \code{>=0}, then it will be rescaled to 
#' \code{[0, 1]}; if the signal \code{<=0}, then it will be rescaled to 
#' \code{[-1, 0]}; and if the signal in \code{(-Inf, +Inf)}, then it will be 
#' rescaled to \code{[-1, 1]}.
#' @param verbose A single logical value specifying to display detailed 
#' messages (when \code{verbose=TRUE}) or not (when \code{verbose=FALSE}).
#' @param decay.fun A signal decay function. Available options include 
#' 'Weibull', 'exponential', and 'linear' (see \code{\link{signalDecay}}).
#' Users may also define a custom decay model with at least two arguments, 
#' e.g., \code{function(x, signal) \{ ... \}}, which should returns a vector of  
#' projected signals of the same length as \code{x}. Additional arguments may  
#' include any variable available as graph vertex attribute.
#' @param aggregate.fun A function used to aggregate the projected signals. 
#' It must be provided as a unary function, e.g., \code{function(x) { ... }}, 
#' which should aggregate a vector of signals to a single scalar value. 
#' Available options include 'mean', 'wmean', 'log.wmean', and 'exp.wmean' 
#' (See \code{\link{signalAggregation}}).
#' @param kns Deprecated from PathwaySpace 1.0.1; use 'k' instead.
#' @return A preprocessed \linkS4class{PathwaySpace} class object.
#' @author Mauro Castro
#' @seealso \code{\link{buildPathwaySpace}}
#' @examples
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
#' # Create a 2D-landscape image
#' ps <- polarProjection(ps)
#'
#' @import methods
#' @docType methods
#' @rdname polarProjection-methods
#' @aliases polarProjection
#' @export
#'
setMethod("polarProjection", "PathwaySpace", function(ps, k = 8, 
  pdist = 0.5, theta = 180, directional = FALSE, rescale = TRUE, 
  verbose = TRUE, decay.fun = signalDecay(), 
  aggregate.fun = signalAggregation(),
  kns = deprecated()) {
  #--- validate the pipeline status
  if (!.checkStatus(ps, "Preprocess")) {
    stop("NOTE: the 'ps' object needs preprocessing!", call. = FALSE)
  }
  ### deprecate
  if (lifecycle::is_present(kns)) {
    deprecate_soft("1.0.1", "polarProjection(kns)", 
      "polarProjection(k)")
    k <- kns
  }
  ###
  if(verbose) message("Validating arguments...")
  .validate.args("singleInteger", "k", k)
  .validate.args("singleNumber", "pdist", pdist)
  .validate.args("singleNumber", "theta", theta)
  .validate.args("singleLogical", "directional", directional)
  .validate.args("singleLogical", "rescale", rescale)
  .validate.args("singleLogical", "verbose", verbose)
  .validate.args("function", "decay.fun", decay.fun)
  .validate.args("function", "aggregate.fun", aggregate.fun)
  if (k < 1) {
    stop("'k' should be >=1", call. = FALSE)
  }
  n <- gs_vcount(ps)
  if (k > n) k <- n
  if (pdist < 0 || pdist > 1) {
    stop("'pdist' should be in [0,1]", call. = FALSE)
  }
  if (theta <= 0 || theta > 360) {
    msg <- paste0("'polar.theta' should be an angle of projection\n",
      "with degrees in (0,360]")
    stop(msg, call. = FALSE)
  }
  
  #--- validate functions
  if(!missing(decay.fun)){
    gs_vertex_attr(ps, "decayFunction") <- decay.fun
  }
  .validate_aggregate_fun(aggregate.fun)
  
  #--- pack args
  pars <- list(k = k, pdist = pdist, theta = theta, rescale = rescale, 
    projection="Polar", aggregate.fun = aggregate.fun,
    directional = directional)
  for (nm in names(pars)) {
    ps@pars$ps[[nm]] <- pars[[nm]]
  }
  ps <- .polarProjection(ps, verbose)
  ps <- .updateStatus(ps, "PolarProjection")
  if (.checkStatus(ps, "CircularProjection")) {
    if(verbose) message("-- circular projection replaced by polar.")
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
#' @param baseline A fraction (in \code{[0,1]}) of the signal scale of a 
#' PathwaySpace image. This term only affects the image baseline projection, 
#' which represents a silhouette of the graph's layout outlined in the 
#' resulting image. When \code{baseline = 0} (i.e. lower level of the signal 
#' scale), the baseline will extend over the entire image space, so no 
#' silhouette will be visible.
#' @param pdist A term (in \code{[0,1]}) determining a distance unit for the
#' signal decay function. This distance will affect the extent over which 
#' the convolution operation projects the image baseline.
#' @param fill.cavity A single logical value specifying to fill cavities 
#' in the silhouette mask (when \code{verbose=TRUE}) or not 
#' (when \code{verbose=FALSE}).
#' @param verbose A single logical value specifying to display detailed 
#' messages (when \code{verbose=TRUE}) or not (when \code{verbose=FALSE}).
#' @return A preprocessed \linkS4class{PathwaySpace} class object.
#' @author Mauro Castro and TCGA Network.
#' @seealso \code{\link{circularProjection}}
#' @examples
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
setMethod("silhouetteMapping", "PathwaySpace", function(ps, baseline = 0.01,
  pdist = 0.05, fill.cavity = TRUE, verbose = TRUE) {
  #--- validate the pipeline status
  if (!.checkStatus(ps, "Preprocess")) {
    stop("NOTE: the 'ps' object needs preprocessing!", call. = FALSE)
  }
  if(verbose) message("Validating arguments...")
  #--- validate argument types
  .validate.args("singleNumber", "baseline", baseline)
  .validate.args("singleNumber", "pdist", pdist)
  .validate.args("singleLogical", "fill.cavity", fill.cavity)
  .validate.args("singleLogical", "verbose", verbose)
  #--- validate argument values
  if (baseline < 0 || baseline > 1) {
    stop("'baseline' should be in [0,1]", call. = FALSE)
  }
  if (pdist < 0 || pdist > 1) {
    stop("'pdist' should be in [0,1]", call. = FALSE)
  }
  #--- pack args (for default projection)
  k <- min(8, gs_vcount(ps))
  pars <- list(baseline = baseline, pdist = pdist, k = k, 
    fill.cavity = fill.cavity, decay.fun = signalDecay())
  for (nm in names(pars)) {
    ps@pars$ps$silh[[nm]] <- pars[[nm]]
  }
  #--- run ps pipeline
  if(verbose) message("Mapping graph silhouette...")
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
#' @param verbose A single logical value specifying to display detailed 
#' messages (when \code{verbose=TRUE}) or not (when \code{verbose=FALSE}).
#' @param segm_fun A segmentation function used to detect summits
#' (see \code{\link{summitWatershed}}).
#' @param ... Additional arguments passed to the segmentation function.
#' @return A preprocessed \linkS4class{PathwaySpace} class object.
#' @author Mauro Castro and TCGA Network.
#' @seealso \code{\link{circularProjection}}
#' @examples
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
    minsize = 30, threshold = 0.5, verbose = TRUE, 
    segm_fun = summitWatershed, ...) {
    #--- validate the pipeline status
    if (!.checkStatus(ps, "Projection")) {
        msg <- paste0("NOTE: the 'ps' object needs to be\n",
            "evaluated by a 'projection' method!")
        stop(msg, call. = FALSE)
    }
    #--- validate argument types
    .validate.args("singleInteger", "maxset", maxset)
    .validate.args("singleInteger", "minsize", minsize)
    .validate.args("singleNumber", "threshold", threshold)
    .validate.args("singleLogical", "verbose", verbose)
    .validate.args("function", "segm_fun", segm_fun)
    #--- validate argument values
    if (maxset < 1) {
        stop("'maxset' should be >=1", call. = FALSE)
    }
    if (minsize < 1) {
        stop("'minsize' should be >=1", call. = FALSE)
    }
    if (threshold < 0 || threshold > 1) {
        stop("'threshold' should be in [0,1]", call. = FALSE)
    }
    #--- pack args
    pars <- list(maxset = maxset, minsize = minsize,
        summit_threshold = threshold, segm_fun = segm_fun,
        segm_arg = list(...=...))
    for (nm in names(pars)) {
        ps@pars$ps$summit[[nm]] <- pars[[nm]]
    }
    #--- run ps pipeline
    ps <- .summitMapping(ps, verbose, ...=...)
    ps <- .updateStatus(ps, "Summits")
    return(ps)
})

#' @title Accessors for Fetching Slots from a PathwaySpace Object
#'
#' @description \code{getPathwaySpace} retrives information from
#' individual slots available in a PathwaySpace object.
#'
#' @param ps A preprocessed \linkS4class{PathwaySpace} class object
#' @param what A single character value specifying which information should 
#' be retrieved from the slots.
#' Options: "nodes", "edges", "graph", "image", "pars", "misc", 
#' "signal","projections", "status", "silhouette", "summits", 
#' "summit_mask", "summit_contour"
#' @return Content from slots in the \linkS4class{PathwaySpace} object.
#' @examples
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
      "projections", "status", "signal", "silhouette", "summits", 
      "summit_mask", "summit_contour")
    if (!what %in% opts) {
        opts <- paste0(opts, collapse = ", ")
        stop("'what' must be one of:\n", opts, call. = FALSE)
    }
    if (what == "nodes") {
        obj <- ps@nodes
    } else if (what == "edges") {
        obj <- ps@edges
    } else if (what == "graph") {
      obj <- ps@graph   
    } else if (what == "image") {
      obj <- ps@image   
    } else if (what == "pars") {
      obj <- ps@pars
    } else if (what == "misc") {
      obj <- ps@misc
    } else if (what == "projections") {
        obj <- ps@projections
    } else if (what == "status") {
        obj <- ps@status
    } else if (what == "signal") {
      obj <- gs_vertex_attr(ps, "signal")
    } else if (what == "silhouette") {
        obj <- ps@projections$xfloor
    } else if (what == "summits") {
        obj <- ps@projections$summits$lset
    } else if (what == "summit_mask") {
        obj <- ps@projections$summits$mset
    } else if (what == "summit_contour") {
        obj <- ps@projections$summits$cset
    }
    return(obj)
})

################################################################################
### Accessors
################################################################################

#-------------------------------------------------------------------------------
# show summary information on screen
setMethod("show", "PathwaySpace", function(object) {
  message("A PathwaySpace-class object for:")
  summary(getGraphSpace(object, what = "graph"))
  cat("+ status:", .summariseStatus(object))
})

#-------------------------------------------------------------------------------
#' @title Accessor Functions for PathwaySpace Objects
#'
#' @description Get or set 'signal' and 'decay' functions in a 
#' \linkS4class{PathwaySpace} class object.
#'
#' @param x A \linkS4class{PathwaySpace} class object.
#' @param value The new value of the attribute.
#' @return Updated \linkS4class{PathwaySpace} object.
#' @examples
#' data('gtoy1', package = 'RGraphSpace')
#' ps <- buildPathwaySpace(gtoy1, nrc = 100)
#' 
#' # Check vertex names
#' names(ps)
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
#' #----
#' 
#' # Access decay function of a specific vertex
#' vertexDecay(ps)[["n3"]]
#' 
#' # Modify decay function of a specific vertex
#' vertexDecay(ps)[["n3"]] <- signalDecay(method = "linear")
#' 
#' # Modify decay functions of two vertices
#' vertexDecay(ps)[c("n1","n3")] <- list( signalDecay() )
#' 
#' # Modify decay functions of all vertices
#' vertexDecay(ps) <- signalDecay(method = "weibull", shape = 2)
#' 
#' @import methods
#' @docType methods
#' @rdname vertexSignal-accessors
#' @aliases vertexSignal
#' @aliases vertexSignal<-
#' @aliases vertexDecay
#' @aliases vertexDecay<-
#' @export
setMethod("vertexSignal", "PathwaySpace", function(x){
  gs_vertex_attr(x, "signal")
})

#' @rdname vertexSignal-accessors
#' @export
setMethod("vertexSignal<-", "PathwaySpace",
  function(x, value) {
    if (!is.numeric(value) || !is.vector(value)){
      stop("'value' must be a numeric vector or scalar.", call. = FALSE)
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
    gs_vertex_attr(x, "decayFunction") <- value
    return(x)
  }
)

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
    
    # Call the GraphSpace method first
    x <- callNextMethod()
    
    x <- .validate_ps_containers(x)
    
    return(x)
  }
)

#' @rdname PathwaySpace-accessors
#' @aliases gs_edge_attr<-
#' @export
setReplaceMethod(
  "gs_edge_attr","PathwaySpace", function(x, name, ..., value) {
    
    # Call the GraphSpace method first
    x <- callNextMethod()
    
    x <- .validate_ps_containers(x)
    
    return(x)
  }
)

#-------------------------------------------------------------------------------
.validate_ps_containers <- function(ps) {
  ps <- .validate_signal(ps)
  ps <- .validate_weights(ps)
  ps <- .validate_decayFunction(ps)
  return(ps)
}

#-------------------------------------------------------------------------------
.validate_signal <- function(ps) {
  signal <- gs_vertex_attr(ps, "signal")
  if(is.null(signal)){
    stop("'signal' vertex attribute must be available.", call. = FALSE)
  }
  if (!is.numeric(signal)){
    stop("vertex 'signal' variable must be numeric.", call. = FALSE)
  }
  signal <- .revise_signal(signal)
  ps@pars$ps$zscale <- .get_signal_scale(signal)
  ps@nodes$signal <- signal
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
      stop("'weight' edge attribute must be available.", call. = FALSE)
    }
    if (!is.numeric(weight)){
      stop("edge 'weight' variable must be numeric.", call. = FALSE)
    }
    ps@edges$weight <- .revise_weights(weight)
  }
  return(ps)
}
.revise_weights <- function(wt){
  if (all(is.na(wt))) wt[] <- 1
  wt[is.na(wt)] <- min(wt, na.rm = TRUE)
  if (sd(wt) > 0) {
    wt <- wt-min(wt)
    wt <- wt/max(wt)
  } else {
    wt[] <- 1
  }
  return(wt)
}

#-------------------------------------------------------------------------------
.validate_decayFunction <- function(ps){
  att <- names(gs_vertex_attr(ps))
  if(! "decayFunction" %in% att){
    stop("Missing a vertex 'decayFunction' attribute.", call. = FALSE)
  }
  decayFunction <- gs_vertex_attr(ps, "decayFunction")
  # check function, args, and vertex attributes
  lg <- unlist(lapply(decayFunction, is.function))
  if(!all(lg)){
    msg1 <- "Each vertex 'decay function' must be a function, e.g.,\n"
    msg2 <- "function(x, signal) { ... }"
    stop(msg1, msg2, call. = FALSE)
  }
  not_used <- lapply(decayFunction, .check_decay_args, nodes=ps@nodes)
  if(.all_equal_fun(decayFunction)){
    ps@pars$ps$decay$fun <- decayFunction[[1]]
    ps@pars$ps$decay$info <- "global-defined-decay"
  } else {
    ps@pars$ps$decay$fun <- "customized"
    ps@pars$ps$decay$info <- "local-defined-decay"
  }
  ps@pars$ps$decay$is_default_args <- .is_default_args(decayFunction)
  return(ps)
}
.check_decay_args <- function(decay_fun, nodes, args = c("x","signal")){
  fargs <- formalArgs(args(decay_fun))
  missing_args <- setdiff(args, fargs)
  if (length(missing_args) > 0) {
    msg <- paste0("Invalid 'decay function':")
    msg <- paste0(msg, " expected arguments 'x' and 'signal' not found.")
    stop(msg, call. = FALSE)
  }
  decay_args <- c(args, setdiff(fargs, args))
  
  if(!all(decay_args %in% colnames(nodes))){
    extra_args <- decay_args[!decay_args %in% colnames(nodes)]
    msg1 <- "Each 'decay function' argument must correspond to a vertex attribute.\n"
    msg2 <- paste0(sQuote(extra_args, q=FALSE), collapse = ", ")
    msg2 <- paste0("The following argument(s) do not match any vertex attribute: ",
      msg2)
    stop(msg1, msg2, call. = FALSE)
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

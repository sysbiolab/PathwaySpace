#' @title  Constructor of PathwaySpace-class objects.
#'
#' @description \code{buildPathwaySpace} is a constructor of
#' PathwaySpace-class objects.
#'
#' @param gs A \code{\link[RGraphSpace]{GraphSpace}} object. Alternatively, 
#' an \code{\link[igraph]{igraph}} object with vertex coordinates assigned 
#' to \code{x} and \code{y} vertex attributes, and vertex labels assigned 
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

#' @title Creating 2D-landscape images from graph objects.
#'
#' @description \code{circularProjection} implements a convolution
#' algorithm to project signals across a 2D-coordinate system.
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
#' @param decay.fun A signal decay function. Available: 'Weibull',
#' 'exponential', and 'linear' functions (see \code{\link{signalDecay}}).
#' @param aggregate.fun A signal aggregation function. Available: 
#' weighted 'linear', 'log', and 'exponential' functions 
#' (see \code{\link{signalAggregation}}).
#' @param kns Deprecated from PathwaySpace 1.0.1; use 'k' instead.
#' @param knn Deprecated from PathwaySpace 1.0.1; use 'k' instead.
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
  kns = deprecated(), knn = deprecated()) {
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
    if (lifecycle::is_present(knn)) {
      deprecate_soft("1.0.1", "circularProjection(knn)", 
        "circularProjection(k)")
      k <- knn
    }
    ###
    if(verbose) message("Validating arguments...")
    #--- validate argument types
    .validate.args("singleInteger", "k", k)
    .validate.args("singleNumber", "pdist", pdist)
    .validate.args("singleLogical", "rescale", rescale)
    .validate.args("singleLogical", "verbose", verbose)
    .validate.args("function", "decay.fun", decay.fun)
    .validate.args("function", "aggregate.fun", aggregate.fun)
    #--- validate argument values
    if (k < 1) {
        stop("'k' should be >=1", call. = FALSE)
    }
    n <- length(getPathwaySpace(ps, "vertex"))
    if (k > n) k <- n
    if (pdist < 0 || pdist > 1) {
        stop("'pdist' should be in [0,1]", call. = FALSE)
    }
    #--- pack args
    pars <- list(k = k, pdist = pdist, rescale = rescale, 
      projection="Circular", decay_fun = decay.fun, 
      aggregate_fun = aggregate.fun)
    for (nm in names(pars)) {
        ps@pars$proj[[nm]] <- pars[[nm]]
    }
    #--- run ps pipeline
    ps <- .circularProjection(ps, verbose)
    ps <- .updateStatus(ps, "CircularProjection")
    if (.checkStatus(ps, "PolarProjection")) {
        if(verbose) message("-- polar projection replaced by circular.")
        ps <- .updateStatus(ps, "PolarProjection", FALSE)
    }
    # ps <- .removeSilhouette(ps, verbose)
    return(ps)
})

#' @title Creating 2D-landscape images from graph objects.
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
#' @param decay.fun A signal decay function. Available: 'Weibull', 
#' 'exponential', and 'linear' functions (see \code{\link{signalDecay}}).
#' @param aggregate.fun A signal aggregation function. Available: 
#' weighted 'linear', 'log', and 'exponential' functions 
#' (see \code{\link{signalAggregation}}).
#' @param kns Deprecated from PathwaySpace 1.0.1; use 'k' instead.
#' @param knn Deprecated from PathwaySpace 1.0.1; use 'k' instead.
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
  pdist = 0.5, theta = 180, directional = FALSE, 
  rescale = TRUE, verbose = TRUE, decay.fun = signalDecay(), 
  aggregate.fun = signalAggregation(), 
  kns = deprecated(), knn = deprecated()) {
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
    if (lifecycle::is_present(knn)) {
      deprecate_soft("1.0.1", "polarProjection(knn)", 
        "polarProjection(k)")
      k <- knn
    }
    ###
    if(verbose) message("Validating arguments...")
    .validate.args("singleInteger", "k", k)
    .validate.args("singleNumber", "pdist", pdist)
    .validate.args("singleNumber", "theta", theta)
    .validate.args("singleLogical", "rescale", rescale)
    .validate.args("singleLogical", "verbose", verbose)
    .validate.args("singleLogical", "directional", directional)
    .validate.args("function", "decay.fun", decay.fun)
    .validate.args("function", "aggregate.fun", aggregate.fun)
    if (k < 1) {
        stop("'k' should be >=1", call. = FALSE)
    }
    n <- length(getPathwaySpace(ps, "vertex"))
    if (k > n) k <- n
    if (pdist < 0 || pdist > 1) {
        stop("'pdist' should be in [0,1]", call. = FALSE)
    }
    if (theta <= 0 || theta > 360) {
        msg <- paste0("'polar.theta' should be an angle of projection\n",
            "with degrees in (0,360]")
        stop(msg, call. = FALSE)
    }
    pars <- list(k = k, pdist = pdist, theta = theta,  
      rescale = rescale, projection="Polar", directional = directional,
      decay_fun = decay.fun, aggregate_fun = aggregate.fun)
    for (nm in names(pars)) {
        ps@pars$proj[[nm]] <- pars[[nm]]
    }
    ps <- .polarProjection(ps, verbose)
    ps <- .updateStatus(ps, "PolarProjection")
    if (.checkStatus(ps, "CircularProjection")) {
        if(verbose) message("-- circular projection replaced by polar.")
        ps <- .updateStatus(ps, "CircularProjection", FALSE)
    }
    # ps <- .removeSilhouette(ps, verbose)
    return(ps)
})

#' @title Decorating PathwaySpace images with graph silhouettes.
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
    n <- length(getPathwaySpace(ps, "vertex"))
    k <- min(8, n)
    pars <- list(baseline = baseline, pdist = pdist, k = k, 
      fill.cavity = fill.cavity, decay_fun = signalDecay())
    for (nm in names(pars)) {
        ps@pars$silh[[nm]] <- pars[[nm]]
    }
    #--- run ps pipeline
    if(verbose) message("Mapping graph silhouette...")
    ps <- .silhouetteCircular(ps, verbose)
    ps <- .updateStatus(ps, "Silhouette")
    return(ps)
})

#' @title Mapping summits on PathwaySpace images.
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
        ps@pars$summit[[nm]] <- pars[[nm]]
    }
    #--- run ps pipeline
    ps <- .summitMapping(ps, verbose, ...=...)
    ps <- .updateStatus(ps, "Summits")
    return(ps)
})

#' @title Accessors for fetching slots from a PathwaySpace object.
#'
#' @description \code{getPathwaySpace} retrives information from
#' individual slots available in a PathwaySpace object.
#'
#' @param ps A preprocessed \linkS4class{PathwaySpace} class object
#' @param what A single character value specifying which information should 
#' be retrieved from the slots.
#' Options: 'graph','gxy','gxyz','pars','misc','status','summits',
#' 'summit_mask', and 'summit_contour'.
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
    opts <- c("gspace","vertex", "vsignal", "nodes", "edges", 
        "gxy", "gxyz", "pars", "misc", "status", "silhouette",
        "summits", "summit_mask", "summit_contour")
    if (!what %in% opts) {
        opts <- paste0(opts, collapse = ", ")
        stop("'what' must be one of:\n", opts, call. = FALSE)
    }
    if (what == "gspace") {
        obj <- ps@gspace
    } else if (what == "vertex") {
        obj <- ps@vertex
    } else if (what == "vsignal") {
        obj <- ps@vsignal
    } else if (what == "nodes") {
        obj <- ps@nodes
    } else if (what == "edges") {
        obj <- ps@edges
    } else if (what == "gxy") {
        obj <- ps@gxy
    } else if (what == "gxyz") {
        obj <- ps@gxyz
    } else if (what == "pars") {
        obj <- ps@pars
    } else if (what == "misc") {
        obj <- ps@misc
    } else if (what == "status") {
        obj <- ps@status
    } else if (what == "silhouette") {
        obj <- ps@misc$xfloor
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
# show summary information on screen
setMethod("show", "PathwaySpace", function(object) {
    message("A PathwaySpace-class object:\n--status:")
    print(getPathwaySpace(object), quote = FALSE)
})

#-------------------------------------------------------------------------------
#' @title Accessor function for PathwaySpace objects.
#' 
#' @description Get length of a PathwaySpace object.
#' 
#' @param x A \linkS4class{PathwaySpace} class object.
#' @return A non-negative integer of length 1.
#' @examples
#' data('gtoy1', package = 'RGraphSpace')
#' ps <- buildPathwaySpace(gtoy1, nrc = 100)
#' length(ps)
#' 
#' @rdname length-methods
#' @aliases length
#' @export
setMethod("length", "PathwaySpace", function(x) length(x@vertex))

#-------------------------------------------------------------------------------
#' @title Accessor functions for PathwaySpace objects.
#'
#' @description Get and set 'vertex' names of a \linkS4class{PathwaySpace}
#' class object.
#'
#' @param x A \linkS4class{PathwaySpace} class object.
#' @return A character vector.
#' @examples
#' data('gtoy1', package = 'RGraphSpace')
#' ps <- buildPathwaySpace(gtoy1, nrc = 100)
#' names(ps)
#'
#' @import methods
#' @docType methods
#' @rdname names-methods
#' @aliases names
#' @export
#' 
setMethod("names", "PathwaySpace", function(x) x@vertex)

#-------------------------------------------------------------------------------
#' @title Accessor functions for fetching slots from a PathwaySpace object.
#'
#' @description Get or set 'signal' for a \linkS4class{PathwaySpace}
#' class object.
#'
#' @param ps A \linkS4class{PathwaySpace} class object.
#' @param value A numeric vector with values representing signal
#' intensities. This vector should be aligned to the "vertex" slot.
#' @return A numeric vector.
#' @examples
#' data('gtoy1', package = 'RGraphSpace')
#' ps <- buildPathwaySpace(gtoy1, nrc = 100)
#' 
#' # Check vertex names
#' names(ps)
#' 
#' # Check vertex signal
#' vertexSignal(ps)
#' 
#' # Set signal to a given vertex
#' vertexSignal(ps)[1] <- 1
#' 
#' # Set signal to a given vertex
#' vertexSignal(ps)["n3"] <- 1
#' 
#' # Set '1s' to all vertices
#' vertexSignal(ps) <- 1
#' 
#' @import methods
#' @docType methods
#' @rdname vertexSignal-methods
#' @aliases vertexSignal
#' @aliases vertexSignal<-
#' @export
setMethod("vertexSignal", "PathwaySpace", function(ps) ps@vsignal)

#' @rdname vertexSignal-methods
#' @export
setMethod("vertexSignal<-", "PathwaySpace",
    function(ps, value) {
        ps@vsignal[] <- value
        validObject(ps)
        ps <- .update.vsignal(ps)
        return(ps)
    }
)
.update.vsignal <- function(ps) {
    vsignal <- vertexSignal(ps)
    vsignal <- .revise.vertex.signal(vsignal)
    ps@gxy[,"vsignal"] <- vsignal
    ps@nodes[,"vsignal"] <- vsignal
    zscale <- .get.signal.scale(vsignal)
    ps@pars$zscale <- zscale
    return(ps)
}
.revise.vertex.signal <- function(vsignal){
    vsignal[is.nan(vsignal)] <- NA
    vsignal[vsignal == Inf] <- NA
    vsignal[vsignal == -Inf] <- NA
    if (all(is.na(vsignal))) vsignal[] <- 0
    return(vsignal)
}

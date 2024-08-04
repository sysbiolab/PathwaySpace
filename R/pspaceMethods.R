#' @title  Constructor of PathwaySpace-class objects.
#'
#' @description \code{buildPathwaySpace} is a constructor of
#' PathwaySpace-class objects.
#'
#' @param g An \code{igraph} object. It must include graph
#' layout information, with vertex coordinates assigned to \code{x} 
#' and \code{y} vertex attributes. It must also include vertex labels 
#' assigned to the \code{name} vertex attribute.
#' @param nrc A single positive integer indicating the number of rows and 
#' columns (in pixels) for a square image matrix. This argument will 
#' affect the resulting image size and resolution.
#' @param mar A single numeric value (in \code{[0,1]}) indicating the size of
#' the outer margins as a fraction of the image matrix.
#' @param verbose A single logical value specifying to display detailed 
#' messages (when \code{verbose=TRUE}) or not (when \code{verbose=FALSE}).
#' @return A preprocessed \linkS4class{PathwaySpace} class object.
#' @author Vinicius Chagas, Victor Apolonio, Mauro Castro,
#' and TCGA Network.
#' @seealso \code{\link[igraph]{undirected_graph}}
#' @examples
#' # Load a demo igraph
#' data('gtoy1', package = 'PathwaySpace')
#'
#' # Create a new PathwaySpace object
#' pts <- buildPathwaySpace(gtoy1, nrc = 100)
#' # note: adjust 'nrc' to increase image resolution
#' 
#' @importFrom igraph degree vcount ecount which_mutual
#' @importFrom igraph as_edgelist as_adjacency_matrix
#' @importFrom igraph simplify V E 'V<-' 'E<-' is_directed
#' @importFrom stats quantile sd
#' @importFrom scales rescale
#' @importFrom RANN nn2
#' @aliases buildPathwaySpace
#' @export
#'
buildPathwaySpace <- function(g, nrc = 500, mar = 0.075, verbose = TRUE) {
    if(verbose) message("Validating argument types and values...")
    #--- validate argument types
    .validate.args("singleInteger", "nrc", nrc)
    .validate.args("singleNumber", "mar", mar)
    .validate.args("singleLogical", "verbose", verbose)
    #--- validate argument values
    if (nrc < 2) {
        stop("'nrc' should be >=2", call. = FALSE)
    }
    if (mar < 0 || mar > 1) {
        stop("'mar' should be in [0,1]", call. = FALSE)
    }
    #--- validate the igraph object
    if(verbose) message("Validating 'g' object...")
    .validate.igraph(g)
    #--- build PathwaySpace-class
    pts <- .buildPathwaySpace(g, nrc, mar, verbose)
    pts <- .updateStatus(pts, "Preprocess")
    return(pts)
}

#' @title Creating 2D-landscape images from graph objects.
#'
#' @description \code{circularProjection} implements a convolution
#' algorithm to project signal across a 2D-coordinate system.
#'
#' @param pts A \linkS4class{PathwaySpace} class object.
#' @param knn A single positive integer determining the k-nearest signal 
#' sources used in the signal convolution operation.
#' @param pdist A term (in \code{[0,1]}) determining a distance unit for the signal
#' convolution related to the image space. This distance will affect the
#' extent over which the convolution operation projects the signal between
#' source- and destination points.
#' @param rescale A single logical value indicating whether to rescale 
#' the signal. If the signal \code{>=0}, then it will be rescaled to 
#' \code{[0,1]}; if the signal \code{<=0}, then it will be rescaled to 
#' \code{[-1,0]}; and if the signal in \code{(-Inf,+Inf)}, then it will be 
#' rescaled to \code{[-1,1]}.
#' @param verbose A single logical value specifying to display detailed 
#' messages (when \code{verbose=TRUE}) or not (when \code{verbose=FALSE}).
#' @param decay_fun A signal decay function. Available: 'Weibull',
#' 'exponential', and 'linear' functions (see \code{\link{weibullDecay}}).
#' @param ... Additional arguments passed to the decay function.
#' @return A preprocessed \linkS4class{PathwaySpace} class object.
#' @author Vinicius Chagas, Victor Apolonio, Mauro Castro, 
#' and TCGA Network.
#' @seealso \code{\link{buildPathwaySpace}}
#' @examples
#' # Load a demo igraph
#' data('gtoy1', package = 'PathwaySpace')
#'
#' # Create a new PathwaySpace object
#' pts <- buildPathwaySpace(gtoy1, nrc = 100)
#' # note: adjust 'nrc' to increase image resolution
#'
#' # Create a 2D-landscape image
#' pts <- circularProjection(pts)
#'
#' @import methods
#' @docType methods
#' @rdname circularProjection-methods
#' @aliases circularProjection
#' @export
#'
setMethod("circularProjection", "PathwaySpace", function(pts, knn = 8,
    pdist = 0.15, rescale = TRUE, verbose = TRUE, 
    decay_fun = weibullDecay, ...) {
    #--- validate the pipeline status
    if (!.checkStatus(pts, "Preprocess")) {
        stop("NOTE: the 'pts' object needs preprocessing!", call. = FALSE)
    }
    if(verbose) message("Validating argument types and values...")
    #--- validate argument types
    .validate.args("singleInteger", "knn", knn)
    .validate.args("singleNumber", "pdist", pdist)
    .validate.args("singleLogical", "rescale", rescale)
    .validate.args("singleLogical", "verbose", verbose)
    .validate.args("function", "decay_fun", decay_fun)
    #--- validate argument values
    if (knn < 1) {
        stop("'knn' should be >=1", call. = FALSE)
    }
    k <- length(getPathwaySpace(pts, "vertex"))
    if (knn > k) {
        if(verbose)
            message("-- provided knn > k vertices; knn will be adjusted to k.")
        knn <- k
    }
    if (pdist < 0 || pdist > 1) {
        stop("'pdist' should be in [0,1]", call. = FALSE)
    }
    #--- pack args
    pars <- list(knn = knn, pdist = pdist, rescale = rescale, 
        projection="Circular", decay_fun = decay_fun, 
        decay_args = list(...=...))
    for (nm in names(pars)) {
        pts@pars[[nm]] <- pars[[nm]]
    }
    #--- run pts pipeline
    pts <- .circularProjection(pts, verbose, ...=...)
    pts <- .updateStatus(pts, "CircularProjection")
    if (.checkStatus(pts, "PolarProjection")) {
        if(verbose) message("-- polar projection replaced by circular.")
        pts <- .updateStatus(pts, "PolarProjection", FALSE)
    }
    # pts <- .removeSummits(pts, verbose)
    # pts <- .removeSilhouette(pts, verbose)
    return(pts)
})

#' @title Creating 2D-landscape images from graph objects.
#'
#' @description \code{polarProjection} implements a convolution algorithm
#' to project a signal across a 2D-coordinate system.
#'
#' @param pts A \linkS4class{PathwaySpace} class object.
#' @param knn A single positive integer determining the k-nearest signal 
#' sources used in the signal convolution operation.
#' @param pdist A term (in \code{[0,1]}) determining a distance unit for the
#' signal convolution related to length between any two connected vertices. 
#' This distance will affect the extent over which the convolution operation 
#' projects the signal between source- and destination points along the 
#' polar coordinates of the edges.
#' @param rescale A single logical value indicating whether to rescale 
#' the signal. If the signal \code{>=0}, then it will be rescaled to 
#' \code{[0,1]}; if the signal \code{<=0}, then it will be rescaled to 
#' \code{[-1,0]}; and if the signal in \code{(-Inf,+Inf)}, then it will be 
#' rescaled to \code{[-1,1]}.
#' @param theta Angle of projection (degrees in \code{(0,360]}).
#' @param directional If directional edges are available, this argument can 
#' be used to orientate the signal projection on directed graphs.
#' @param verbose A single logical value specifying to display detailed 
#' messages (when \code{verbose=TRUE}) or not (when \code{verbose=FALSE}).
#' @param decay_fun A signal decay function. Available: 'Weibull', 
#' 'exponential', and 'linear' functions (see \code{\link{weibullDecay}}).
#' @param ... Additional arguments passed to the decay function.
#' @return A preprocessed \linkS4class{PathwaySpace} class object.
#' @author Vinicius Chagas, Victor Apolonio, Mauro Castro,
#' and TCGA Network.
#' @seealso \code{\link{buildPathwaySpace}}
#' @examples
#' # Load a demo igraph
#' data('gtoy1', package = 'PathwaySpace')
#'
#' # Create a new PathwaySpace object
#' pts <- buildPathwaySpace(gtoy1, nrc = 100)
#' # note: adjust 'nrc' to increase image resolution
#'
#' # Create a 2D-landscape image
#' pts <- polarProjection(pts)
#'
#' @import methods
#' @docType methods
#' @rdname polarProjection-methods
#' @aliases polarProjection
#' @export
#'
setMethod("polarProjection", "PathwaySpace", function(pts, knn = 8, 
    pdist = 0.5, rescale = TRUE, theta = 180, directional = FALSE, 
    verbose = TRUE, decay_fun = weibullDecay, ...) {
    #--- validate the pipeline status
    if (!.checkStatus(pts, "Preprocess")) {
        stop("NOTE: the 'pts' object needs preprocessing!", call. = FALSE)
    }
    if(verbose) message("Validating argument types and values...")
    .validate.args("singleInteger", "knn", knn)
    .validate.args("singleNumber", "pdist", pdist)
    .validate.args("singleNumber", "theta", theta)
    .validate.args("singleLogical", "rescale", rescale)
    .validate.args("singleLogical", "verbose", verbose)
    .validate.args("singleLogical", "directional", directional)
    .validate.args("function", "decay_fun", decay_fun)
    if (knn < 1) {
        stop("'knn' should be >=1", call. = FALSE)
    }
    k <- length(getPathwaySpace(pts, "vertex"))
    if (knn > k) {
        if(verbose)
            message("-- provided knn > k vertices; knn will be adjusted to k.")
        knn <- k
    }
    if (pdist < 0 || pdist > 1) {
        stop("'pdist' should be in [0,1]", call. = FALSE)
    }
    if (theta <= 0 || theta > 360) {
        msg <- paste0("'polar.theta' should be an angle of projection\n",
            "with degrees in (0,360]")
        stop(msg, call. = FALSE)
    }
    pars <- list(knn = knn, pdist = pdist, theta = theta,  
        rescale = rescale, projection="Polar", directional = directional,
        decay_fun = decay_fun, decay_args = list(...=...))
    for (nm in names(pars)) {
        pts@pars[[nm]] <- pars[[nm]]
    }
    pts <- .polarProjection(pts, verbose, ...=...)
    pts <- .updateStatus(pts, "PolarProjection")
    if (.checkStatus(pts, "CircularProjection")) {
        if(verbose) message("-- circular projection replaced by polar.")
        pts <- .updateStatus(pts, "CircularProjection", FALSE)
    }
    # pts <- .removeSummits(pts, verbose)
    # pts <- .removeSilhouette(pts, verbose)s
    return(pts)
})

#' @title Decorating PathwaySpace images with graph silhouettes.
#'
#' @description \code{silhouetteMapping} constructs an image baseline used
#' to outline the graph layout in a PathwaySpace image.
#'
#' @param pts A \linkS4class{PathwaySpace} class object.
#' @param baseline A fraction (in \code{[0,1]}) of the signal scale of a 
#' PathwaySpace image. This term only affects the image baseline projection, 
#' which represents a silhouette of the graph's layout outlined in the 
#' resulting image. When \code{baseline = 0} (i.e. lower level of the signal 
#' scale), the baseline will extend over the entire image space, so no 
#' silhouette will be visible.
#' @param pdist A term (in \code{0,1}) determining a distance unit for
#' the signal convolution related to the image space. This distance
#' will affect the extent over which the convolution operation
#' projects the image baseline.
#' @param verbose A single logical value specifying to display detailed 
#' messages (when \code{verbose=TRUE}) or not (when \code{verbose=FALSE}).
#' @return A preprocessed \linkS4class{PathwaySpace} class object.
#' @author Vinicius Chagas, Victor Apolonio, Mauro Castro,
#' and TCGA Network.
#' @seealso \code{\link{circularProjection}}
#' @examples
#' # Load a demo igraph
#' data('gtoy1', package = 'PathwaySpace')
#'
#' # Create a new PathwaySpace object
#' pts <- buildPathwaySpace(gtoy1, nrc = 100)
#' # note: adjust 'nrc' to increase image resolution
#'
#' # Create a 2D-landscape image
#' pts <- circularProjection(pts)
#'
#' # Map graph silhouette
#' pts <- silhouetteMapping(pts)
#'
#' @import methods
#' @docType methods
#' @rdname silhouetteMapping-methods
#' @aliases silhouetteMapping
#' @export
#'
setMethod("silhouetteMapping", "PathwaySpace", function(pts, baseline = 0.01,
    pdist = 0.05, verbose = TRUE) {
    #--- validate the pipeline status
    if (!.checkStatus(pts, "Projection")) {
        msg <- paste0("NOTE: the 'pts' object needs to be evaluated\n",
            "by a 'projection' method!")
        stop(msg, call. = FALSE)
    }
    if (.checkStatus(pts, "PolarProjection")) {
        msg <- paste0("NOTE: silhouette decoration not available\n",
            "for polar projection.")
        stop(msg, call. = FALSE)
    }
    if(verbose) message("Validating argument types and values...")
    #--- validate argument types
    .validate.args("singleNumber", "baseline", baseline)
    .validate.args("singleNumber", "pdist", pdist)
    .validate.args("singleLogical", "verbose", verbose)
    #--- validate argument values
    if (baseline < 0 || baseline > 1) {
        stop("'baseline' should be in [0,1]", call. = FALSE)
    }
    if (pdist < 0 || pdist > 1) {
        stop("'pdist' should be in [0,1]", call. = FALSE)
    }
    #--- pack args
    pars <- list(baseline = baseline, pdist = pdist)
    for (nm in names(pars)) {
        pts@pars[[nm]] <- pars[[nm]]
    }
    #--- run pts pipeline
    if(verbose) message("Mapping graph silhouette...")
    pts <- .silhouetteCircular(pts, verbose)
    pts <- .updateStatus(pts, "Silhouette")
    # pts <- .removeSummits(pts, verbose)
    return(pts)
})

#' @title Mapping summits on PathwaySpace images.
#'
#' @description The \code{summitMapping} method implements a segmentation
#' strategy to identify summits on a 2D-landscape image 
#' (see \code{\link{summitWatershed}}).
#'
#' @param pts A \linkS4class{PathwaySpace} class object.
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
#' @author Vinicius Chagas, Victor Apolonio, Mauro Castro, 
#' and TCGA Network.
#' @seealso \code{\link{circularProjection}}
#' @examples
#' # Load a demo igraph
#' data('gtoy1', package = 'PathwaySpace')
#'
#' # Create a new PathwaySpace object
#' pts <- buildPathwaySpace(gtoy1, nrc = 100)
#' # note: adjust 'nrc' to increase image resolution
#'
#' # Create a 2D-landscape image
#' pts <- circularProjection(pts)
#'
#' # Map summits in a 2D-landscape image
#' pts <- summitMapping(pts)
#'
#' @import methods
#' @docType methods
#' @rdname summitMapping-methods
#' @aliases summitMapping
#' @export
#'
setMethod("summitMapping", "PathwaySpace", function(pts, maxset = 30, 
    minsize = 30, threshold = 0.5, verbose = TRUE, 
    segm_fun = summitWatershed, ...) {
    #--- validate the pipeline status
    if (!.checkStatus(pts, "Projection")) {
        msg <- paste0("NOTE: the 'pts' object needs to be\n",
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
        pts@pars[[nm]] <- pars[[nm]]
    }
    #--- run pts pipeline
    pts <- .summitMapping(pts, verbose, ...=...)
    pts <- .updateStatus(pts, "Summits")
    return(pts)
})

#' @title Accessors for fetching slots from a PathwaySpace object.
#'
#' @description \code{getPathwaySpace} retrives information from
#' individual slots available in a PathwaySpace object.
#'
#' @param pts A preprocessed \linkS4class{PathwaySpace} class object
#' @param what A single character value specifying which information should 
#' be retrieved from the slots.
#' Options: 'graph','gxy','gxyz','pars','misc','status','summits',
#' 'summit_mask', and 'summit_contour'.
#' @return Content from slots in the \linkS4class{PathwaySpace} object.
#' @examples
#' # Load a demo igraph
#' data('gtoy1', package = 'PathwaySpace')
#'
#' # Create a new PathwaySpace object
#' pts <- buildPathwaySpace(gtoy1, nrc = 100)
#' # note: adjust 'nrc' to increase image resolution
#'
#' # Get the 'status' slot in pts
#' status <- getPathwaySpace(pts, what = 'status')
#'
#' @import methods
#' @docType methods
#' @rdname getPathwaySpace-methods
#' @aliases getPathwaySpace
#' @export
setMethod("getPathwaySpace", "PathwaySpace", function(pts, what = "status") {
    opts <- c("vertex", "vsignal", "vweight", "edges", 
        "gxy", "gxyz", "pars", "misc", "status", "silhouette",
        "summits", "summit_mask", "summit_contour")
    if (!what %in% opts) {
        opts <- paste0(opts, collapse = ", ")
        stop("'what' must be one of:\n", opts, call. = FALSE)
    }
    if (what == "vertex") {
        obj <- pts@vertex
    } else if (what == "vsignal") {
        obj <- pts@vsignal
    } else if (what == "vweight") {
        obj <- pts@vweight
    } else if (what == "edges") {
        obj <- pts@edges
    } else if (what == "gxy") {
        obj <- pts@gxy
    } else if (what == "gxyz") {
        obj <- pts@gxyz
    } else if (what == "pars") {
        obj <- pts@pars
    } else if (what == "misc") {
        obj <- pts@misc
    } else if (what == "status") {
        obj <- pts@status
    } else if (what == "silhouette") {
        obj <- pts@misc$silhouette
    } else if (what == "summits") {
        obj <- pts@misc$summits$lset
    } else if (what == "summit_mask") {
        obj <- pts@misc$summits$mset
    } else if (what == "summit_contour") {
        obj <- pts@misc$summits$cset
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
#' data('gtoy1', package = 'PathwaySpace')
#' pts <- buildPathwaySpace(gtoy1, nrc = 100)
#' length(pts)
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
#' data('gtoy1', package = 'PathwaySpace')
#' pts <- buildPathwaySpace(gtoy1, nrc = 100)
#' names(pts)
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
#' @param pts A \linkS4class{PathwaySpace} class object.
#' @param value A numeric vector with values representing signal
#' intensities. This vector should be aligned to the "vertex" slot.
#' @return A numeric vector.
#' @examples
#' data('gtoy1', package = 'PathwaySpace')
#' pts <- buildPathwaySpace(gtoy1, nrc = 100)
#' vertexSignal(pts)
#'
#' @import methods
#' @docType methods
#' @rdname vertexSignal-methods
#' @aliases vertexSignal
#' @aliases vertexSignal<-
#' @export
setMethod("vertexSignal", "PathwaySpace", function(pts) pts@vsignal)

#' @rdname vertexSignal-methods
#' @export
setMethod("vertexSignal<-", "PathwaySpace",
    function(pts, value) {
        pts@vsignal[] <- value
        validObject(pts)
        pts <- .update.vsignal(pts)
        return(pts)
    }
)
.update.vsignal <- function(pts) {
    vsignal <- vertexSignal(pts)
    vsignal <- .revise.vertex.signal(vsignal)
    pts@gxy[,"vsignal"] <- vsignal
    zscale <- .getSignalScale(vsignal)
    pts@pars$zscale <- zscale
    return(pts)
}
.revise.vertex.signal <- function(vsignal){
    vsignal[is.nan(vsignal)] <- NA
    vsignal[vsignal == Inf] <- NA
    vsignal[vsignal == -Inf] <- NA
    if (all(is.na(vsignal))) vsignal[] <- 0
    return(vsignal)
}

#-------------------------------------------------------------------------------
#' @title Accessor functions for fetching slots from a PathwaySpace object.
#'
#' @description Get or set 'weights' for a \linkS4class{PathwaySpace}
#' class object.
#'
#' @param pts A \linkS4class{PathwaySpace} class object.
#' @param value A numeric vector with values representing vertex weights.
#' This vector should be aligned to the "vertex" slot.
#' @return A numeric vector.
#' @examples
#' data('gtoy1', package = 'PathwaySpace')
#' pts <- buildPathwaySpace(gtoy1, nrc = 100)
#' vertexWeight(pts)
#'
#' @import methods
#' @docType methods
#' @rdname vertexWeight-methods
#' @aliases vertexWeight
#' @aliases vertexWeight<-
#' @export
setMethod("vertexWeight", "PathwaySpace", function(pts) pts@vweight)

#' @rdname vertexWeight-methods
#' @export
setMethod("vertexWeight<-", "PathwaySpace",
    function(pts, value) {
        pts@vweight[] <- value
        validObject(pts)
        pts <- .update.vweight(pts)
        return(pts)
    }
)
.update.vweight <- function(pts) {
    vweight <- vertexWeight(pts)
    vweight <- .revise.vertex.weight(vweight)
    pts@pars$vwscale <- .get.vwscale(vweight)
    pts@gxy[,"vweight"] <- vweight
    return(pts)
}
.revise.vertex.weight <- function(vweight){
    if (all(is.na(vweight))) vweight[] <- 1
    vweight[is.na(vweight)] <- min(vweight, na.rm = TRUE)
    if (sd(vweight) > 0) {
        vweight <- scales::rescale(vweight, to = c(1, 2))
    } else {
        vweight[] <- 1
    }
    return(vweight)
}





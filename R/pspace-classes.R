
#-------------------------------------------------------------------------------
#' @title Projection Data Container for PathwaySpace Objects
#'
#' @description
#' \code{SpaceProjection} is an S4 class that stores the intermediate and
#' final matrices produced during signal projection in a
#' \linkS4class{PathwaySpace} object. It is created internally by
#' \code{\link{circularProjection}} and is not intended to be constructed
#' directly by the user.
#'
#' @slot coordinates A numeric matrix with one row per graph node and four
#'   columns (\code{X}, \code{Y}, \code{Xint}, \code{Yint}), storing
#'   continuous and integer grid coordinates for each node.
#' @slot floor A numeric matrix of dimensions \code{nrc x nrc} representing
#'   the projection floor, used to mask regions outside the graph silhouette.
#' @slot signal A numeric matrix of dimensions \code{nrc x nrc} storing the
#'   smoothed signal values before final scaling.
#' @slot result A numeric matrix of dimensions \code{nrc x nrc} containing
#'   the final projected signal, scaled to \code{[0, 1]} (or
#'   \code{[-1, 1]} for \code{"negpos"} scale types).
#'
#' @return A \code{SpaceProjection} object.
#'
#' @seealso
#' \linkS4class{PathwaySpace}
#'
#' @aliases SpaceProjection-class
#' @exportClass SpaceProjection
setClass("SpaceProjection",
  slots = c(
    coordinates = "matrix", # nNodes x 4 -- X, Y, Xint, Yint
    floor = "matrix",       # nrc x nrc -- grid floor matrix
    signal = "matrix",      # nrc x nrc -- smoothed signal matrix
    result = "matrix"       # nrc x nrc -- final projection matrix
  ),
  prototype = list(
    coordinates = matrix(nrow = 0, ncol = 4, 
      dimnames = list(NULL, c("X", "Y", "Xint", "Yint"))),
    floor = matrix(nrow = 0, ncol = 0),
    signal = matrix(nrow = 0, ncol = 0),
    result = matrix(nrow = 0, ncol = 0)
  )
)

setValidity("SpaceProjection", function(object) {
  errors <- character()
  
  # floor, signal and projection must share dimensions
  grid_slots <- list(object@floor, object@signal, object@result)
  grid_dims  <- lapply(grid_slots, dim)
  non_empty  <- grid_dims[sapply(grid_dims, function(d) prod(d) > 0)]
  if (length(non_empty) > 1 && length(unique(non_empty)) > 1) {
    errors <- c(errors,
      "Slots '@floor', '@signal', and '@result' must have the same dimensions.")
  }
  
  # coordinates must have 4 columns when non-empty
  if (nrow(object@coordinates) > 0) {
    expected <- c("X", "Y", "Xint", "Yint")
    if (!identical(colnames(object@coordinates), expected)) {
      errors <- c(errors,
        "Slot '@coordinates' must have columns: X, Y, Xint, Yint.")
    }
  }
  
  if (length(errors) == 0) TRUE else errors
})

#-------------------------------------------------------------------------------
#' @title PathwaySpace: An S4 class for signal projection on image spaces
#'
#' @description
#' \code{PathwaySpace} extends the \linkS4class[RGraphSpace]{GraphSpace} class with
#' signal projection slots. It stores projected signal matrices,
#' projection parameters, and workflow status, and is the main object used
#' by the \pkg{PathwaySpace} package.
#' 
#' @slot nodes Inherited from \linkS4class[RGraphSpace]{GraphSpace}.
#' @slot edges Inherited from \linkS4class[RGraphSpace]{GraphSpace}.
#' @slot graph Inherited from \linkS4class[RGraphSpace]{GraphSpace}.
#' @slot image Inherited from \linkS4class[RGraphSpace]{GraphSpace}.
#' @slot fdata Inherited from \linkS4class[RGraphSpace]{GraphSpace}.
#' @slot pars Inherited from \linkS4class[RGraphSpace]{GraphSpace}.
#' @slot misc Inherited from \linkS4class[RGraphSpace]{GraphSpace}.
#' @slot uuid Inherited from \linkS4class[RGraphSpace]{GraphSpace}.
#' @slot projection A \linkS4class{SpaceProjection} object storing the 
#' intermediate and final matrices produced by a projection method.
#' @slot pars_ps A list with PathwaySpace parameters.
#' @slot status A vector containing the processing status of the
#' PathwaySpace object.
#' 
#' @return A \code{PathwaySpace} object.
#' 
#' @section Constructor:
#' See \code{\link{buildPathwaySpace}}.
#'
#' @seealso
#' \linkS4class[RGraphSpace]{GraphSpace},
#' \linkS4class{SpaceProjection},
#' \code{\link{buildPathwaySpace}},
#' \code{\link{circularProjection}}
#' 
#' @author Sysbiolab Team, Mauro Castro (\email{mauro.castro@@ufpr.br})
#' @aliases PathwaySpace-class
#' @importClassesFrom RGraphSpace GraphSpace
#' @exportClass PathwaySpace
setClass("PathwaySpace",
  contains = "GraphSpace",
  slots = c(
    projection = "SpaceProjection",
    pars_ps = "list",
    status = "character"
  ),
  prototype = list(
    projection = new("SpaceProjection"),
    pars_ps = list(),
    status = character()
  )
)
setValidity("PathwaySpace", function(object) {
    
    errors <- character()
    
    ## Check projection slot
    if (!is(object@projection, "SpaceProjection")) {
        errors <- c(errors, "'@projection' slot must be a 'SpaceProjection' object.")
    }
    
    ## Check status slot
    if (!is.character(object@status)) {
        errors <- c(errors, "'@status' slot must be a character vector.")
    }
    
    ## Check pars_ps slot
    if (!is.list(object@pars_ps)) {
        errors <- c(errors, "'@pars_ps' slot must be a list.")
    }
    
    ## Return
    if (length(errors) == 0) TRUE else errors
    
})

#-------------------------------------------------------------------------------
.migrate_ps_pars <- function(ps) {
  if (!inherits(ps, "PathwaySpace")) {
    rlang::abort("'ps' must be a PathwaySpace object.")
  }
  if( !.hasSlot(ps, "projection") || !is(ps@projection, "SpaceProjection") ){
    rlang::warn("Outdated 'PathwaySpace' object: updating on the fly.")
    nrc <- ps@pars$ps$nrc %||% 500
    ps@pars$is.normalized <- FALSE
    gs <- new("GraphSpace",
      nodes = ps@nodes,
      edges = ps@edges,
      graph = ps@graph,
      image = ps@image,
      fdata = ps@fdata,
      pars = ps@pars,
      misc = ps@misc,
      uuid = ps@uuid
    )
    ps <- buildPathwaySpace(gs, nrc = nrc, verbose = FALSE)
  }
  ps
}

#-------------------------------------------------------------------------------
# Ensures all GraphSpace slots (including future additions) are
# carried over correctly during as(gs, "PathwaySpace") coercion.
setReplaceMethod("coerce", c("PathwaySpace", "GraphSpace"), 
    function(from, to = "GraphSpace", value) {
    for (what in slotNames("GraphSpace")) {
        methods::slot(from, what) <- methods::slot(value, what)
    }
    from
})

#-------------------------------------------------------------------------------
# show summary information on screen
setMethod("show", "PathwaySpace", function(object) {
    message("A PathwaySpace-class object for:")
    withCallingHandlers(
        callNextMethod(),
        message = function(m) invokeRestart("muffleMessage")
    )
    cat("+ status:", .summariseStatus(object))
})

#-------------------------------------------------------------------------------
.summariseStatus <- function(ps){
  sts <- c(
    .checkStatus(ps, "Preprocess"),
    .checkStatus(ps, "Projection"),
    .checkStatus(ps, "Silhouette"),
    .checkStatus(ps, "Summits"))
  sts[] <- ifelse(sts, "[x]", "[ ]")
  sts <- paste(names(sts), sts, collapse = "  ", sep="")
  sts
}
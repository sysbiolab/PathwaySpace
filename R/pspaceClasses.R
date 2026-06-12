#' @title PathwaySpace: An S4 class for signal propagation on image spaces
#'
#' @slot nodes A data frame with xy-vertex coordinates.
#' @slot edges A data frame with edges.
#' @slot graph An igraph object.
#' @slot image A raster background image matrix.
#' @slot pars A list inherited GraphSpace parameters.
#' @slot misc A list with intermediate objects for downstream methods.
#' @slot projections A list with processed objects for downstream methods.
#' @slot pars_ps A list with PathwaySpace parameters.
#' @slot status A vector containing the processing status of the
#' PathwaySpace object.
#' @author Sysbiolab Team, Mauro Castro (\email{mauro.castro@@ufpr.br})
#'
#' @method circularProjection \code{\link{circularProjection}}
#' @method polarProjection \code{\link{polarProjection}}
#' @method silhouetteMapping \code{\link{silhouetteMapping}}
#' @method summitMapping \code{\link{summitMapping}}
#' @method getPathwaySpace \code{\link{getPathwaySpace}}
#' @method plotPathwaySpace \code{\link{plotPathwaySpace}}
#' @aliases PathwaySpace-class
#' @importClassesFrom RGraphSpace GraphSpace
#' @return An S4 class object.
#' @section Constructor:
#' see \code{\link{buildPathwaySpace}} constructor.
#' @exportClass PathwaySpace
#' @importClassesFrom RGraphSpace GraphSpace
#' @export
setClass("PathwaySpace",
    contains = "GraphSpace",
    slots = c(
        projections = "list",
        pars_ps = "list",
        status = "character"
    ),
    prototype = list(
        projections = list(),
        pars_ps = list(),
        status = character()
    )
)
setValidity("PathwaySpace", function(object) {
    
    errors <- character()
    
    ## Check projections slot
    if (!is.list(object@projections)) {
        errors <- c(errors, "'projections' must be a list.")
    }
    
    ## Check status slot
    if (!is.character(object@status)) {
        errors <- c(errors, "'status' must be a character vector.")
    }
    
    ## Check pars_ps slot
    if (!is.list(object@pars_ps)) {
        errors <- c(errors, "'pars_ps' must be a list.")
    }
    
    ## Return
    if (length(errors) == 0) TRUE else errors
    
})

#-------------------------------------------------------------------------------
.migrate_ps_pars <- function(ps) {
    if (!inherits(ps, "PathwaySpace")) {
        rlang::abort("'ps' must be a PathwaySpace object.")
    }
    # guard against old objects without pars_ps slot
    if (!.hasSlot(ps, "pars_ps")) {
        rlang::warn("Outdated 'PathwaySpace' object: updating on the fly.")
        ps <- new("PathwaySpace",
            nodes       = ps@nodes,
            edges       = ps@edges,
            graph       = ps@graph,
            image       = ps@image,
            fdata       = ps@fdata,
            pars        = ps@pars,
            misc        = ps@misc,
            uuid        = ps@uuid,
            projections = ps@projections,
            status      = ps@status,
            pars_ps     = if (!is.null(ps@pars$ps)) ps@pars$ps else list()
        )
        ps@pars$ps <- NULL
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


#' @title PathwaySpace: An S4 class for signal propagation on image spaces
#'
#' @slot nodes A data frame with xy-vertex coordinates.
#' @slot edges A data frame with edges.
#' @slot graph An igraph object.
#' @slot image A raster background image matrix.
#' @slot pars A list with parameters.
#' @slot misc A list with intermediate objects for downstream methods.
#' @slot projections A list with processed objects for downstream methods.
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
        status = "character"
    ),
    prototype = list(
        projections = list(),
        status = character()
    )
)
setValidity("PathwaySpace", function(object) {
    errors <- character()
    
    ## Inherit and check GraphSpace validity
    gs_valid <- getValidity(getClass("GraphSpace"))(object)
    if (!isTRUE(gs_valid)) {
        errors <- c(errors, paste0("GraphSpace validity: ", gs_valid))
    }
    
    ## Check projections slot
    if (!is.list(object@projections)) {
        errors <- c(errors, "'projections' must be a list.")
    }
    
    ## Check status slot
    if (!is.character(object@status)) {
        errors <- c(errors, "'status' must be a character vector.")
    }
    
    ## Return
    if (length(errors) == 0) TRUE else errors
    
})


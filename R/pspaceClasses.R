#' @title PathwaySpace: An S4 class for signal propagation on image spaces.
#'
#' @slot gspace A \code{\link[RGraphSpace]{GraphSpace-class}} object.
#' @slot vertex A character vector with vertex names.
#' @slot vsignal A numerical vector with vertex signals.
#' @slot nodes A data frame listing graph nodes.
#' @slot edges A data frame listing graph edges.
#' @slot gxy A data frame with xy-vertex coordinates (numerical).
#' @slot gxyz A numerical matrix with x-cols and y-rows coordinates,
#' and a z-signal.
#' @slot pars A list with parameters.
#' @slot misc A list with intermediate objects for downstream methods.
#' @slot status A vector containing the processing status of the
#' PathwaySpace object.
#' @author Mauro Castro, \email{mauro.castro@@ufpr.br}
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
#'
## Class PathwaySpace
setClass("PathwaySpace",
    slot = c(
        gspace = "GraphSpace",
        vertex = "character",
        vsignal = "numeric",
        nodes = "data.frame",
        edges = "data.frame",
        gxy = "matrix",
        gxyz = "matrix",
        pars = "list",
        misc = "list",
        status = "character"
    ),
    prototype = list(
        gspace = new(Class = "GraphSpace"),
        vertex = character(),
        vsignal = numeric(),
        nodes = data.frame(),
        edges = data.frame(),
        gxy = matrix(nrow = 0, ncol = 2),
        gxyz = matrix(nrow = 0, ncol = 0),
        pars = list(),
        misc = list(),
        status = character()
    )
)
setValidity("PathwaySpace", function(object) {
    vertex <- object@vertex
    vsignal <- object@vsignal
    slot_lengths <- c(length(vertex), length(vsignal))
    if (length(unique(slot_lengths)) != 1){
        msg <- paste0("lengths of slots 'vertex' and 'vsignal' differ.")
        return(msg)
    }
    if (anyDuplicated(vertex) > 0){
        return("slot 'vertex' should contain unique names.")
    }
    if (length(vertex) != nrow(object@gxy)){
        return("names in slots 'vertex' and 'gxy' differ.")
    }
    if ( any(vertex != rownames(object@gxy)) ){
        return("names in slots 'vertex' and 'gxy' differ.")
    }
    TRUE
})


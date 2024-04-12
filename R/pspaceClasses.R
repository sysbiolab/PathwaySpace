#' @title PathwaySpace: An S4 class for signal propagation on image spaces.
#'
#' @slot vertex A character vector with vertex names.
#' @slot vsignal A numerical vector with vertex signals.
#' @slot vweight A numerical vector with vertex weights.
#' @slot edges  A data frame with edges. 
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
#' @return An S4 class object.
#' @section Constructor:
#' see \code{\link{buildPathwaySpace}} constructor.
#' @exportClass PathwaySpace
#'
## Class PathwaySpace
setClass("PathwaySpace",
    slot = c(
        vertex = "character",
        vsignal = "numeric",
        vweight = "numeric",
        edges = "data.frame",
        gxy = "matrix",
        gxyz = "matrix",
        pars = "list",
        misc = "list",
        status = "character"
    ),
    prototype = list(
        vertex = character(),
        vsignal = numeric(),
        vweight = numeric(),
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
    slot_lengths <- c(length(vertex), 
        length(object@vsignal), 
        length(object@vweight))
    if (length(unique(slot_lengths)) != 1){
        msg <- paste0("lengths of slots 'vertex', 'vsignal', ", 
            "and 'vweight' differ.")
        return(msg)
    }
    if (anyDuplicated(vertex)>0){
        return("slot 'vertex' should contain unique names.")
    }
    if (length(vertex) != nrow(object@gxy)){
        return("names in slots 'vertex' and 'gxy' differ.")
    }
    TRUE
})


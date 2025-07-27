#-------------------------------------------------------------------------------
#' @title Weibull decay function
#'
#' @description The \code{weibullDecay} function is used by
#' PathwaySpace's methods for signal convolution and projection.
#'
#' @param x A numeric vector of distances (in [0,1]).
#' @param signal A single numeric value representing a signal.
#' @param decay A decay factor (in [0,1]). This term indicates how much the
#' \code{signal} decreases as a function of distance in pathway space. 
#' For example, at a specific distance defined by the \code{pdist} parameter 
#' (see \code{\link{circularProjection}}), the signal intensity will
#' be the initial signal multiplied by \code{decay}.
#' @param shape A parameter (>=1) of a Weibull function. When \code{shape=1}
#' the Weibull decay follows an exponential decay. When \code{shape>1}
#' the function is first convex, then concave with an inflection point.
#' @return A numeric vector; if missing 'x', it will return decay function.
#' @author Mauro Castro.
#' @seealso \code{\link{expDecay}}, \code{\link{linearDecay}}
#' @examples
#' x <- seq(0, 2, 0.01)
#' y <- weibullDecay(x, signal = 1, decay = 0.5, shape = 4)
#' plot(x, y)
#'
#' @aliases weibullDecay
#' @export
#'
weibullDecay <- function(x, signal = 1, decay = 0.001, shape = 1.05) {
    y <- signal * decay^(x^shape)
    return(y)
}
attributes(weibullDecay)$name <- "weibullDecay"

#-------------------------------------------------------------------------------
#' @title Exponential decay function
#'
#' @description The \code{expDecay} function is used by PathwaySpace's methods
#' for signal convolution and projection.
#'
#' @param x A numeric vector of distances (in [0,1]).
#' @param signal A single numeric value representing a signal.
#' @param decay A decay factor (in [0,1]). This term indicates how much the
#' \code{signal} decreases as a function of distance in pathway space. 
#' For example, at a specific distance defined by the \code{pdist} parameter 
#' (see \code{\link{circularProjection}}), the signal intensity will
#' be the initial signal multiplied by \code{decay}.
#' the \code{\link{weibullDecay}} function.
#' @return A numeric vector; if missing 'x', it will return decay function.
#' @author Mauro Castro.
#' @seealso \code{\link{weibullDecay}}, \code{\link{linearDecay}}
#' @examples
#' x <- seq(0, 2, 0.01)
#' y <- expDecay(x, signal = 1, decay = 0.5)
#' plot(x, y)
#'
#' @aliases expDecay
#' @export
#'
expDecay <- function(x, signal = 1, decay = 0.001) {
    y <- signal * decay^x
    return(y)
}
attributes(expDecay)$name <- "expDecay"

#-------------------------------------------------------------------------------
#' @title A simple linear decay function
#'
#' @description The \code{linearDecay} function is used by PathwaySpace's
#' methods for signal convolution and projection.
#'
#' @param x A numeric vector of distances (in [0,1]).
#' @param signal A single numeric value representing a signal.
#' \code{\link{weibullDecay}} and \code{\link{expDecay}} functions.
#' @return A numeric vector; if missing 'x', it will return decay function.
#' @author Mauro Castro.
#' @seealso \code{\link{weibullDecay}}, \code{\link{expDecay}}
#' @examples
#' x <- seq(0, 2, 0.01)
#' y <- linearDecay(x, signal = 1)
#' plot(x, y)
#'
#' @aliases linearDecay
#' @export
#'
linearDecay <- function(x, signal = 1) {
    y <- signal * (1 - x)
    y[(y * sign(signal)) < 0] <- 0
    return(y)
}
attributes(linearDecay)$name <- "linearDecay"


#-------------------------------------------------------------------------------
#' @title Signal decay functions
#'
#' @description Signal decay functions for \code{\link{PathwaySpace}}
#' internal calls.
#' 
#' @param method A character string specifying a method for
#' signal decay (any of \code{weibull}, \code{exp}, or \code{linear}), 
#' returning \code{\link{weibullDecay}}, \code{\link{expDecay}}, 
#' or \code{\link{linearDecay}} functions, respectively.
#' @param decay A decay factor (in (0,1]) passed to the 
#' \code{\link{weibullDecay}} and \code{\link{expDecay}} functions.
#' @param shape A parameter (>=1) passed to the \code{\link{weibullDecay}} 
#' function.
#' @return An the function.
#' @author Mauro Castro.
#' @seealso \code{\link{circularProjection}} and \code{\link{polarProjection}}
#' @examples
#' signalDecay()
#' 
#' @rdname signalDecay
#' @export
#'
signalDecay <- function(method = c("weibull", "exp", "linear"), 
    decay = 0.001, shape = 1.05){
    method <- match.arg(method)
    .validate.args("singleNumber", "decay", decay)
    .validate.args("singleNumber", "shape", shape)
    if(decay < 0 || decay > 1){
        stop("'decay' must be in [0,1]")
    }
    if(shape < 1){
        stop("'shape' must be >=1")
    }
    if(method=="weibull"){
        if(decay==0) decay <- .Machine$double.xmin
        f <- function(x, signal){
            y <- signal * decay^(x^shape)
            return(y)
        }
        body(f) <- do.call("substitute", list(body(f), 
            list(decay = decay, shape = shape)))
        attributes(f)$name <- "weibullDecay"
        return(f)
    } else if(method=="exp"){
        if(!missing(shape)){
            warning("'shape' ignored unless method is 'weibull'.")
        }
        if(decay==0) decay <- .Machine$double.xmin
        f <- function(x, signal){
            y <- signal * decay^x
            return(y)
        }
        body(f) <- do.call("substitute", list(body(f), list(decay = decay)))
        attributes(f)$name <- "expDecay"
        return(f)
    } else if(method=="linear"){
        if(!missing(decay)){
            warning("'decay' ignored unless method is 'weibull' or 'exp'.")
        }
        if(!missing(shape)){
            warning("'shape' ignored unless method is 'weibull'.")
        }
        f <- function(x, signal){
            y <- signal * (1 - x)
            y[(y * sign(signal)) < 0] <- 0
            return(y)
        }
        attributes(f)$name <- "linearDecay"
        return(f)
    }
}

#-------------------------------------------------------------------------------
#' @title Signal aggregation functions
#'
#' @description Signal aggregation functions for \code{\link{PathwaySpace}}
#' internal calls. The aggregation should be symmetric with respect to signal 
#' polarity, ensuring that opposite signals produce corresponding outputs.
#' 
#' @param method A character string specifying the method for
#' signal aggregation, returning either a customized \code{\link{mean}} or 
#' \code{\link{weighted.mean}} function.
#' @return An aggregation function.
#' @author Mauro Castro.
#' @seealso \code{\link{circularProjection}}, \code{\link{polarProjection}}, 
#' \code{\link{weighted.mean}}
#' @examples
#' signalAggregation()
#' 
#' @importFrom stats weighted.mean
#' @rdname signalAggregation
#' @export
#' 
signalAggregation <- function(method = c("mean", "wmean", "log.wmean", 
    "exp.wmean")){
    method <- match.arg(method)
    if(method=="mean"){
        f <- function(x){
            y <- mean(x, na.rm = TRUE)
            return(y)
        }
        attributes(f)$name <- "meanSignal"
        return(f)
    } else if(method=="wmean"){
        f <- function(x){
            y <- weighted.mean(x, abs(x), na.rm = TRUE)
            return(y)
        }
        attributes(f)$name <- "weightedSignal"
        return(f)
    } else if(method=="log.wmean"){
        f <- function(x){
            y <- weighted.mean(x, log1p(abs(x)), na.rm = TRUE)
            return(y)
        }
        attributes(f)$name <- "logWeightedSignal"
        return(f)
    } else if(method=="exp.wmean"){
        f <- function(x){
            w <- abs(x)
            w <- (w / sum(w, na.rm = TRUE))^2
            y <- weighted.mean(x, w, na.rm = TRUE)
            return(y)
        }
        attributes(f)$name <- "expWeightedSignal"
        return(f)
    }

}
attributes(signalAggregation)$name <- "signalAggregation"


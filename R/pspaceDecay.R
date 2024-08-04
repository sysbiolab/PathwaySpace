#' @title Weibull decay function.
#'
#' @description The \code{weibullDecay} function is used by
#' PathwaySpace's methods for signal convolution and projection.
#'
#' @param x A numeric vector of distances (in [0,1]).
#' @param signal A single numeric value representing a signal.
#' @param decay The rate (in [0,1]) at which the signal decays.
#' This term indicates how much the \code{signal} decreases at a certain
#' distance in \code{x}. At the distance defined by the \code{pdist} term
#' (see \code{\link{circularProjection}}), the signal's value will
#' correspond to the initial signal multiplied by \code{1 - decay}.
#' @param shape A parameter (>=1) of a Weibull function. When \code{shape=1}
#' the Weibull decay follows an exponential decay. When \code{shape>1}
#' the function is first convex, then concave with an inflection point.
#' @return A numeric vector.
#' @author Vinicius Chagas, Victor Apolonio, and
#' Mauro Castro (\email{mauro.castro@@ufpr.br})
#' @seealso \code{\link{expDecay}}, \code{\link{linearDecay}}
#' @examples
#' x <- c(1:100) / 100
#' y <- weibullDecay(x, 1)
#' plot(x, y)
#'
#' @aliases weibullDecay
#' @export
#'
weibullDecay <- function(x, signal, decay = 0.999, shape = 1.05) {
    y <- signal * (1 - decay)^(x^shape)
    return(y)
}
attributes(weibullDecay)$name <- "weibullDecay"

#' @title Exponential decay function.
#'
#' @description The \code{expDecay} function is used by PathwaySpace's methods
#' for signal convolution and projection.
#'
#' @param x A numeric vector of distances (in [0,1]).
#' @param signal A single numeric value representing a signal.
#' @param decay The rate (in [0,1]) at which the signal decays.
#' This term indicates how much the \code{signal} decreases at a certain
#' distance in \code{x}. At the distance defined by the \code{pdist} term
#' (see \code{\link{circularProjection}}), the signal's value will
#' correspond to the initial signal multiplied by \code{1 - decay}.
#' @param ... Not used; argument implemented for call compatibility with
#' the \code{\link{weibullDecay}} function.
#' @return A numeric vector.
#' @author Vinicius Chagas, Victor Apolonio, and
#' Mauro Castro (\email{mauro.castro@@ufpr.br})
#' @seealso \code{\link{weibullDecay}}, \code{\link{linearDecay}}
#' @examples
#' x <- c(1:100) / 100
#' y <- expDecay(x, 1)
#' plot(x, y)
#'
#' @aliases expDecay
#' @export
#'
expDecay <- function(x, signal, decay = 0.999, ...) {
    y <- signal * (1 - decay)^x
    return(y)
}
attributes(expDecay)$name <- "expDecay"

#' @title A simple linear decay function.
#'
#' @description The \code{linearDecay} function is used by PathwaySpace's
#' methods for signal convolution and projection.
#'
#' @param x A numeric vector of distances (in [0,1]).
#' @param signal A single numeric value representing a signal.
#' @param ... Not used; argument implemented for call compatibility with
#' \code{\link{weibullDecay}} and \code{\link{expDecay}} functions.
#' @return A numeric vector.
#' @author Vinicius Chagas, Victor Apolonio, and
#' Mauro Castro (\email{mauro.castro@@ufpr.br})
#' @seealso \code{\link{weibullDecay}}, \code{\link{expDecay}}
#' @examples
#' x <- c(1:100) / 100
#' y <- linearDecay(x, 1)
#' plot(x, y)
#'
#' @aliases linearDecay
#' @export
#'
linearDecay <- function(x, signal, ...) {
    y <- signal * (1 - x)
    y[(y * sign(signal)) < 0] <- 0
    return(y)
}
attributes(linearDecay)$name <- "linearDecay"

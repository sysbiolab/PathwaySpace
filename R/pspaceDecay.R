#-------------------------------------------------------------------------------
#' @title Constructor of Weibull decay functions
#'
#' @description
#' The `weibullDecay()` constructor either creates a decay function or 
#' returns a `ggplot` object for visualizing the decay model. It is a utility 
#' function used internally by \code{\link{circularProjection}} and 
#' \code{\link{polarProjection}}.
#' 
#' @param decay A decay factor (in [0,1]). This term indicates how much a
#' \code{signal} decreases as a function of distance in pathway space. 
#' For example, at a specific distance defined by the \code{pdist} parameter, 
#' the signal intensity will be the initial signal multiplied by \code{decay}.
#' @param pdist A distance normalization term (in (0, 1]) at which the signal 
#' reaches `signal * decay`. This parameter is used to anchor the decay to 
#' a meaningful distance (see `details`). Also, when \code{pdist = 1}, it will
#' represent the diameter of the inscribed circle within the coordinate space
#' of a `PathwaySpace` object.
#' @param shape A parameter (>=1) of a Weibull function. When \code{shape=1}
#' the Weibull decay follows an exponential decay. When \code{shape>1}
#' the function is first convex, then concave with an inflection point.
#' @param plot A logical value indicating whether to return a `ggplot` object.
#' @param demo.signal A numeric value in `[-Inf, Inf]`, only passed when 
#' \code{plot = TRUE} to visualize the decay curve with a specific signal 
#' intensity. The value is ignored by the function constructor, as the decay 
#' function itself is returned without using an initial signal.
#' @return Returns either a function of the form 
#' \code{function(x, signal) \{ ... \}} or, if \code{plot = TRUE}, a `ggplot`
#' object illustrating the decay model.
#' @author Sysbiolab Team
#' @seealso \code{\link{linearDecay}}, \code{\link{expDecay}}
#' @details 
#' The `weibullDecay()` constructor creates a decay model based on the Weibull 
#' distribution. It describes how a signal decreases as a function of distance, 
#' controlled by both a decay rate and a shape parameter.
#' 
#' The decay function is defined as:
#' 
#' \deqn{y = signal \times decay^{\left(\frac{x}{pdist}\right)^{shape}}}
#' 
#' where \eqn{signal} represents the initial intensity, \eqn{decay} controls 
#' the rate of attenuation, \eqn{x} is a vector of normalized distances,
#' and \eqn{shape} adjusts the curvature of the decay. When \eqn{shape = 1}, 
#' the function follows an exponential decay. For \eqn{shape > 1}, the curve 
#' transitions from convex to concave, exhibiting an inflection point. The 
#' \eqn{pdist} parameter anchors the model such that:
#' \itemize{
#'   \item \eqn{y = signal} when \eqn{x = 0}
#'   \item \eqn{y = signal \times decay} when \eqn{x = pdist}
#' }
#' 
#' @examples
#' # Return a decay function
#' decay_fun <- weibullDecay(decay = 0.5, pdist = 0.4, shape = 2)
#'
#' # Plot decay model parameters
#' # weibullDecay(decay = 0.5, pdist = 0.4, shape = 2, plot = TRUE)
#' 
#' @importFrom ggplot2 geom_point geom_hline rel
#' @aliases weibullDecay
#' @export
#' 
weibullDecay <- function(decay = 0.001, pdist = 0.15, shape = 1.05,
  plot = FALSE, demo.signal = 1) {
  
  .validate.ps.args("singleNumber", "decay", decay)
  .validate.ps.args("singleNumber", "pdist", pdist)
  .validate.ps.args("singleNumber", "shape", shape)
  .validate.ps.args("singleNumber", "demo.signal", demo.signal)

  if(decay < 0 || decay > 1){
    stop("'decay' must be in [0,1]", call. = FALSE)
  }
  if(shape < 1){
    stop("'shape' must be >=1", call. = FALSE)
  }
  if(pdist <= 0 || pdist > 1){
    stop("'pdist' must be in (0,1]", call. = FALSE)
  }
  
  if(decay==0) decay <- .Machine$double.xmin
  if(decay==1) decay <- 1 - (1/.Machine$longdouble.max.exp)
  f <- function(x, signal){
    y <- signal * decay^( (x/pdist)^shape )
    return(y)
  }
  body(f) <- do.call("substitute", list(body(f),
    list(decay = decay, shape = shape, pdist = pdist)))
  attributes(f)$name <- "weibullDecay"
  
  if(plot){
    x <- seq(0, 1, 0.02)
    y <- f(x, demo.signal)
    p <- .plot_decay(x, y, demo.signal, decay, pdist, name = "weibullDecay", shape)
    return(p)
  } else {
    if(!missing(demo.signal)){
      warning("The value of 'demo.signal' is ignored by the function constructor.")
    }
    return(f)
  }
  
}

#-------------------------------------------------------------------------------
#' @title Constructor of exponential decay functions
#'
#' @description
#' The `expDecay()` constructor either creates a decay function or 
#' returns a `ggplot` object for visualizing the decay model. It is a utility 
#' function used internally by \code{\link{circularProjection}} and 
#' \code{\link{polarProjection}}.
#' 
#' @param decay A decay factor (in [0,1]). This term indicates how much a
#' \code{signal} decreases as a function of distance in pathway space. 
#' For example, at a specific distance defined by the \code{pdist} parameter, 
#' the signal intensity will be the initial signal multiplied by \code{decay}.
#' @param pdist A distance normalization term (in (0, 1]) at which the signal 
#' reaches `signal * decay`. This parameter is used to anchor the decay to 
#' a meaningful distance (see `details`). Also, when \code{pdist = 1}, it will
#' represent the diameter of the inscribed circle within the coordinate space
#' of a `PathwaySpace` object.
#' @param plot A logical value indicating whether to return a `ggplot` object.
#' @param demo.signal A numeric value in `[-Inf, Inf]`, only passed when 
#' \code{plot = TRUE} to visualize the decay curve with a specific signal 
#' intensity. The value is ignored by the function constructor, as the decay 
#' function itself is returned without using an initial signal.
#' @return Returns either a function of the form 
#' \code{function(x, signal) \{ ... \}} or, if \code{plot = TRUE}, a `ggplot`
#' object illustrating the decay model.
#' @author Sysbiolab Team
#' @seealso \code{\link{linearDecay}}, \code{\link{weibullDecay}}
#' @details 
#' The `expDecay()` constructor creates an exponential decay model. It describes 
#' how a signal decreases as a function of distance, controlled by a decay 
#' rate parameter.
#' 
#' The decay function is defined as:
#' 
#' \deqn{y = signal \times decay^{\left(\frac{x}{pdist}\right)}}
#' 
#' where \eqn{signal} represents the initial intensity, \eqn{decay} controls 
#' the rate of attenuation, and \eqn{x} is a vector of normalized distances.
#' The \eqn{pdist} parameter anchors the model such that:
#' \itemize{
#'   \item \eqn{y = signal} when \eqn{x = 0}
#'   \item \eqn{y = signal \times decay} when \eqn{x = pdist}
#' }
#' @examples
#' # Return a decay function
#' decay_fun <- expDecay(decay = 0.25, pdist = 0.5)
#'
#' # Plot decay model parameters
#' # expDecay(decay = 0.25, pdist = 0.5, plot = TRUE)
#' 
#' @aliases expDecay
#' @export
#'
expDecay <- function(decay = 0.001, pdist = 0.15, plot = FALSE, demo.signal = 1) {
  
  .validate.ps.args("singleNumber", "decay", decay)
  .validate.ps.args("singleNumber", "pdist", pdist)
  .validate.ps.args("singleLogical", "plot", plot)
  .validate.ps.args("singleNumber", "demo.signal", demo.signal)
  
  if(decay < 0 || decay > 1){
    stop("'decay' must be in [0,1]", call. = FALSE)
  }
  if(pdist <= 0 || pdist > 1){
    stop("'pdist' must be in [0,1]", call. = FALSE)
  }
  
  if(decay==0) decay <- .Machine$double.xmin
  if(decay==1) decay <- 1 - (1/.Machine$longdouble.max.exp)
  f <- function(x, signal){
    y <- signal * decay^(x/pdist)
    return(y)
  }
  body(f) <- do.call("substitute", list(body(f), list(decay = decay, 
    pdist = pdist)))
  attributes(f)$name <- "expDecay"
  
  if(plot){
    x <- seq(0, 1, 0.02)
    y <- f(x, demo.signal)
    p <- .plot_decay(x, y, demo.signal, decay, pdist, name = "expDecay")
    return(p) 
  } else {
    if(!missing(demo.signal)){
      warning("The value of 'demo.signal' is ignored by the function constructor.")
    }
    return(f)
  }
  
}

#-------------------------------------------------------------------------------
#' @title Constructor of linear decay functions
#'
#' @description
#' The `linearDecay()` constructor either creates a decay function or 
#' returns a `ggplot` object for visualizing the decay model. It is a utility 
#' function used internally by \code{\link{circularProjection}} and 
#' \code{\link{polarProjection}}.
#' 
#' @param decay A decay factor (in [0,1]). This term indicates how much a
#' \code{signal} decreases as a function of distance in pathway space. 
#' For example, at a specific distance defined by the \code{pdist} parameter, 
#' the signal intensity will be the initial signal multiplied by \code{decay}.
#' @param pdist A distance normalization term (in (0, 1]) at which the signal 
#' reaches `signal * decay`. This parameter is used to anchor the decay to 
#' a meaningful distance (see `details`). Also, when \code{pdist = 1}, it will
#' represent the diameter of the inscribed circle within the coordinate space
#' of a `PathwaySpace` object.
#' @param plot A logical value indicating whether to return a `ggplot` object.
#' @param demo.signal A numeric value in `[-Inf, Inf]`, only passed when 
#' \code{plot = TRUE} to visualize the decay curve with a specific signal 
#' intensity. The value is ignored by the function constructor, as the decay 
#' function itself is returned without using an initial signal.
#' @return Returns either a function of the form 
#' \code{function(x, signal) \{ ... \}} or, if \code{plot = TRUE}, a `ggplot`
#' object illustrating the decay model.
#' @author Sysbiolab Team
#' @seealso \code{\link{expDecay}}, \code{\link{weibullDecay}}
#' @details
#' The `linearDecay()` constructor creates a simple linear decay model. It 
#' describes how a signal decreases proportionally with distance.
#'
#' The decay function is defined as:
#' \deqn{y = signal \times \left(1 - (1 - decay) \times \frac{x}{pdist}\right)}
#'
#' where \eqn{signal} represents the initial intensity, \eqn{decay} defines 
#' the relative signal level at \eqn{pdist}, and \eqn{x} is a vector of 
#' normalized distances. The signal decreases uniformly from its initial 
#' value to \eqn{pdist}, which is a reference distance that anchors 
#' the model such that:
#' \itemize{
#'   \item \eqn{y = signal} when \eqn{x = 0}
#'   \item \eqn{y = signal \times decay} when \eqn{x = pdist}
#' }
#' 
#' This makes the linear form consistent with the exponential and Weibull decay
#' functions, both of which also reach \eqn{signal \times decay} at the
#' reference distance.
#' 
#' @examples
#' # Return a decay function
#' decay_fun <- linearDecay(decay = 0.5, pdist = 0.25)
#'
#' # Plot decay model parameters
#' # linearDecay(decay = 0.5, pdist = 0.25, plot = TRUE)
#'
#' @aliases linearDecay
#' @export
#'
linearDecay <- function(decay = 0.001, pdist = 0.15, plot = FALSE, 
  demo.signal = 1) {
  
  .validate.ps.args("singleNumber", "decay", decay)
  .validate.ps.args("singleNumber", "pdist", pdist)
  .validate.ps.args("singleLogical", "plot", plot)
  .validate.ps.args("singleNumber", "demo.signal", demo.signal)
  
  if(decay < 0 || decay > 1){
    stop("'decay' must be in [0,1]", call. = FALSE)
  }
  if(pdist <= 0 || pdist > 1){
    stop("'pdist' must be in [0,1]", call. = FALSE)
  }
  
  f <- function(x, signal){
    # linear decay
    y <- signal * ( 1 - ( (1 - decay) * x/pdist ) )
    # clip 'y' to prevent flipping sign
    y <- pmax(y * sign(signal), 0) * sign(signal)
    return(y)
  }
  body(f) <- do.call("substitute", list(body(f), 
    list(decay = decay, pdist = pdist)))
  attributes(f)$name <- "linearDecay"
  
  if(plot){
    x <- seq(0, 1, 0.02)
    y <- f(x, demo.signal)
    p <- .plot_decay(x, y, demo.signal, decay, pdist, name = "linearDecay")
    return(p)
  } else {
    if(!missing(demo.signal)){
      warning("The value of 'demo.signal' is ignored by the function constructor.")
    }
    return(f)
  }
  
}

#-------------------------------------------------------------------------------
.plot_decay <- function(x, y, signal, decay, pdist, name, shape){
  if(missing(shape)){
    call <- paste0(name,"(...)\n", 
      "decay = ", format(decay, digits=2), "\n",
      "pdist = ", format(pdist, digits=2)
    )
  } else {
    call <- paste0(name,"(...)\n", 
      "decay = ", format(decay, digits=2), "\n", 
      "pdist = ", format(pdist, digits=2), "\n", 
      "shape = ", format(shape, digits=2)
    )
  }
  si <- ifelse(signal>=0, TRUE, FALSE)
  pi <- ifelse(pdist<0.5, TRUE, FALSE)
  data <- data.frame(x = x, y = y)
  if(signal>=0 && signal <= 1){
    limits <- c(0, 1.025)
    breaks <- pretty(c(0, 1))
  } else if(signal >= -1 && signal < 0){
    limits <- c(-1.025, 0)
    breaks <- pretty(c(-1, 0))
  } else {
    mx <- signal*1.025
    limits <- range(c(0, 1.025 * sign(mx),mx))
    breaks <- pretty(limits)
    breaks <- breaks[abs(breaks)<abs(mx)]
  }
  et1 <- element_text(size=12)
  et2 <- element_text(size=13)
  p <- ggplot(aes(x = x, y = y), data = data) +
    scale_y_continuous(limits = limits, breaks = breaks) + 
    scale_x_continuous(breaks = seq(0, 1, 0.2)) + 
    labs(y="Signal Intensity", x="Normalized Distance") + 
    # geom_vline(xintercept=0, color="grey40") +
    # geom_hline(yintercept=0, color="grey40") +
    annotate("segment", x = 0, xend = pdist, 
      y = decay*signal, yend = decay*signal, 
      color = "green4", linetype = "21", linewidth = 0.75) +
    annotate("segment", x = pdist, xend = pdist, 
      y = decay*signal*ifelse(si,0,1), yend = signal*ifelse(si,decay,1), 
      color = "green4", linetype = "21", linewidth = 0.75) +
    annotate("text", 
      x = 0, y = decay*signal, 
      label = "italic(S)[0] %*% italic(decay)", parse = TRUE,
      color = "green4", vjust = -0.3, hjust = 0, size=4) + 
    annotate("text", 
      x = pdist*ifelse(pi,1.05,0.95), y = signal * ifelse(si,0.05,0.95), 
      label="italic(pdist)", parse = TRUE, color = "green4", vjust = 0,
      hjust = ifelse(pi,0,1), size=4) + 
    annotate("text", 
      x = max(x), y = signal, label=call, 
      vjust = ifelse(si,1,0), hjust = 1, size=4) + 
    geom_point(color = "grey50") + 
    theme(aspect.ratio = 0.85, axis.text = et1, axis.title = et2,
      panel.background = element_rect(fill = "grey95"), 
      panel.grid.major = element_line(linewidth = rel(1.2), colour = "white"),
      panel.grid.minor = element_line(linewidth = rel(0.6), colour = "white")
      )
  
  # Add the new point to the plot
  p <- p + 
    geom_point(data = data.frame(x = 0, y = signal), 
      aes(x = x, y = y), color = "red", size = 6, shape = 18) +
    geom_point(data = data.frame(x = pdist, y = signal * decay),
      aes(x = x, y = y), color = "green3", size = 4, shape = 19) +
    annotate("text", x = 0.02, y = signal, 
      label = "italic(S)[0] * ' (initial signal)'", parse = TRUE,
      color = "red", hjust = 0, vjust = ifelse(si,0,1), size = 4)
  
  return(p)
}

#-------------------------------------------------------------------------------
#' @title Signal aggregation functions
#'
#' @description Signal aggregation functions for \code{\link{circularProjection}}
#' and \code{\link{polarProjection}} internal calls. The aggregation should be 
#' symmetric with respect to signal polarity, ensuring that opposite signals 
#' produce corresponding outputs.
#' 
#' @param method A character string specifying the method for
#' signal aggregation, returning either a customized \code{\link{mean}} or 
#' \code{\link{weighted.mean}} function.
#' @return Returns a function of the form: \code{function(x) { ... }}
#' @author Sysbiolab Team
#' @seealso \code{\link{circularProjection}}, \code{\link{polarProjection}}, 
#' \code{\link{weighted.mean}}
#' @examples
#' aggregate.fun <- signalAggregation()
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


#-------------------------------------------------------------------------------
#' @title Polar transformation functions
#'
#' @description Creates polar transformation functions for 
#' \code{\link{polarProjection}} internal calls. These functions are used to 
#' adjusts signal decay according to point-to-edge angular distances, 
#' with options to attenuate angular shapes.
#' 
#' @param method String indicating the transformation to apply. 
#' Must be one of: "power", "gaussian", or "logistic".
#' @param s Single numeric value in \code{[0, 1]}. Controls the spread around 
#' the \code{x} mean of the Gaussian function.
#' @param k Single numeric value \code{>=1}. Controls the steepness of 
#' the logistic function.
#' @param m Single numeric value in \code{[0, 1]}. Specifies the midpoint of 
#' the logistic function.
#' @return Returns a function of the form: \code{function(x, beta) { ... }}, 
#' that applies the specified shape-based transformation.
#' @author Sysbiolab Team
#' @seealso \code{\link{polarProjection}}
#' @details 
#' The polar transformation controls how much the projected signal decays as 
#' a function of the angular distance between a point in pathway space and 
#' a reference edge axis. The function returned by \code{polarDecay()} expects 
#' two arguments, with the following signature: 
#' \code{function(x, beta) { ... }}.
#' 
#' **Power:**
#' \deqn{x^{\beta}}
#' where \eqn{x} is a vector of normalized angular distances (in \code{[0, 1]})
#' and \eqn{beta} is a non-negative exponent that controls the rate of signal 
#' decay. Increasing \eqn{beta} results in a steeper decay rate, modulating 
#' the angular span of the projection.
#' 
#' **Gaussian:**
#' \deqn{\exp\left(-\frac{(1-x)^2}{2\sigma^2}\right)^{\beta}}
#' where \eqn{sigma} controls the spread around the mean, creating 
#' fuzzier effect on projections.
#' 
#' **Logistic:**
#' \deqn{(1 / (1 + \exp(k (x - m))))^{\beta}}
#' where \eqn{k} is the steepness and \eqn{m} is the function's midpoint, 
#' making more gradual transitions.
#'   
#' These transformations are intended to be plugged into the higher-level 
#' \code{\link{polarProjection}} function, allowing user control over the 
#' polar projection profiles. 
#' 
#' @examples
#' polar.fun <- polarDecay("power")
#' 
#' @md
#' @rdname polarDecay
#' @export
#' 
polarDecay <- function(method = c("power", "gaussian", "logistic"), 
    s = 0.5, k = 10, m = 0.5) {
    
    .validate.ps.args("singleNumber", "s", s)
    .validate.ps.args("singleNumber", "k", k)
    .validate.ps.args("singleNumber", "m", m)
    method <- match.arg(method)
    
    if(s < 0 || s>1)stop("'s' must be in [0,1]", call. = FALSE)
    if(k < 1)stop("'k' must be >=1", call. = FALSE)
    if(m < 0 || m>1)stop("'m' must be in [0,1]", call. = FALSE)
    
    if (method == "power") {
        if(!missing(s)){
            warning("'s' ignored unless method is 'gaussian'.")
        }
        if(!missing(k)){
            warning("'k' ignored unless method is 'logistic'.")
        }
        if(!missing(m)){
            warning("'m' ignored unless method is 'logistic'.")
        }
        f <- function(x, beta){
            y <- x ^ beta
            return(y)
        }
        attributes(f)$name <- "power"
        return(f)
        
    } else if (method == "gaussian") {
        if(!missing(k)){
            warning("'k' ignored unless method is 'logistic'.")
        }
        if(!missing(m)){
            warning("'m' ignored unless method is 'logistic'.")
        }
        f <- function(x, beta){
            y <- exp( - ((1 - x)^2) / (2*s^2) ) ^ beta
            return(y)
        }
        body(f) <- do.call("substitute", list(body(f), list(s = s)))
        attributes(f)$name <- "gaussian"
        return(f)
        
    } else if (method == "logistic") {
        if(!missing(s)){
            warning("'s' ignored unless method is 'gaussian'.")
        }
        f <- function(x, beta){
            y <- ( 1 / (1 + exp(-k * (x - m))) ) ^ beta
            return(y)
        }
        body(f) <- do.call("substitute", list(body(f), list(k = k, m = m)))
        attributes(f)$name <- "logistic"
        return(f)
        
    }
}

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#' @title Deprecated function
#'
#' @description 
#' Use \code{\link{weibullDecay}}, \code{\link{expDecay}}, 
#' and \code{\link{linearDecay}}.
#' 
#' @param ... Deprecated arguments
#' @return Stop unconditionally
#' @author Sysbiolab Team
#' @examples
#' decay.fun <- weibullDecay()
#' 
#' @rdname signalDecay
#' @export
#'
signalDecay <- function(...){
  lifecycle::deprecate_stop("1.0.3", "signalDecay()", "weibullDecay()")
}


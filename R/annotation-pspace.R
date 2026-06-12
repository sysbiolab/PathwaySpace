
#-------------------------------------------------------------------------------
#' @title Annotation Functions for PathwaySpace Plots
#'
#' @description
#' \code{annotation_pspace_signal()} annotates a \code{ggplot}-based
#' \code{PathwaySpace} plot with the projected signal layer, rendered as a
#' raster heatmap bounded to the normalized unit space \code{[0, 1]}.
#'
#' @param ps A \linkS4class{PathwaySpace} object containing a valid signal
#'   projection. Run \code{\link{circularProjection}} before calling this
#'   function.
#' @param colors A character vector of colors used to build the signal
#'   palette. Defaults to \code{\link{pspace.cols}()}.
#' @param bg.color A string specifying the background color, used for
#'   zero-signal or masked regions. Defaults to \code{"grey95"}.
#' @param si.color A single color for silhouette.
#' (see \code{\link{silhouetteMapping}}).
#' @param si.alpha A transparency level in `[0, 1]`, used to adjust the 
#' opacity of the silhouette. This parameter is useful for improving the 
#' perception of a background image, when one is available.
#' @param zlab The title for the 'z' axis of the image signal.
#' @param zlim The 'z' limits of the plot (a numeric vector with two numbers).
#' If NULL, limits are determined from the range of the input values.
#' @param slices An integer specifying the number of discrete color levels
#'   used to quantize the signal. For \code{"negpos"} scale types, this is
#'   rounded up to the nearest even number.
#' @param interpolate A logical value indicating whether to apply linear
#'   interpolation when rendering the raster at a different resolution than
#'   its native size. Defaults to \code{FALSE}.
#'
#' @return A list of \code{ggplot2} layer objects that can be added to a
#'   \code{ggplot()} call with \code{+}.
#'
#' @seealso
#' \code{\link{circularProjection}}
#'
#' @examples
#' data("gtoy1", package = "RGraphSpace")
#' ps <- buildPathwaySpace(gtoy1, nrc = 100)
#' vertexSignal(ps) <- 1
#' ps <- circularProjection(ps)
#' 
#' \dontrun{
#' ggplot(ps) +
#'   annotation_pspace_signal(ps, si.alpha = 0.8) +
#'   theme_gspace_coords(is_norm = TRUE)
#' }
#'
#' @importFrom ggplot2 geom_raster geom_point scale_fill_gradientn
#' @importFrom grDevices as.raster colorRampPalette
#' @importFrom ggnewscale new_scale_fill
#' @rdname annotation_pspace
#' @export
annotation_pspace_signal <- function(ps,
  colors = pspace.cols(), 
  bg.color = "grey95", 
  si.color = "grey85",
  si.alpha = 1,
  zlab = "Density", 
  zlim = NULL, 
  slices = 25,
  interpolate = FALSE) {
  
  if (missing(ps) || !inherits(ps, "PathwaySpace")) {
    rlang::abort("'ps' must be a PathwaySpace object.")
  }
  .validate.colors("allColors","colors", colors)
  .validate.colors("singleColor", "si.color", si.color)
  .validate.ps.args("singleNumber", "si.alpha", si.alpha)
  .validate.ps.args("singleString", "zlab", zlab)
  .validate.ps.args("singleInteger", "slices", slices)
  .validate.ps.args("singleLogical", "interpolate", interpolate)
  if(!is.na(bg.color) ){
    .validate.colors("singleColor", "bg.color", bg.color)
  }
  if(!is.null(zlim)) {
    .validate.ps.args("numeric_vec", "zlim", zlim)
    if(length(zlim)!=2) 
      stop("'zlim' should be a numeric vector of lenght 2.", call. = FALSE)
  }
  
  #--- get slots from ps
  gxyz <- getPathwaySpace(ps, "projections")$gxyz
  pars_ps <- getPathwaySpace(ps, "projections")$pars_ps
  
  if (is.null(gxyz)) {
    rlang::abort(c(
      "x" = "No projections found in 'ps'.",
      "i" = "Run `circularProjection()` before plotting."
    ))
  }
  
  #--- set colors
  if(pars_ps$configs$scale.type=="negpos"){
    slices <- ceiling(slices/2) * 2
  }
  colors <- colorRampPalette(colors)(slices)
  
  # set zlim
  if(is.null(zlim)){
    zlim <- pars_ps$configs$zlim
  } else {
    gxyz[gxyz < zlim[1]] <- zlim[1]
    gxyz[gxyz > zlim[2]] <- zlim[2]
  }
  if (all(zlim == 0)) zlim[2] <- 1
  
  #--- trim colors and set zlim args
  cl <- .trimcols(colors, bg.color, zlim, pars_ps)
  cl <- .set_zlim_args(cl, zlim)
  
  #--- slice the signal matrix
  bks <- seq(zlim[1], zlim[2], length.out = slices)
  gxyz[,] <- bks[cut(as.numeric(gxyz), breaks = sort(unique(bks)),
    include.lowest = TRUE)]
  
  #--- melt the signal matrix
  gxyz <- data.frame(arrayInd(seq_along(gxyz), .dim = dim(gxyz)), 
    as.numeric(gxyz))
  colnames(gxyz) <- c("Y", "X", "Z")
  
  #--- rescale matrix indices to plot coordinates
  gxyz$X <- scales::rescale(gxyz$X)
  gxyz$Y <- scales::rescale(gxyz$Y)
  
  #--- apply opacity effect, scaling alpha to z
  if(si.alpha < 1){
    si.color <- grDevices::adjustcolor(si.color, si.alpha)
    si.alpha <- .scale_alpha(si.alpha, gxyz, zlim, pars_ps)
  } else {
    si.alpha <- 1
  }
  
  #--- main projection
  X <- Y <- Z <- NULL
  list(
  ggplot2::geom_raster(data = gxyz, 
    mapping = ggplot2::aes(X, Y, fill = Z), interpolate = interpolate, 
    na.rm=TRUE, alpha = si.alpha),
  ggplot2::scale_fill_gradientn(name = zlab, limits = cl$zlim,
    breaks = cl$breaks, labels = names(cl$breaks), colours = cl$pal, 
    aesthetics = "fill", na.value = si.color),
    ggnewscale::new_scale_fill()
  )
  
}

#-------------------------------------------------------------------------------
.scale_alpha <- function(si.alpha, gxyz, zlim, pars_ps){
  az <- gxyz$Z
  mxz <- max(abs(zlim))
  if(pars_ps$configs$scale.type == "negpos") {
    slim <- 0.5 * (1 - si.alpha)
    slim <- slim * mxz
    az[az < 0 & az < -slim] <- mxz
    az[az > 0 & az > slim] <- mxz
    az <- abs(az)
  } else if(pars_ps$configs$scale.type == "neg") {
    az[az < zlim[2]] <- zlim[1]
  } else {
    az[az > zlim[1]] <- zlim[2]
  }
  az <- az/mxz
  alpha <- az^10 + si.alpha
  return(alpha)
}

#-------------------------------------------------------------------------------
.set_zlim_args <- function(cl, zlim){
  # adjust labels for z-axis middle and tips
  bks_names <- cl$breaks
  bks_names <- format(bks_names,  trim = TRUE)
  n <- length(bks_names)
  bks_names[!seq_len(n) %in% c(1, ceiling(n/2), n)] <- ""
  names(cl$breaks) <- bks_names
  # expand 'zlim' and palette tips
  expand <- TRUE
  if(expand){
    tips <- (zlim[2] - zlim[1]) * 0.1
    cl$zlim <- c(zlim[1] - tips, zlim[2] + tips)
    cl$pal <- c(cl$pal[1], cl$pal, cl$pal[length(cl$pal)])
  } else {
    cl$zlim <- zlim
  }
  return(cl)
}

#-------------------------------------------------------------------------------
.trimcols <- function(colors, bg.color, zlim, pars_ps) {
  if(pars_ps$configs$scale.type == "negpos") {
    
    if (is.na(bg.color)) {
      if (length(colors) %% 2 == 1){
        bg.color <- colors[(length(colors)+1)/2]
      } else {
        bg.color <- colors[(length(colors)/2)+1]
      } 
    } else {
      if (length(colors) %% 2 == 1){
        colors[(length(colors)+1)/2] <- bg.color
      } else {
        n <- length(colors)/2
        colors <- c(colors[seq_len(n)], bg.color, 
          colors[(n+1):(n*2)])
      } 
    }
    cols <- colorRampPalette(colors)(17)
    bg <- cols[9]
    
  } else if(pars_ps$configs$scale.type=="neg") {
    
    if (is.na(bg.color)) {
      bg.color <- colors[length(colors)]
      colors <- colors[-length(colors)]
    }
    cols <- colorRampPalette(c(colors,bg.color))(16)
    bg <- cols[length(cols)]
    cols <- cols[-length(cols)]
    
  } else {
    
    if (is.na(bg.color)) {
      bg.color <- colors[1]
      colors <- colors[-1]
    }
    cols <- colorRampPalette(c(bg.color, colors))(16)
    bg <- cols[1]
    cols <- cols[-1]
    
  }
  
  bkIn <- seq(zlim[1], zlim[2], length.out = length(cols))
  bkOut <- pretty(zlim, n = 11)
  pal <- .trimRamp(bkIn, bkOut, cols)
  
  return(list(breaks = bkOut, pal = pal, bg = bg))
}

#-------------------------------------------------------------------------------
.trimRamp <- function(bkIn, bkOut, cols) {
  cols <- t(col2rgb(cols)/255)
  bkOut <- ifelse(bkOut < bkIn[1], bkIn[1], bkOut)
  bkOut <- ifelse(bkOut > bkIn[length(bkIn)], bkIn[length(bkIn)], bkOut)
  bnout <- .bincode(bkOut, bkIn, right = TRUE, include.lowest = TRUE)
  rcol <- lapply(unique(bnout), function(i){
    j <- bnout == i
    .getcolor(bkOut[j], bkIn[i], bkIn[i+1], cols[i, ], cols[i+1, ])
  })
  rcol <- unlist(rcol)
  return(rcol)
}
.getcolor <- function (x, bk1, bk2, c1, c2) {
  c1 <- grDevices::convertColor(c1, "sRGB", "Lab")
  c2 <- grDevices::convertColor(c2, "sRGB", "Lab")
  rcol <- matrix(ncol = 3, nrow = length(x))
  for (i in seq_len(3)) {
    xx <- (x - bk2) * (c2[i] - c1[i]) / (bk2 - bk1) + c2[i]
    rcol[, i] <- xx
  }
  rcol <- grDevices::convertColor(rcol, "Lab", "sRGB")
  .xlim <- function(x) {
    x[x < 0] <- 0; x[x > 1] <- 1
    return(x)
  }
  rcol[, ] <- .xlim(as.numeric(rcol))
  .rgb2hex <- function(r, g, b){
    grDevices::rgb(r, g, b, maxColorValue = 1)
  }
  rcol <- .rgb2hex(rcol[, 1], rcol[, 2], rcol[, 3])
  return(rcol)
}

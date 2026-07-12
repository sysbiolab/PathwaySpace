#' @title Plotting 2D-landscape images for the PathwaySpace package
#'
#' @description \code{plotPathwaySpace} is a wrapper function to 
#' create dedicated ggplot graphics for PathwaySpace-class objects.
#'
#' @param ps A \linkS4class{PathwaySpace} class object.
#' @param colors A vector of colors.
#' @param bg.color A single color for background.
#' @param si.color A single color for silhouette.
#' (see \code{\link{silhouetteMapping}}).
#' @param si.alpha A transparency level in `[0, 1]`, used to adjust the 
#' opacity of the silhouette. This parameter is useful for improving the 
#' perception of a background image, when one is available.
#' @param theme Name of a custom PathwaySpace theme. These themes 
#' (from 'th0' to 'th3') consist mainly of preconfigured ggplot settings, 
#' which the user can subsequently refine using \code{\link[ggplot2]{ggplot2}}.
#' @param title A string for the title.
#' @param xlab The title for the 'x' axis of a 2D-image space.
#' @param ylab The title for the 'y' axis of a 2D-image space.
#' @param zlab The title for the 'z' axis of the image signal.
#' @param font.size A single numeric value passed to plot annotations.
#' @param font.color A single color passed to plot annotations.
#' @param zlim The 'z' limits of the plot (a numeric vector with two numbers).
#' If NULL, limits are determined from the range of the input values.
#' @param slices A single positive integer value used to split 
#' the image signal into equally-spaced intervals.
#' @param add.image A logical value indicating whether to add a background 
#' image, when one is available (see \code{\link[RGraphSpace]{GraphSpace}}).
#' @param add.grid A logical value indicating whether to add gridlines to 
#' the image space. However, gridlines will only appear when the image 
#' is decorated with graph silhouettes (see \code{\link{silhouetteMapping}}).
#' @param grid.color A color passed to \code{\link[ggplot2]{geom_point}}.
#' @param marks A vector of vertex names to be highlighted in the 
#' image space. This argument overrides 'add.labels'.
#' @param mark.size A size argument passed to \code{\link[ggplot2]{geom_text}}.
#' @param mark.color A color passed to \code{\link[ggplot2]{geom_text}}.
#' @param mark.padding A box padding argument passed to 
#' \code{\link[ggrepel]{geom_text_repel}}.
#' @param mark.line.width A line width argument passed to 
#' \code{\link[ggrepel]{geom_text_repel}}.
#' @param use.dotmark A logical value indicating whether "marks" should be 
#' represented as dots.
#' @param add.summits A logical value indicating whether to add contour 
#' lines to 'summits' (when summits are available; 
#' see \code{\link{summitMapping}}).
#' @param label.summits A logical value indicating whether to label summits.
#' @param summit.color A color passed to 'summits'.
#' @param add.marks Deprecated. Use \code{marks} instead.
#' @return A ggplot-class object.
#' @author Sysbiolab Team, Mauro Castro.
#' @seealso \code{\link{circularProjection}}, \code{\link{polarProjection}}
#' @examples
#' # Load a demo igraph
#' data('gtoy1', package = 'RGraphSpace')
#'
#' # # Check graph validity
#' gs <- GraphSpace(gtoy1)
#' 
#' gs <- normalizeGraphSpace(gs)
#' 
#' # Create a PathwaySpace object
#' ps <- buildPathwaySpace(gs, nrc = 300)
#' # note: adjust 'nrc' to increase image resolution
#'
#' # Set '1s' as vertex signal
#' vertexSignal(ps) <- 1
#'
#' # Create a 2D-landscape image
#' ps <- circularProjection(ps, k = 2,
#'    decay.fun = weibullDecay(pdist = 0.4))
#'
#' # Plot a 2D-landscape image
#' plotPathwaySpace(ps)
#' 
#' @import methods
#' @docType methods
#' @importFrom ggplot2 ggplot annotate element_text theme theme_bw theme_gray
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous expansion
#' @importFrom ggplot2 aes scale_fill_gradientn annotation_raster
#' @importFrom ggplot2 coord_fixed geom_raster labs
#' @importFrom ggplot2 margin element_blank element_rect element_line
#' @importFrom grid unit
#' @importFrom ggrepel geom_text_repel
#' @importFrom grDevices convertColor col2rgb rgb as.raster
#' @importFrom grDevices adjustcolor colorRampPalette
#' @importFrom stats runif
#' @importFrom RGraphSpace plotGraphSpace theme_gspace_coords gs_image
#' @rdname plotPathwaySpace-methods
#' @aliases plotPathwaySpace
#' @export
#'
setMethod("plotPathwaySpace", "PathwaySpace", 
  function(ps, 
    title = activeFeature(ps), 
    colors = pspace.cols(), bg.color = "grey95", 
    si.color = "grey85", si.alpha = 1,
    theme = "th0",
    xlab = "Graph coordinates 1", 
    ylab = "Graph coordinates 2", 
    zlab = "Density", 
    font.size = 1, font.color = "white",
    zlim = NULL, slices = 25, add.image = FALSE,
    add.grid = TRUE, grid.color = "white", 
    marks = FALSE, mark.size = 3, mark.color = "white", 
    mark.padding = 0.5, mark.line.width = 0.5, use.dotmark = FALSE, 
    add.summits = TRUE, label.summits = TRUE, summit.color = "white",
    add.marks = deprecated()) {
    
    if (lifecycle::is_present(add.marks)) {
      lifecycle::deprecate_soft("1.4.2", "plotPathwaySpace(add.marks)",
        details = "Use `plotPathwaySpace(marks)` instead.")
      marks <- add.marks
    }
    
    .check_updated_ps(ps)
    
    #--- validate the ps object and args
    if (!.checkStatus(ps, "Projection") && !.checkStatus(ps, "Silhouette")) {
      rlang::abort(c(
        "The 'ps' object has not been evaluated by a 'projection' method.",
        "i" = "Run a projection method on 'ps' before calling this function."
      ))
    }
    .validate.ps.args("singleNumber", "si.alpha", si.alpha)
    .validate.ps.args("singleString", "xlab", xlab)
    .validate.ps.args("singleString", "ylab", ylab)
    .validate.ps.args("singleString", "zlab", zlab)
    .validate.ps.args("singleNumber", "font.size", font.size)
    .validate.ps.args("singleInteger", "slices", slices)
    .validate.ps.args("singleLogical", "add.grid", add.grid)
    .validate.ps.args("singleLogical", "add.summits", add.summits)
    .validate.ps.args("singleLogical", "label.summits", label.summits)
    .validate.ps.args("singleNumber", "mark.size", mark.size)
    .validate.ps.args("singleNumber", "mark.padding", mark.padding)
    .validate.ps.args("singleNumber","mark.line.width", mark.line.width)
    .validate.ps.args("singleLogical", "use.dotmark", use.dotmark)
    .validate.ps.args("singleLogical", "add.image", add.image)
    .validate.colors("singleColor", "mark.color", mark.color)
    .validate.colors("allColors","colors", colors)
    .validate.colors("singleColor", "si.color", si.color)
    .validate.colors("singleColor", "font.color", font.color)
    .validate.colors("singleColor", "grid.color", grid.color)
    .validate.colors("singleColor", "summit.color", summit.color)
    if(!is.na(bg.color) ){
      .validate.colors("singleColor", "bg.color", bg.color)
    }
    if(!is.null(title)){
      .validate.ps.args("singleString", "title", title)
    }
    theme <- match.arg(theme, choices = c("th0", "th1", "th2", "th3"))
    if(!is.null(zlim)) {
      .validate.ps.args("numeric_vec", "zlim", zlim)
      if(length(zlim)!=2) 
        rlang::abort("'zlim' should be a numeric vector of lenght 2.")
      if (zlim[1] == zlim[2])
        rlang::abort("'zlim' must have two distinct values.")
      zlim <- sort(zlim)
    }
    if (si.alpha < 0 || si.alpha > 1) {
      rlang::abort("'si.alpha' should be in [0,1]")
    }
    
    #--- get slots from ps
    silstatus <- .checkStatus(ps, "Silhouette")
    summits <- getPathwaySpace(ps, "summits")
    cset <- getPathwaySpace(ps, "summit_contour")
    projection <- getPathwaySpace(ps, "projection")
    pars_gs <- getGraphSpace(ps, "pars")
    pars_ps <- getPathwaySpace(ps, "pars")
    gxy <- projection@coordinates
    gxyz <- projection@result

    #--- set colors
    if(pars_ps$configs$scale.type=="negpos"){
      slices <- ceiling(slices/2) * 2
    }
    colors <- colorRampPalette(colors)(slices)
    
    # set scales
    if(is.null(zlim)){
      zlim <- pars_ps$configs$zlim %||% pars_ps$configs$scaling
    } else {
      gxyz[gxyz < zlim[1]] <- zlim[1]
      gxyz[gxyz > zlim[2]] <- zlim[2]
    }
    
    #--- set gspace theme
    gs_theme <- theme_gspace_coords(theme = theme, 
      is_norm = TRUE, xlab = xlab, ylab = ylab, 
      txt_size = font.size, leg_size = font.size, 
      bg_color = "grey95")
    
    gs_pars <- attributes(gs_theme)$gspace_pars
    
    #--- trim colors and set zlim args
    cl <- .trimcols(colors, bg.color, zlim, pars_ps)
    cl <- .set_zlim_args(cl, zlim)
    bks <- seq(zlim[1], zlim[2], length.out = slices)
    gxyz[, ] <- bks[cut(as.numeric(gxyz), breaks = sort(unique(bks)),
      include.lowest = TRUE)]
    
    #--- get grid
    gridln <- .getGrid(gxyz, gs_pars$axis.ticks)
    gridln <- as.numeric(gridln)
    gridln <- data.frame(
      arrayInd(seq_along(gridln), .dim = dim(gxyz)), gridln)
    colnames(gridln) <- c("Y", "X", "L")
    
    #--- set main input for ggplot
    gxyz <- data.frame(arrayInd(seq_along(gxyz), .dim = dim(gxyz)),
      as.numeric(gxyz))
    colnames(gxyz) <- c("Y", "X", "Z")
    
    #--- scale coordinates to plot space
    gxyz$X <- scales::rescale(gxyz$X)
    gxyz$Y <- scales::rescale(gxyz$Y)
    gxyz$L <- gridln$L
    
    #--- set a bg color effect, scaling alpha to z
    if(si.alpha < 1){
      si.color <- adjustcolor(si.color, si.alpha)
      si.alpha <- .scale_alpha(si.alpha, gxyz, zlim, pars_ps)
    }
    
    #--- initialize ggplot
    ggp <- .set_pspace(gxyz, zlab, cl, si.color) + gs_theme
    
    #--- add image
    if(pars_gs$image.space %||% FALSE){
      img <- gs_image(ps)
      if(add.image){
        ggp <- .add_image(ggp, img)
      } else {
        ggi <- .add_image(ggp, img)
        if(add.grid) ggi <- .add_grid(ggi, gxyz, grid.color)
      }
    }
    ggp <- ggp + ggplot2::labs(fill = zlab)
      
    #--- add main projection
    ggp <- ggp + ggplot2::geom_raster(interpolate = FALSE, 
      na.rm=TRUE, alpha = si.alpha)
    
    #--- add a grid
    if(add.grid) ggp <- .add_grid(ggp, gxyz, grid.color)
    
    #--- add contour lines if available
    bl <- add.summits || label.summits
    if (bl && !is.null(summits) && length(summits) > 0 && sum(cset) > 0) {
      ggp <- .add_contour(ggp, gxyz, summits, cset, 
        summit.color, mark.size, add.summits, label.summits)
    }
    
    #--- add marks if available
    if(isTRUE(marks)){
      ggp <- .add_marks(ggp, gxy, pars_ps, marks = rownames(gxy), 
        mark.size, mark.color, mark.padding, mark.line.width, use.dotmark)
    } else if(is.character(marks)){
      .validate.ps.args("allCharacter", "marks", marks)
      ggp <- .add_marks(ggp, gxy, pars_ps, marks = marks, mark.size,
        mark.color, mark.padding, mark.line.width, use.dotmark)
    }
    
    #--- add annotations
    if(.checkStatus(ps, "Projection")){
      ggp <- .custom_annotations(ggp, title, pars_ps, font.size, 
        font.color, silstatus, si.color)
    }
    
    return(ggp)
    
  }
)

#-------------------------------------------------------------------------------
.custom_annotations <- function(ggp, title, pars_ps, font.size, 
  font.color, silstatus, si.color){
  if(silstatus){
    if(si.color=="grey85"){
      fcol <- "grey20"
    } else {
      fcol <- font.color
    }
    sep <- "; "
    xlab <- 0.01
    hjust <- 0
  } else {
    fcol <- font.color
    sep <- "\n"
    xlab <- 0.99
    hjust <- 1
  }
  if(!is.null(title)){
    ggp <- ggp + ggplot2::annotate("text", label = title,
      colour = fcol, size = font.size*4, x = 0, y = 0.99, 
      hjust = 0, vjust = 1)
  }
  dfun <- pars_ps$decay$fun
  if(!is.null(dfun)){
    if(dfun == "weibullDecay"){
      dfun <- "Weibull decay"
    } else if(dfun == "expDecay"){
      dfun <- "Exponential decay"
    } else if(dfun == "linearDecay"){
      dfun <- "Linear decay"
    } else {
      dfun <- "Custom decay"
    }
    pars_ps$dfun <- dfun
  } else {
    pars_ps$dfun <- "Custom decay"
  }
  if(pars_ps$projection=="Polar"){
    annot <- pars_ps[c("projection", "dfun", "k", "beta")]
    annot$k <- paste0("k = ", annot$k, "; ")
    annot$beta <- paste0("beta = ", pars_ps$beta)
  } else {
    annot <- pars_ps[c("projection", "dfun", "k")]
    annot$k <- paste0("k = ", annot$k)
  }
  annot$projection <- paste0(annot$projection, " projection", sep)
  annot$dfun <- paste0(annot$dfun, sep)
  annot <- paste(unlist(annot), collapse = "")
  ggp <- ggp + ggplot2::annotate("text", label = annot,
    colour = fcol, size = font.size*3, x = xlab, 
    y = 0.01, hjust = hjust, vjust = 0)
  return(ggp)
}

#-------------------------------------------------------------------------------
.set_pspace <- function(gxyz, zlab, cl, si.color){
  X <- Y <- Z <- NULL
  ggp <- ggplot2::ggplot(gxyz, ggplot2::aes(X, Y, fill = Z)) +
    ggplot2::scale_fill_gradientn(limits = cl$zlim,
      breaks = cl$breaks, labels = names(cl$breaks),
      colours = cl$pal, aesthetics = "fill", na.value = si.color) +
    ggplot2::coord_fixed()
  return(ggp)
}

#-------------------------------------------------------------------------------
.add_grid <- function(ggp, gxyz, grid.color){
  dt <- gxyz[gxyz$L == 1, c("X", "Y")]
  ggp <- ggp + ggplot2::annotate(geom = "point", x = dt$X, 
    y = dt$Y, color = grid.color, size = 0.2, pch = 15)
  return(ggp)
}

#-------------------------------------------------------------------------------
.add_image <- function(ggp, image){
  ggp <- ggp + annotation_raster(image, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
  return(ggp)
}

#-------------------------------------------------------------------------------
.add_contour <- function(ggp, gxyz, summits, cset, summit.color, 
  mark.size, add.summits, label.summits) {
  setnames <- names(summits)
  if (is.null(setnames)) {
    setnames <- as.character(seq_along(summits))
  }
  cset <- data.frame(which(!is.na(cset) | is.na(cset), arr.ind = TRUE),
    as.numeric(cset))
  colnames(cset) <- c("Y", "X", "C")
  gxyz <- cbind(gxyz, C = cset$C)
  xy.tx <- NULL
  concav <- sort(unique(gxyz$C))[-1]
  for (id in concav) {
    xy.cv <- gxyz[gxyz$C == id, c("X", "Y")]
    xy.tx <- rbind(xy.tx, colMeans(xy.cv))
    if(add.summits){
      ggp <- ggp + ggplot2::annotate(geom = "tile", x = xy.cv[, 1], 
        y = xy.cv[, 2], color = summit.color,
        fill = summit.color, linewidth = 0.2)
    }
  }
  rownames(xy.tx) <- setnames
  xy.tx <- as.data.frame(xy.tx)
  if(label.summits){
    if(add.summits){
      ggp <- ggp + ggplot2::annotate(geom = "text", x = xy.tx[, 1],
        y = xy.tx[, 2], label = rownames(xy.tx), vjust = 0.5, hjust = 0.5,
        color = summit.color, size = mark.size,
        fontface = "bold")
    } else {
      nudgex <- sign(xy.tx$X - 0.5) * xy.tx$X
      nudgey <- sign(xy.tx$Y - 0.5) * xy.tx$Y
      ggp <- ggp + ggrepel::geom_text_repel(
        mapping = aes(label = rownames(xy.tx),segment.size = 0.4), 
        xlim = c(0.1, 0.9), ylim = c(0.1, 0.9),
        data = xy.tx, min.segment.length = 0.1,
        fontface = "bold", force = 3, segment.linetype = 1,
        max.overlaps = nrow(xy.tx) + 5, point.padding = 0, seed = 123,
        max.iter = 20000, max.time = 30, nudge_x = nudgex * 0.075,
        nudge_y = nudgey * 0.075, size = mark.size, colour = summit.color,
        segment.colour = summit.color)
    }
  }
  return(ggp)
}

#-------------------------------------------------------------------------------
.add_marks <- function(ggp, gxy, pars_ps, marks, mark.size,
  mark.color, mark.padding, mark.line.width, use.dotmark) {
  
  # scale coordinates to plot space
  gxy_df <- as.data.frame(gxy)
  gxy_df <- gxy_df[, c("X", "Y")]
  gxy_df$X <- scales::rescale(gxy_df$X, from = c(1, pars_ps$nrc))
  gxy_df$Y <- scales::rescale(gxy_df$Y, from = c(1, pars_ps$nrc))
  
  # match  marks
  marks <- marks[!duplicated(marks)]
  if (is.null(names(marks))) names(marks) <- marks
  names(marks) <- ifelse(is.na(names(marks)), marks, names(marks))
  names(marks) <- ifelse(names(marks) == "", marks, names(marks))
  idx_df <- .get_mark_idx(marks, gxy_df)
  if(any(is.na(idx_df$idx))){
    rlang::abort("All 'marks' should be annotated in the 'PathwaySpace' object.")
  }
  
  # set df
  gxy_df <- gxy_df[idx_df$idx, , drop = FALSE]
  gxy_df$ID <- idx_df$label
  ID <- NULL
  
  # call to ggplot
  if (use.dotmark) {
    ggp <- ggp + ggplot2::annotate("point", shape = 5, 
      colour = mark.color, size = mark.size * 0.4, stroke = 0.4, 
      x = gxy_df$X, y = gxy_df$Y)
  }
  
  nudgex <- sign(gxy_df$X - 0.5) * gxy_df$X * 0.1
  nudgey <- sign(gxy_df$Y - 0.5) * gxy_df$Y * 0.1
  nudgex <- pmax(nudgex,-0.1)
  nudgey <- pmax(nudgey,-0.1)
  nudgex <- pmin(nudgex, 0.1)
  nudgey <- pmin(nudgey, 0.1)
  
  ggp <- ggp + ggrepel::geom_text_repel(
    mapping = aes(label = ID, segment.size = mark.line.width), 
    xlim = c(0.1, 0.9), ylim = c(0.1, 0.9),
    data = gxy_df, min.segment.length = 0.1,
    fontface = "bold", force = 3, segment.linetype = "2121", 
    max.overlaps = nrow(gxy_df) + 5, point.padding = 0, seed = 123, 
    max.iter = 20000, max.time = 30, nudge_x = nudgex, 
    nudge_y = nudgey, size = mark.size, colour = mark.color, 
    segment.colour = mark.color, box.padding = mark.padding)
  
  return(ggp)
}
.get_mark_idx <- function(marks, gxy_df){
  idx_df <- data.frame(name=marks, label=names(marks))
  idx_df$idx1 <- match(marks, rownames(gxy_df))
  idx_df$idx2 <- match(names(marks), rownames(gxy_df))
  idx_df$idx <- sapply(seq_len(nrow(idx_df)), function(i){
    idx <- idx_df[i, c(3,4)]
    idx[!is.na(idx)][1]
  })
  idx_df <- idx_df[,c(1,2,5)]
  return(idx_df)
}

#-------------------------------------------------------------------------------
#--- grid lines
.getGrid <- function(gxyz, ticks = c(0.2, 0.4, 0.6, 0.8), ndots = 100) {
  ticks <- ticks[ticks>0 & ticks <1]
  nc <- ncol(gxyz)
  nr <- nrow(gxyz)
  by <- ceiling(nc/ndots)
  ic <- ceiling(ticks * nc)
  ir <- ceiling(ticks * nr)
  grid1 <- grid2 <- array(0, dim = dim(gxyz))
  grid1[ir, seq(1, ncol(grid1), by = by)] <- 1L
  grid2[seq(1, nrow(grid2), by = by), ic] <- 1L
  grid1 <- (grid1 + grid2) > 0
  grid1[grid1] <- 1
  grid1[!is.na(gxyz)] <- 0
  return(grid1)
}


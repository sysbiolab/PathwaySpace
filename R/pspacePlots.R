#' @title Plotting 2D-landscape images for the PathwaySpace package.
#'
#' @description \code{plotPathwaySpace} is a wrapper function to 
#' create dedicated ggplot graphics for PathwaySpace-class objects.
#'
#' @param pts A \linkS4class{PathwaySpace} class object.
#' @param colors A vector of colors. Each color is a specific tone used to 
#' create a customized color palette that matches the signal type projected 
#' in the image space, such as binary or numeric values in (-Inf, +Inf). The 
#' palette will be created with 5-color tones, interpolated according to the
#' provided sequence in the vector of colors. The proportion of each color 
#' tone can be adjusted by the 'trim.colors' argument. The primary objective 
#' of this color palette construction is to enhance the visibility of summits 
#' and valleys in the image space. If you wish to bypass the automatic palette 
#' generation and use the 'colors' input as-is, simply set 'trim.colors' 
#' to NULL.
#' @param trim.colors An vector with 5 positive integer numbers. This argument
#' can be used to adjust the proportion of each color tone in the palette.
#' @param bg.color A single color for background.
#' @param theme.name Name of a custom PathwaySpace theme. These themes 
#' (from 'th0' to 'th3') consist mainly of preconfigured ggplot settings, 
#' which the user can subsequently fine-tune within the resulting 
#' ggplot object.
#' @param title A string for the title.
#' @param font.size A single numeric value passed to ggplot themes.
#' @param font.color A single color passed to ggplot themes.
#' @param xlab The title for the 'x' axis of a 2D-image space.
#' @param ylab The title for the 'y' axis of a 2D-image space.
#' @param zlab The title for the 'z' axis of the image signal.
#' @param zlim The 'z' limits of the plot (a numeric vector with two numbers).
#' If NULL, limits are determined from the range of the input values.
#' @param slices A single positive integer value used to split 
#' the image signal into equally-spaced intervals.
#' @param add.grid A logical value indicating whether to add gridlines to 
#' the image space. However, gridlines will only appear when the image 
#' is decorated with graph silhouettes (see \code{\link{silhouetteMapping}}).
#' @param grid.color A color passed to \code{\link[ggplot2]{geom_point}}.
#' @param add.contour A logical value indicating whether to add contour 
#' lines to 'summits' (see \code{\link{summitMapping}}).
#' @param contour.color A color passed to \code{\link[ggplot2]{geom_tile}}.
#' @param label.summits A logical value indicating whether to label summits,
#' (when summits are available).
#' @param marks A logical value indicating whether to add 'marks' to vertex 
#' positions. Alternatively, this could be a vector listing vertex names.
#' @param mark.size A font size argument passed to 
#' \code{\link[ggplot2]{geom_text}}.
#' @param mark.color A color passed to \code{\link[ggrepel]{geom_text_repel}}.
#' @param mark.padding A box padding argument passed to 
#' \code{\link[ggrepel]{geom_text_repel}}.
#' @param mark.line.width A line width argument passed to 
#' \code{\link[ggrepel]{geom_text_repel}}.
#' @param use.dotmark A logical value indicating whether "marks" should be 
#' represented as dots.
#' @return A ggplot-class object.
#' @author Vinicius Chagas, Victor Apolonio, Mauro Castro, 
#' and TCGA Network.
#' @seealso \code{\link{circularProjection}}
#' @examples
#' # Load a demo igraph
#' data('gtoy1', package = 'PathwaySpace')
#'
#' # Create a new PathwaySpace object
#' pts <- buildPathwaySpace(gtoy1, nrc = 100)
#' # note: adjust 'nrc' to increase image resolution
#'
#' # Create a 2D-landscape image
#' pts <- circularProjection(pts)
#'
#' # Map graph silhouette
#' pts <- silhouetteMapping(pts)
#'
#' # Plot a 2D-landscape image
#' plotPathwaySpace(pts)
#' 
#' @import methods
#' @docType methods
#' @importFrom ggplot2 ggplot annotate element_text theme theme_bw theme_gray
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous expansion
#' @importFrom ggplot2 aes scale_fill_gradientn
#' @importFrom ggplot2 coord_fixed geom_raster
#' @importFrom ggplot2 margin element_blank element_rect element_line
#' @importFrom ggrepel geom_text_repel
#' @importFrom grDevices convertColor col2rgb rgb
#' @importFrom grDevices adjustcolor colorRampPalette
#' @importFrom RGraphSpace plotGraphSpace
#' @rdname plotPathwaySpace-methods
#' @aliases plotPathwaySpace
#' @export
#'
setMethod("plotPathwaySpace", "PathwaySpace", 
    function(pts, colors = pspace.cols(), trim.colors = c(3, 2, 1, 2, 3), 
        bg.color = "grey85", theme.name = c("th0", "th1", "th2", "th3"),
        title = "PathwaySpace", font.size = 1, font.color = "white",
        xlab = "Pathway coordinates 1", ylab = "Pathway coordinates 2", 
        zlab = "Density", zlim = NULL, slices = 25, add.grid = TRUE, 
        grid.color = "white", add.contour = TRUE, contour.color = "white",
        label.summits = TRUE, marks = FALSE, mark.size = 3, 
        mark.color = "white", mark.padding = 0.5, mark.line.width = 0.5, 
        use.dotmark = FALSE) {
        #--- validate the pts object and args
        if (!.checkStatus(pts, "Projection")) {
            stop("NOTE: 'pts' needs to be evaluated by a 'projection' method!",
                call. = FALSE)
        }
        theme.name <- match.arg(theme.name)
        .validate.args("singleString", "title", title)
        .validate.args("singleNumber", "font.size", font.size)
        .validate.args("singleString", "xlab", xlab)
        .validate.args("singleString", "ylab", ylab)
        .validate.args("singleString", "zlab", zlab)
        .validate.args("singleInteger", "slices", slices)
        .validate.args("singleLogical", "add.grid", add.grid)
        .validate.args("singleLogical", "add.contour", add.contour)
        .validate.args("singleLogical", "label.summits", label.summits)
        .validate.args("singleNumber", "mark.padding", mark.padding)
        .validate.args("singleLogical", "use.dotmark", use.dotmark)
        .validate.args("numeric_vec","mark.size", mark.size)
        .validate.args("singleNumber","mark.line.width", mark.line.width)
        .validate.plot.args("trim.colors", trim.colors)
        .validate.plot.args("bg.color", bg.color)
        .validate.plot.args("marks", marks)
        .validate.colors("singleColor", "font.color", font.color)
        .validate.colors("singleColor", "grid.color", grid.color)
        .validate.colors("singleColor", "contour.color", contour.color)
        .validate.colors("allColors","mark.color", mark.color)
        .validate.colors("allColors","colors", colors)
        if(!is.null(zlim)) {
          .validate.args("numeric_vec", "zlim", zlim)
          if(length(zlim)!=2) 
            stop("'zlim' should be a numeric vector of lenght 2.")
        }
        #--- get slots from pts
        gxy <- getPathwaySpace(pts, "gxy")
        gxyz <- getPathwaySpace(pts, "gxyz")
        pars <- getPathwaySpace(pts, "pars")
        summits <- getPathwaySpace(pts, "summits")
        cset <- getPathwaySpace(pts, "summit_contour")
        silstatus <- .checkStatus(pts, "Silhouette")
        ldim <- dim(gxyz)
        #--- set colors
        if(!is.null(trim.colors)){
            colors <- .pspacePalette(colors, trim.colors)
        }
        if(pars$zscale$scale.type=="negpos"){
          slices <- ceiling(slices/2) * 2
        }
        # set scale
        if(is.null(zlim)){
            if(pars$rescale){
              zlim <- pars$zscale$scaling
            } else {
              # zlim <- pars$zscale$range
              mx <- pars$zscale$maxsig
              if(pars$zscale$scale.type=="negpos"){
                zlim <- c(-mx, mx)
              } else if(pars$zscale$scale.type=="neg"){
                zlim <- c(-mx, 0)
              } else {
                zlim <- c(0, mx)
              }
            }
        }
        gxyz[gxyz < zlim[1]] <- zlim[1]
        gxyz[gxyz > zlim[2]] <- zlim[2]
        if (all(zlim == 0)) zlim[2] <- 1
        #--- trim colors and set theme args
        cl <- .trimcols(colors, bg.color, zlim, pars)
        cl <- .set.theme.bks(theme.name, cl)
        cl <- .set.theme.zlim(cl, zlim)
        # slice gxyz image
        bks <- seq(zlim[1], zlim[2], length.out = slices)
        gxyz[, ] <- bks[cut(as.numeric(gxyz), breaks = sort(unique(bks)),
            include.lowest = TRUE)]
        #--- get grid
        gridln <- .getGrid(gxyz, cl$axis.ticks)
        # gridln <- gridln[nrow(gridln):1, ]
        gridln <- as.numeric(gridln)
        gridln <- data.frame(
            arrayInd(seq_along(gridln), .dim = ldim), gridln)
        colnames(gridln) <- c("Y", "X", "L")
        #--- set main input for ggplot
        gxyz <- as.numeric(gxyz)
        gxyz <- data.frame(arrayInd(seq_along(gxyz), .dim = ldim), gxyz)
        colnames(gxyz) <- c("Y", "X", "Z")
        #--- set xy coordinates
        rx <- range(gxyz$X)
        ry <- range(gxyz$Y)
        gxyz$X <- scales::rescale(gxyz$X)
        gxyz$Y <- scales::rescale(gxyz$Y)
        gxyz$L <- gridln$L
        #--- get ggplot object
        ggp <- .get.ggplot(gxyz, xlab, ylab, zlab, cl, add.grid, grid.color)
        #--- add contour lines if available
        bl <- add.contour || label.summits
        if (bl && !is.null(summits) && length(summits) > 0 && sum(cset) > 0) {
            ggp <- .add.contour(ggp, gxyz, summits, cset, 
              contour.color, mark.size, add.contour, label.summits)
        }
        #--- add marks if available
        bl <- is.logical(marks) && marks
        if (bl || is.character(marks)) {
            if(bl) marks <- rownames(gxy)
            ggp <- .add.marks(ggp, gxy, rx, ry, marks, mark.color,
                mark.size, mark.padding, mark.line.width, use.dotmark)
        }
        #--- add annotations
        ggp <- .custom.annotations(ggp, title, pars, font.size, 
            font.color, silstatus)
        #--- apply custom theme
        ggp <- .custom.themes(ggp, theme.name, font.size, bg.color)
        return(ggp)
    }
)

#-------------------------------------------------------------------------------
.custom.annotations <- function(ggp, title, pars, font.size, 
    font.color, silstatus){
    ggp <- ggp + ggplot2::annotate("text", label = title,
        colour = font.color, size = font.size*4.5, x = 0, y = 0.99, 
        hjust = 0, vjust = 1)
    dfun <- attributes(pars$decay_fun)$name
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
        pars$dfun <- dfun
    }
    if(silstatus){
        sep <- "; "
        fcol <- "grey20"
        xlab <- 0.01
        hjust <- 0
    } else {
        sep <- "\n"
        fcol <- font.color
        xlab <- 0.99
        hjust <- 1
    }
    if(pars$projection=="Polar"){
        annot <- pars[c("projection", "dfun", "kns", "theta")]
        annot$kns <- paste0("kns = ", annot$kns, "; ")
        annot$theta <- paste0("theta = ", pars$theta)
    } else {
        annot <- pars[c("projection", "dfun", "kns")]
        annot$kns <- paste0("kns = ", annot$kns)
    }
    annot$projection <- paste0(annot$projection, " projection", sep)
    annot$dfun <- paste0(annot$dfun, sep)
    annot <- paste(unlist(annot), collapse = "")
    ggp <- ggp + ggplot2::annotate("text", label = annot,
        colour = fcol, size = font.size*3.5, x = xlab, 
        y = 0.01, hjust = hjust, vjust = 0)
    return(ggp)
}

#-------------------------------------------------------------------------------
.get.ggplot <- function(gxyz, xlab, ylab, zlab, cl, add.grid, grid.color){
    X <- Y <- Z <- NULL
    ggp <- ggplot2::ggplot(gxyz, ggplot2::aes(X, Y, fill = Z)) +
        ggplot2::scale_x_continuous(name = xlab, breaks = cl$axis.ticks,
            labels = format(cl$axis.ticks), position = cl$x.position,
            limits = cl$xylim, expand = ggplot2::expansion(mult = 0)) +
        ggplot2::scale_y_continuous(name = ylab, breaks = cl$axis.ticks,
            labels = format(cl$axis.ticks), limits = cl$xylim,
            expand = ggplot2::expansion(mult = 0)) +
        ggplot2::scale_fill_gradientn(name = zlab, limits = cl$zlim,
            breaks = cl$breaks, labels = names(cl$breaks),
            colours = cl$pal, aesthetics = "fill", na.value = cl$bg) +
        ggplot2::coord_fixed() + ggplot2::geom_raster(interpolate = FALSE)
    if (add.grid) {
        dt <- gxyz[gxyz$L == 1, c("X", "Y")]
        ggp <- ggp + ggplot2::annotate(geom = "point", x = dt$X, 
            y = dt$Y, color = grid.color, size = 0.2, pch = 15)
    }
    return(ggp)
}

#-------------------------------------------------------------------------------
.add.contour <- function(ggp, gxyz, summits, cset, contour.color, 
    mark.size, add.contour, label.summits) {
    setnames <- names(summits)
    if (is.null(setnames)) {
        setnames <- seq_along(setnames)
    }
    cset <- data.frame(which(!is.na(cset) | is.na(cset), arr.ind = TRUE),
        as.numeric(cset))
    colnames(cset) <- c("Y", "X", "C")
    gxyz <- cbind(gxyz, C = cset$C)
    xy.tx <- NULL
    concav <- sort(unique(gxyz$C))[-1]
    for (i in seq_along(concav)) {
        xy.cv <- gxyz[gxyz$C == i, c("X", "Y")]
        xy.tx <- rbind(xy.tx, colMeans(xy.cv))
        if(add.contour){
          ggp <- ggp + ggplot2::annotate(geom = "tile", x = xy.cv[, 1], 
            y = xy.cv[, 2], color = contour.color,
            fill = contour.color, linewidth = 0.2)
        }
    }
    rownames(xy.tx) <- setnames
    xy.tx <- as.data.frame(xy.tx)
    if(label.summits){
      if(add.contour){
        ggp <- ggp + ggplot2::annotate(geom = "text", x = xy.tx[, 1],
          y = xy.tx[, 2], label = rownames(xy.tx), vjust = 0.5, hjust = 0.5,
          color = contour.color, size = mark.size,
          fontface = "bold")
      } else {
        nudgex <- sign(xy.tx$X - 0.5) * xy.tx$X
        nudgey <- sign(xy.tx$Y - 0.5) * xy.tx$Y
        ggp <- ggp + ggrepel::geom_text_repel(mapping = aes(label = rownames(xy.tx),
          segment.size = 0.4), data = xy.tx, min.segment.length = 0.1,
          fontface = "bold", force = 3, segment.linetype = 1,
          max.overlaps = nrow(xy.tx) + 5, point.padding = 0, seed = 123,
          max.iter = 20000, max.time = 30, nudge_x = nudgex * 0.075,
          nudge_y = nudgey * 0.075, size = mark.size, colour = contour.color,
          segment.colour = contour.color)
      }
    }
    return(ggp)
}

#-------------------------------------------------------------------------------
.add.marks <- function(ggp, gxy, rx, ry, marks, mark.color,
    mark.size, mark.padding, mark.line.width, use.dotmark) {
    gxy <- as.data.frame(gxy)
    if (is.null(names(marks))) names(marks) <- marks
    names(marks) <- ifelse(names(marks) == "", marks, names(marks))
    marks <- marks[marks %in% rownames(gxy)]
    if (length(mark.color) > 1) {
        if (is.null(names(mark.color))) {
            stop("'mark.color' should be named for this call!", call. = FALSE)
        }
        if (!all(marks %in% names(mark.color))) {
            stop("All 'marks' should be listed in 'mark.color' names!",
                call. = FALSE)
        }
        mark.color <- mark.color[marks]
    }
    if (length(mark.size) > 1) {
        if (is.null(names(mark.size))) {
            stop("'mark.size' should be named for this call!", call. = FALSE)
        }
        if (!all(marks %in% names(mark.size))) {
            stop("All 'marks' should be listed in 'mark.size' names!",
                call. = FALSE)
        }
        mark.size <- mark.size[marks]
    }
    gxy$X <- scales::rescale(gxy$X, from = rx)
    gxy$Y <- scales::rescale(gxy$Y, from = ry)
    gxy <- gxy[, c("X", "Y")]
    gxy$mark.color <- mark.color
    gxy$mark.size <- mark.size
    gxy <- gxy[marks, , drop = FALSE]
    gxy$ID <- names(marks); ID <- NULL
    if (use.dotmark) {
        ggp <- ggp + ggplot2::annotate("point", shape = 5, 
            colour = mark.color, size = mark.size * 0.4, stroke = 0.4, 
            x = gxy$X, y = gxy$Y)
    }
    nudgex <- sign(gxy$X - 0.5) * gxy$X
    nudgey <- sign(gxy$Y - 0.5) * gxy$Y
    ggp <- ggp + ggrepel::geom_text_repel(mapping = aes(label = ID,
        segment.size = mark.line.width), data = gxy, min.segment.length = 0.1,
        fontface = "bold", force = 3, segment.linetype = "2121", 
        max.overlaps = nrow(gxy) + 5, point.padding = 0, seed = 123, 
        max.iter = 20000, max.time = 30, nudge_x = nudgex * 0.15, 
        nudge_y = nudgey * 0.15, size = mark.size, colour = mark.color, 
        segment.colour = mark.color, box.padding = mark.padding)
    return(ggp)
}

#-------------------------------------------------------------------------------
.custom.themes <- function(gg, theme.name, font.size, bg.color) {
    et1 <- ggplot2::element_text(size = 14 * font.size)
    et2 <- ggplot2::element_text(size = 12 * font.size)
    if (theme.name == "th0") {
        gg <- .custom.th0(gg, font.size, bg.color)
    } else if (theme.name == "th1") {
        gg <- .custom.th1(gg, font.size, bg.color)
    } else if (theme.name == "th2") {
        gg <- .custom.th2(gg, font.size, bg.color)
    } else {
        gg <- .custom.th3(gg, font.size, bg.color)
    }
    return(gg)
}
.custom.th0 <- function(gg, font.size, bg.color) {
    et1 <- ggplot2::element_text(size = 14 * font.size)
    et2 <- ggplot2::element_text(size = 12 * font.size)
    gg <- gg + ggplot2::theme(axis.title = et1, axis.text = et2,
        legend.title = et2, legend.text = et2,
        panel.background = element_rect(fill = bg.color))
    return(gg)
}
.custom.th1 <- function(gg, font.size,
    bg.color) {
    et1 <- ggplot2::element_text(size = 14 * font.size)
    et2 <- ggplot2::element_text(size = 12 * font.size)
    gg <- gg + ggplot2::theme_bw() +
        ggplot2::theme(axis.title = et1,
            axis.text = et2, legend.title = et2,
            legend.text = et2, legend.margin = margin(0, 0, 0, 0), 
            plot.margin = margin(1, 1, 1, 1), 
            legend.background = element_blank(),
            legend.box.background = element_blank(),
            plot.background = element_blank(),
            panel.grid.minor = element_line(linewidth = 0.7, 
                colour = bg.color),
            panel.grid.major = element_line(linewidth = 0.7,
                colour = bg.color),
            axis.ticks = element_line(linewidth = 0.7),
            axis.line = element_blank(),
            panel.border = element_rect(linewidth = 1.2))
    return(gg)
}
.custom.th2 <- function(gg, font.size, bg.color) {
    et1 <- ggplot2::element_text(size = 14 * font.size)
    et2 <- ggplot2::element_text(size = 12 * font.size)
    gg <- gg + ggplot2::theme_gray() + ggplot2::theme(axis.title = et1,
        axis.text = et2, legend.title = et2,
        legend.text = et2, legend.margin = margin(0, 0, 0, 0), 
        plot.margin = margin(5, 10, 0, 10), 
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        plot.background = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks = element_line(linewidth = 0.7),
        axis.line = element_blank(), panel.border = element_blank(),
        panel.background = element_rect(fill = bg.color))
    return(gg)
}
.custom.th3 <- function(gg, font.size, bg.color) {
    et1 <- ggplot2::element_text(size = 14 * font.size)
    et2 <- ggplot2::element_text(size = 12 * font.size, hjust=0.5)
    gg <- gg + ggplot2::theme_gray() + 
      ggplot2::theme(axis.title = et1, axis.text = et2, 
        legend.title = element_text(size = 12 * font.size, vjust = 1), 
        legend.text = et2,
        legend.margin = margin(0, 0, 0, 0),
        legend.position = "bottom", plot.margin = margin(5, 5, 5, 5), 
        legend.box.margin = margin(0, 0, 0, 0), 
        # legend.text.align = 0.5,
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        plot.background = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks = element_line(linewidth = 0.5),
        axis.line = element_blank(), panel.border = element_blank(),
        panel.background = element_rect(fill = bg.color))
    return(gg)
}
.set.theme.bks <- function(theme.name, cl=list()){
    if (theme.name %in% c("th2")) {
        cl$axis.ticks <- seq(0.1, 0.9, 0.2)
        cl$xylim <- c(-0.01, 1.01)
        cl$x.position <- "bottom"
        cl$justify <- "right"
    } else if (theme.name %in% c("th3")) {
        cl$axis.ticks <- c(0.25, 0.5, 0.75)
        cl$xylim <- c(-0.01, 1.01)
        cl$x.position <- "top"
        cl$justify <- "centre"
    } else {
        cl$axis.ticks <- seq(0, 1, 0.2)
        cl$xylim <- c(-0.05, 1.05)
        cl$x.position <- "bottom"
        cl$justify <- "right"
    }
    return(cl)
}
.set.theme.zlim <- function(cl, zlim){
    # adjust labels for z-axis midle and tips
    bks_names <- cl$breaks
    bks_names <- format(bks_names,  trim = TRUE)
    n <- length(bks_names)
    bks_names[!seq_len(n) %in% c(1, ceiling(n/2), n)] <- ""
    # bks_names <- format(bks_names, justify=cl$justify)
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
#--- get grid lines
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

#-------------------------------------------------------------------------------
.pspacePalette <- function(colors, trim.colors) {
    if(length(colors) != 5){
        colors <- colorRampPalette(colors)(5)
    } 
    tms <- trim.colors * 3
    offset <- list()
    offset[[1]] <- c(0.03, 0.09, 0.07, 0)
    offset[[2]] <- c(0.12, 0.26, 0.08, 0)
    offset[[3]] <- c(0.25, 0.11, 0.09, 0)
    offset[[4]] <- c(0.00, 0.31, 0.10, 0)
    offset[[5]] <- c(0.31, 0.30, 0.09, 0)
    cols <- lapply(seq_along(colors), function(i){
        cl <- adjustcolor(colors[i], offset = offset[[i]])
        colorRampPalette(c(colors[i], cl))(tms[i])
    })
    cols[[4]] <- rev(cols[[4]])
    cols[[5]] <- rev(cols[[5]])
    cols <- unlist(cols)
    cols <- colorRampPalette(cols)(25)
    return(cols)
}

#-------------------------------------------------------------------------------
.trimcols <- function(colors, bg.color, zlim, pars) {
    if(pars$zscale$scale.type == "neg") {
        if (is.null(bg.color)) {
            bg.color <- colors[length(colors)]
            colors <- colors[-length(colors)]
        }
        cols <- colorRampPalette(c(colors,bg.color))(16)
        bg <- cols[length(cols)]
        cols <- cols[-length(cols)]
    } else if(pars$zscale$scale.type=="negpos") {
        if (is.null(bg.color)) {
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
    } else {
        if (is.null(bg.color)) {
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

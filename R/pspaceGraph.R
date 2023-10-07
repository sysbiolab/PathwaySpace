#' @title Plotting graphs for the PathwaySpace package.
#'
#' @description \code{plotGraphSpace} is a wrapper function to 
#' create dedicated ggplot graphics for PathwaySpace-class objects.
#'
#' @param pts A \linkS4class{PathwaySpace} class object.
#' @param node.size  A single numeric value for node size. 
#' @param edge.width  A single numeric value for edge width. 
#' @param node.color A single color for node color.
#' @param edge.color A single color for edge color.
#' @param bg.color A single color for background.
#' @param font.size A single numeric value passed to ggplot themes.
#' @param theme.name Name of a custom PathwaySpace theme. These themes 
#' (from 'th0' to 'th3') consist mainly of preconfigured ggplot settings, 
#' which the user can subsequently fine-tune within the resulting 
#' ggplot object.
#' @param xlab The title for the 'x' axis of a 2D-image space.
#' @param ylab The title for the 'y' axis of a 2D-image space.
#' @param marks A logical value indicating whether to add 'marks' to vertex 
#' positions. Alternatively, this could be a vector listing vertex names.
#' @param mark.size A font size argument passed to \code{\link{geom_text}}.
#' @param mark.color A color passed to \code{\link{geom_text}}.
#' @return A ggplot-class object.
#' @seealso \code{\link{circularProjection}}
#' @examples
#' # Load a demo igraph
#' data('gtoy1', package = 'PathwaySpace')
#'
#' # Create a new PathwaySpace object
#' pts <- buildPathwaySpace(gtoy1)
#'
#' plotGraphSpace(pts)
#' 
#' @import methods
#' @importFrom ggplot2 geom_point geom_segment
#' @importFrom ggforce geom_circle
#' @importFrom grid arrow unit
#' @docType methods
#' @rdname plotGraphSpace-methods
#' @aliases plotGraphSpace
#' @export
#'
setMethod("plotGraphSpace", "PathwaySpace", 
    function(pts, node.size = 1, edge.width = 0.5, 
        node.color = "grey50", edge.color = "grey80", bg.color = "grey95",
        font.size = 1, theme.name = c("th0", "th1", "th2", "th3"),
        xlab = "Pathway coordinates 1", ylab = "Pathway coordinates 2",
        marks = FALSE, mark.size = 3, mark.color = "grey20") {
        #--- validate the pts object and args
        theme.name <- match.arg(theme.name)
        .validate.args("singleNumber", "node.size", node.size)
        .validate.args("singleNumber", "edge.width", edge.width)
        .validate.colors("singleColor", "node.color", node.color)
        .validate.colors("singleColor", "edge.color", edge.color)
        .validate.plot.args("bg.color", bg.color)
        .validate.args("singleNumber", "font.size", font.size)
        .validate.args("singleString", "xlab", xlab)
        .validate.args("singleString", "ylab", ylab)
        .validate.plot.args("marks", marks)
        .validate.args("numeric_vec","mark.size", mark.size)
        .validate.colors("singleColor","mark.color", mark.color)
        #--- get slots from pts
        gxy <- getPathwaySpace(pts, "gxy")
        pars <- getPathwaySpace(pts, "pars")
        edges <- getPathwaySpace(pts, "edges")
        #--- set main input for ggplot
        gxy <- as.data.frame(gxy)
        from <- c(1, pars$nrc)
        gxy$X <- scales::rescale(gxy$X, from = from)
        gxy$Y <- scales::rescale(gxy$Y, from = from)
        #--- get segments
        exy <- .get.exy(gxy, edges)
        #--- get ggplot object
        cl <- .set.theme.bks(theme.name)
        if(pars$is.directed){
            ggp <- .get.dir.graph(gxy, exy, xlab, ylab, cl, 
                node.size, edge.width, node.color, edge.color)
        } else {
            ggp <- .get.und.graph(gxy, exy, xlab, ylab, cl,
                node.size, edge.width, node.color, edge.color)
        }
        #--- add marks if available
        bl <- is.logical(marks) && marks
        if (bl || is.character(marks)) {
            if(bl) marks <- rownames(gxy)
            ggp <- .add.node.marks(ggp, gxy, marks, mark.color, 
                mark.size, node.size)
        }
        #--- apply custom theme
        ggp <- .custom.themes(ggp, theme.name, 
            font.size=font.size, bg.color=bg.color)
        return(ggp)
    }
)

#-------------------------------------------------------------------------------
.get.exy <- function(gxy, edges){
    exy <- data.frame(
        x1 = gxy[edges$vertex1,"X"], 
        x2 = gxy[edges$vertex2,"X"], 
        y1 = gxy[edges$vertex1,"Y"], 
        y2 = gxy[edges$vertex2,"Y"])
    exy$emode <- edges$emode
    return(exy)
}

#-------------------------------------------------------------------------------
.get.und.graph <- function(gxy, exy, xlab, ylab, cl, node.size, 
    edge.width, node.color, edge.color){
    X <- Y <- x1 <- x2 <- y1 <- y2 <- NULL
    ggp <- ggplot2::ggplot(gxy, ggplot2::aes(X, Y)) +
        ggplot2::scale_x_continuous(name = xlab, breaks = cl$axis.ticks,
            labels = format(cl$axis.ticks), position = cl$x.position,
            limits = cl$xylim, expand = ggplot2::expansion(mult = 0)) +
        ggplot2::scale_y_continuous(name = ylab, breaks = cl$axis.ticks,
            labels = format(cl$axis.ticks), limits = cl$xylim,
            expand = ggplot2::expansion(mult = 0)) +
        ggplot2::coord_fixed() + 
        geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), 
            data = exy, size=edge.width, colour = edge.color)
    ggp <- ggp + ggplot2::geom_point(colour=node.color, 
        size=node.size, stroke=ifelse(node.size<0.5, node.size, 0.5))
    return(ggp)
}
.get.dir.graph <- function(gxy, exy, xlab, ylab, cl, node.size, 
    edge.width, node.color, edge.color){
    node.size <- node.size * 0.004
    offset <- node.size + 0.02
    ends <- c("last","both")[exy$emode]
    arr <- grid::arrow(angle=30, length = grid::unit(0.02, "npc"), 
        type = "open", ends=ends)
    exy$offset.start <- c(0, offset)[exy$emode]
    exy$offset.end <- c(offset, offset)[exy$emode]
    exy <- .offset.exy(exy)
    X <- Y <- R <- x1 <- x2 <- y1 <- y2 <- NULL
    gxy$R <- node.size
    ggp <- ggplot2::ggplot(gxy, ggplot2::aes(X, Y)) +
        ggplot2::scale_x_continuous(name = xlab, breaks = cl$axis.ticks,
            labels = format(cl$axis.ticks), position = cl$x.position,
            limits = cl$xylim, expand = ggplot2::expansion(mult = 0)) +
        ggplot2::scale_y_continuous(name = ylab, breaks = cl$axis.ticks,
            labels = format(cl$axis.ticks), limits = cl$xylim,
            expand = ggplot2::expansion(mult = 0)) +
        ggplot2::coord_fixed() + 
        ggplot2::geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), 
            arrow = arr, data = exy,  size=edge.width, colour = edge.color)
    ggp <- ggp + ggforce::geom_circle(
        aes(x0 = X, y0 = Y, r = R), fill = node.color,
        colour=node.color, show.legend = FALSE)
    return(ggp)
}

#-------------------------------------------------------------------------------
.offset.exy <- function(exy){
    exy$dx <- exy$x2 - exy$x1
    exy$dy <- exy$y2 - exy$y1
    exy$dist <- sqrt( exy$dx^2 + exy$dy^2 )
    exy$px <- exy$dx/exy$dist
    exy$py <- exy$dy/exy$dist
    exy$x1 <- exy$x1 + (exy$px * exy$offset.start)
    exy$y1 <- exy$y1 + (exy$py * exy$offset.start)
    exy$x2 <- exy$x2 - (exy$px * exy$offset.end)
    exy$y2 <- exy$y2 - (exy$py * exy$offset.end)
    return(exy)
}


#-------------------------------------------------------------------------------
.add.node.marks <- function(ggp, gxy, marks, mark.color, 
    mark.size, node.size) {
    node.size <- node.size * 0.005
    vjust <- 0
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
    gxy$mark.color <- mark.color
    gxy$mark.size <- mark.size
    gxy <- gxy[marks, , drop = FALSE]
    gxy$ID <- names(marks); ID <- NULL
    ggp <- ggp + ggplot2::geom_text(mapping = aes(label = ID), 
        data = gxy, fontface = "bold",
        nudge_y = node.size,
        size = mark.size,
        vjust=vjust,
        colour = mark.color)
    return(ggp)
}


################################################################################
### Package documentation
################################################################################
#' @keywords internal
#' @author
#' The Cancer Genome Atlas Analysis Network.
#'
#' @description
#' For a given graph containing vertices, edges, and a signal associated with
#' the vertices, the PathwaySpace package performs a convolution operation,
#' which involves a weighted combination of neighboring vertices and their
#' associated signals. The package then uses a decay function to propagate
#' these signals, creating geodesic paths on a 2D-image space.
#'
#' @details
#'
#' \tabular{ll}{
#' Package: \tab PathwaySpace\cr
#' Type: \tab Software\cr
#' License: \tab Artistic-2.0\cr
#' Maintainer: \tab Mauro Castro \email{mauro.a.castro@@gmail.com}\cr
#' }
#'
#' @section Index:
#' \tabular{ll}{
#' \link{PathwaySpace-class}: 
#' \tab An S4 class for signal propagation on pathway spaces.\cr
#' \link{buildPathwaySpace}: 
#' \tab Constructor of PathwaySpace-class objects.\cr
#' \link{circularProjection}: 
#' \tab Creating 2D-landscape images from graph objects.\cr
#' \link{polarProjection}: 
#' \tab Creating 2D-landscape images from graph objects.\cr
#' \link{silhouetteMapping}: 
#' \tab Mapping graph silhouettes on PathwaySpace images.\cr
#' \link{summitMapping}: 
#' \tab Mapping summits on a 2D-landscape image.\cr
#' \link{getPathwaySpace}: 
#' \tab Accessory method for fetching slots from a PathwaySpace object.\cr
#' \link{plotImageSpace}: 
#' \tab Plotting 2D-landscape images for the PathwaySpace package.\cr
#' }
#' Further information is available in the vignettes by typing
#' \code{vignette('PathwaySpace')}. Documented topics are also available in
#' HTML by typing \code{help.start()} and selecting the PathwaySpace package
#' from the menu.
#'
#' @references
#' The Cancer Genome Atlas Analysis Network (2023). PathwaySpace: Spatial propagation 
#' of network signals along geodesic paths. R package version 0.99.
#'
"_PACKAGE"
#> [1] '_PACKAGE'

################################################################################
### Documentation for some 'toy' datasets
################################################################################
#' @title Toy 'igraph' objects.
#'
#' @description Small 'igraph' objects used for workflow demonstrations.
#' All graphs include 'x', 'y', and 'name' vertex attributes.
#'
#' @format igraph
#'
#' @usage data(gtoy1)
#'
#' @source This package.
#'
#' @docType data
#' @keywords gtoys
#' @name gtoys
#' @aliases gtoy1 gtoy2
#' @return A pre-processed igraph object.
#' @examples
#' data(gtoy1)
#' data(gtoy2)
NULL

################################################################################
### Documentation for the 'gimage' dataset
################################################################################
#' @title An image matrix.
#'
#' @description An image matrix used for workflow demonstrations.
#'
#' @format matrix
#'
#' @usage data(gimage)
#'
#' @source This package.
#'
#' @docType data
#' @keywords gimage
#' @name gimage
#' @aliases gimage
#' @return An image matrix.
#' @examples
#' data(gimage)
NULL

################################################################################
### Documentation for the 'PCv12_pruned_igraph' dataset
################################################################################
#' @title A pruned and laid out igraph object from Pathway Commons V12.
#'
#' @description
#' This igraph object was created from a 'sif' file available from 
#' the Pathway Commons V12 (Rodchenkov et al., 2020), which was filtered to 
#' keep interactions from the following sources: CTD, Recon, HumanCyc, 
#' DrugBank, MSigDB, DIP, BioGRID, IntAct, BIND, and PhosphoSite. The igraph 
#' was additionally pruned and laid out by a force-directed algorithm aiming 
#' signal projection on PathwaySpace's images. Edges with the smallest 
#' betweenness centrality were pruned using 'backward elimination' and
#' 'forward selection' strategies. The resulting graph represents the 
#' main connected component with the minimum number of edges.
#'
#' @format igraph
#'
#' @usage data(PCv12_pruned_igraph)
#'
#' @source Pathway Commons V12.
#'
#' @author Chris Wong, Mauro Castro, and TCGA Network.
#'
#' @references
#' Rodchenkov et al. Pathway Commons 2019 Update: integration, analysis and
#' exploration of pathway data. Nucleic Acids Research 48(D1):D489–D497, 2020.
#' \doi{10.1093/nar/gkz946}
#' \href{https://www.pathwaycommons.org/archives/PC2/v12/}{Pathway Commons V12}
#'
#' @docType data
#' @keywords PCv12_pruned_igraph
#' @name PCv12_pruned_igraph
#' @aliases PCv12_pruned_igraph
#' @return An igraph object.
#' @examples
#' data(PCv12_pruned_igraph)
#' ## Suggestion to vizualize this igraph in R:
#' library(RedeR)
#' library(igraph)
#' rdp <- RedPort()
#' ## calld(rdp)
#' ## resetd(rdp)
#' ## addGraph(rdp, PCv12_pruned_igraph)
NULL

################################################################################
### Documentation for the 'CGC_20211118' dataset
################################################################################
#' @title COSMIC-CGC genes mapped to PathwaySpace images.
#'
#' @description
#' A data frame listing 'GeneSymbol' and 'Entrez' IDs from the
#' COSMIC-CGC database (Sondka et al., 2020). These genes are used to
#' demonstrate the PathwaySpace's summit mapping pipeline, which assigns 
#' summits to an image space.
#'
#' @format data.frame
#'
#' @usage data(CGC_20211118)
#'
#' @source COSMIC-CGC database (release v95, tier 1 collection).
#'
#' @references
#' Sondka et al. The COSMIC Cancer Gene Census: describing genetic
#' dysfunction across all human cancers. Nat Rev Cancer 18, 696-705, 2018.
#' Doi: 10.1038/s41568-018-0060-1.
#'
#' @docType data
#' @keywords CGC_20211118
#' @name CGC_20211118
#' @aliases CGC_20211118
#' @return A data.frame object.
#' @examples
#' data(CGC_20211118)
NULL


################################################################################
### Documentation for the 'Hallmarks_v2023_1_Hs_symbols' dataset
################################################################################
#' @title A list with Hallmark gene sets (v2023.1).
#'
#' @description
#' A list with Human gene symbols from the MSigDB's Hallmark gene set 
#' collection (Liberzon et al., 2015). These gene sets are used to
#' demonstrate the PathwaySpace's summit mapping pipeline, which assigns 
#' summits to an image space.
#'
#' @format list
#'
#' @usage data(Hallmarks_v2023_1_Hs_symbols)
#'
#' @source MSigDB database (v2023.1).
#'
#' @references
#' Liberzon et al. The Molecular Signatures Database (MSigDB) hallmark 
#' gene set collection. Cell Systems 1(5):417-425, 2015
#' Doi: 10.1016/j.cels.2015.12.004.
#'
#' @docType data
#' @keywords Hallmarks_v2023_1_Hs_symbols
#' @name Hallmarks_v2023_1_Hs_symbols
#' @aliases Hallmarks_v2023_1_Hs_symbols
#' @aliases hallmarks
#' @return A list object.
#' @examples
#' data(Hallmarks_v2023_1_Hs_symbols)
NULL


#' @title A simple vector of colors for PathwaySpace images.
#'
#' @param n Number of colors.
#' @return A vector with hexadecimal color codes.
#' @seealso \code{\link{plotImageSpace}}
#' @examples
#' pspace.cols()
#'
#' @aliases pspace.cols
#' @export
#'
pspace.cols <- function(n=5) {
    colors <- c("#303f9d","#578edb","#63b946","#f3930c","#a60d0d")
    colorRampPalette(colors)(n)
}



################################################################################
##########################  misc. accessory functions  #########################
################################################################################


################################################################################
### Count the number of local neighbors
################################################################################
.countPxNeighbors <- function(x) {
    x <- xx <- .expandMask(x)
    nr <- nrow(x)
    nc <- ncol(x)
    ii <- jj <- c(-1, 0, 1)
    for (i in 2:(nr - 1)) {
        for (j in which(x[i, ] == 1)) {
            xx[i, j] <- sum(x[ii + i, jj + j] == 1) - 1
        }
    }
    xx <- xx[-c(1, nr), -c(1, nc)]
    return(xx)
}
.expandMask <- function(x, val = 0) {
    xx <- array(val, dim = dim(x) + 2)
    nr <- nrow(xx)
    nc <- ncol(xx)
    xx[-c(1, nr), -c(1, nc)] <- x
    return(xx)
}
.reduceMask <- function(x) {
    nr <- nrow(x)
    nc <- ncol(x)
    x <- x[-c(1, nr), -c(1, nc)]
    return(x)
}

################################################################################
### Find outlines in a labeled image mask (labels near bg's '0s')
################################################################################
.findOutlines <- function(mask) {
    ids <- sort(unique(as.numeric(mask)), decreasing = TRUE)
    ids <- ids[ids > 0]
    x <- mask
    x[, ] <- 0
    for (id in ids) {
        tp <- .findPxEdges(mask == id)
        x[tp == 1] <- id
    }
    x
}

################################################################################
### Erode and dilatate pixels of a mask
################################################################################
.openPxEdges <- function(mask, kernel = .kernelMask(2)) {
    mask <- .erodePxEdges(mask, kernel)
    mask <- .dilatePxEdges(mask, kernel)
    return(mask)
}

################################################################################
### Dilatate pixels of a mask
################################################################################
.dilatePxEdges <- function(mask, kernel = .kernelMask(2)) {
    ids <- as.numeric(names(table(mask)))
    ids <- ids[ids > 0]
    for (id in ids) {
        tp <- .dledge(mask == id, kernel = kernel)
        mask[tp == 1] <- id
    }
    mask
}
.dledge <- function(dx, kernel) {
    nr <- nrow(dx)
    nc <- ncol(dx)
    n <- (dim(kernel)[1] - 1)/2
    ed <- .findPxEdges(dx)
    idx <- which(ed == 1, arr.ind = TRUE)
    if (nrow(idx) > 0) {
        for (id in seq_len(nrow(idx))) {
            i <- idx[id, 1]
            j <- idx[id, 2]
            ii <- (i - n):(i + n)
            jj <- (j - n):(j + n)
            bi <- ii > 0 & ii <= nr
            bj <- jj > 0 & jj <= nc
            dx[ii[bi], jj[bj]] <- dx[ii[bi], jj[bj]] + kernel[bi, bj]
        }
        dx[dx > 0] <- 1
    }
    dx
}

################################################################################
### Erode pixels of a mask
################################################################################
.erodePxEdges <- function(mask, kernel = .kernelMask(2)) {
    x <- mask
    x[x != 0] <- 1
    tp <- .eledge(x, kernel = kernel)
    mask[tp == 0] <- 0
    mask
}
.eledge <- function(ex, kernel) {
    nr <- nrow(ex)
    nc <- ncol(ex)
    n <- (dim(kernel)[1] - 1)/2
    ed <- .findBgEdges(ex)
    idx <- which(ed == 1, arr.ind = TRUE)
    if (nrow(idx) > 0) {
        for (id in seq_len(nrow(idx))) {
            i <- idx[id, 1]
            j <- idx[id, 2]
            ii <- (i - n):(i + n)
            jj <- (j - n):(j + n)
            bi <- ii > 0 & ii <= nr
            bj <- jj > 0 & jj <= nc
            ex[ii[bi], jj[bj]] <- ex[ii[bi], jj[bj]] - kernel[bi, bj]
        }
        ex[ex < 0] <- 0
    } else {
        ex[] <- 0
    }
    ex
}

################################################################################
### Remove tips with less than n neighbors
################################################################################
.removeTips <- function(mask, n = 3) {
    tp <- .countPxNeighbors(mask == 1)
    array(as.numeric(tp > n), dim = dim(tp))
}

################################################################################
### Find neighbors in an image mask
################################################################################
# Find bg (0s) near to pixels (1s)
.findBgEdges <- function(mask) {
    ed <- .expandPxEdges(mask)
    ed[mask == 1] <- 0
    ed
}
# Find pixels (1s) near bg (0s)
.findPxEdges <- function(mask) {
    ed <- .countPxNeighbors(mask)
    ed <- ed > 0 & ed < 8
    ed
}

################################################################################
### Expand pixels (in one unit)
################################################################################
.expandPxEdges <- function(mask) {
    nrc <- dim(mask)
    idx <- which(mask == 1, arr.ind = TRUE)
    ii <- idx[, 1]
    jj <- idx[, 2]
    ii <- c(ii - 1, ii, ii, ii + 1)
    jj <- c(jj, jj - 1, jj + 1, jj)
    ix <- (ii > 0 & ii <= nrc[1]) & (jj > 0 & jj <= nrc[2])
    ii <- ii[ix]
    jj <- jj[ix]
    idx <- mask[cbind(ii, jj)]
    mask[cbind(ii, jj)][idx == 0] <- 1
    return(mask)
}

################################################################################
### Fill cavity of a labeled mask
################################################################################
.fillCavity <- function(mask) {
    xm <- mask
    xm[is.na(xm)] <- 0
    ids <- sort(unique(as.numeric(xm)), decreasing = TRUE)
    ids <- ids[ids > 0]
    for (id in ids) {
        idx <- which(xm == id, arr.ind = TRUE)
        rr <- range(idx[, 1])
        rr <- rr[1]:rr[2]
        rc <- range(idx[, 2])
        rc <- rc[1]:rc[2]
        if (length(rr) > 3 && length(rc) > 3) {
            x1 <- x2 <- xm[rr, rc, drop = FALSE]
            x2[x2 != id] <- 0
            x2 <- x2 == 0
            x2 <- .expandMask(x2, val = 1)
            #--- fill large cavities
            x3 <- .labelMask(x2)
            x3 <- .reduceMask(x3)
            x1[x1 == 0 & x3 > 1] <- id
            #--- fill small spots
            x3 <- .removeTips(x2, n = 3)
            x3 <- .reduceMask(x3)
            x4 <- .reduceMask(x2)
            x1[x3 != x4] <- id
            mask[rr, rc] <- x1
        }
    }
    mask
}
.countCavity <- function(mask) {
    xm <- mask
    xm[is.na(xm)] <- 0
    ids <- sort(unique(as.numeric(xm)))
    ids <- ids[ids > 0]
    nc <- vapply(ids, function(id) {
        idx <- which(xm == id, arr.ind = TRUE)
        rr <- range(idx[, 1])
        rr <- rr[1]:rr[2]
        rc <- range(idx[, 2])
        rc <- rc[1]:rc[2]
        if (length(rr) > 3 && length(rc) > 3) {
            x1 <- x2 <- xm[rr, rc, drop = FALSE]
            x2[x2 != id] <- 0
            x2 <- x2 == 0
            x2 <- .expandMask(x2, val = 1)
            x2 <- .labelMask(x2)
            n <- length(table(x2)) - 2
        } else {
            n <- 0
        }
        n
    }, numeric(1))
    return(nc)
}

################################################################################
### Make a kernel mask
################################################################################
.kernelMask <- function(r = 5, shape = c("round", "diamond", "square")) {
    shape <- match.arg(shape)
    if (shape == "round") {
        kmask <- .rdKernel(r)
    } else if (shape == "diamond") {
        kmask <- .dmKernel(r)
    } else {
        kmask <- .sqKernel(r)
    }
    kmask
}
.rdKernel <- function(r = 5) {
    len <- r * 2 + 1
    x <- matrix(0, nrow = len, ncol = len)
    cnt <- (len + 1)/2
    r <- cnt - 1
    a <- (cnt - r):(cnt + r)
    b <- (cnt - r):(cnt + r)
    for (i in a) {
        for (j in b) {
            if (sqrt((i - cnt)^2 + (j - cnt)^2) <= r) {
                x[i, j] <- 1
            }
        }
    }
    x
}
.dmKernel <- function(r = 5) {
    len <- r * 2 + 1
    x <- matrix(0, nrow = len, ncol = len)
    cnt <- (len + 1)/2
    x[cnt, cnt] <- 1
    x <- .expandKernel(x)
    x
}
.sqKernel <- function(r = 5) {
    len <- r * 2 + 1
    x <- matrix(0, nrow = len, ncol = len)
    x[, ] <- 1
    x
}
.expandKernel <- function(x) {
    while (!.anyOuterPx(x)) {
        x <- .expandPxEdges(x)
    }
    x
}
.anyOuterPx <- function(x) {
    d <- dim(x)
    b <- c(x[1, ], x[d[1], ], x[, 1], x[, d[2]])
    ifelse(any(b == 1), TRUE, FALSE)
}

################################################################################
### Label contiguous pixels in an image mask
################################################################################
.labelMask <- function(mask) {
    mask[is.na(mask)] <- 0
    mask <- mask > 0
    idxmask <- data.frame(which(mask | !mask, arr.ind = TRUE))
    idxmask$key <- seq_len(nrow(idxmask))
    idxmask$val <- as.numeric(mask)
    idxmask$label <- NA
    label <- 1
    seed <- idxmask[which(idxmask$val == 1)[1], , drop = FALSE]
    while (sum(idxmask$val) > 0) {
        idxmask$label[seed$key] <- label
        idxmask$val[seed$key] <- 0
        seed <- .expandSeed(seed, idxmask, nrc = dim(mask))
        if (nrow(seed) == 0) {
            seed <- idxmask[which(idxmask$val == 1)[1], , drop = FALSE]
            label <- label + 1
        }
    }
    idxmask <- idxmask[!is.na(idxmask$label), , drop = FALSE]
    mask[cbind(idxmask$row, idxmask$col)] <- idxmask$label
    mask
}
.expandSeed <- function(seed, idx, nrc) {
    ii <- seed[, "row"]
    jj <- seed[, "col"]
    ii <- c(ii - 1, ii - 1, ii - 1, ii, ii, ii + 1, ii + 1, ii + 1)
    jj <- c(jj - 1, jj, jj + 1, jj - 1, jj + 1, jj - 1, jj, jj + 1)
    ix <- (ii > 0 & ii <= nrc[1]) & (jj > 0 & jj <= nrc[2])
    ii <- ii[ix]
    jj <- jj[ix]
    ij <- ii + nrc[1] * (jj - 1)
    seed <- idx[unique(ij), , drop = FALSE]
    seed <- seed[seed$val > 0, , drop = FALSE]
    return(seed)
}

################################################################################
### Relabel a mask with labeled objects
################################################################################
.relabel <- function(mask) {
    lbs <- sort(unique(as.numeric(mask)))
    lbs <- lbs[lbs > 0]
    rmask <- mask
    for (i in seq_along(lbs)) {
        rmask[mask == lbs[i]] <- i
    }
    return(rmask)
}
.relabelBySize <- function(mask) {
    objs <- table(mask)
    lbs <- as.integer(names(objs))
    objs <- objs[lbs > 0]
    lbs <- lbs[lbs > 0]
    lbs <- lbs[order(objs, decreasing = TRUE)]
    rmask <- mask
    for (i in seq_along(lbs)) {
        rmask[mask == lbs[i]] <- i
    }
    return(rmask)
}
.relabelBySignal <- function(x, mask, decreasing = TRUE) {
    mask[is.na(mask)] <- 0
    lbs <- sort(unique(as.numeric(mask)))
    lbs <- lbs[lbs > 0]
    sig <- vapply(lbs, function(lb) {
        quantile(x[mask == lb], 0.95, na.rm = TRUE)
    }, numeric(1))
    lbs <- lbs[order(sig, decreasing = decreasing)]
    rmask <- mask
    for (i in seq_along(lbs)) {
        rmask[mask == lbs[i]] <- i
    }
    return(rmask)
}

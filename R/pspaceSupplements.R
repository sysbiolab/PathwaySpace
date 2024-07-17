################################################################################
### Main constructor of PathwaySpace-class objects
################################################################################
.buildPathwaySpace <- function(g, nrc = 500, mar = 0.075, verbose = TRUE) {
    #--- get graph data
    if(verbose) message("Extracting 'x', 'y', and 'name' vertex attributes...")
    g <- igraph::simplify(g)
    X <- igraph::V(g)$x
    Y <- igraph::V(g)$y
    vertex <- igraph::V(g)$name
    if(igraph::is_directed(g)){
        if(verbose) message("Extracting directed edges...")
        edges <- .get.directed.edges(g, vertex)
    } else {
        if(verbose) message("Extracting undirected edges...")
        edges <- .get.undirected.edges(g, vertex)
    }
    g <- .validate.vertex.signal(g, verbose)
    vsignal <- igraph::V(g)$signal
    vweight <- igraph::V(g)$weight
    names(vsignal) <- names(vweight) <- vertex
    vsignal <- .revise.vertex.signal(vsignal)
    vweight <- .revise.vertex.weight(vweight)
    #--- get signal scale and put pars in a list
    zscale <- .getSignalScale(vsignal)
    vwscale <- .get.vwscale(vweight)
    pars <- list(nrc = nrc, mar = mar, 
        zscale = zscale, vwscale=vwscale,
        is.directed = igraph::is_directed(g))
    #--- combine xy and center nodes
    gxy <- cbind(X = X, Y = Y)
    rownames(gxy) <- vertex
    gxy <- .centerNodes(gxy, nrc, mar)
    # get dists not related to vsignal
    if(verbose) message("Computing distances between vertices...")
    edges <- .get.edge.dist(gxy, edges)
    # add vsignal and vweight to gxy
    gxy <- cbind(gxy, vsignal = vsignal, vweight = vweight)
    #---create a PathwaySpace object
    plnames <- c("Preprocess", "CircularProjection", "PolarProjection",
        "Silhouette", "Summits")
    status <- rep("[ ]", length(plnames))
    names(status) <- plnames
    pts <- new(Class = "PathwaySpace",  vertex = vertex,
        vsignal = vsignal, vweight = vweight, 
        edges = edges, gxy = gxy, pars = pars, 
        misc=list(g=g), status = status)
    return(pts)
}

#-------------------------------------------------------------------------------
.get.directed.edges <- function(g, vertex){
    idxmutual <- igraph::which_mutual(g)
    E(g)$emode <- 1
    E(g)$emode[idxmutual] <- 2
    emode <- igraph::as_adjacency_matrix(g, sparse = FALSE, attr = "emode")
    bl <- lower.tri(emode) & emode==2
    emode[bl] <- 0
    edges <- arrayInd(seq_len(prod(dim(emode))), dim(emode), useNames = TRUE)
    edges <- as.data.frame(edges)
    colnames(edges) <- c("vertex1", "vertex2")
    edges$emode <- as.numeric(emode)
    edges <- edges[edges$emode>0,]
    edges$name1 <- vertex[edges$vertex1]
    edges$name2 <- vertex[edges$vertex2]
    edges <- edges[order(edges$vertex1,edges$vertex2), ]
    rownames(edges) <- NULL
    edges$eweight <- 1
    return(edges)
}
.get.undirected.edges <- function(g, vertex){
    edges <- as_edgelist(g, names = FALSE)
    rownames(edges) <- colnames(edges) <- NULL
    edges <- as.data.frame(edges)
    colnames(edges) <- c("vertex1", "vertex2")
    edges$name1 <- vertex[edges$vertex1]
    edges$name2 <- vertex[edges$vertex2]
    edges <- edges[order(edges$vertex1,edges$vertex2), ]
    edges$eweight <- 1
    edges$emode <- 0
    return(edges)
}

#-------------------------------------------------------------------------------
#--- get point-to-vertice dists for circular and polar projections
.get.knn.dists <- function(lpts, gxy, knn, bg.rm = TRUE){
    # message("Computing distance between k-nearest neighbors...")
    if(bg.rm){
        bg <- is.na(gxy[,"vsignal"]) | gxy[,"vsignal"]==0
        if(all(bg)){
            bg.rm <- FALSE
        } else {
            gxy <- cbind(gxy, key=seq_len(nrow(gxy)))
            gxy <- gxy[!bg, , drop=FALSE]
        }
    }
    if (knn == 1) {
        min.ksearch <- ceiling(nrow(gxy) * 0.1)
        min.ksearch <- min(max(min.ksearch, 10), nrow(gxy))
    } else {
        min.ksearch <- min(max(knn, 10), nrow(gxy))
    }
    nnpg <- RANN::nn2(gxy[, c("X", "Y"), drop=FALSE], 
        lpts[, c("X", "Y"), drop=FALSE], k = min.ksearch)
    names(nnpg) <- c("nn", "dist")
    if(bg.rm){
        nnpg$nn[,] <- as.numeric(gxy[as.numeric(nnpg$nn),"key"])
    }
    return(nnpg)
}

#-------------------------------------------------------------------------------
#--- get point-to-vertice dists for silhouettes
.get.silhouette.dists <- function(lpts, gxy, knn){
    min.ksearch <- min(max(knn, 10), nrow(gxy))
    nnbg <- RANN::nn2(gxy[, c("X", "Y")], lpts[, c("X", "Y")],
        k = min.ksearch)
    names(nnbg) <- c("nn", "dist")
    return(nnbg)
}

################################################################################
### Accessors of the constructor
################################################################################
.updateStatus <- function(pts, name, check = TRUE) {
    pts@status[name] <- ifelse(check, "[x]", "[ ]")
    return(pts)
}
.checkStatus <- function(pts, name) {
    if(name=="Projection"){
        projections <- c("PolarProjection","CircularProjection")
        sts <- any(pts@status[projections] == "[x]")
    } else {
        sts <- pts@status[name] == "[x]"
    }
    sts
}
.removeSummits <- function(pts, verbose=TRUE){
    if (.checkStatus(pts, "Summits")) {
        if(verbose) message("-- available summits removed.")
        i <- which(names(pts@misc) == "summits")
        pts@misc <- pts@misc[-i]
        pts <- .updateStatus(pts, "Summits", FALSE)
    } 
    return(pts)
}
.removeSilhouette <- function(pts, verbose=TRUE){
    if (.checkStatus(pts, "Silhouette")) {
        if(verbose) message("-- available silhouette removed.")
        i <- which(names(pts@misc) == "xfloor")
        pts@misc <- pts@misc[-i]
        pts <- .updateStatus(pts, "Silhouette", FALSE)
    } 
    return(pts)
}
.getSignalScale <- function(zsignal) {
    zscale <- list(range = range(zsignal, na.rm = TRUE))
    zscale$maxsig <- max(abs(zscale$range))
    if (.all_binaryValues(zsignal)) {
        # message("binary")
        zscale$signal.type <- "binary"
        zscale$scale.type <- "pos"
        zscale$range <- c(0, 1)
        zscale$scaling <- c(0, 1)
    } else if (zscale$range[1] == zscale$range[2]) {
        # message("found a single value; it will be treated as binary")
        zscale$signal.type <- "binary"
        if (zscale$range[1] > 0) {
            zscale$range[1] <- 0
            zscale$scaling <- c(0, 1)
            zscale$scale.type <- "pos"
        } else {
            zscale$range[2] <- 0
            zscale$scaling <- c(-1, 0)
            zscale$scale.type <- "neg"
        }
    } else if (sum(sign(zscale$range)) == 0) {
        # message("continuous in (-Inf,Inf)")
        zscale$signal.type <- "continuous"
        zscale$scale.type <- "negpos"
        zscale$scaling <- c(-1, 1)
    } else if (all(zscale$range >= 0)) {
        # message("continuous >=0")
        zscale$signal.type <- "continuous"
        zscale$scale.type <- "pos"
        zscale$scaling <- c(0, 1)
    } else if (all(zscale$range <= 0)) {
        # message("continuous <=0")
        zscale$signal.type <- "continuous"
        zscale$scale.type <- "neg"
        zscale$scaling <- c(-1, 0)
    }
    return(zscale)
}
.get.vwscale <- function(vweight){
    if (sd(vweight, na.rm = TRUE) > 0) {
        vwscale <- TRUE
    } else {
        vwscale <- FALSE
    }
    return(vwscale)
}
.get.edge.dist <- function(gxy, edges) {
    dts <- sqrt((gxy[edges$name1, "X"] - gxy[edges$name2, "X"])^2 +
        (gxy[edges$name1, "Y"] - gxy[edges$name2, "Y"])^2)
    edges$edist <- dts
    return(edges)
}
.centerNodes <- function(gxy, nrc, mar = 0.1) {
    #-- get scaling factors and margins
    nmar <- floor(nrc * mar)
    nrc.mar <- nrc - (nmar * 2)
    #-- rescale in [0,1], keep aspect ratio
    gxy[, 1] <- gxy[, 1] - min(gxy[, 1])
    gxy[, 2] <- gxy[, 2] - min(gxy[, 2])
    maxij <- max(gxy[, c(1, 2)])
    if(maxij==0){
        gxy[,] <- nrc/2
    } else {
        gxy[, 1] <- gxy[, 1] / maxij
        gxy[, 2] <- gxy[, 2] / maxij
        #-- set center and aspect ratio
        ci <- 1 - max(gxy[, 1]) - min(gxy[, 1])
        cj <- 1 - max(gxy[, 2]) - min(gxy[, 2])
        gxy[, 1] <- gxy[, 1] + ci / 2
        gxy[, 2] <- gxy[, 2] + cj / 2
        #-- rescale in nrc.mar
        gxy[, 1] <- gxy[, 1] * (nrc.mar - 1)
        gxy[, 2] <- gxy[, 2] * (nrc.mar - 1)
        #-- set center and coordinates
        gxy[, 1] <- gxy[, 1] + 1
        gxy[, 2] <- gxy[, 2] + 1
        gxy[, 1] <- gxy[, 1] + nmar
        gxy[, 2] <- gxy[, 2] + nmar
    }
    gxy <- cbind(gxy, Xint = round(gxy[, 1]), Yint = round(gxy[, 2]))
    return(gxy)
}

################################################################################
### Main PathwaySpace calls
################################################################################
.circularProjection <- function(pts, verbose = TRUE, ...) {
    if(verbose) message("Using circular projection...")
    pars <- getPathwaySpace(pts, "pars")
    gxy <- getPathwaySpace(pts, "gxy")
    lpts <- .pointsInMatrix(pars$nrc)
    nnpg <- .get.knn.dists(lpts, gxy, pars$knn)
    if (pars$vwscale) {
        if(verbose) message("Scaling projection to vertex weight...")
    }
    #--- compute landscape signal
    if(verbose) message("Running signal processing and convolution...")
    xsig <- array(0, c(pars$nrc, pars$nrc))
    if (nrow(gxy) > 0) {
        xsig[lpts[, c("Y", "X")]] <- .get.ldsig(gxy,
            pars, nnpg, ... = ...)
    }
    pts@misc$xsig <- xsig
    # image(.transpose.and.flip(xsig))
    #--- rescale xsig to the original signal
    # use silhouette if available
    if (.checkStatus(pts, "Silhouette")){
      xfloor <- pts@misc$xfloor
      gxyz <- .rescale.landscape(xsig = xsig, pars = pars, xfloor = xfloor)
    } else {
      gxyz <- .rescale.landscape(xsig = xsig, pars = pars)
    }
    pts@gxyz <- gxyz
    # image(.transpose.and.flip(gxyz))
    #---return the gxyz landscape
    x <- gxyz
    x[x == 0] <- NA
    sz <- round(sum(abs(x), na.rm = TRUE) / sum(!is.na(x)), 4)
    sz <- ifelse(is.nan(sz), 0, sz)
    if(verbose) message("Density: ", sz, " per pixel!")
    return(pts)
}

#-------------------------------------------------------------------------------
.polarProjection <- function(pts, verbose = TRUE, ...) {
    edges <- getPathwaySpace(pts, "edges")
    if(nrow(edges)==0){
        msg <- paste0("Note, this 'PathwaySpace' object does not ",
            "contain edges.\nThe 'polarProjection' method requires ",
            "at least one edge for projection.")
        stop(msg)
    }
    pars <- getPathwaySpace(pts, "pars")
    if(pars$directional){
        if(pars$is.directed){
            message("Using directional polar projection...")
        } else {
            stop("'directional' used with undirected graphs.")
        }
    } else {
        if(verbose) message("Using undirectional polar projection...")
    }
    gxy <- getPathwaySpace(pts, "gxy")
    lpts <- .pointsInMatrix(pars$nrc)
    nnpg <- .get.knn.dists(lpts, gxy, pars$knn)
    if (pars$vwscale) {
        if(verbose) message("Scaling projection to vertex weight...")
    }
    #--- for polar, 'pdist' is applied at edge dist
    rg <- range(edges$edist)
    rg <- rg[1] + ( (rg - rg[1]) * pars$pdist )
    edges$edistNorm <- scales::rescale(edges$edist, to=rg)
    #--- compute landscape signal
    if(verbose) message("Running signal processing and convolution...")
    xsig <- array(0, c(pars$nrc, pars$nrc))
    if (nrow(gxy) > 0) {
        xsig[lpts[, c("Y", "X")]] <- .get.ldsig.polar(lpts,
            gxy, edges, pars, nnpg, ...)
    }
    # image(.transpose.and.flip(xsig))
    #--- rescale xsig to the original signal
    gxyz <- .rescale.landscape(xsig = xsig, pars = pars)
    pts@misc$xsig <- xsig
    pts@gxyz <- gxyz
    #---return the gxyz landscape
    x <- gxyz
    x[x == 0] <- NA
    sz <- round(sum(abs(x), na.rm = TRUE) / sum(!is.na(x)), 4)
    sz <- ifelse(is.nan(sz), 0, sz)
    if(verbose) message("Density: ", sz, " per pixel!")
    return(pts)
}

#-------------------------------------------------------------------------------
# If 'pars$rescale = F', it will rescale 'xsig' to the original range
.rescale.landscape <- function(xsig, pars, xfloor = NULL) {
    if (!is.null(xfloor)) xsig[xfloor == 0] <- NA
    bl <- all(range(xsig, na.rm = TRUE) == 0)
    if(bl){
        x <- xsig
    } else {
        if (!pars$rescale) {
            x <- scales::rescale(xsig, to = pars$zscale$range)
        } else if (pars$zscale$signal.type == "binary") {
            x <- scales::rescale(xsig, to = pars$zscale$scaling)
        } else {
            endpoints <- c(0.01, 1e-04)
            if (pars$zscale$scale.type == "negpos") {
                to <- pars$zscale$range/pars$zscale$maxsig
                x <- .rescale.negpos(xsig, pars, to, endpoints = endpoints)
            } else if (pars$zscale$scale.type == "pos") {
                to <- pars$zscale$scaling
                x <- .rescale.pos(xsig, pars, to, endpoints = endpoints)
            } else if (pars$zscale$scale.type == "neg") {
                to <- pars$zscale$scaling
                x <- .rescale.neg(xsig, pars, to, endpoints = endpoints)
            }
            #- maybe export 'endpoints' in future versions;
            #- it will cut noise and outliers at the endpoints;
            #- outliers outside that range are moved to the nearby endpoint;
            #- noise are set to background (i.e. zero).
        }
    }
    return(x)
}
.rescale.negpos <- function(xsig, pars, to, endpoints) {
    pr <- c(endpoints[1], 1 - endpoints[2])
    # pos
    xp <- xsig
    xp[xp <= 0] <- NA
    qt <- quantile(as.numeric(xp), na.rm = TRUE, probs = pr)
    xp[xp < qt[1]] <- NA
    xp[xp > qt[2]] <- qt[2]
    xp <- scales::rescale(xp, to = c(0, to[2]))
    # neg
    pr <- c(endpoints[2], 1 - endpoints[1])
    xn <- xsig
    xn[xn >= 0] <- NA
    qt <- quantile(as.numeric(xn), na.rm = TRUE, probs = pr)
    xn[xn < qt[1]] <- qt[1]
    xn[xn > qt[2]] <- NA
    xn <- scales::rescale(xn, to = c(to[1], 0))
    # update xsig
    xsig[!is.na(xn)] <- xn[!is.na(xn)]
    xsig[!is.na(xp)] <- xp[!is.na(xp)]
    return(xsig)
}
.rescale.pos <- function(xsig, pars, to, endpoints) {
    pr <- c(endpoints[1], 1 - endpoints[2])
    xp <- xsig
    xp[xp <= 0] <- NA
    qt <- quantile(as.numeric(xp), na.rm = TRUE, probs = pr)
    xp[xp < qt[1]] <- NA
    xp[xp > qt[2]] <- qt[2]
    xp <- scales::rescale(xp, to = to)
    xsig[!is.na(xp)] <- xp[!is.na(xp)]
    return(xsig)
}
.rescale.neg <- function(xsig, pars, to, endpoints) {
    pr <- c(endpoints[2], 1 - endpoints[1])
    xn <- xsig
    xn[xn >= 0] <- NA
    qt <- quantile(as.numeric(xn), na.rm = TRUE, probs = pr)
    xn[xn < qt[1]] <- qt[1]
    xn[xn > qt[2]] <- NA
    xn <- scales::rescale(xn, to = to)
    xsig[!is.na(xn)] <- xn[!is.na(xn)]
    return(xsig)
}

#-------------------------------------------------------------------------------
#--- main circular projection calls
.get.ldsig <- function(gxy, pars, nnpg, ...) {
    if (pars$zscale$maxsig == 0) {
        return(0)
    }
    decay_fun <- pars$decay_fun
    nn <- ncol(nnpg$nn)
    dsig <- vapply(seq_len(nrow(nnpg$nn)), function(ii) {
        p_idx <- nnpg$nn[ii, ]
        p_dst <- nnpg$dist[ii, ]
        p_sig <- gxy[p_idx, "vsignal"]
        #--- scaling projection on pdist
        pdist <- pars$pdist
        if (pars$vwscale) {
            p_wgt <- gxy[p_idx, "vweight"]
            pdist <- pdist * p_wgt
        }
        #--- project observed signal
        x <- p_dst / (pars$nrc * pdist)
        s <- p_sig / pars$zscale$maxsig
        dsig <- decay_fun(x, s, ... = ...)
        return(dsig)
    }, numeric(nn))
    if(is.matrix(dsig)){
        dsig <- t(dsig)
        Z <- .summ.dsig(dsig, pars)
    } else {
        Z <- dsig
    }
    return(Z)
}

#-------------------------------------------------------------------------------
#--- main polar projection calls
.get.ldsig.polar <- function(lpts, gxy, edges, pars, nnpg, ...) {
    if (pars$zscale$maxsig == 0) {
        return(0)
    }
    decay_fun <- pars$decay_fun
    #--- get polar coords for edges
    edlist <- .edge.list(edges, gxy)
    etheta <- .edge.list.theta(edlist, gxy)
    edleng <- .edge.attr(edges, gxy, "edistNorm")
    if(pars$directional){
        emode <- .edge.mode(edges, gxy)
        for(i in seq_along(emode)){
            idx <- emode[[i]] == 1
            edlist[[i]] <- edlist[[i]][idx]
            etheta[[i]] <- etheta[[i]][idx]
            edleng[[i]] <- edleng[[i]][idx]
        }
    }
    idx <- which(names(etheta) %in% rownames(gxy))
    etheta <- etheta[idx]
    edlist <- edlist[idx]
    edleng <- edleng[idx]
    minlen <- min(unlist(edleng))
    #--- start polar signal projection
    nn <- ncol(nnpg$nn)
    dsig <- vapply(seq_len(nrow(nnpg$nn)), function(ii) {
        p_idx <- nnpg$nn[ii, ]
        p_dst <- nnpg$dist[ii, ]
        p_sig <- gxy[p_idx, "vsignal"]
        #--- scaling projection on pdist
        pdist <- 1 #controlled at edge dist
        if (pars$vwscale) {
            p_wgt <- gxy[p_idx, "vweight"]
            pdist <- pdist * p_wgt
        }
        #--- get edge's theta and dist
        p_et <- etheta[p_idx]
        p_el <- edleng[p_idx]
        #--- compute theta for p1 vs. knn
        dx <- lpts[ii, "X"] - gxy[p_idx, "X"]
        dy <- lpts[ii, "Y"] - gxy[p_idx, "Y"]
        p_theta <- atan2(dy, dx)
        #--- get delta theta (multi directional)
        mdr <- .get.delta.theta(p_et, p_el, p_theta, pars, minlen)
        dth <- mdr[1, ]
        wln <- mdr[2, ]
        #--- update pdist with wln (the signal tends to zero at wln)
        wln <- wln / pars$nrc
        pdist <- pdist * wln
        #--- approx. project. to the area of a sector for isolated nodes
        if(!pars$directional){
            rs <- 0.5 * .deg2rad(pars$theta)
            rs <- sqrt(rs / pi)
            dth[is.na(dth)] <- rs
        }
        #--- update pdist with dtheta
        pdist <- pdist * dth
        #--- project observed signal
        x <- p_dst / (pars$nrc * pdist)
        s <- p_sig / pars$zscale$maxsig
        dsig <- decay_fun(x, s, ... = ...)
        return(dsig)
    }, numeric(nn))
    if(is.matrix(dsig)){
        dsig <- t(dsig)
        Z <- .summ.dsig(dsig, pars)
    } else {
        Z <- dsig
    }
    return(Z)
}
.get.delta.theta <- function(p_et, p_el, p_theta, pars, minlen) {
    nth <- log(0.25, base = (360 - pars$theta) / 360)
    mdr <- vapply(seq_along(p_et), function(i) {
        eth <- p_et[[i]]
        eln <- p_el[[i]]
        if (length(eth) > 0) {
            dth <- abs(p_theta[i] - eth)
            dth <- pi - abs(dth - pi)
            dth <- (dth / pi)^nth
            if (length(eln) > 1) {
                wln <- sum(eln * (dth / sum(dth)))
            } else {
                wln <- eln
            }
            dth <- max(dth)
            res <- c(dth, wln)
        } else {
            res <- c(NA, minlen)
        }
        res
    }, numeric(2))
    return(mdr)
}
.edge.list <- function(edges, gxy) {
    nms <- rownames(gxy)
    el <- lapply(nms, function(nm) {
        e <- c(edges$name2[edges$name1==nm],
            edges$name1[edges$name2==nm])
        which(nms %in% e)
    })
    names(el) <- nms
    return(el)
}
.edge.attr <- function(edges, gxy, eattr = "edistNorm") {
    nms <- rownames(gxy)
    el <- lapply(nms, function(nm) {
        idx1 <- edges$name1 == nm
        att1 <- edges[[eattr]][idx1]
        names(att1) <- edges$name2[idx1]
        idx2 <- edges$name2 == nm
        att2 <- edges[[eattr]][idx2]
        names(att2) <- edges$name1[idx2]
        e <- c(att1, att2)
        e[order(match(names(e), nms))]
    })
    names(el) <- nms
    return(el)
}
.edge.mode <- function(edges, gxy) {
    nms <- rownames(gxy)
    el <- lapply(nms, function(nm) {
        idx1 <- edges$name1 == nm
        att1 <- edges$emode[idx1]
        names(att1) <- edges$name2[idx1]
        att1[] <- as.numeric(att1>0)
        idx2 <- edges$name2 == nm
        att2 <- edges$emode[idx2]
        names(att2) <- edges$name1[idx2]
        att2[] <- as.numeric(att2==2)
        e <- c(att1, att2)
        e[order(match(names(e), nms))]
    })
    names(el) <- nms
    return(el)
}
#-------------------------------------------------------------------------------
#--- summarize projected signals
.summ.dsig <- function(dsig, pars) {
    # sort projected signals
    dsig <- matrix(dsig[order(row(dsig), -abs(dsig))],
        nrow = nrow(dsig), byrow = TRUE)
    if (pars$knn > 1) {
        knn <- min(pars$knn, ncol(dsig))
        dsig <- dsig[, seq_len(knn), drop = FALSE]
        if (pars$zscale$signal.type == "continuous") {
            # stack all signals projected at a given point
            stfun <- function(sig) {
                sw <- abs(sig)
                sum(sig * (sw / sum(sw, na.rm = TRUE)), na.rm = TRUE)
            }
            Z <- apply(dsig, 1, stfun)
        } else {
            # get mean signal
            Z <- apply(dsig, 1, mean)
        }
    } else {
        Z <- dsig[, 1]
    }
    return(Z)
}

################################################################################
### Polar transformation functions
################################################################################

#-------------------------------------------------------------------------------
#--- compute theta between vertices
.edge.list.theta <- function(edlist, gxy) {
    theta <- lapply(seq_along(edlist), function(e1) {
        e2 <- edlist[[e1]]
        if (length(e2) > 0) {
            res <- .edges.theta(e1, e2, gxy)
        } else {
            res <- numeric()
        }
        res
    })
    names(theta) <- rownames(gxy)
    return(theta)
}
.edges.theta <- function(e1, e2, gxy) {
    p1 <- gxy[e1, c("X", "Y")]
    p2 <- gxy[e2, c("X", "Y"), drop = FALSE]
    dx <- p1["X"] - p2[,"X"]
    dy <- p1["Y"] - p2[,"Y"]
    dist <- sqrt(dx^2 + dy^2)
    theta <- ifelse(dist < 1e-07, NA, atan2(dy, dx))
    names(theta) <- rownames(p2)
    return(theta)
}
#-------------------------------------------------------------------------------
.deg2rad <- function(x) {
    x * pi / 180
}
.rad2deg <- function(x) {
    x * 180 / pi
}

################################################################################
### Map graph's silhouette
################################################################################
.silhouetteCircular <- function(pts, verbose = TRUE) {
    pars <- getPathwaySpace(pts, "pars")
    gxy <- getPathwaySpace(pts, "gxy")
    lpts <- .pointsInMatrix(pars$nrc)
    #--- get point-to-vertice dists for silhouettes
    nnbg <- .get.silhouette.dists(lpts, gxy, pars$knn)
    #--- get landscape floor
    xfloor <- array(0, c(pars$nrc, pars$nrc))
    if (pars$baseline > 0 && pars$pdist > 0) {
        xfloor[lpts[, c("Y", "X")]] <- .get.ldfloor(gxy, pars, nnbg)
    }
    xfloor <- .cutfloor(xfloor, pars)
    xfloor[gxy[, c("Yint", "Xint")]] <- 1
    #---return silhouette
    nbg <- sum(xfloor == 1)
    sz <- round(nbg / prod(dim(xfloor)) * 100, 2)
    if(verbose) message("Silhouette: ", sz, "% of the landscape area!")
    #--- update 'pts' with a silhouette matrix
    pts@misc$xfloor <- xfloor
    #--- add xfloor and rescale xsig to the original signal
    xsig <- pts@misc$xsig
    gxyz <- .rescale.landscape(xsig = xsig, pars = pars, xfloor = xfloor)
    pts@gxyz <- gxyz
    return(pts)
}

#-------------------------------------------------------------------------------
.get.ldfloor <- function(gxy, pars, nnbg) {
    if (pars$zscale$maxsig == 0) {
        return(0)
    }
    decay_fun <- pars$decay_fun
    decay_args <- pars$decay_args
    nn <- ncol(nnbg$nn)
    sfloor <- vapply(seq_len(nrow(nnbg$nn)), function(ii) {
        p_dst <- nnbg$dist[ii, ]
        #--- scaling projection on pdist
        pdist <- pars$pdist
        if (pars$vwscale) {
            p_idx <- nnbg$nn[ii, ]
            p_wgt <- gxy[p_idx, "vweight"]
            pdist <- pdist * p_wgt
        }
        #--- get floor of an ideal (max) signal
        x <- p_dst / (pars$nrc * pdist)
        s <- 1 / pars$zscale$maxsig
        sfloor <- base::do.call(
            decay_fun, c(list(x = x, signal = s), decay_args))
        return(sfloor)
    }, numeric(nn))
    sfloor <- t(sfloor)
    sfloor <- matrix(sfloor[order(row(sfloor), -sfloor)],
        nrow = nrow(sfloor), byrow = TRUE)
    sfloor <- sfloor[, seq_len(pars$knn), drop = FALSE]
    Z <- apply(sfloor, 1, mean)
    return(Z)
}

#-------------------------------------------------------------------------------
.cutfloor <- function(xfloor, pars) {
    rg <- range(xfloor, na.rm = TRUE)
    if (rg[1] < rg[2]) {
        bln <- pars$baseline
        xfloor <- xfloor - rg[1]
        xfloor <- xfloor / max(xfloor, na.rm = TRUE)
        x <- xfloor
        x[x < bln] <- 0
        x[x > 0] <- 1
        x <- .fillCavity(x > 0)
        xfloor[x == 0] <- 0
        xfloor[x >= bln] <- 1
    } else {
        if (pars$baseline == 0) {
            xfloor[, ] <- 1
        } else {
            xfloor[, ] <- 0
        }
    }
    xfloor[1, ] <- 0
    xfloor[, 1] <- 0
    xfloor[, ncol(xfloor)] <- 0
    xfloor[nrow(xfloor), ] <- 0
    return(xfloor)
}

################################################################################
### Other accessor functions
################################################################################

#-------------------------------------------------------------------------------
.pointsInMatrix <- function(nrc) {
    d <- c(nrc,nrc)
    pts <- arrayInd(seq_len(prod(d)), d, useNames = TRUE)
    colnames(pts) <- c("Y", "X")
    return(pts)
}

#-------------------------------------------------------------------------------
.remove.duplicated.nodes <- function(gxy) {
    tp1 <- paste0(gxy[, c("X")], ":", gxy[, c("Y")])
    tp2 <- paste0(gxy[, c("Y")], ":", gxy[, c("X")])
    idx <- duplicated(tp1) | duplicated(tp2) | duplicated(rownames(gxy))
    gxy <- gxy[!idx, ]
    return(gxy)
}

#-------------------------------------------------------------------------------
# pairwise vertice-to-vertice dists 
.get.pairwise.dists <- function(gxy){
    dx <- outer(gxy[, "X"], gxy[, "X"], FUN = "-")
    dy <- outer(gxy[, "Y"], gxy[, "Y"], FUN = "-")
    nngg <- sqrt(dx^2 + dy^2)
    return(nngg)
}

#-------------------------------------------------------------------------------
.transpose.and.flip <- function(mtx) {
    mtx <- t(mtx[nrow(mtx):1, ])
    return(mtx)
}

################################################################################
### Map summits for enrichment analysis
################################################################################
.summitMapping <- function(pts, verbose = TRUE, ...) {
    pars <- getPathwaySpace(pts, "pars")
    gxyz <- getPathwaySpace(pts, "gxyz")
    gxy <- getPathwaySpace(pts, "gxy")
    pts@misc$summits <- .find.summits(gxyz = gxyz,
        gxy = gxy, maxset = pars$maxset, minsize = pars$minsize,
        threshold = pars$summit_threshold,
        segm_fun = pars$segm_fun, ... = ...)
    return(pts)
}

#-------------------------------------------------------------------------------
.find.summits <- function(gxyz, gxy, maxset, minsize, threshold,
    segm_fun, ...) {
    #-- get coords
    lpts <- as.matrix(gxy[, c("X", "Y")])
    lpts[, "X"] <- round(lpts[, "X"])
    lpts[, "Y"] <- round(lpts[, "Y"])
    #-- apply th
    smt <- gxyz
    smt[smt < threshold] <- 0
    #-- run watershed
    # smt <- summitWatershed(smt, tolerance=tolerance, ext=1)
    # smt <- base::do.call(segm_fun, c(list(x = smt), pars$segm_arg))
    smt <- segm_fun(smt, ... = ...)
    xx <- .openPxEdges(smt > 0)
    smt[xx == 0] <- 0
    
    #--- retrive itens
    lset <- .retrive_items(smt, lpts)
    
    #--- apply minsize
    len <- unlist(lapply(lset, length))
    len <- which(len < minsize)
    nset <- length(lset)
    if (length(len) > 0) {
        lset <- lset[-len]
        smt[smt %in% len] <- 0
        nset <- nset - length(len)
        smt <- .relabel(smt)
        if (length(lset) > 0)
            names(lset) <- seq_along(lset)
    }
    
    #--- relabel by size and apply maxset
    # smt <- .relabelBySize(smt)
    # smt[smt > maxset] <- 0
    
    #--- relabel by signal and apply maxset
    smt <- .relabelBySignal(gxyz, smt)
    smt[smt > maxset] <- 0
    
    #--- retrive itens
    lset <- .retrive_items(smt, lpts)
    nset <- length(lset)
    
    #--- get summit outlines
    cset <- .findOutlines(smt)
    
    return(list(lset = lset, cset = cset, mset = smt, nset = nset))
}
.retrive_items <- function(smt, lpts){
  nset <- length(table(as.numeric(smt))) - 1
  lset <- lapply(seq_len(nset), function(i) {
    xi <- smt == i
    idx <- xi[lpts[, c("Y", "X")]]
    rownames(lpts)[idx]
  })
  names(lset) <- seq_along(lset)
  return(lset)
}


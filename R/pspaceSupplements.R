################################################################################
### Main constructor of PathwaySpace-class objects
################################################################################
.buildPathwaySpace <- function(gs, nrc = 500, verbose = TRUE) {
  
    #--- get GraphSpace data
    # get pars
    pars <- getGraphSpace(gs, "pars")
    # get nodes
    nodes <- getGraphSpace(gs, "nodes")
    nodes <- nodes[, c("X", "Y", "name", "weight")]
    # get edges
    edges <- getGraphSpace(gs, "edges")
    edges <- edges[, c("vertex1", "vertex2", "name1",
      "name2", "emode", "weight")]
    
    #--- initialize v-signal
    nodes$vsignal <- 0
    vsignal <- nodes$vsignal
    vertex <- nodes$name
    names(vsignal) <- vertex
    
    #--- update pars
    pars$nrc <- nrc
    pars$zscale <- .get.signal.scale(vsignal)
    pars <- .update.zlim(pars)
    
    #--- check weights
    nodes$weight <- .revise.weights(nodes$weight)
    pars$vweight <- .get.wscale(nodes$weight)
    if(nrow(edges)>0){
      edges$weight <- .revise.weights(edges$weight)
      pars$eweight <- .get.wscale(edges$weight)
    } else {
      pars$eweight <- FALSE
    }
    
    #--- rescale coordinates to nrc dimension
    if(verbose) message("Rescaling 'x' and 'y' coordinates...")
    gxy <- .rescale.coord(nodes, nrc)
    
    #--- get distances -not- related to v-signal
    if(nrow(edges) > 0){
      if(verbose) message("Computing distances between connected vertices...")
      edges <- .get.edge.dist(gxy, edges)
    } else {
      edges$edist <- numeric()
    }
    
    #---create a PathwaySpace object
    if(verbose) message("Creating a 'PathwaySpace' object...")
    pnames <- c("Preprocess", "CircularProjection", "PolarProjection",
        "Silhouette", "Summits")
    status <- rep("[ ]", length(pnames))
    names(status) <- pnames
    ps <- new(Class = "PathwaySpace",  
      gspace = gs, vertex = vertex, vsignal = vsignal, 
      nodes = nodes, edges = edges, gxy = gxy, 
      pars = pars, status = status)
    return(ps)
}

################################################################################
### Accessors of the constructor
################################################################################
.updateStatus <- function(ps, name, check = TRUE) {
    ps@status[name] <- ifelse(check, "[x]", "[ ]")
    return(ps)
}
.checkStatus <- function(ps, name) {
    if(name=="Projection"){
        projections <- c("PolarProjection","CircularProjection")
        sts <- any(ps@status[projections] == "[x]")
    } else {
        sts <- ps@status[name] == "[x]"
    }
    sts
}
.removeSummits <- function(ps, verbose=TRUE){
    if (.checkStatus(ps, "Summits")) {
        if(verbose) message("-- available summits removed.")
        i <- which(names(ps@misc) == "summits")
        ps@misc <- ps@misc[-i]
        ps <- .updateStatus(ps, "Summits", FALSE)
    } 
    return(ps)
}
.removeSilhouette <- function(ps, verbose=TRUE){
    if (.checkStatus(ps, "Silhouette")) {
        if(verbose) message("-- available silhouette removed.")
        i <- which(names(ps@misc) == "xfloor")
        ps@misc <- ps@misc[-i]
        ps <- .updateStatus(ps, "Silhouette", FALSE)
        ps@gxyz <- .rescale.landscape(xsig = ps@misc$xsig, pars = ps@pars)
    } 
    return(ps)
}
.get.signal.scale <- function(vsignal) {
    zscale <- list(range = range(vsignal, na.rm = TRUE))
    zscale$maxsig <- max(abs(zscale$range))
    if (.all_binaryValues(vsignal)) {
        # message("binary")
        zscale$signal.type <- "binary"
        zscale$scale.type <- "pos"
        zscale$range <- c(0, 1)
        zscale$scaling <- c(0, 1)
    } else if (zscale$range[1] == zscale$range[2]) {
        # message("found a single value; it will be treated as binary")
        zscale$signal.type <- "binary"
        if (zscale$range[1] >= 0) {
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
.update.zlim <- function(pars, rescale = TRUE){
  if(rescale){
    zlim <- pars$zscale$scaling
  } else {
    mx <- pars$zscale$maxsig
    if(mx==0) mx <- 1
    if(pars$zscale$scale.type=="negpos"){
      zlim <- c(-mx, mx)
    } else if(pars$zscale$scale.type=="neg"){
      zlim <- c(-mx, 0)
    } else {
      zlim <- c(0, mx)
    }
  }
  pars$zlim <- zlim
  return(pars)
}
.get.edge.dist <- function(gxy, edges) {
    dts <- sqrt((gxy[edges$name1, "X"] - gxy[edges$name2, "X"])^2 +
        (gxy[edges$name1, "Y"] - gxy[edges$name2, "Y"])^2)
    edges$edist <- dts
    return(edges)
}
.rescale.coord <- function(nodes, nrc, from = c(0, 1)){
  nodes$X <- scales::rescale(nodes$X, to = c(1, nrc), from = from)
  nodes$Y <- scales::rescale(nodes$Y, to = c(1, nrc), from = from)
  nodes$Xint <- round(nodes$X)
  nodes$Yint <- round(nodes$Y)
  nodes <- nodes[,c("X", "Y", "Xint", "Yint", "vsignal", "weight")]
  nodes <- as.matrix(nodes)
  return(nodes)
}
.revise.weights <- function(weight){
  if (all(is.na(weight))) weight[] <- 1
  weight[is.na(weight)] <- min(weight, na.rm = TRUE)
  if (sd(weight) > 0) {
    weight <- scales::rescale(weight, to = c(0, 1))
  } else {
    weight[] <- 1
  }
  return(weight)
}
.get.wscale <- function(weight){
  if (sd(weight, na.rm = TRUE) > 0) {
    wscale <- TRUE
  } else {
    wscale <- FALSE
  }
  return(wscale)
}

################################################################################
### Main PathwaySpace calls
################################################################################
.circularProjection <- function(ps, verbose = TRUE) {
    if(verbose) message("Using circular projection...")
    pars <- getPathwaySpace(ps, "pars")
    gxy <- getPathwaySpace(ps, "gxy")
    lpts <- .get.points.in.matrix(pars$nrc)
    nnpg <- .get.point.distances(lpts, gxy, pars$proj$k)
    
    # if (pars$wscale) {
    #     if(verbose) message("Scaling projection to vertex weight...")
    # }
    
    #--- compute landscape signal
    if(verbose) message("Running signal convolution...")
    xsig <- array(0, c(pars$nrc, pars$nrc))
    if (nrow(gxy) > 0) {
        xsig[lpts[, c("Y", "X")]] <- .get.ldsig(gxy, pars, nnpg)
    }
    ps@misc$xsig <- xsig
    # image(.transpose.and.flip(xsig))
    #--- rescale xsig to the original signal
    # use silhouette if available
    if (.checkStatus(ps, "Silhouette")){
      xfloor <- ps@misc$xfloor
      gxyz <- .rescale.landscape(xsig = xsig, pars = pars, xfloor = xfloor)
    } else {
      gxyz <- .rescale.landscape(xsig = xsig, pars = pars)
    }
    ps@gxyz <- gxyz
    ps@pars <- .update.zlim(pars, pars$proj$rescale)
    return(ps)
}

#-------------------------------------------------------------------------------
.polarProjection <- function(ps, verbose = TRUE) {
    edges <- getPathwaySpace(ps, "edges")
    if(nrow(edges)==0){
        msg <- paste0("Note, this 'PathwaySpace' object does not ",
            "contain edges.\nThe 'polarProjection' method requires ",
            "at least one edge for projection.")
        stop(msg)
    }
    pars <- getPathwaySpace(ps, "pars")
    if(pars$proj$directional){
        if(pars$is.directed){
            message("Using polar projection on directed graph...")
        } else {
            stop("'directional' used with undirected graph.")
        }
    } else {
        if(verbose) message("Using polar projection on undirected graph...")
    }
    gxy <- getPathwaySpace(ps, "gxy")
    lpts <- .get.points.in.matrix(pars$nrc)
    nnpg <- .get.point.distances(lpts, gxy, pars$proj$k)
    if (pars$vweight) {
        if(verbose) message("Scaling projection to vertex weight...")
    }
    #--- for polar, 'pdist' is applied at edge dist
    rg <- range(edges$edist)
    rg <- rg[1] + ( (rg - rg[1]) * pars$proj$pdist )
    edges$edistNorm <- scales::rescale(edges$edist, to=rg)
    #--- compute landscape signal
    if(verbose) message("Running signal convolution...")
    xsig <- array(0, c(pars$nrc, pars$nrc))
    if (nrow(gxy) > 0) {
        xsig[lpts[, c("Y", "X")]] <- .get.ldsig.polar(lpts,
            gxy, edges, pars, nnpg)
    }
    ps@misc$xsig <- xsig
    # image(.transpose.and.flip(xsig))
    #--- rescale xsig to the original signal
    gxyz <- .rescale.landscape(xsig = xsig, pars = pars)
    ps@gxyz <- gxyz
    ps@pars <- .update.zlim(pars, pars$proj$rescale)
    return(ps)
}

#-------------------------------------------------------------------------------
#--- get point-to-vertices distances for circular and polar projections
.get.point.distances <- function(lpts, gxy, k, min.ksearch = 30){
  gxy <- cbind(gxy, key=seq_len(nrow(gxy)))
  bg <- is.na(gxy[,"vsignal"]) | gxy[,"vsignal"]==0
  if(!all(bg)) gxy <- gxy[!bg, , drop=FALSE]
  if (k == 1) {
    ksearch <- ceiling(nrow(gxy) * 0.01)
    ksearch <- min(max(ksearch, min.ksearch), nrow(gxy))
  } else {
    ksearch <- min(max(k, min.ksearch), nrow(gxy))
  }
  nnpg <- RANN::nn2(gxy[, c("X", "Y")], 
    lpts[, c("X", "Y")], k = ksearch)
  names(nnpg) <- c("nn", "dist")
  nnpg$nn[,] <- as.numeric(gxy[as.numeric(nnpg$nn),"key"])
  return(nnpg)
}

#-------------------------------------------------------------------------------
# If 'rescale = F', it will rescale 'xsig' to the original range
.rescale.landscape <- function(xsig, pars, xfloor = NULL) {
    if (!is.null(xfloor)) xsig[xfloor == 0] <- NA
    bl <- all(range(xsig, na.rm = TRUE) == 0)
    if(bl){
        x <- xsig
    } else {
        if (!pars$proj$rescale) {
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
            #- export 'endpoints' in future versions;
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
    xp[xp < qt[1]] <- qt[1]
    xp[xp > qt[2]] <- qt[2]
    xp <- scales::rescale(xp, to = c(0, to[2]))
    # neg
    pr <- c(endpoints[2], 1 - endpoints[1])
    xn <- xsig
    xn[xn >= 0] <- NA
    qt <- quantile(as.numeric(xn), na.rm = TRUE, probs = pr)
    xn[xn < qt[1]] <- qt[1]
    xn[xn > qt[2]] <- qt[2]
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
    xp[xp < qt[1]] <- qt[1]
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
    xn[xn > qt[2]] <- qt[2]
    xn <- scales::rescale(xn, to = to)
    xsig[!is.na(xn)] <- xn[!is.na(xn)]
    return(xsig)
}

#-------------------------------------------------------------------------------
#--- main circular projection calls
.get.ldsig <- function(gxy, pars, nnpg) {
    if (pars$zscale$maxsig == 0) {
        return(0)
    }
    decay_fun <- pars$proj$decay_fun
    nn <- ncol(nnpg$nn)
    dsig <- vapply(seq_len(nrow(nnpg$nn)), function(ii) {
        p_idx <- nnpg$nn[ii, ]
        p_dst <- nnpg$dist[ii, ]
        p_sig <- gxy[p_idx, "vsignal"]
        #--- scaling projection on pdist
        pdist <- pars$proj$pdist
        if (pars$vweight) {
            p_wgt <- gxy[p_idx, "weight"]
            pdist <- pdist * p_wgt
        }
        #--- project observed signal
        x <- p_dst / (pars$nrc * pdist)
        s <- p_sig / pars$zscale$maxsig
        dsig <- decay_fun(x, s)
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
.get.ldsig.polar <- function(lpts, gxy, edges, pars, nnpg) {
    if (pars$zscale$maxsig == 0) {
        return(0)
    }
    decay_fun <- pars$proj$decay_fun
    #--- get polar coords for edges
    edlist <- .edge.list(edges, gxy)
    etheta <- .edge.list.theta(edlist, gxy)
    edleng <- .edge.attr(edges, gxy, "edistNorm")
    if(pars$proj$directional){
        emode <- .edge.mode(edges, gxy)
        for(i in seq_along(emode)){
            idx <- emode[[i]] == 1
            edlist[[i]] <- edlist[[i]][idx]
            etheta[[i]] <- etheta[[i]][idx]
            edleng[[i]] <- edleng[[i]][idx]
        }
    }
    #--- start polar signal projection
    minlen <- min(unlist(edleng))
    nn <- ncol(nnpg$nn)
    dsig <- vapply(seq_len(nrow(nnpg$nn)), function(ii) {
        p_idx <- nnpg$nn[ii, ]
        p_dst <- nnpg$dist[ii, ]
        p_sig <- gxy[p_idx, "vsignal"]
        #--- scaling projection on pdist
        pdist <- 1 #controlled at edge dist
        if (pars$vweight) {
            p_wgt <- gxy[p_idx, "weight"]
            pdist <- pdist * p_wgt
        }
        #--- get edge's theta and dist
        p_et <- etheta[p_idx]
        p_el <- edleng[p_idx]
        #--- compute theta for p1 vs. k
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
        if(!pars$proj$directional){
            rs <- 0.5 * .deg2rad(pars$proj$theta)
            rs <- sqrt(rs / pi)
            dth[is.na(dth)] <- rs
        } else {
          dth[is.na(dth)] <- 0
        }
        #--- update pdist with dtheta
        pdist <- pdist * dth
        #--- project observed signal
        x <- p_dst / (pars$nrc * pdist)
        s <- p_sig / pars$zscale$maxsig
        dsig <- decay_fun(x, s)
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
    nth <- log(0.25, base = (360 - pars$proj$theta) / 360)
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
  # aggregate signals
  if (pars$proj$k > 1) {
    k <- min(pars$proj$k, ncol(dsig))
    dsig <- dsig[, seq_len(k), drop = FALSE]
    Z <- apply(dsig, 1, pars$proj$aggregate_fun)
  } else {
    Z <- dsig[, 1]
  }
  #Note: NaNs may result from some aggregate functions when signal 
  # values are only 0s
  Z[is.na(Z)] <- 0
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
.silhouetteCircular <- function(ps, verbose = TRUE) {
    pars <- getPathwaySpace(ps, "pars")
    gxy <- getPathwaySpace(ps, "gxy")
    lpts <- .get.points.in.matrix(pars$nrc)
    #--- get point-to-vertices distances for silhouettes
    nnbg <- .get.silhouette.dists(lpts, gxy, pars$silh$k)
    #--- get landscape floor
    xfloor <- array(0, c(pars$nrc, pars$nrc))
    if (pars$silh$pdist > 0) {
        xfloor[lpts[, c("Y", "X")]] <- .get.ldfloor(gxy, pars, nnbg)
    }
    xfloor <- .cutfloor(xfloor, pars)
    xfloor[gxy[, c("Yint", "Xint")]] <- 1
    #---return silhouette
    nbg <- sum(xfloor == 1)
    sz <- round(nbg / prod(dim(xfloor)) * 100, 2)
    if(verbose) message("Silhouette: ", sz, "% of the landscape area!")
    #--- update 'ps' with xfloor
    ps@misc$xfloor <- xfloor
    #--- rescale gxyz
    if (.checkStatus(ps, "Projection")){
      xsig <- ps@misc$xsig
    } else {
      xsig <- array(0, c(pars$nrc, pars$nrc))
    }
    ps@gxyz <- .rescale.landscape(xsig = xsig, pars = pars, xfloor = xfloor)
    return(ps)
}

#-------------------------------------------------------------------------------
#--- get point-to-vertices dists for silhouettes
.get.silhouette.dists <- function(lpts, gxy, k){
  min.ksearch <- min(max(k, 10), nrow(gxy))
  nnbg <- RANN::nn2(gxy[, c("X", "Y")], lpts[, c("X", "Y")],
    k = min.ksearch)
  names(nnbg) <- c("nn", "dist")
  return(nnbg)
}

#-------------------------------------------------------------------------------
.get.ldfloor <- function(gxy, pars, nnbg) {
    decay_fun <- pars$silh$decay_fun
    nn <- ncol(nnbg$nn)
    sfloor <- vapply(seq_len(nrow(nnbg$nn)), function(ii) {
        p_dst <- nnbg$dist[ii, ]
        #--- scaling projection on pdist
        pdist <- pars$silh$pdist
        #--- get floor of an ideal (max) signal
        x <- p_dst / (pars$nrc * pdist)
        s <- 1
        sfloor <- decay_fun(x, s)
        return(sfloor)
    }, numeric(nn))
    sfloor <- t(sfloor)
    sfloor <- matrix(sfloor[order(row(sfloor), -sfloor)],
        nrow = nrow(sfloor), byrow = TRUE)
    sfloor <- sfloor[, seq_len(pars$silh$k), drop = FALSE]
    Z <- apply(sfloor, 1, mean)
    return(Z)
}

#-------------------------------------------------------------------------------
.cutfloor <- function(xfloor, pars) {
    rg <- range(xfloor, na.rm = TRUE)
    if (rg[1] != rg[2]) {
        xfloor <- xfloor - rg[1]
        xfloor <- xfloor / max(xfloor, na.rm = TRUE)
        mask <- xfloor
        mask[mask < pars$silh$baseline] <- 0
        mask[mask > 0] <- 1
        if(pars$silh$fill.cavity) mask <- .fillCavity(mask)
        xfloor[mask == 0] <- 0
        xfloor[mask > 0] <- 1
    } else {
      xfloor[, ] <- 0
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
.get.points.in.matrix <- function(nrc) {
    d <- c(nrc,nrc)
    lpts <- arrayInd(seq_len(prod(d)), d, useNames = TRUE)
    colnames(lpts) <- c("Y", "X")
    return(lpts)
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
.summitMapping <- function(ps, verbose = TRUE, ...) {
    pars <- getPathwaySpace(ps, "pars")
    gxyz <- getPathwaySpace(ps, "gxyz")
    gxy <- getPathwaySpace(ps, "gxy")
    ps@misc$summits <- .find.summits(gxyz = gxyz,
        gxy = gxy, maxset = pars$summit$maxset, 
        minsize = pars$summit$minsize,
        threshold = pars$summit$summit_threshold,
        segm_fun = pars$summit$segm_fun, ...=...)
    return(ps)
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
    # smt <- base::do.call(segm_fun, c(list(x = smt), pars$summit$segm_arg))
    smt <- segm_fun(smt, ...=...)
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


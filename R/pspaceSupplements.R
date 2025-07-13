################################################################################
### Main constructor of PathwaySpace-class objects
################################################################################
.buildPathwaySpace <- function(gs, nrc = 500, verbose = TRUE) {
  
    #--- update pars
    gs@pars$ps$nrc <- nrc
    
    #--- initialize vertex signal
    att <- names(gs_vertex_attr(gs))
    if( ! "signal" %in% att ){
      gs_vertex_attr(gs, "signal") <- 0
    }
    
    #--- initialize vertex function
    att <- names(gs_vertex_attr(gs))
    if( ! "decayFunction" %in% att ){
      gs_vertex_attr(gs, "decayFunction") <- signalDecay()
    }
    
    #--- initialize edge weight
    att <- names(gs_edge_attr(gs))
    if( ! "weight" %in% att ){
      gs_edge_attr(gs, "weight") <- 1
    }
    
    #--- create a PathwaySpace object
    if(verbose) message("Creating a 'PathwaySpace' object...")
    ps <- as(gs, "PathwaySpace")
        
    #--- initialize status
    pnames <- c("Preprocess", "CircularProjection", "PolarProjection",
      "Silhouette", "Summits")
    status <- rep("[ ]", length(pnames))
    names(status) <- pnames
    ps@status <- status
    
    .validate_ps_containers(ps)
    
    return(ps)
}

################################################################################
### Main circularProjection calls
################################################################################
.circularProjection <- function(ps, verbose = TRUE) {
  
    if(verbose) message("Using circular projection...")
    pars <- getPathwaySpace(ps, "pars")
    nodes <- getPathwaySpace(ps, "nodes")
    
    if(verbose) message("Mapping 'x' and 'y' coordinates...")
    gxy <- .rescale_coord(nodes, pars$ps$nrc)
    lpts <- .get_points_in_matrix(pars$ps$nrc)
    nnpg <- .get_near_neighbours(lpts, gxy, nodes$signal, pars$ps$k)
    
    # project signal
    if(verbose) message("Running signal convolution...")
    xsig <- array(0, c(pars$ps$nrc, pars$ps$nrc))
    if (nrow(nodes) > 0) {
        xsig[lpts[, c("Y", "X")]] <- .get_ldsig(nodes, pars, nnpg)
    }
    # image(.transpose_and_flip(xsig))
    
    # add/update projections
    ps@projections$gxy <- gxy
    ps@projections$xsig <- xsig
    ps <- .update_projections(ps)
    return(ps)
}

#-------------------------------------------------------------------------------
.get_ldsig <- function(nodes, pars, nnpg) {
  if (pars$ps$zscale$maxsig == 0) {
    return(0)
  }
  dist <- as.numeric(nnpg$dist)
  nn <- as.numeric(nnpg$nn)
  nnu <- unique(nn)
  nnl <- .find_nn_positions(nn, nnu)
  dsig <- rep(NA, length(dist))
  if(pars$ps$decay$is_default_args){
    for(i in seq_along(nnu)){
      nnli <- nnl[[i]]
      vertex <- nnu[i]
      decay_fun <- nodes$decayFunction[[vertex]]
      signal <- nodes$signal[vertex]
      x <- dist[nnli] / (pars$ps$nrc * pars$ps$pdist)
      signal <- signal / pars$ps$zscale$maxsig
      dsig[nnli] <- decay_fun(x=x, signal=signal)
    }
  } else {
    # Use do.call (slower, but flexible)
    for(i in seq_along(nnu)){
      nnli <- nnl[[i]]
      vertex <- nnu[i]
      decay_fun <- nodes$decayFunction[[vertex]]
      decay_args <- formalArgs(decay_fun)
      args_list <- as.list(nodes[vertex, decay_args, drop=FALSE])
      args_list$x <- dist[nnli] / (pars$ps$nrc * pars$ps$pdist)
      args_list$signal <- args_list$signal / pars$ps$zscale$maxsig
      dsig[nnli] <- do.call(decay_fun, args_list)
    }
  }
  if(ncol(nnpg$nn)>1){
    dsig <- matrix(dsig, nrow = nrow(nnpg$nn), ncol = ncol(nnpg$nn))
    Z <- .summ_dsig(dsig, pars)
  } else {
    Z <- dsig
  }
  return(Z)
}

#-------------------------------------------------------------------------------
# Find positions of unique 'nnu' ids in 'nn'
.find_nn_positions <- function(nn, nnu) {
  .find_positions <- function(i, x, ord, x_sorted){
    int <- findInterval(c(i - 0.5, i + 0.5), x_sorted, 
      left.open = TRUE, rightmost.closed=TRUE, checkSorted=FALSE)
    ord[seq(int[1]+1,int[2])]
  }
  nno <- order(nn)
  nns <- nn[nno]
  lapply(nnu, function(i) .find_positions(i, nn, nno, nns) )
}

#-------------------------------------------------------------------------------
#--- get point-to-vertices distances for circular and polar projections
.get_near_neighbours <- function(lpts, gxy, signal, k, min.ksearch = 30){
  gxy <- cbind(gxy, key=seq_len(nrow(gxy)))
  bg <- is.na(signal) | signal==0
  if(!all(bg)) gxy <- gxy[!bg, , drop=FALSE]
  if (k == 1) {
    ksearch <- ceiling(nrow(gxy) * 0.01)
    ksearch <- min(max(ksearch, min.ksearch), nrow(gxy))
  } else {
    ksearch <- min(max(k, min.ksearch), nrow(gxy))
  }
  nnpg <- RANN::nn2(gxy[, c("X", "Y"), drop=FALSE], 
    lpts[, c("X", "Y"), drop=FALSE], k = ksearch)
  names(nnpg) <- c("nn", "dist")
  nnpg$nn[,] <- as.numeric(gxy[as.numeric(nnpg$nn),"key"])
  return(nnpg)
}

#-------------------------------------------------------------------------------
#--- summarize projected signals
.summ_dsig <- function(dsig, pars) {
  # sort projected signals
  dsig <- matrix(dsig[order(row(dsig), -abs(dsig))],
    nrow = nrow(dsig), byrow = TRUE)
  # aggregate signals
  if (pars$ps$k > 1) {
    k <- min(pars$ps$k, ncol(dsig))
    dsig <- dsig[, seq_len(k), drop = FALSE]
    Z <- apply(dsig, 1, pars$ps$aggregate.fun)
  } else {
    Z <- dsig[, 1]
  }
  # Note: NaNs may result from some aggregate functions when signal 
  # values are only 0s
  Z[is.na(Z)] <- 0
  return(Z)
}

################################################################################
### Main polarProjection calls
################################################################################
.polarProjection <- function(ps, verbose = TRUE) {
  
    nodes <- getPathwaySpace(ps, "nodes")
    edges <- getPathwaySpace(ps, "edges")
    if(nrow(edges)==0){
        msg <- paste0("Note, this 'PathwaySpace' object does not ",
            "contain edges.\nThe 'polarProjection' method requires ",
            "at least one edge for projection.")
        stop(msg)
    }
    pars <- getPathwaySpace(ps, "pars")
    if(pars$ps$directional){
        if(pars$is.directed){
            message("Using polar projection on directed graph...")
        } else {
            stop("'directional' used with undirected graph.")
        }
    } else {
        if(verbose) message("Using polar projection on undirected graph...")
    }

    if(verbose) message("Mapping 'x' and 'y' coordinates...")
    gxy <- .rescale_coord(nodes, pars$ps$nrc)
    lpts <- .get_points_in_matrix(pars$ps$nrc)
    nnpg <- .get_near_neighbours(lpts, gxy, nodes$signal, pars$ps$k)
    
    if(verbose){
      message("Computing distances between connected vertices...")
      if (.is_variable(edges$weight)) {
        # 'weight' is mapped to theta when variable
        if(verbose) message("Scaling projection to edge weight...")
      }
    }
    edges$edist <- .get_edge_dist(edges, gxy)
    # for polar, 'pdist' is scaled to edge dist and polar coordinates
    nnpg <- .add_scaled_pdist(nnpg, lpts, gxy, edges, pars)
    
    # project signal
    if(verbose) message("Running signal convolution...")
    xsig <- array(0, c(pars$ps$nrc, pars$ps$nrc))
    if (nrow(gxy) > 0) {
        xsig[lpts[, c("Y", "X")]] <- .get_ldsig_polar(nodes, pars, nnpg)
    }
    # image(.transpose_and_flip(xsig))
    
    # add/update projections
    ps@projections$gxy <- gxy
    ps@projections$xsig <- xsig
    ps <- .update_projections(ps)
    return(ps)
}

#-------------------------------------------------------------------------------
.get_ldsig_polar <- function(nodes, pars, nnpg) {
  if (pars$ps$zscale$maxsig == 0) {
    return(0)
  }
  pdist <- as.numeric(nnpg$pdist)
  dist <- as.numeric(nnpg$dist)
  nn <- as.numeric(nnpg$nn)
  nnu <- unique(nn)
  nnl <- .find_nn_positions(nn, nnu)
  dsig <- rep(NA, length(dist))
  if(pars$ps$decay$is_default_args){
    for(i in seq_along(nnu)){
      nnli <- nnl[[i]]
      vertex <- nnu[i]
      decay_fun <- nodes$decayFunction[[vertex]]
      signal <- nodes$signal[vertex]
      x <- dist[nnli] / (pars$ps$nrc * pdist[nnli])
      signal <- signal / pars$ps$zscale$maxsig
      dsig[nnli] <- decay_fun(x=x, signal=signal)
    }
  } else {
    for(i in seq_along(nnu)){
      nnli <- nnl[[i]]
      vertex <- nnu[i]
      decay_fun <- nodes$decayFunction[[vertex]]
      decay_args <- formalArgs(decay_fun)
      args_list <- as.list(nodes[vertex, decay_args, drop=FALSE])
      args_list$x <- dist[nnli] / (pars$ps$nrc * pdist[nnli])
      args_list$signal <- args_list$signal / pars$ps$zscale$maxsig
      dsig[nnli] <- do.call(decay_fun, args_list)
    }
  }
  if(ncol(nnpg$nn)>1){
    dsig <- matrix(dsig, nrow = nrow(nnpg$nn), ncol = ncol(nnpg$nn))
    Z <- .summ_dsig(dsig, pars)
  } else {
    Z <- dsig
  }
  return(Z)
}

#-------------------------------------------------------------------------------
# add scaled pdist to polar coordinates
.add_scaled_pdist <- function(nnpg, lpts, gxy, edges, pars){
  # scale edist
  rg <- range(edges$edist)
  rg <- rg[1] + ( (rg - rg[1]) * pars$ps$pdist )
  edges$scaled_eleng <- scales::rescale(edges$edist, to=rg)
  # get polar coords for edges
  edlist <- .edge_list(edges, gxy)
  etheta <- .edge_list_theta(edlist, gxy)
  edleng <- .edge_attr(edges, gxy, "scaled_eleng")
  edwght <- .edge_attr(edges, gxy, "weight")
  if(pars$ps$directional){
    emode <- .edge_mode(edges, gxy)
    for(i in seq_along(emode)){
      idx <- emode[[i]] == 1
      edlist[[i]] <- edlist[[i]][idx]
      etheta[[i]] <- etheta[[i]][idx]
      edleng[[i]] <- edleng[[i]][idx]
      edwght[[i]] <- edwght[[i]][idx]
    }
  }
  # scale pdist to polar coordinates
  minlen <- min(unlist(edleng))
  nc <- ncol(nnpg$nn)
  pdist <- vapply(seq_len(nrow(nnpg$nn)), function(ii) {
    p_idx <- nnpg$nn[ii, ]
    p_dst <- nnpg$dist[ii, ]
    pdist <- 1 #pdist was set to scaled edge length
    # get edge's theta and dist
    p_et <- etheta[p_idx]
    p_el <- edleng[p_idx]
    p_wt <- edwght[p_idx]
    # compute theta for p1 vs. k
    dx <- lpts[ii, "X"] - gxy[p_idx, "X"]
    dy <- lpts[ii, "Y"] - gxy[p_idx, "Y"]
    p_theta <- atan2(dy, dx)
    # get delta theta (multi directional)
    mdr <- .get_delta_theta(p_et, p_el, p_wt, p_theta, pars, minlen)
    dth <- mdr[1, ]
    wln <- mdr[2, ]
    #--- update pdist with wln (the signal tends to zero at wln)
    wln <- wln / pars$ps$nrc
    pdist <- pdist * wln
    # approx. to the area of a sector for isolated nodes
    if(!pars$ps$directional){
      rs <- 0.5 * .deg2rad(pars$ps$theta)
      rs <- sqrt(rs / pi)
      dth[is.na(dth)] <- rs
    } else {
      dth[is.na(dth)] <- 0
    }
    # update pdist with dtheta
    pdist * dth
  }, numeric(nc))
  nnpg$pdist <- t(pdist)
  return(nnpg)
}
.get_delta_theta <- function(p_et, p_el, p_wt, p_theta, pars, minlen) {
    nth <- log(0.25, base = (360 - pars$ps$theta) / 360)
    mdr <- vapply(seq_along(p_et), function(i) {
        eth <- p_et[[i]]
        eln <- p_el[[i]]
        ewt <- p_wt[[i]]
        if (length(eth) > 0) {
            dth <- abs(p_theta[i] - eth)
            dth <- pi - abs(dth - pi)
            dth <- (dth / pi)^nth
            if (length(eln) > 1) {
                wln <- sum(eln * (dth / sum(dth)))
            } else {
                wln <- eln
            }
            dth <- max(dth*ewt)
            res <- c(dth, wln)
        } else {
            res <- c(NA, minlen)
        }
        res
    }, numeric(2))
    return(mdr)
}
.edge_list <- function(edges, gxy) {
    nms <- rownames(gxy)
    el <- lapply(nms, function(nm) {
        e <- c(edges$name2[edges$name1==nm],
            edges$name1[edges$name2==nm])
        which(nms %in% e)
    })
    names(el) <- nms
    return(el)
}
.edge_attr <- function(edges, gxy, eattr = "scaled_eleng") {
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
.edge_mode <- function(edges, gxy) {
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
#--- compute theta between vertices
.edge_list_theta <- function(edlist, gxy) {
  theta <- lapply(seq_along(edlist), function(e1) {
    e2 <- edlist[[e1]]
    if (length(e2) > 0) {
      res <- .edges_theta(e1, e2, gxy)
    } else {
      res <- numeric()
    }
    res
  })
  names(theta) <- rownames(gxy)
  return(theta)
}
.edges_theta <- function(e1, e2, gxy) {
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
### Rescale projections
################################################################################
# If 'rescale = F', it will rescale 'xsig' to the original range
.update_projections <- function(ps) {
  xfloor <- ps@projections$xfloor
  xsig <- ps@projections$xsig
  pars <- ps@pars
  if(is.null(pars$ps$rescale)) pars$ps$rescale <- TRUE
  if (!is.null(xfloor)) xsig[xfloor == 0] <- NA
  bl <- all(range(xsig, na.rm = TRUE) == 0)
  if(bl){
    gxyz <- xsig
  } else {
    if (!pars$ps$rescale) {
      gxyz <- scales::rescale(xsig, to = pars$ps$zscale$range)
    } else if (pars$ps$zscale$signal.type == "binary") {
      gxyz <- scales::rescale(xsig, to = pars$ps$zscale$scaling)
    } else {
      endpoints <- c(0.01, 1e-04)
      if (pars$ps$zscale$scale.type == "negpos") {
        to <- pars$ps$zscale$range/pars$ps$zscale$maxsig
        gxyz <- .rescale_negpos(xsig, pars, to, endpoints = endpoints)
      } else if (pars$ps$zscale$scale.type == "pos") {
        to <- pars$ps$zscale$scaling
        gxyz <- .rescale_pos(xsig, pars, to, endpoints = endpoints)
      } else if (pars$ps$zscale$scale.type == "neg") {
        to <- pars$ps$zscale$scaling
        gxyz <- .rescale_neg(xsig, pars, to, endpoints = endpoints)
      }
      #- export 'endpoints' in future versions;
      #- it will cut noise and outliers at the endpoints;
      #- outliers outside that range are moved to the nearby endpoint;
      #- noise are set to background (i.e. zero).
    }
  }
  pars$ps$configs <- .get_configs(pars)
  ps@projections$pars <- pars
  ps@projections$gxyz <- gxyz
  return(ps)
}
.get_configs <- function(pars){
  zconfig <- pars$ps$zscale
  if(is.null(pars$ps$rescale)){
    zconfig$rescaled <- TRUE
  } else {
    zconfig$rescaled <- pars$ps$rescale
  }
  if(zconfig$rescaled){
    zconfig$zlim <- zconfig$scaling
  } else {
    mx <- zconfig$maxsig
    if(mx==0) mx <- 1
    if(zconfig$scale.type=="negpos"){
      zlim <- c(-mx, mx)
    } else if(zconfig$scale.type=="neg"){
      zlim <- c(-mx, 0)
    } else {
      zlim <- c(0, mx)
    }
    zconfig$zlim <- zlim
  }
  return(zconfig)
}
.rescale_negpos <- function(xsig, pars, to, endpoints) {
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
.rescale_pos <- function(xsig, pars, to, endpoints) {
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
.rescale_neg <- function(xsig, pars, to, endpoints) {
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

################################################################################
### Map silhouette
################################################################################
.silhouetteCircular <- function(ps, verbose = TRUE) {
  
    nodes <- getPathwaySpace(ps, "nodes")
    pars <- getPathwaySpace(ps, "pars")
    
    if(verbose) message("Mapping 'x' and 'y' coordinates...")
    gxy <- .rescale_coord(nodes, pars$ps$nrc)
    lpts <- .get_points_in_matrix(pars$ps$nrc)
    
    #--- get point-to-vertices distances for silhouettes
    nnbg <- .get_silhouette_dists(lpts, gxy, pars$ps$silh$k)
    
    #--- project floor
    xfloor <- array(0, c(pars$ps$nrc, pars$ps$nrc))
    if (nrow(nodes)>0 & pars$ps$silh$pdist > 0) {
        xfloor[lpts[, c("Y", "X")]] <- .get_ldfloor(pars, nnbg)
    }
    xfloor <- .cutfloor(xfloor, pars)
    xfloor[gxy[, c("Yint", "Xint")]] <- 1
    
    #---return silhouette
    nbg <- sum(xfloor == 1)
    sz <- round(nbg / prod(dim(xfloor)) * 100, 2)
    if(verbose) message("Silhouette: ", sz, "% of the landscape area!")
    
    #--- add/update projections
    ps@projections$gxy <- gxy
    ps@projections$xfloor <- xfloor
    if (!.checkStatus(ps, "Projection")){
      ps@projections$xsig <- array(0, c(pars$ps$nrc, pars$ps$nrc))
    }
    ps <- .update_projections(ps)
    return(ps)
}

#-------------------------------------------------------------------------------
#--- get point-to-vertices dists for silhouettes
.get_silhouette_dists <- function(lpts, gxy, k){
  min.ksearch <- min(max(k, 10), nrow(gxy))
  nnbg <- RANN::nn2(gxy[, c("X", "Y")], lpts[, c("X", "Y")],
    k = min.ksearch)
  names(nnbg) <- c("nn", "dist")
  return(nnbg)
}

#-------------------------------------------------------------------------------
.get_ldfloor <- function(pars, nnbg) {
    decay.fun <- pars$ps$silh$decay.fun
    nn <- ncol(nnbg$nn)
    sfloor <- vapply(seq_len(nrow(nnbg$nn)), function(ii) {
        p_dst <- nnbg$dist[ii, ]
        #--- scaling projection on pdist
        pdist <- pars$ps$silh$pdist
        #--- get floor of an ideal (max) signal
        x <- p_dst / (pars$ps$nrc * pdist)
        s <- 1
        sfloor <- decay.fun(x, s)
        return(sfloor)
    }, numeric(nn))
    sfloor <- t(sfloor)
    sfloor <- matrix(sfloor[order(row(sfloor), -sfloor)],
        nrow = nrow(sfloor), byrow = TRUE)
    sfloor <- sfloor[, seq_len(pars$ps$silh$k), drop = FALSE]
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
        mask[mask < pars$ps$silh$baseline] <- 0
        mask[mask > 0] <- 1
        if(pars$ps$silh$fill.cavity) mask <- .fillCavity(mask)
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
### Map summits for enrichment analysis
################################################################################
.summitMapping <- function(ps, verbose = TRUE, ...) {
    pars <- getPathwaySpace(ps, "pars")
    gxy <- ps@projections$gxy
    gxyz <- ps@projections$gxyz
    ps@projections$summits <- .find_summits(
        gxy = gxy, gxyz = gxyz, 
        maxset = pars$ps$summit$maxset, 
        minsize = pars$ps$summit$minsize,
        threshold = pars$ps$summit$summit_threshold,
        segm_fun = pars$ps$summit$segm_fun, ...=...)
    return(ps)
}

#-------------------------------------------------------------------------------
.find_summits <- function(gxy, gxyz, maxset, minsize, threshold,
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
    # smt <- base::do.call(segm_fun, c(list(x = smt), pars$ps$summit$segm_arg))
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

################################################################################
### Diverse internal accessors
################################################################################
.updateStatus <- function(ps, name, check = TRUE) {
  ps@status[name] <- ifelse(check, "[x]", "[ ]")
  return(ps)
}
.checkStatus <- function(ps, name) {
  if(name=="Projection"){
    projections <- c("PolarProjection","CircularProjection")
    sts <- any(ps@status[projections] == "[x]")
    names(sts) <- "Projection"
  } else {
    sts <- ps@status[name] == "[x]"
    names(sts) <- name
  }
  sts[is.na(sts)] <- FALSE
  sts
}
.summariseStatus <- function(ps){
  sts <- c(.checkStatus(ps, "Preprocess"),
    .checkStatus(ps, "Projection"),
    .checkStatus(ps, "Silhouette"),
    .checkStatus(ps, "Summits"))
  sts[] <- ifelse(sts, "[x]", "[ ]")
  sts <- paste(names(sts), sts, collapse = "  ", sep="")
  sts
}
.removeSummits <- function(ps, verbose=TRUE){
  if (.checkStatus(ps, "Summits")) {
    if(verbose) message("-- available summits removed.")
    i <- which(names(ps@projections) == "summits")
    ps@projections <- ps@projections[-i]
    ps <- .updateStatus(ps, "Summits", FALSE)
  } 
  return(ps)
}
.removeSilhouette <- function(ps, verbose=TRUE){
  if (.checkStatus(ps, "Silhouette")) {
    if(verbose) message("-- available silhouette removed.")
    i <- which(names(ps@projections) == "xfloor")
    ps@projections <- ps@projections[-i]
    ps <- .updateStatus(ps, "Silhouette", FALSE)
    ps <- .update_projections(ps)
  }
  return(ps)
}
.get_signal_scale <- function(signal) {
  zscale <- list(range = range(signal, na.rm = TRUE))
  zscale$maxsig <- max(abs(zscale$range))
  if (.all_binaryValues(signal)) {
    # "binary"
    zscale$signal.type <- "binary"
    zscale$scale.type <- "pos"
    zscale$range <- c(0, 1)
    zscale$scaling <- c(0, 1)
  } else if (zscale$range[1] == zscale$range[2]) {
    # when one single value; it will be treated as binary
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
    # continuous in (-Inf,Inf)
    zscale$signal.type <- "continuous"
    zscale$scale.type <- "negpos"
    zscale$scaling <- c(-1, 1)
  } else if (all(zscale$range >= 0)) {
    # continuous >=0
    zscale$signal.type <- "continuous"
    zscale$scale.type <- "pos"
    zscale$scaling <- c(0, 1)
  } else if (all(zscale$range <= 0)) {
    # continuous <=0
    zscale$signal.type <- "continuous"
    zscale$scale.type <- "neg"
    zscale$scaling <- c(-1, 0)
  }
  return(zscale)
}
.rescale_coord <- function(nodes, nrc, from = c(0, 1)){
  gxy <- nodes[, c("x", "y")]
  colnames(gxy) <- c("X", "Y")
  rownames(gxy) <- nodes$name
  gxy$X <- scales::rescale(gxy$X, to = c(1, nrc), from = from)
  gxy$Y <- scales::rescale(gxy$Y, to = c(1, nrc), from = from)
  gxy$Xint <- round(gxy$X)
  gxy$Yint <- round(gxy$Y)
  gxy <- as.matrix(gxy)
  return(gxy)
}
.is_variable <- function(weight){
  if (length(weight) > 1 && sd(weight, na.rm = TRUE) > 0) {
    wscale <- TRUE
  } else {
    wscale <- FALSE
  }
  return(wscale)
}
.get_edge_dist <- function(edges, gxy) {
  dts <- sqrt((gxy[edges$name1, "X"] - gxy[edges$name2, "X"])^2 +
      (gxy[edges$name1, "Y"] - gxy[edges$name2, "Y"])^2)
  return(dts)
}

#-------------------------------------------------------------------------------
.validate_aggregate_fun <- function(aggregate.fun){
  fargs <- formalArgs(args(aggregate.fun))
  fargs <- fargs[ !fargs %in% c("...", "na.rm")]
  if(length(fargs)>1){
    stop("'aggregate.fun' must be a unary function, e.g, function(x) {...}",
      call. = FALSE)
  }
  TRUE
}

#-------------------------------------------------------------------------------
.get_points_in_matrix <- function(nrc) {
  d <- c(nrc,nrc)
  lpts <- arrayInd(seq_len(prod(d)), d, useNames = TRUE)
  colnames(lpts) <- c("Y", "X")
  return(lpts)
}

################################################################################
### Other accessors
################################################################################

#-------------------------------------------------------------------------------
# .remove_duplicated_nodes <- function(gxy) {
#     tp1 <- paste0(gxy[, c("X")], ":", gxy[, c("Y")])
#     tp2 <- paste0(gxy[, c("Y")], ":", gxy[, c("X")])
#     idx <- duplicated(tp1) | duplicated(tp2) | duplicated(rownames(gxy))
#     gxy <- gxy[!idx, ]
#     return(gxy)
# }

#-------------------------------------------------------------------------------
# pairwise vertice-to-vertice dists 
# .get_pairwise_dists <- function(gxy){
#     dx <- outer(gxy[, "X"], gxy[, "X"], FUN = "-")
#     dy <- outer(gxy[, "Y"], gxy[, "Y"], FUN = "-")
#     nngg <- sqrt(dx^2 + dy^2)
#     return(nngg)
# }

#-------------------------------------------------------------------------------
# .transpose_and_flip <- function(mtx) {
#     mtx <- t(mtx[nrow(mtx):1, ])
#     return(mtx)
# }

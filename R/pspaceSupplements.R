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
    
    #--- create a PathwaySpace object
    if(verbose) message("Creating a 'PathwaySpace' object...")
    ps <- as(gs, "PathwaySpace")
        
    #--- initialize status
    pnames <- c("Preprocess", "CircularProjection", "PolarProjection",
      "Silhouette", "Summits")
    status <- rep("[ ]", length(pnames))
    names(status) <- pnames
    ps@status <- status
    
    ps <- .validate_ps_containers(ps)
    
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
  nnpg <- .get_near_points(lpts, gxy, nodes, pars$ps$nrc, pars$ps$pdist)
  
  # project signal
  if(verbose) message("Running signal convolution...")
  xsig <- array(0, c(pars$ps$nrc, pars$ps$nrc))
  if (nrow(nodes) > 0) {
    xsig[lpts[, c("Y", "X")]] <- .get_ldsig(nodes, pars, nnpg, lpts)
  }
  # image(.transpose_and_flip(xsig))
  
  # add/update projections
  ps@projections$gxy <- gxy
  ps@projections$xsig <- xsig
  ps <- .update_projections(ps)
  return(ps)
}

#-------------------------------------------------------------------------------
.get_ldsig <- function(nodes, pars, nnpg, lpts) {
  if (pars$ps$zscale$maxsig == 0) {
    return(0)
  }
  nsig <- list()
  if(pars$ps$decay$is_default_args){
    for(i in seq_len(nrow(nodes))){
      decay_fun <- nodes$decayFunction[[i]]
      signal <- nodes$signal[i]
      x <- nnpg$dist[[i]] / (pars$ps$nrc * pars$ps$pdist)
      signal <- signal / pars$ps$zscale$maxsig
      nsig[[i]] <- decay_fun(x=x, signal=signal)
    }
  } else {
    # Use do.call (slower, but flexible)
    for(i in seq_len(nrow(nodes))){
      decay_fun <- nodes$decayFunction[[i]]
      decay_args <- formalArgs(decay_fun)
      args_list <- as.list(nodes[i, decay_args, drop=FALSE])
      args_list$x <- nnpg$dist[[i]] / (pars$ps$nrc * pars$ps$pdist)
      args_list$signal <- args_list$signal / pars$ps$zscale$maxsig
      nsig[[i]] <- do.call(decay_fun, args_list)
    }
  }
  nsig <- unlist(nsig)
  nn <- unlist(nnpg$nn)
  nnu <- unique(nn)
  nnl <- .find_nnu_positions(nnu, nn)
  n <- max(sapply(nnl, length))
  dsig_mt <- matrix(0, nrow = nrow(lpts), ncol = n)
  for(i in seq_along(nnl)){
    s <- nsig[nnl[[i]]]
    dsig_mt[nnu[i], seq_along(s) ] <- s
  }
  Z <- .summ_dsig_mt(dsig_mt, pars)
  return(Z)
}

#-------------------------------------------------------------------------------
# Find positions of unique 'nnu' in 'nn'
.find_nnu_positions <- function(nnu, nn) {
  .find_positions <- function(i, nn_idx, starts, ends){
    nn_idx[seq(starts[i], ends[i])]
  }
  nn_idx <- order(nn)
  nn_ord <- nn[nn_idx]
  nn_rev <- rev(nn_ord)
  starts <- match(nnu, nn_ord)
  ends <- length(nn_rev) - match(nnu, nn_rev) + 1
  lapply(seq_along(nnu), function(i) {
    .find_positions(i, nn_idx, starts, ends)
  } )
}

#-------------------------------------------------------------------------------
#--- summarize projected signals
.summ_dsig_mt <- function(dsig, pars) {
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

#-------------------------------------------------------------------------------
#--- get vertex-to-point distances for circular projections
# .get_near_neighbours <- function(lpts, gxy, nodes, nrc, pdist){
#   eradius <- .estimate_radius(nodes, pdist)
#   eradius <- nrc * eradius
#   ksearch <- ceiling(pi * max(eradius)^2)
#   nnpg <- RANN::nn2(lpts[, c("X", "Y")], gxy[, c("X", "Y")], k = ksearch)
#   names(nnpg) <- c("nn", "dist")
#   nnpg$nn <- lapply(seq_len(nrow(nnpg$nn)), function(i){nnpg$nn[i,]} )
#   nnpg$dist <- lapply(seq_len(nrow(nnpg$dist)), function(i){nnpg$dist[i,]} )
#   return(nnpg)
# }

#-------------------------------------------------------------------------------
#--- get vertex-to-point distances for circular and polar projections
.get_near_points <- function(lpts, gxy, nodes, nrc, pdist){
  eradius <- .estimate_radius(nodes, pdist)
  eradius <- nrc * eradius
  nnpg <- list(nn = list(), dist = list())
  lpts <- as.data.frame(lpts)
  blocks <- unique(lpts$X)
  first <- match(blocks, lpts$X)
  last  <- length(lpts$X) - match(blocks, rev(lpts$X)) + 1
  blocks <- data.frame(value = blocks, first, last)
  for(i in seq_len(nrow(gxy))){
    r <- eradius[i]
    if(r>0){
      p <- gxy[i, ]
      x_low <- max(ceiling(p["X"] - r), 1)
      x_high <- min(floor(p["X"] + r), nrow(blocks))
      lpts_filt <- lpts[seq(blocks[x_low, "first"], blocks[x_high, "last"]), ]
      d2 <- (lpts_filt$X - p["X"])^2 + (lpts_filt$Y - p["Y"])^2
      within_radius <- which(d2 <= r^2)
      if(length(within_radius)>0){
        nn_dists <- sqrt(d2[within_radius])
        nearby_idx <- lpts_filt[within_radius, ]
        nn_idx <- with(nearby_idx, (X - 1) * nrc + Y)
      } else {
        nn_idx <- numeric()
        nn_dists <- numeric()
      }
    } else {
      nn_idx <- numeric()
      nn_dists <- numeric()
    }
    nnpg$nn[[i]] <- nn_idx
    nnpg$dist[[i]] <- nn_dists
  }
  return(nnpg)
}
.estimate_radius <- function(nodes, pdist, decay_end = 0.01){
  maxlen <- 1/pdist
  x <- seq(maxlen, 0, length.out=100)
  fun <- nodes$decayFunction
  signal <- abs(nodes$signal)
  signal[is.na(signal)] <- 0
  signal <- signal/max(signal)
  eradius <- vapply(seq_along(fun), function(i){
    f <- fun[[i]]
    s <- signal[i]
    if(s==0){
      r <- 0
    } else {
      xs <- vapply(x, function(l){ f(l, s) }, numeric(1))
      idx <- findInterval(decay_end, xs)
      r <- ifelse(length(idx)>0, x[idx[1]], x[1])
    }
    return(r)
  }, numeric(1))
  eradius <- eradius/maxlen
  return(eradius)
}

################################################################################
### Main polarProjection calls
################################################################################
.polarProjection <- function(ps, verbose = TRUE) {
  
    nodes <- getPathwaySpace(ps, "nodes")
    edges <- getPathwaySpace(ps, "edges")
    # node_degree <- degree(ps@graph)
    # node_degree <- node_degree[rownames(nodes)]
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
    
    if(verbose){
      message("Computing distances between connected vertices...")
    }
    # for polar, 'pdist' is scaled to edge dist and polar coordinates
    gxy <- .rescale_coord(nodes, pars$ps$nrc)
    edges$edist <- .get_edge_dist(edges, gxy)
    
    # adjust pdist to search a larger projection area
    if(pars$ps$pdist.anchor=="frame"){
      pars$ps$psearch <- min(1, pars$ps$pdist * 1.1)
    } else {
      # get psearch within edge bounds
      pars$ps$psearch <- max(edges$edist)/pars$ps$nrc
    }
    
    if(verbose) message("Computing point-to-vertex linear distances...")
    lpts <- .get_points_in_matrix(pars$ps$nrc)
    nnpg <- .get_near_points(lpts, gxy, nodes, pars$ps$nrc, pars$ps$pdist)
    if(verbose) message("Computing point-to-edge angular distances...")
    nnpg <- .get_angular_dist(nnpg, lpts, gxy, edges, pars)
    nnpg <- .scale_pdist_polar(nnpg, pars)
    
    # project signal
    if(verbose) message("Running signal convolution...")
    xsig <- array(0, c(pars$ps$nrc, pars$ps$nrc))
    if (nrow(gxy) > 0) {
      xsig[lpts[, c("Y", "X")]] <- .get_ldsig_polar(nodes, pars, nnpg, lpts)
    }
    # image(.transpose_and_flip(xsig))
    
    # add/update projections
    ps@projections$gxy <- gxy
    ps@projections$xsig <- xsig
    ps <- .update_projections(ps)
    return(ps)
}

#-------------------------------------------------------------------------------
.get_ldsig_polar <- function(nodes, pars, nnpg, lpts) {
  if (pars$ps$zscale$maxsig == 0) {
    return(0)
  }
  nsig <- list()
  if(pars$ps$decay$is_default_args){
    for(i in seq_len(nrow(nodes))){
      decay_fun <- nodes$decayFunction[[i]]
      signal <- nodes$signal[i]
      x <- nnpg$dist[[i]] / (nnpg$pdist_dth[[i]] * pars$ps$nrc)
      signal <- signal / pars$ps$zscale$maxsig
      nsig[[i]] <- decay_fun(x=x, signal=signal)
    }
  } else {
    # Use do.call (slower, but flexible)
    for(i in seq_len(nrow(nodes))){
      decay_fun <- nodes$decayFunction[[i]]
      decay_args <- formalArgs(decay_fun)
      args_list <- as.list(nodes[i, decay_args, drop=FALSE])
      args_list$x <- nnpg$dist[[i]] / (nnpg$pdist_dth[[i]] * pars$ps$nrc)
      args_list$signal <- args_list$signal / pars$ps$zscale$maxsig
      nsig[[i]] <- do.call(decay_fun, args_list)
    }
  }
  nsig <- unlist(nsig)
  nn <- unlist(nnpg$nn)
  nnu <- unique(nn)
  nnl <- .find_nnu_positions(nnu, nn)
  n <- max(sapply(nnl, length))
  dsig_mt <- matrix(0, nrow = nrow(lpts), ncol = n)
  for(i in seq_along(nnl)){
    s <- nsig[nnl[[i]]]
    dsig_mt[nnu[i], seq_along(s) ] <- s
  }
  Z <- .summ_dsig_mt(dsig_mt, pars)
  return(Z)
}

#-------------------------------------------------------------------------------
.get_angular_dist <- function(nnpg, lpts, gxy, edges, pars){
  # scale edist ..only used when pdist.anchor = edge
  rg <- range(edges$edist)
  rg <- rg[1] + ( (rg - rg[1]) * pars$ps$pdist )
  edges$scaled_eleng <- scales::rescale(edges$edist, to=rg)
  # get polar coords
  edlist <- .edge_list(edges, gxy)
  etheta <- .edge_list_theta(edlist, gxy)
  edleng <- .edge_attr(edges, gxy, "scaled_eleng")
  if(pars$ps$directional){
    emode <- .edge_mode(edges, gxy)
    for(i in seq_along(emode)){
      idx <- emode[[i]] == 1
      edlist[[i]] <- edlist[[i]][idx]
      etheta[[i]] <- etheta[[i]][idx]
      edleng[[i]] <- edleng[[i]][idx]
    }
  }
  # get normalized point-to-edge angular distances
  minlen <- min(unlist(edleng))
  for(i in seq_len(length(nnpg$nn))){
    p_idx <- nnpg$nn[[i]]
    p_dst <- nnpg$dist[[i]]
    # get edge's theta and dist
    p_et <- etheta[[i]]
    p_el <- edleng[[i]]
    # compute theta for p vs. i
    dx <- lpts[p_idx, "X"] - gxy[i, "X"]
    dy <- lpts[p_idx, "Y"] - gxy[i, "Y"]
    p_theta <- atan2(dy, dx)
    # get delta theta (multi directional)
    if(length(p_et)>0){
      mdr <- .get_delta_theta(p_theta, p_et, p_el, pars)
      dth <- mdr[, 1]
      wln <- mdr[, 2]
    } else {
      if(pars$ps$directional){
        dth <- 0
      } else {
        # for an isolated node 'a', dth will aim the area of a sector
        # pdist correction by dth: pdist_dth = pdist * dth^beta
        # beta = log(0.25) / log((360 - theta) / 360)
        # pdist * sqrt( theta/360 ) ==> pdist_dth
        theta <- 360 - 360 * exp( log(0.25)/pars$ps$beta )
        dth <- (theta/360) ^ (1/(2*pars$ps$beta))
      }
      dth <- rep(dth, length(p_dst))
      wln <- rep(minlen, length(p_dst))
    }
    nnpg$dtheta[[i]] <- dth
    nnpg$edleng[[i]] <- wln
  }
  return(nnpg)
}
.get_delta_theta <- function(p_theta, p_et, p_el, pars) {
  # dth: normalized point-to-edge angular distances
  # wln: aggregate p_el by dth
  mdr <- t(vapply(seq_along(p_theta), function(i) {
    dth <- abs(p_theta[i] - p_et)
    dth <- abs(dth - pi)/pi
    if (length(p_el) > 1) {
      wln <- sum(p_el * (dth / sum(dth)))
    } else {
      wln <- p_el
    }
    dth <- max(dth)
    c(dth, wln)
  }, numeric(2)))
  return(mdr)
}

#-------------------------------------------------------------------------------
.scale_pdist_polar <- function(nnpg, pars) {
  pfun <- pars$ps$polar.fun
  for(i in seq_len(length(nnpg$nn))){
    if(pars$ps$pdist.anchor=="frame"){
      pdist <- pars$ps$pdist
    } else {
      pdist <- nnpg$edleng[[i]] / pars$ps$nrc
    }
    nnpg$pdist_dth[[i]] <- pdist * pfun(nnpg$dtheta[[i]], pars$ps$beta)
  }
  return(nnpg)
}

#-------------------------------------------------------------------------------
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
  dx <- p2[,"X"] - p1["X"]
  dy <- p2[,"Y"] - p1["Y"]
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
    
    #--- get point-to-vertex distances for silhouettes
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
#--- get point-to-vertex distances for silhouettes
.get_silhouette_dists <- function(lpts, gxy, k = 10){
  min.ksearch <- min(max(k, 10), nrow(gxy))
  nnbg <- RANN::nn2(gxy[, c("X", "Y")], lpts[, c("X", "Y")],
    k = min.ksearch)
  names(nnbg) <- c("nn", "dist")
  return(nnbg)
}

#-------------------------------------------------------------------------------
.get_ldfloor <- function(pars, nnbg) {
    decay_fun <- pars$ps$silh$decay.fun
    nn <- ncol(nnbg$nn)
    sfloor <- vapply(seq_len(nrow(nnbg$nn)), function(ii) {
        p_dst <- nnbg$dist[ii, ]
        #--- scaling projection on pdist
        pdist <- pars$ps$silh$pdist
        #--- get floor of an ideal (max) signal
        x <- p_dst / (pars$ps$nrc * pdist)
        sfloor <- decay_fun(x=x, signal=1)
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
.validate_polar_fun <- function(polar.fun){
  fargs <- formalArgs(args(polar.fun))
  fargs <- fargs[ !fargs %in% c("x", "beta")]
  if(length(fargs)>1){
    stop("the 'polar.fun' must have the signature: function(x, beta) {...}",
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

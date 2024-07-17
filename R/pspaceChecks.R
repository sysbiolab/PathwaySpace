
#-------------------------------------------------------------------------------
# Validate igraph for PathwaySpace
.validate.igraph <- function(g, full.check=TRUE) {
    if (!is(g, "igraph")) {
        stop("'g' should be an 'igraph' class object.", call. = FALSE)
    }
    #--- checks
    if (igraph::vcount(g)<2) {
        stop("'g' should have at least two vertices.", call. = FALSE
        ) 
    }
    if(full.check){
      X <- igraph::V(g)$x
      Y <- igraph::V(g)$y
      if (is.null(X) || is.null(Y)) {
        stop("'x' and 'y' vertex attributes should be available in 'g'.",
          call. = FALSE
        ) 
      }
      if (!.is_numericVector(X) || !.is_numericVector(Y)) {
        msg <- paste0("vertex attributes 'x' and 'y' should be ",
          "numeric vectors.")
        stop(msg, call. = FALSE
        )
      }
    }
    vertex <- igraph::V(g)$name
    if (is.null(vertex)) {
        stop("'name' vertex attribute should be available in 'g'.",
            call. = FALSE
        )
    }
    if (!.all_characterValues(vertex)) {
        msg <- paste0("vertex attribute 'name' should be a ",
            "character vector.")
        stop(msg, call. = FALSE
        )
    }
    if (anyDuplicated(vertex)>0) {
        msg <- paste0("vertex attribute 'name' should not contain ",
            "duplicated names.")
        stop(msg, call. = FALSE
        )
    }
    if (igraph::ecount(g)==0) {
        warning("The input 'g' graph does not have any edges.", call. = FALSE)
    }
}
.validate.vertex.signal <- function(g, verbose = TRUE){
    vsignal <- igraph::V(g)$signal
    vweight <- igraph::V(g)$weight
    if (!is.null(vsignal)) {
        if(verbose) message("Extracting 'signal' from 'g'...", appendLF=FALSE)
        if(.is_numericVector(vsignal)){
            if(verbose) message("OK!")
        } else {
            if(verbose) message(" not valid (need a numeric vector)")
            igraph::V(g)$signal <- 0
        }
    } else {
        igraph::V(g)$signal <- 0
    }
    if (!is.null(vweight)) {
        if(verbose) message("Extracting 'weight' from 'g'...", appendLF=FALSE)
        if(.is_numericVector(vweight)){
            if(verbose) message("OK!")
        } else {
            if(verbose) message(" not valid (need a numeric vector)")
            igraph::V(g)$weight <- 1
        }
    } else {
        igraph::V(g)$weight <- 1
    }
    return(g)
}
#-------------------------------------------------------------------------------
.validate.args <- function(check, name, para) {
    if (check == "numeric_vec") {
        msg <- paste0("'", name, "' should be a numeric vector.")
        if (!is.numeric(para) || !is.vector(para)) stop(msg, call. = FALSE)
    } else if (check == "numeric_mtx") {
        msg <- paste0("'", name, "' should be a numeric matrix.")
        if (!is.numeric(para) || !is.matrix(para)) 
          stop(msg, call. = FALSE)
    } else if (check == "square_numeric_mtx") {
      msg <- paste0("'", name, "' should be a named numeric square matrix.")
      if (!is.numeric(para) || !is.matrix(para)) 
        stop(msg, call. = FALSE)
      if(ncol(para) != nrow(para)) 
        stop(msg, call. = FALSE)
      if(is.null(colnames(para)) || is.null(rownames(para))) 
        stop(msg, call. = FALSE)
      if(any(colnames(para) != rownames(para)))
        stop(msg, call. = FALSE)
    } else if (check == "allCharacter") {
        msg <- paste0("'", name, "' should be a vector with strings.")
        if (!.all_characterValues(para)) stop(msg, call. = FALSE)
    } else if (check == "allBinary") {
        msg <- paste0("'", name, "' should be a vector with binary values.")
        if (!.all_binaryValues(para)) stop(msg, call. = FALSE)
    } else if (check == "allInteger") {
        msg <- paste0("'", name, "' should be a vector with integer values.")
        if (!.all_integerValues(para)) stop(msg, call. = FALSE)
    } else if (check == "singleString") {
        msg <- paste0("'", name, "' should be a single string.")
        if (!.is_singleString(para)) stop(msg, call. = FALSE)
    } else if (check == "singleInteger") {
        msg <- paste0("'", name, "' should be a single integer value.")
        if (!.is_singleInteger(para)) stop(msg, call. = FALSE)
    } else if (check == "singleNumber") {
        msg <- paste0("'", name, "' should be a single numeric value.")
        if (!.is_singleNumber(para)) stop(msg, call. = FALSE)
    } else if (check == "function") {
        msg <- paste0("'", name, "' should be a function.")
        if (!is.function(para)) stop(msg, call. = FALSE)
    } else if (check == "singleLogical") {
        msg <- paste0("'", name, "' should be a single logical value.")
        if (!.is_singleLogical(para)) stop(msg, call. = FALSE)
    }
}

#-------------------------------------------------------------------------------
.validate.colors <- function(check, name, para) {
    if (name == "singleColor") {
        if (!.is_singleColor(para)) {
            msg <- paste0("'", name, "' should be a single color.")
            stop(msg, call. = FALSE)
        }
    } else if (name == "allColors") {
        if (!.is_color(para)) {
            msg <- paste0("'", name, "' should be a vector with colors.")
            stop(msg, call. = FALSE)
        }
    }
}

#-------------------------------------------------------------------------------
# Validate custom plot args
.validate.plot.args <- function(name, para) {
    if (name == "colors") {
        if (length(para) > 1) {
            if (!.is_color(para)) {
                msg <- paste0(
                    "'", name,
                    "' arg has invalid color specification."
                )
                stop(msg, call. = FALSE)
            }
        } else {
            opt <- c("pspaceRedBluPal")
            if (!para %in% opt) {
                stop("invalid color argument.", call. = FALSE)
            }
        }
    } else if (name == "trim.colors") {
        if (!is.null(para)){
            if (!.all_integerValues(para)) {
                msg <- paste0(
                    "'", name,
                    "' should be a vector with integer values."
                )
                stop(msg, call. = FALSE)
            } else if (length(para) != 5) {
                msg <- paste0("'", name, "' should be a vector of length 5.")
                stop(msg, call. = FALSE)
            }
        }
    } else if (name == "bg.color") {
        if (!is.null(para)) {
            if (!.is_singleColor(para)) {
                msg <- paste0("'", name, "' should be a single color.")
                stop(msg, call. = FALSE)
            }
        }
    } else if (name == "marks") {
        if (!.is_singleLogical(para)) {
            if (!.all_characterValues(para)) {
                msg <- paste0("'", name, "' should be a vector with strings.")
                stop(msg, call. = FALSE)
            }
        }
    }
}

#-------------------------------------------------------------------
.is_singleNumber <- function(para) {
    (is.integer(para) || is.numeric(para)) &&
        length(para) == 1L && !is.na(para)
}
.is_singleInteger <- function(para) {
    lg <- (is.integer(para) || is.numeric(para)) &&
        length(para) == 1L && !is.na(para)
    if (lg) {
        lg <- ((para + 1) / (ceiling(para) + 1)) == 1
    }
    return(lg)
}
.is_singleString <- function(para) {
    is.character(para) && length(para) == 1L && !is.na(para)
}
.is_singleLogical <- function(para) {
    is.logical(para) && length(para) == 1L && !is.na(para)
}
.all_binaryValues <- function(para) {
    all(para %in% c(0, 1, NA))
}
.all_integerValues <- function(para) {
    para <- abs(para)
    lg <- (all(is.integer(para)) || all(is.numeric(para))) &&
        !any(is.na(para))
    if (lg) {
        lg <- all(((para + 1) / (ceiling(para) + 1)) == 1)
    }
    return(lg)
}
.all_characterValues <- function(para) {
    all(is.character(para)) && !any(is.na(para))
}
.is_numericVector <- function(para){
    is.vector(para) && .all_numericValues(para)
}
.all_numericValues <- function(para) {
    (all(is.numeric(para)) || all(is.integer(para))) && 
        !any(is.na(para))
}
.is_color <- function(x) {
    res <- try(col2rgb(x), silent = TRUE)
    return(!"try-error" %in% class(res))
}
.is_singleColor <- function(para) {
    .is_color(para) && length(para) == 1L && !is.na(para)
}

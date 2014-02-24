

index2xy.CdfEnvAffy <- function(object, i) {
  indices2xy(i, nc = object@nrow) - getOption("BioC")$affy$xy.offset
}
xy2index.CdfEnvAffy <- function(object, x, y) {
  o <- getOption("BioC")$affy$xy.offset
  xy2indices(x+o, y+o, nc = object@nrow)
}

setClass("CdfEnvAffy",
         representation(envir = "environment",
                        envName = "character",
                        index2xy = "function",
                        xy2index = "function",
                        nrow = "integer",
                        ncol = "integer",
                        probeTypes = "character",
                        chipType = "character"),
         prototype = list(index2xy = index2xy.CdfEnvAffy,
           xy2index = xy2index.CdfEnvAffy))


## ---
setAs("CdfEnvAffy", "environment", function(from) from@envir )

setAs("CdfEnvAffy", "Cdf",
      function(from) {
        m <- matrix(as.numeric(NA), from@nrow, from@ncol)
        l <- ls(from@envir)
        for (i in seq(along=l)) {
          tmp <- indexProbes(from, "pm", l[i])[[1]]
          m[tmp] <- i
          tmp <- indexProbes(from, "mm", l[i])[[1]]
          m[tmp] <- i
        }
        cdf <- new("Cdf", cdfName=from@chipType, name=m, name.levels=l)
        return(cdf)
      }
      )


geneNames.CdfEnvAffy <- function(object) {
  ls(as(object, "environment"))
}
setMethod("geneNames", "CdfEnvAffy", geneNames.CdfEnvAffy)

## ---
setMethod("[", signature(x="CdfEnvAffy", i="character",
                         j="missing", drop="missing"),
function(x, i, j, drop=FALSE) {
  if( !missing(j)) {
    stop("Improper subsetting. Only one vector of IDs should be given.\n")
  }
  
  if (is.matrix(i)) {
    if (! is.integer(i)) {
      stop("not implemented")
         ##stop("when a matrix, 'i' should be of mode 'integer'")
    }
    y <- x
    y@envName <- paste(x@envName, "-subsetXYcoords", sep="")
    y@envir <- new.env(hash=TRUE, parent = emptyenv())
    
    ## make a Cdf (faster lookup for XY or indexes).
    cdfenv <- get(x@envName)
    cdf <- as(cdfenv, "Cdf")
    idx <- xy2index(x, i)
    for (i in idx) {
      id <- cdf@names.level[cdf@names[i]]
      tmp <- indexProbes(y, y@probeTypes, id)
      ## implementation not complete
      ## pm or mm to be sorted
      ## and idx appended to tmp
      assign(id, tmp, envir=y@envir)
    }
    
  } else {
    if (! is.character(i)) {
      stop("when not a matrix, 'i' should be of mode 'character'")
    }
    y <- x
    y@envName <- paste(x@envName, "-subsetProbeSets", sep="")
    y@envir <- new.env(hash=TRUE, parent = emptyenv())
    

    for (id in i) {
      ##tmp <- indexProbes(x, x@probeTypes, id)
      tmp <- do.call(cbind, lapply(x@probeTypes, function(pt) indexProbes(x, pt, id)[[1]]))
      assign(id, tmp, envir=y@envir)
    }
  }
  return(y)
})



## ---
indexProbes.CdfEnvAffy <- function(object, which, probeSetNames=NULL) {

  probeTypes <- object@probeTypes

  ##FIXME: hack for compatibility with 'affy'
  if (identical(which, "both") || identical(which, c("pm", "mm", "both"))) {
    which <- probeTypes
    warning("The use of \"both\" in 'which' is deprecated.")
  }
  ##

  if ( ! all(which %in% probeTypes))
    stop(paste("'which' can only take values from:", paste(probeTypes, collapse=", ")))

  i.probes <- match(which, probeTypes)

  envir <- as(object, "environment")

  if(is.null(probeSetNames))
    probeSetNames <- ls(envir)

  ans <-  mget(probeSetNames, envir=envir, ifnotfound=list(NA))

  ## this kind of thing could be included in 'multiget' as
  ## an extra feature. A function could be specified to
  ## process what is 'multi'-get on the fly
  for (i in seq(along=ans)) {

    if ( is.na(ans[[i]][1]) )
      next

    ##as.vector cause it might be a matrix if all probe types
    tmp <- as.vector(ans[[i]][, i.probes])

    ans[[i]] <- tmp
  }

  return(ans)
}

setMethod("indexProbes", signature("CdfEnvAffy", which = "character"),
         indexProbes.CdfEnvAffy)

## ---
setGeneric("index2xy", def = function(object, ...) standardGeneric("index2xy"))
setMethod("index2xy", signature(object="CdfEnvAffy"),
         function(object, ...) object@index2xy(object, ...))

## ---
setGeneric("xy2index", def = function(object, ...) standardGeneric("xy2index"), useAsDefault = FALSE)
setMethod("xy2index", signature(object="CdfEnvAffy"),
          function(object, ...) object@xy2index(object, ...))

## ---

plot.CdfEnvAffy <- function(x, xlab = "", ylab = "", main = x@chipType, ...) {
  plot(0, 0, xlim = range(0, x@nrow), ylim = range(0, x@ncol), type="n", xlab = xlab, ylab = ylab, main = main, ...)
}

setMethod("plot", signature(x="CdfEnvAffy", y="missing"), plot.CdfEnvAffy)

## ---

setMethod("show", signature("CdfEnvAffy"),
          function(object) {
            cat("Instance of class CdfEnvAffy:\n")
            cat(" name     :", object@envName, "\n")
            cat(" chip-type:", object@chipType, "\n")
            cat(" size     :", object@nrow, "x", object@ncol, "\n")
            cat("", length(ls(as(object, "environment"))), "probe set(s) defined.\n")
          })

## ---

validCdfEnvAffy <- function(cdfenv, verbose=TRUE) {

  if (verbose)
    cat("Validating CdfEnvAffy:\n")
  envir <- as(cdfenv, "environment")
  keys <- ls(envir)

  ## probe types
  if (verbose)
    cat("  Checking probe types.\n")
  n <- length(cdfenv@probeTypes)
  tmp <- rep(FALSE, n)
  for (i in seq(along=keys)) {
    if (ncol(get(keys[i], envir = envir)) != n)
      tmp[i] <- TRUE
  }
  if (n > 0 && sum(tmp) != 0)
    valid <- FALSE
  else
    valid <- TRUE
  r.probeTypes <- list(valid=valid, invalid.ones=keys[which(tmp)])
  if (verbose)
    cat(sum(tmp), "invalid ones.\n")
  ## XY
  if (verbose)
    cat("  Checking XY coordinates.\n")
  tmp <- rep(FALSE, n)
  for (i in seq(along=keys)) {
    ip <- indexProbes(cdfenv, which = cdfenv@probeTypes, probeSetNames = keys[i])[[1]]
    xy <- index2xy(cdfenv, ip)
    if (any(xy[, 1] > cdfenv@nrow, na.rm = TRUE) || any(xy[, 2] > cdfenv@ncol, na.rm = TRUE))
      tmp[i] <- TRUE
  }
  if (n > 0 && sum(tmp) != 0)
    valid <- FALSE
  else
    valid <- TRUE
  r.xy <- list(valid=valid, invalid.ones=keys[which(tmp)])
  if (verbose)
    cat(sum(tmp), "invalid ones\n")
  
  r.details <- list(probeTypes=r.probeTypes, xy=r.xy)

  r <- all( unlist(lapply(r.details, function(x) x$valid)) )
  attr(r, "details") <- r.details
  return(r)
}

## ---

printValidCdfEnvAffy <- function(x) {
  printDetails <- function(y) {
    if (y$valid) {
      cat("  valid.\n")
      return()
    }
    n <- length(y$invalid.ones)
    if (n == 1)
      cat(paste(" ", n, "invalid probe set\n"))
    else
      cat(paste(" ", n, "invalid probe sets\n"))
    if (n <= 5)
      cat(paste(y$invalid.ones, collapse=" "), "\n")
    else
      cat(paste(paste(y$invalid.ones[1:5], collapse=" "), "...\n"))
  }

  r.details <- attr(x, "details")

  cat("Probe types:\n")
  printDetails(r.details$probeTypes)

  cat("XY coordinates:\n")
  printDetails(r.details$xy)

}

## ---

validAffyBatch <- function(abatch, cdfenv) {
  stopifnot(is(abatch, "AffyBatch"),
            is(cdfenv, "CdfEnvAffy"))

  if ( (abatch@nrow != cdfenv@nrow) || (abatch@ncol != cdfenv@ncol))
    valid <- FALSE
  else
    valid <- TRUE

  return(valid)
    #r.dim <- list(valid = valid)

  #return(r.dim)
}


## ---

# setMethod("initialize", "CdfEnvAffy",
#           function(.Object) {

#           })

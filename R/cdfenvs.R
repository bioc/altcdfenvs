wrapCdfEnvAffy <- function(cdfenv, nrow.chip, ncol.chip, chiptype) {
  object <- new("CdfEnvAffy", envir = cdfenv,
                nrow = as.integer(nrow.chip), ncol = as.integer(ncol.chip),
                probeTypes = c("pm", "mm"),
                chipType = chiptype)
  valid <- validCdfEnvAffy(object, verbose=FALSE)
  if ( ! valid ) {
    printValidCdfEnvAffy(valid)
    stop("invalid CdfEnvAffy")
  }
  return(object)
#   class(cdfenv) <- c("environment", "CdfEnv")
#   attr(cdfenv, "nrow") <- nrow.chip
#   attr(cdfenv, "ncol") <- ncol.chip
#   attr(cdfenv, "chipType") <- chiptype
#   return(cdfenv)
}

getCdfEnvAffy <- function(abatch) {
  if (class(abatch) != "AffyBatch")
    stop("arg must be of class 'AffyBatch'.")

  cdfenv <- getCdfInfo(abatch)
  
  cdfenv <- wrapCdfEnvAffy(cdfenv, abatch@nrow, abatch@ncol, abatch@cdfName)
  return(cdfenv)
}

buildCdfEnv.matchprobes <- function(matches, ids, probes.pack, abatch=NULL, nrow.chip=NULL, ncol.chip=NULL, chiptype=NULL, mm=NA, simplify = TRUE, x.colname = "x", y.colname = "y") {

  if (! (is.list(matches) && length(matches) > 0) && length(matches[[1]]) < 3)
    stop("arg 'matches' should be a list a returned by matchprobes.")

  if (length(matches[[1]]) != length(ids))
    stop("'matches' and 'ids' must have the same length.")

  if ( ! is.null(abatch)) {
    if (class(abatch) != "AffyBatch")
      stop("abatch must be of class 'AffyBatch'.")
    nrow.chip <- abatch@nrow
    ncol.chip <- abatch@ncol
    chiptype <- abatch@cdfName
  }

  if (is.null(nrow.chip) || is.null(ncol.chip) || is.null(chiptype))
    stop("nrow.chip, ncol.chip or chiptype not defined.")
  
  do.call("library", list(probes.pack))
  probe.tab <- get(probes.pack, envir=as.environment(paste("package:", probes.pack, sep="")))

  cdfenv <- new.env(hash=TRUE)

  for (i in seq(along=matches$match)) {
    xy <- getxy.probeseq(probeseq=probe.tab, i.row=matches$match[[i]],
                         x.colname = x.colname, y.colname = y.colname)
    if (nrow(xy) == 0 && simplify) {      
      next
    }
    assign(ids[i],
           cbind(xy2indices(xy[, 1], xy[, 2], nr=nrow.chip), mm),
           envir=cdfenv)
  }
  
  cdfenv <- wrapCdfEnvAffy(cdfenv, nrow.chip, ncol.chip, chiptype)  
  return(cdfenv)  
}

getxy.probeseq <- function(ppset.id=NULL, probeseq=NULL, i.row=NULL, offset.one=TRUE,
                           x.colname = "x", y.colname = "y") {
  if (sum(c(is.null(ppset.id), is.null(i.row))) != 1)
    stop("specify one and only one of 'ppset.id', 'i.row'")

  if (is.null(i.row))
    i.row <- probeseq$Probe.Set.Name %in% ppset.id

  mm.offset <- rep(0, length=length(i.row))
  mm.offset[i.row < 0] <- 1
  i.row <- abs(i.row)
  
  xy <- cbind(probeseq[[x.colname]][i.row], probeseq[[y.colname]][i.row] + mm.offset) + 1
  colnames(xy) <- c("x", "y")
  return(xy)
}


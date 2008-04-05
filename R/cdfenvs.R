wrapCdfEnvAffy <- function(cdfenv, nrow.chip, ncol.chip, chiptype, check = TRUE, verbose = FALSE) {
  object <- new("CdfEnvAffy", envir = cdfenv,
                envName = chiptype,
                nrow = as.integer(nrow.chip), ncol = as.integer(ncol.chip),
                probeTypes = c("pm", "mm"),
                chipType = chiptype)
  if (check) {
    valid <- validCdfEnvAffy(object, verbose=verbose)
    if ( ! valid ) {
      printValidCdfEnvAffy(valid)
      stop("invalid CdfEnvAffy")
    }
  }
  return(object)
}

getCdfEnvAffy <- function(abatch) {
  if (! is(abatch, "AffyBatch"))
    stop("arg must be of class 'AffyBatch'.")

  cdfenv <- getCdfInfo(abatch)

  cdfenv <- wrapCdfEnvAffy(cdfenv, abatch@nrow, abatch@ncol, abatch@cdfName)
  return(cdfenv)
}


getxy.probeseq <- function(ppset.id=NULL, probeseq=NULL, i.row=NULL,
                           xy.offset=NULL,
                           x.colname = "x", y.colname = "y") {

  if ( is.null(xy.offset) ) {
    xy.offset <- getOption("BioC")$affy$xy.offset
  }
  
  if (sum(c(is.null(ppset.id), is.null(i.row))) != 1)
    stop("specify one and only one of 'ppset.id', 'i.row'")

  if (is.null(probeseq))
    stop("the argument 'probeseq' must be specified !")
  
  if (is.null(i.row))
    i.row <- probeseq$Probe.Set.Name %in% ppset.id

    
  mm.offset <- rep(0, length=length(i.row))
  mm.offset[i.row < 0] <- 1
  i.row <- abs(i.row)

  xy <- cbind(probeseq[[x.colname]][i.row],
              probeseq[[y.colname]][i.row] + mm.offset) + xy.offset
  
  colnames(xy) <- c("x", "y")
  return(xy)
}


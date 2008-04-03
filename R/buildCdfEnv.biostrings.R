
mmProbes <- function(probes)
{  
  len_probe <- unique(nchar(probes$sequence))
  if (length(len_probe) > 1)
    stop(paste("Different length for probes",
               "(and the handling of that case is not implemented)."))
  if (len_probe != 25)
    stop("The expected probe length is 25 bp.")
  mmpos <- 13


  ## First write. Using DNAStringSet, short and elegant...
  ## but unfortunately unbearably slow
  ##  mmseq <-
  ##     lapply(as.list(probes$stringset),
  ##            function(x) {
  ##              replaceLetterAtLoc(x, mmpos,
  ##                                 as.character(complement(x[mmpos])))
  ##            })
  pmprobe <- substr(probes$sequence, mmpos, mmpos)
  mmprobe <- rep(as.character(NA), length=length(pmprobe))
  mmprobe[grep("[Aa]", pmprobe)] <- "T"
  mmprobe[grep("[Tt]", pmprobe)] <- "A"
  mmprobe[grep("[Gg]", pmprobe)] <- "G"
  mmprobe[grep("[Cc]", pmprobe)] <- "C"
  mmseq <- paste(substr(probes$sequence, 1, mmpos-1),
                 mmprobe,
                 substr(probes$sequence, mmpos+1, nchar(probes$sequence)),
                 sep = "")
                   
  return(mmseq)
}


setClass("AffyProbesMatch",
         representation(pm = "list", mm = "list",
                        labels = "character", chip_type = "character",
                        probes = "ANY"), # should be class "probetable" - S4 don't seem to cope with it
         validity = function(x) {
           if (length(x@pm) != length(x@mm))
             return("mm and pm should have identical lengths")
           if (length(x@pm) != length(x@labels))
             return("labels and pm should have identical lengths")
           if (any(duplicated(x@labels)))
             return("labels should be unique.")
           if (length(x@chip_type) != 1)
             return("chip_type should be *one* chip type name")
           if (! all(lapply(x@pm, function(y) inherits(y, "ByPos_MIndex"))))
             return("all pm should inherit from ByPos_MIndex")
           if (! all(lapply(x@mm, function(y) inherits(y, "ByPos_MIndex"))))
             return("all mm should inherit from ByPos_MIndex")
           return(TRUE)
         })

## setMethod("merge", signaturex="AffyProbesMatch", y="AffyProbesMatch",
##           function() {
##
##          })


matchAffyProbes <- function(probes, targets, chip_type)
{

  if (! inherits(probes, "probetable")) {
    stop(paste("'probes' should inherit from class 'probetable'."))
  }

  stringset <- DNAStringSet(probes$sequence)
  
  if (inherits(targets, "character")) {
    targets <- as.list(targets)
    for (ii in seq(along = targets)) {
      if (is.na(targets[[ii]])) {
        stop(paste("Target", ii, "is NA."))
      }
      targets[[ii]] <- DNAString(targets[[ii]])
    }
  } else if (inherits(targets, "list")) {
    for (ii in seq(along = targets)) {
      if (! inherits(targets[[ii]], "DNAString")) {
        stop("Invalid 'targets'.")
      }
    }
  } else if (! inherits(targets, "DNAString")) {
    stop("Invalid 'targets'.")
  }
  
  labels <- names(targets)
  if (is.null(labels)) {
    labels <- as.character(seq(along=targets))
  }
  
  pmdict <- PDict(stringset)
  mindex_pm <- vector("list", length = length(targets))
  for (ii in seq(along = targets)) {
    mindex_pm[[ii]] <- matchPDict(pmdict, targets[[ii]])
  }
  
  mmseq <- mmProbes(probes)
  mmdict <- PDict(mmseq)
  mindex_mm <- vector("list", length = length(targets))
  for (ii in seq(along = targets)) {
    mindex_mm[[ii]] <- matchPDict(mmdict, targets[[ii]])
  }
  
  apm <- new("AffyProbesMatch", pm = mindex_pm, mm = mindex_mm,
             labels = labels, chip_type = chip_type, probes=probes)
  return(apm)  
}



buildCdfEnv.biostrings <- function(apm, probes.pack,
                                   abatch=NULL,
                                   nrow.chip=NULL, ncol.chip=NULL,
                                   chiptype=NULL, mm=NA, simplify = TRUE,
                                   x.colname = "x", y.colname = "y",
                                   verbose = FALSE) {

  validObject(apm)

  if ( ! is.null(abatch)) {
    if (! is(abatch, "AffyBatch"))
      stop("abatch must be of class 'AffyBatch'.")
    nrow.chip <- abatch@nrow
    ncol.chip <- abatch@ncol
    chiptype <- abatch@cdfName
  }

  if (is.null(nrow.chip) || is.null(ncol.chip) || is.null(chiptype))
    stop("nrow.chip, ncol.chip or chiptype not defined.")

  probetab <- apm@probetab
  
  cdfenv <- new.env(hash=TRUE)

  if (verbose) {
    cat("Processing the matches:\n")
    pbt <- new("ProgressBarText", length(apm@pm),
               barsteps = as.integer(20))
    open(pbt)
  }
  for (i in seq(along = apm@pm)) {
    if (verbose)
      update(pbt)

    xy <- getxy.probeseq(probeseq=probetab,
                         i.row=countIndex(apm@pm[[i]]) > 0,
                         x.colname = x.colname, y.colname = y.colname)
    if (nrow(xy) == 0 && simplify) {
      next
    }
    assign(ids[i],
           cbind(xy2indices(xy[, 1], xy[, 2], nr=nrow.chip), mm),
           envir=cdfenv)
  }
  if (verbose)
    close(pbt)

  cdfenv <- wrapCdfEnvAffy(cdfenv, nrow.chip, ncol.chip, chiptype)
  return(cdfenv)
}

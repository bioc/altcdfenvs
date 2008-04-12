
mmProbes <- function(probes)
{  
  len_probe <- unique(nchar(probes$sequence))
  if (length(len_probe) > 1)
    stop(paste("Different length for probes",
               "(and the handling of that case is not implemented)."))
  if (len_probe != 25)
    stop(paste("The expected probe length is 25 bp, not ",
               len_probe, ".", sep=""))
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
  mmprobe[grep("[Gg]", pmprobe)] <- "C"
  mmprobe[grep("[Cc]", pmprobe)] <- "G"
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
 validity = function(obj) {
   if (length(obj@pm) != length(obj@mm))
     return("mm and pm should have identical lengths")
   if (length(obj@pm) != length(obj@labels))
     return("labels and pm should have identical lengths")
   if (any(duplicated(obj@labels)))
     return("labels should be unique.")
   if (length(obj@chip_type) != 1)
     return("chip_type should be *one* chip type name")
   if (! all(unlist(lapply(obj@pm,
                           function(y) inherits(y, "integer"))))) {
     return("all pm should inherit from numeric")
   }
   if (! all(unlist(lapply(obj@mm,
                           function(y) inherits(y, "integer"))))) {
     return("all mm should inherit from numeric")
   }
   return(TRUE)
 })

setMethod("show",
          signature = c("AffyProbesMatch"),
          function(object) {
            cat("AffyProbesMatch:\n")
            cat(paste(length(object@pm),
                      "target(s) sequences matched",
                      "against", nrow(object@probes),
                      "probes of chip type", object@chip_type,
                      ".\n"))
          }
          )

setMethod("combine",
          signature = c("AffyProbesMatch", "AffyProbesMatch"),
          function(x, y, ...) {
            if (x@chip_type != y@chip_type)
              stop("Both 'chip_type' must be identical.")
            if (! identical(x@probes, y@probes))
              stop("Both probe tables must be identical.")
            pm <- c(x@pm, y@pm)
            mm <- c(x@mm, y@mm)
            labels <- c(x@labels, y@labels)
            chip_type <- x@chip_type
            probetable <- x@probes
            res <- new("AffyProbesMatch",
                       pm = pm, mm = mm,
                       labels = labels,
                       chip_type = chip_type,
                       probes = probetable)
            return(res)
          })

toHypergraph <-
  function(object, ...) {
    stop("Not available for the given signature.")
  }
setGeneric("toHypergraph")
setMethod("toHypergraph",
          signature = c("AffyProbesMatch"),
          function(object, simplify=TRUE, ...) {
            if (simplify) {
              target_match <-
                unlist(lapply(object@pm, function(x) length(x) > 0))
              probe_match <- rep(FALSE, length=nrow(object@probes))
               for (i in which(target_match)) {
                 probe_match[object@pm[[i]]] <- TRUE
               }
            } else {
              target_match <- rep(TRUE, length=length(object@pm))
              probe_match <- rep(TRUE, length=nrow(object@probes))
            }

            i_match <- rep(as.integer(NA), nrow(object@probes))
            i_match[probe_match] <- seq(along=which(probe_match))
            
            nodes <-
              paste(as.character(object@probes[["x"]][probe_match]),
                    as.character(object@probes[["y"]][probe_match]),
                    sep = "-")

            hEdges <- lapply(object@pm[target_match],
                             function(x) Hyperedge(nodes[i_match[x]]))
            names(hEdges) <- object@labels[target_match]
            
            hg <- new("Hypergraph",
                      nodes = nodes,
                      hyperedges = hEdges)
            return(hg)
          }
          )

setMethod("toHypergraph",
          signature = c("CdfEnvAffy"),
          function(object, ...)
          {
            targets <- ls(object@envir)
            nodesEnv <- new.env(hash=TRUE, parent=emptyenv())
            for (n in targets) {
              m <- object@envir[[n]]
              labels <- apply(index2xy(object, m[, 1]), 1,
                              function(x) paste(x, collapse="-"))   
              nodesEnv[[n]] <- labels
            }
            nodes <- unlist(as.list(nodesEnv), use.names=FALSE)
            nodes <- unique(nodes)
            hEdges <- lapply(nodesEnv,
                             function(x) Hyperedge(x))
            
            hg <- new("Hypergraph",
                      nodes = nodes,
                      hyperedges = hEdges)
            return(hg)
          })

matchAffyProbes <-
  function(probes, targets, chip_type,
           matchmm = TRUE,
           selectMatches = function(x) which(countIndex(x) > 0),
           ...)
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
    md <- matchPDict(pmdict, targets[[ii]], ...)
    mindex_pm[[ii]] <- selectMatches(md)
  }

  mindex_mm <- vector("list", length = length(targets))
  if (matchmm) {
    mmseq <- mmProbes(probes)
    mmdict <- PDict(mmseq)
    for (ii in seq(along = targets)) {
      md <- matchPDict(mmdict, targets[[ii]], ...)
      mindex_mm[[ii]] <- selectMatches(md)
    }
  }
  
  apm <- new("AffyProbesMatch", pm = mindex_pm, mm = mindex_mm,
             labels = labels, chip_type = chip_type, probes=probes)
  return(apm)  
}




buildCdfEnv.biostrings <- function(apm,
                                   abatch=NULL,
                                   nrow.chip=NULL, ncol.chip=NULL,
                                   simplify = TRUE,
                                   x.colname = "x", y.colname = "y",
                                   verbose = FALSE)
{

  if (verbose)
    cat("validating 'apm'...")
  validObject(apm)
  if (verbose)
    cat("done.\n")
  
  if ( ! is.null(abatch)) {
    if (! is(abatch, "AffyBatch"))
      stop("abatch must be of class 'AffyBatch'.")
    nrow.chip <- abatch@nrow
    ncol.chip <- abatch@ncol
    chip_type <- abatch@cdfName
  }

  if (is.null(nrow.chip) || is.null(ncol.chip))
    stop("nrow.chip, ncol.chip")

  probetab <- apm@probes
  
  cdfenv <- new.env(hash=TRUE)

  if (verbose) {
    cat("Processing the matches:\n")
    pbt <- new("ProgressBarText", length(apm@pm),
               barsteps = as.integer(20))
    open(pbt)
  }

  ##FIXME:
  warning("Check index for MM probes.")
  for (i in seq(along = apm@pm)) {
    if (verbose)
      update(pbt)

    xy <- getxy.probeseq(probeseq = probetab,
                         i.row = apm@pm[[i]],
                         x.colname = x.colname, y.colname = y.colname)
    if (nrow(xy) == 0 && simplify) {
      next
    }
    assign(apm@labels[i],
           cbind(xy2indices(xy[, 1], xy[, 2], nr=nrow.chip),
                 xy2indices(xy[, 1]+1, xy[, 2], nr=nrow.chip)),
           envir=cdfenv)
  }
  if (verbose)
    close(pbt)

  cdfenv <- wrapCdfEnvAffy(cdfenv, nrow.chip, ncol.chip, apm@chip_type)
  return(cdfenv)
}

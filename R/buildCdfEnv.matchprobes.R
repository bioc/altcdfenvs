buildCdfEnv.matchprobes <- function(matches, ids, probes.pack,
                                    abatch=NULL,
                                    nrow.chip=NULL, ncol.chip=NULL,
                                    chiptype=NULL, mm=NA, simplify = TRUE,
                                    x.colname = "x", y.colname = "y",
                                    verbose = FALSE) {

  .Deprecated("buildCdfEnv.biostrings", package="altcdfenvs")
  
  if (! (is.list(matches) && length(matches) > 0)) #&& length(matches[[1]]) < 3)
    stop("arg 'matches' should be a list as returned by 'matchprobes'.")

  if (length(matches[[1]]) != length(ids))
    stop("'matches' and 'ids' must have the same length.")

  if (length(matches[[1]]) != length(unique(ids)))
    stop("Some elements in 'ids' are not unique. You probably do not want this.")
  
  if ( ! is.null(abatch)) {
    if (! is(abatch, "AffyBatch"))
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

  if (verbose) {
    cat("Processing the matches:\n")
    pbt <- new("ProgressBarText", length(matches$match), barsteps = as.integer(20))
    open(pbt)
  }
  for (i in seq(along=matches$match)) {
    if (verbose)
      update(pbt)

    xy <- getxy.probeseq(probeseq=probe.tab, i.row=matches$match[[i]],
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

appendCdfEnvAffy <- function(acdfenv, id, i, nocopy=TRUE) {

  stopifnot(is.character(id), is.matrix(i), is.integer(i))
  tmp.m <- match(colnames(i), acdfenv@probeTypes)

  if (! all(colnames(i) %in% acdfenv@probeTypes)) {
    stop("the CdfEnv only knows about probe types:\n", paste(acdfenv@probeTypes, collapse=","))
  }

  if (exists(id, envir=acdfenv@envir)) {
    stop("'id' already in 'adfenv'")
  }

  if (! nocopy) {
    acdfenv <- copyCdfEnvAffy(acdfenv)
  }

  m <- matrix(as.integer(NA), nrow=nrow(i), ncol=length(acdfenv@probeTypes))
  m[, tmp.m] <- i
  
  ##DEBUG: a consistency check (i.e., the index make sense given the
  ##       geometry of the chip) would be nice...

  assign(id, m, envir=acdfenv@envir)

  return(acdfenv)
    
}

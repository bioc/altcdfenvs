##setGeneric("unique", )

unique.CdfEnvAffy <- function(x, incomparable=FALSE, simplify=TRUE, verbose=FALSE, ...) {
  counts <- countduplicated(x, incomparable=incomparable, verbose=verbose)
  tmp.env<- as(x, "environment")
  ids <- ls(tmp.env)
  ## copy env
  y <- new.env(hash=TRUE)
  if (verbose)
    cat("removing duplicated elements...")
  for (i in ids) {
    tmp.count <- get(i, envir=counts)
    tmp.i <- get(i, envir=tmp.env)
    tmp.ok <- tmp.count == 1
    tmp.i[!tmp.ok] <- NA
    tmp.new <- tmp.i[!apply(tmp.i, 1, function(x) all(is.na(x))), , drop=FALSE]
    if (length(tmp.new) == 0 && simplify) {
      ##if (verbose) {
      ##  cat(paste("removing ", i, " (does not have anymore elements).\n"))
      ##}
      next
    }
    assign(i, tmp.new, envir=y)
  }
  if (verbose)
    cat("done.\n")
  r <- x
  r@envir <- y
  r@envName <- paste(r@envName, "-unique", sep="")
  return(r)
}


countduplicated <- function(x, incomparable=FALSE, verbose=FALSE) {
  if (!is(x, "CdfEnvAffy"))
    stop("x must inherit from 'CdfEnvAffy'")

  if (incomparable != FALSE)
    warning("'incomparable' not yet implemented !")

  if (verbose)
    cat("Initialize...")

  tmp.env <- as(x, "environment")
  ids <- ls(tmp.env)
  p.type <- x@probeTypes

  tmp.count <- rep(as.integer(0),
                   length=sum(unlist(lapply(indexProbes(x, p.type), length))))
  r <- new.env(hash = TRUE)
  if(verbose)
    cat("done.\nCounting probes...")

  for (i in seq(along = ids)) {
    p.i <- get(ids[i], envir = tmp.env)
    tmp.count[p.i] <- tmp.count[p.i] + 1
  }
  if (verbose)
    cat("done.\nAssigning counts...")
  for (i in seq(along = ids)) {
    p.i <- get(ids[i], envir = tmp.env)
    p.i[] <- tmp.count[p.i]
    assign(ids[i], p.i, envir=r)
  }
  if (verbose)
    cat("done.\n")
  return(r)
}

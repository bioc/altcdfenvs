removeIndex <- function(x, i, simplify = TRUE, verbose=FALSE) {
  if (! is.integer(i)) {
    stop("'i' must be of mode integer.")
  }
  remove.me <- rep(FALSE, length=max(i))
  remove.me[i] <- TRUE
  tmp.env<- as(x, "environment")
  ids <- ls(tmp.env)
  ## copy env
  y <- new.env(hash=TRUE)
  if (verbose)
    cat("removing duplicated elements...")
  for (i in ids) {
    tmp.i <- get(i, envir=tmp.env)
    tmp.ok <- (c(tmp.i) > length(remove.me)) | (! remove.me[c(tmp.i)])
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
  r@envName <- paste(r@envName, "-removed", sep="")
  return(r)
}

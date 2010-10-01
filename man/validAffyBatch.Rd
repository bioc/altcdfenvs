\name{validAffyBatch}
\alias{validAffyBatch}
\alias{validCdfEnvAffy}
\alias{printValidCdfEnvAffy}
\title{Check validity of a CdfEnvAffy.}
\description{
  Tries to see if a CdfEnvAffy, or a pair of AffyBatch / CdfEnvAffy is valid.
}
\usage{
validAffyBatch(abatch, cdfenv)
validCdfEnvAffy(cdfenv, verbose=TRUE)
printValidCdfEnvAffy(x)
}
\arguments{
  \item{abatch}{ instance of \code{\link[affy]{AffyBatch-class}} }
  \item{cdfenv}{ instance of \code{\link{CdfEnvAffy-class}} }
  \item{verbose}{ verbose or not}
  \item{x}{object returned by \code{validCdfEnvAffy}}
}
\details{
  The function \code{validAffyBatch} calls in turn
  \code{validCdfEnvAffy}.
}
\seealso{ \code{\link[affy]{AffyBatch-class}}, \code{\link{CdfEnvAffy-class}} }
\examples{
## To be done...
}
\keyword{ manip }


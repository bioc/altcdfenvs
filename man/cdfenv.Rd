\name{cdfenvs}
\alias{wrapCdfEnv}
\alias{getCdfEnv}
\alias{buildCdfEnv.matchprobes}
\alias{CdfEnv}
\title{ functions related to cdfenvs }
\description{
  A set of functions to handle cdfenvs
}
\usage{
wrapCdfEnv(cdfenv, nrow.chip, ncol.chip, chiptype)
getCdfEnv(abatch)
buildCdfEnv.matchprobes(matches, ids, probes.pack, abatch=NULL,
                        nrow.chip=NULL, ncol.chip=NULL, chiptype=NULL, mm=NA)
}
\arguments{
  \item{abatch}{an \code{AffyBatch}}
  \item{cdfenv}{A cdfenv environment}
  \item{chiptype}{A name for the chip type}
  \item{ids}{a vector of probe set identifiers for the matches}
  \item{matches}{a list as returned by the function
    \code{\link[matchprobes]{matchprobes}}}
  \item{mm}{The value to store for MMs}
  \item{ncol.chip}{The number of columns for the chip type}
  \item{nrow.chip}{The number of rows for the chip type}
  \item{probes.pack}{The name of the probe package}
}
\details{
  The class \code{CdfEnv} is defined as S3 to ensure compatibility with
  what is existing while changes are made.
  Migration to S4 is likely to happen after that\ldots
}
\value{
  An instance of class \code{CdfEnv}.
}
\examples{
## To appear
}
\keyword{manip}

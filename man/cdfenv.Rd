\name{cdfenvs}
\alias{wrapCdfEnvAffy}
\alias{getCdfEnvAffy}
\alias{buildCdfEnv.matchprobes}
\alias{CdfEnv}
\title{ functions related to cdfenvs }
\description{
  A set of functions to handle cdfenvs
}
\usage{
wrapCdfEnvAffy(cdfenv, nrow.chip, ncol.chip, chiptype)
getCdfEnvAffy(abatch)
buildCdfEnv.matchprobes(matches, ids, probes.pack, abatch=NULL,
                        nrow.chip=NULL, ncol.chip=NULL, chiptype=NULL,
                        mm=NA, simplify = TRUE,
                        x.colname = "x", y.colname = "y", verbose=FALSE)
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
  \item{simplify}{Simplify the environment created by removing the ids
    without any matching probe}
  \item{x.colname, y.colname}{see the \code{getxy.probeseq}}
  \item{verbose}{verbosity (\code{TRUE} or \code{FALSE})}
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

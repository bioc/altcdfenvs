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
wrapCdfEnvAffy(cdfenv, nrow.chip, ncol.chip, chiptype, check = TRUE,
               verbose = FALSE)
getCdfEnvAffy(abatch)
buildCdfEnv.matchprobes(matches, ids, probes.pack, abatch=NULL,
                        nrow.chip=NULL, ncol.chip=NULL, chiptype=NULL,
                        mm=NA, simplify = TRUE,
                        x.colname = "x", y.colname = "y", verbose=FALSE)
}
\arguments{
  \item{abatch}{an \code{AffyBatch}}
  \item{cdfenv}{A cdfenv environment}
  \item{check}{perform consistency check or not}
  \item{chiptype}{A name for the chip type}
  \item{ids}{a vector of probe set identifiers for the matches}
  \item{matches}{a list as returned by the function
    \code{combineAffyBatch}}
  \item{mm}{The value to store for MMs}
  \item{ncol.chip}{The number of columns for the chip type}
  \item{nrow.chip}{The number of rows for the chip type}
  \item{probes.pack}{The name of the probe package}
  \item{simplify}{Simplify the environment created by removing the ids
    without any matching probe}
  \item{x.colname, y.colname}{see the \code{getxy.probeseq}}
  \item{verbose}{verbosity (\code{TRUE} or \code{FALSE})}
}
\value{
  An instance of class \code{CdfEnvAffy}.
}
\examples{
## See the main vignette
}
\keyword{manip}

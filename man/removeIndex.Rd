\name{removeIndex}
\alias{removeIndex}
\title{ A function to remove probes in an environment }
\description{
  A function to remove probes in an environment, given their index.
}
\usage{
removeIndex(x, i, simplify = TRUE, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ An instance of \code{CdfEnvAffy-class} }
  \item{i}{ A vector of indexes (integers !). }
  \item{simplify}{ Simply the resulting \code{CdfEnvAffy} (see details). }
  \item{verbose}{ verbose output or not. }
}
\details{
  The probes to be removed are set to \code{NA} in the CdfEnvAffy.
  When \code{simplify} is set to \code{TRUE} the probe sets are
  simplified whenever possible. For example, if both pm and mm for the
  same probe pair are set to \code{NA}, then the probe pair is removed
  from the probe set. 
}
\value{
  An instance of \code{CdfEnvAffy-class} is returned.
}
\author{ Laurent Gautier }
\seealso{ \code{\link{CdfEnvAffy-class}} }
\examples{
## to be done
}
\keyword{ manip }

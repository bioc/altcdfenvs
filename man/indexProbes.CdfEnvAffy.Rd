\name{indexProbes.CdfEnvAffy}
\alias{indexProbes.CdfEnvAffy}
\title{ indexes for probes }
\description{
  A function to get the index for probes
}
\usage{
indexProbes.CdfEnvAffy(object, which, probeSetNames = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ \code{CdfEnvAffy} }
  \item{which}{ which kind of probe are of interest (see details). }
  \item{probeSetNames}{ names of the probe sets of interest. If
    \code{NULL}, all the probe sets are considered. }
}
\details{
  The parameter \code{which} let one specify which category of probes
  are of interest. In the case of Affymetrix chips, probes can be \code{"pm"}
  probes or \code{"mm"} probes. It the parameter is set to \code{c("pm",
    "mm")}, both are returned. Should other categories be defined, they
  can be handled as well.
  
  One can note that there is a method \code{indexProbes} defined for the
  class \code{AffyBatch}. This can be considered as a convenience
  extension, since logically the indexes are computed from the cdfenv
  information.
}
\value{
A \code{list} of indexes.
}
\seealso{\code{\link{CdfEnvAffy-class}}, \code{\link[affy]{AffyBatch-class} }}
\examples{
## see "CdfEnvAffy-class"
}
\keyword{ manip }

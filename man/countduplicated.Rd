\name{countduplicated}
\alias{countduplicated}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Count the number of times probes are used }
\description{
  This function counts the number of times the probes in a CdfEnvAffy
  are found in this object.
}
\usage{
countduplicated(x, incomparable = FALSE, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An instance of \code{CdfEnvAffy-class} }
  \item{incomparable}{ (not implemented yet, keep away) }
  \item{verbose}{ verbose or not }
}
\details{
 
}
\value{
  An \code{environment} is returned. Each element in this
  \code{environment} has the same identifier than its corresponding
  probe set in the \code{CdfEnvAffy-class} and contains the number of
  times a probe is in use in the environment (instead of an index number
  in the \code{CdfEnvAffy-class}).
}
\author{ Laurent }
\seealso{ \code{\link{CdfEnvAffy-class}}}
\examples{
}
\keyword{ manip }

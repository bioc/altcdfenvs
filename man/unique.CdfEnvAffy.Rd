\name{unique.CdfEnvAffy}
\alias{unique.CdfEnvAffy}
\title{ Remove duplicated elements from a CdfEnvAffy }
\description{
  Remove duplicated elements from a CdfEnvAffy
}
\usage{
unique.CdfEnvAffy(x, incomparables = FALSE, simplify = TRUE, verbose = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ An instance of \code{CdfEnvAffy-class} }
  \item{incomparables}{ (not yet implemented) }
  \item{simplify}{ simplify the result }
  \item{verbose}{ verbose or not }
  \item{\dots}{ (here for compatibility with the generic \code{unique}) }
}
\details{
  The parameter \code{simplify} has the same function as the one with
  the same name in \code{countduplicated}.
}
\value{
  An instance of \code{CdfEnvAffy-class} in which probes used several
  times are removed.
}
\author{ Laurent }
\section{Warning }{The function differs slightly from the generic
  \code{unique}. Here the elements found in several place a merely removed.}

\seealso{\code{\link{countduplicated}} }
\examples{
##not yet here...
}
\keyword{ manip }

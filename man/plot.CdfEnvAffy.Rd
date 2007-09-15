\name{plot.CdfEnvAffy}
\alias{plot.CdfEnvAffy}
\title{ A function to `plot' a CdfEnvAffy }
\description{
  A function to set the axis and plot the outline for a CdfEnvAffy
}
\usage{
\method{plot}{CdfEnvAffy}(x, xlab = "", ylab = "", main = x@chipType, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a \code{CdfEnvAffy} }
  \item{xlab}{ label for the rows }
  \item{ylab}{ label for the columns }
  \item{main}{ label for the plot. The chip-type by default. }
  \item{\dots}{ optional parameters to be passed to the underlying
    function \code{plot}}
}
\details{
  This function does not `plot' much, but sets the coordinates for
  further plotting (see the examples).
}
\author{ Laurent }
\seealso{\code{\link{CdfEnvAffy-class}} }
\examples{
## See "CdfEnvAffy-class"
}
\keyword{ hplot }


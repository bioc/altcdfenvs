\name{index2xy}
\alias{index2xy}
\alias{index2xy.CdfEnvAffy}
%\alias{index2xy.AltAffyBatch}
\alias{xy2index}
\alias{xy2index.CdfEnvAffy}
%\alias{xy2index.AltAffyBatch}
\title{ Functions to shuttle from indexes to XY coordinates }
\description{
  Functions to shuttle from indexes to XY coordinates.
}
\usage{
index2xy(object, ...)
xy2index(object, ...)
index2xy.CdfEnvAffy(object, i)
%index2xy.AltAffyBatch(object, i)
xy2index.CdfEnvAffy(object, x, y)
%xy2index.AltAffyBatch(object, x, y)
}
\arguments{
  \item{object}{An object of class \code{CdfEnvAffy}.}
  \item{i}{A vector of indexes.}
  \item{x, y}{Vectors of X and Y coordinates.}
  \item{\dots}{Optional parameters (not used).}
}
\details{
  
}
\value{
  A vector of integers (for \code{xy2index} methods), or a matrix of two
  columns (for \code{index2xy} methods).
}
\section{Warning }{}

\seealso{ \code{\link{CdfEnvAffy-class}}}
\examples{
## To be done...
}
\keyword{ manip }


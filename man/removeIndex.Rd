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

## use plasmodiumanopheles chip as an example
library(plasmodiumanophelescdf)

## wrap in a (convenient) CdfEnvAffy object
planocdf <- wrapCdfEnvAffy(plasmodiumanophelescdf, 712, 712, "plasmodiumanophelescdf")
print(planocdf)

## ask for the probe indexed '10759' to be removed
## (note: if one wishes to remove from X/Y coordinates,
## the function xy2index can be of help).
planocdfCustom <- removeIndex(planocdf, as.integer(10759))

## let see what happened (we made this example knowing in which
## probe set the probe indexed '10759' is found).
indexProbes(planocdf, "pm", "200000_s_at")
indexProbes(planocdfCustom, "pm", "200000_s_at")
## The 'second' pm probe (indexed '10579') in the probe set is now set
## to NA.

}
\keyword{ manip }

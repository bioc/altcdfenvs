\name{CdfEnvAffy-class}
\docType{class}
\alias{CdfEnvAffy-class}
\alias{geneNames,CdfEnvAffy-method}
\alias{index2xy,CdfEnvAffy-method}
\alias{indexProbes,CdfEnvAffy,character-method}
\alias{plot,CdfEnvAffy,missing-method}
\alias{show,CdfEnvAffy-method}
\alias{xy2index,CdfEnvAffy-method}
\alias{coerce,CdfEnvAffy,environment-method}
\alias{coerce,CdfEnvAffy,Cdf-method}
\alias{[,CdfEnvAffy-method}

\title{Class "CdfEnvAffy"}
\description{A class to hold the information necessary to handle the
  grouping of probes in set of probes, and to find XY coordinates
 of probes on a chip}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CdfEnvAffy", ...)}.
Typically, there is an instance of the class for each type of chip
(e.g. Hu6800, HG-U95A, etc\ldots).
}
\section{Slots}{
  \describe{
    \item{\code{envir}:}{Object of class \code{"environment"}. It has to
      be tought of as a hashtable: the keys are probe set identifiers,
      or gene names, and the values are indexes.}
    \item{\code{envName}:}{Object of class \code{"character"}. A name
      for the environment.}
    \item{\code{index2xy}:}{Object of class \code{"function"}. The
      function used to resolve index into xy coordinates. Unless you are
    an advanced user, you probably want to ignore this (and rely on the
    default provided with the package).}
    \item{\code{xy2index}:}{Object of class \code{"function"}. The
      function used to resolve xy coordinates into index. Unless you are
    an advanced user, you probably want to ignore this (and rely on the
    default provided with the package).}
    \item{\code{nrow}:}{Object of class \code{"integer"}. The number of
      rows of probes for the chip type.}
    \item{\code{ncol}:}{Object of class \code{"integer"}. The number of
      columns of probes for the chip type.}
    \item{\code{probeTypes}:}{Object of class \code{"character"}. The
      different types of probes stored for each probe set. In the case
      of {\it Affymetrix} chips, the probes are typically perfect match
      (pm) probes or mismatch probes (mm).}
    \item{\code{chipType}:}{Object of class \code{"character"}. The name
      of the chip type the instance is associated with. This is useful
      when one starts to create alternative mappings of the probes on a
      chip (see associated vignette).}
  }
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(object = "CdfEnvAffy")}: subset a cdf, that
      is return a new cdf containing only a subset of the probe
      sets. The subset of probe sets to take is identified as a vector
      of identifiers (mode "character").}
    \item{coerce}{\code{signature(object = "CdfEnvAffy", "environment")}: coerce an
      instance of the class to an \code{enviroment}.}
    \item{coerce}{\code{signature(object = "CdfEnvAffy", "Cdf")}: coerce an
      instance of the class to a \code{Cdf}.}
    \item{geneNames}{\code{signature(object="CdfEnvAffy")}: Return the
    names of the known probe sets (of course, it depends on the associated CDF).}
    \item{index2xy}{\code{signature(object = "CdfEnvAffy", i="integer")}: convert
      index values into XY coordinates.}
    \item{indexProbes}{\code{signature(object = "CdfEnvAffy", which =
	"character", probeSetNames = NULL)}: obtain the indexes for the probes associated wit
      the probe set name \code{probeSetNames}. When \code{probeSetNames}
      is set to \code{NULL} (default), the indexes are returned for the
      probe sets defined on the chip. See \code{indexProbes.CdfEnvAffy}}
    \item{plot}{\code{signature(x = "CdfEnvAffy", y = "missing")}: Plot
      the chip. It mainly sets coordinates for further plotting (see
      examples). See \code{plot.CdfEnvAffy}}
    \item{show}{\code{signature(object = "CdfEnvAffy")}: Print method.}
    \item{xy2index}{\code{signature(object = "CdfEnvAffy", x="integer", y="integer")}: convert XY
      coordinates into index values.}
  }
}
%\references{ ~put references to the literature/web site here ~ }
\author{Laurent Gautier}
%\note{ ~~further notes~~ }

% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{indexProbes.CdfEnvAffy}}, \code{\link{plot.CdfEnvAffy}}
  %  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%  or \code{\link{CLASSNAME-class}} for links to other classes
}
\examples{
## build an instance
library(hgu95acdf)
cdfenv.hgu95a <- wrapCdfEnvAffy(hgu95acdf, 640, 640, "HG-U95A")

show(cdfenv.hgu95a)

## find the indexes for a probe set (pm only)
ip <- indexProbes(cdfenv.hgu95a, "pm", "1000_at")[[1]]
## get the XY coordinates for the probe set
xy <- index2xy(cdfenv.hgu95a, ip)

## plot the chip
plot(cdfenv.hgu95a)

## plot the coordinates
plotLocation(xy)

}
\keyword{classes}

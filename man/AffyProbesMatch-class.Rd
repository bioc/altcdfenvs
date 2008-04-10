\name{AffyProbesMatch-class}
\docType{class}
\alias{AffyProbesMatch-class}
\alias{combine,AffyProbesMatch,AffyProbesMatch-method}
\alias{show,AffyProbesMatch-method}
\alias{toHypergraph,AffyProbesMatch-method}
\title{Class "AffyProbesMatch"}
\description{
Store the results of a call to \code{matchAffyProbes}.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("AffyProbesMatch", ...)}.

An object will store the result of
matching probe sequences against target sequences.
}
\section{Slots}{
  \describe{
    \item{\code{pm}:}{Object of class \code{"list"}:
      each element is vector of index values}
    \item{\code{mm}:}{Object of class \code{"list"}:
      each element is vector of index values}
    \item{\code{labels}:}{Object of class \code{"character"}}
    \item{\code{chip_type}:}{Object of class \code{"character"}
      and of length 1.}
    \item{\code{probes}:}{Object of class \code{"ANY"}:
      the \code{probetable} object used to perform the matches.
    }
  }
}
\section{Methods}{
  \describe{
    \item{combine}{\code{signature(x = "AffyProbesMatch", y =
	"AffyProbesMatch")}: combine two instances. This is can
      be useful when splitting the list of target sequences to
      parallelized the job.}
    \item{show}{\code{signature(x = "AffyProbesMatch")}:Show the
      instance.}
    \item{toHypergraph}{\code{signature(object = "AffyProbesMatch")}:
    build an \code{\link[hypegraph]{Hypergraph}} from the matches.}
  }
}
%\references{ ~put references to the literature/web site here ~ }
%\author{ ~~who you are~~ }
%\note{ ~~further notes~~ }
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%
%\seealso{
%	~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%	or \code{\linkS4class{CLASSNAME}} for links to other classes
%}
\examples{
showClass("AffyProbesMatch")
}
\keyword{classes}

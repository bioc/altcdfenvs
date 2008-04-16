\name{buildCdfEnv.biostrings}
\alias{buildCdfEnv.biostrings}
\title{ Build CDF environments}
\description{
  Build CDF environment from Biostrings matchPDict results
}
\usage{
buildCdfEnv.biostrings(apm, abatch = NULL,
                       nrow.chip = NULL, ncol.chip = NULL,
                       simplify = TRUE,
                       x.colname = "x", y.colname = "y",
                       verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{apm}{ \code{AffyProbesMatch}}
  \item{abatch}{ \code{AffyBatch} }
  \item{nrow.chip}{number of rows for the chip type (see details)}
  \item{ncol.chip}{number of columns for the chip type (see details)}
  \item{simplify}{simplify the environment built (removing
  target names when there is no matching probe)}
  \item{x.colname}{ column name }
  \item{y.colname}{ column name}
  \item{verbose}{ verbose \code{TRUE}/\code{FALSE}}
}
\details{
  %  ~~ If necessary, more details than the description above ~~
  Whenever an \code{abatch} is specified, \code{nrow.chip} and
  \code{ncol.chip} are not needed. Specifying the an \code{AffyBatch}
  in \code{abatch} is the easiest way to specify information about the
  geometry of a chip type.
}
\value{
  An instance of class \code{CdfEnvAffy}.
}
%\references{ ~put references to the literature/web site here ~ }
%\author{ ~~who you are~~ }
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }


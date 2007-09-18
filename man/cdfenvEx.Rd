\name{cdfenvEx}
\alias{cdfenvEx}
\docType{data}
\title{ CdfEnvAffy}
\description{
An example of \code{CdfEnvAffy}
}
\usage{data(cdfenvEx)}
\format{
  The format is:
Formal class 'CdfEnvAffy' [package "altcdfenvs"] with 8 slots
  ..@ index2xy  :function (object, i)  
  ..@ xy2index  :function (object, x, y)  
  ..@ envir     :length 2 <environment>
  ..@ envName   : chr "ZG-DU33"
  ..@ nrow      : int 100
  ..@ ncol      : int 100
  ..@ probeTypes: chr [1:2] "pm" "mm"
  ..@ chipType  : chr "ZG-DU33"
}
\examples{
data(cdfenvEx)

print(cdfenvEx)
}
\keyword{datasets}

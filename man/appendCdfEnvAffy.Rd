\name{appendCdfEnvAffy}
\alias{appendCdfEnvAffy}
\title{ append probe sets to a CdfEnvAffy}
\description{
  append probe sets to a CdfEnvAffy
}
\usage{
appendCdfEnvAffy(acdfenv, id, i, nocopy=FALSE)
}
\arguments{
  \item{acdfenv}{ instance of class \code{CdfEnvAffy}}
  \item{id}{ identifier for the probe set to add }
  \item{i}{ a \code{matrix} of indexes (see details)}
  \item{nocopy}{whether to make a copy of the environment or not (see details)}
}
\details{
  The \code{matrix} \code{i} must have one column per probe type. For
  typical Affymetrix chip types, there are two probe types: \code{"pm"}
  and \code{"mm"}.

  \code{nocopy} set to \code{TRUE} means that the environment is added
  the probe set 'in-situ' (this can boost execution speed if you add a
  lot of probe sets).
}
\value{
  An \code{CdfEnvAffy} is returned
}
\examples{
data(cdfenvEx)

## pm and mm probe set
m <- matrix(1:10, ncol = 2)
colnames(m) <- c("pm", "mm")

appendCdfEnvAffy(cdfenvEx, "blabla", m)

indexProbes(cdfenvEx, c("pm", "mm"), "blabla")

## pm only probe set
m <- matrix(6:9, ncol = 1)
colnames(m) <- c("pm")
appendCdfEnvAffy(cdfenvEx, "blabla2", m)
## note that the unspecified "mm" were set to NA
indexProbes(cdfenvEx, c("pm", "mm"), "blabla2")


}
\keyword{ manip }

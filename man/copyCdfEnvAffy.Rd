\name{copyCdfEnvAffy}
\alias{copyCdfEnvAffy}
\title{ make a copy of a CdfEnvAffy }
\description{
  make a copy of a CdfEnvAffy
}
\usage{
copyCdfEnvAffy(acdfenv)
}
\arguments{
  \item{acdfenv}{ instance of class \code{CdfEnvAffy} }
}
\details{
  Make a copy can be needed since a \code{CdfEnvAffy} contains an \code{environment}
}
\value{
  A \code{CdfEnvAffy}
}
\seealso{\code{\link{CdfEnvAffy-class}}, \code{\link[Biobase]{copyEnv}}}
\keyword{ manip }% at least one, from doc/KEYWORDS

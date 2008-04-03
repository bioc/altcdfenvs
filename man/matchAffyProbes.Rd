\name{matchAffyProbes}
\alias{matchAffyProbes}

\title{ Match the probes on an Affymetrix array }
\description{
  Match the individual probes on an Affymetrix array to
  arbitrary targets.
}
\usage{
matchAffyProbes(probes, targets, chip_type)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{probes}{ a \code{probetable} object }
  \item{targets}{ a vector of references }
  \item{chip_type}{a name for the chip type.}
}
\details{
  
}
\value{
  A two-items list:
  \item{pm}{list of matches for the PM probes (of class \code{MIndex})}
  \item{mm}{list of matches for the MM probes (of class \code{MIndex})}
}
\author{ Laurent Gautier }
\note{


}
\seealso{ \code{\link[Biostrings]{MIndex}}
  and \code{\link{buildCdfEnv.biostrings}}}
\examples{

library(hgu133aprobe)

filename <- system.file("exampleData", "sample.fasta",
                        package="altcdfenvs")

fasta.seq <- readFASTA(filename)

m <- matchAffyProbes(hgu133aprobe,
                     unlist(lapply(fasta.seq, "[[", "seq")),
                     "HG-U133A")

# see the vignette for building a new environment from the begining
# to the end

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }


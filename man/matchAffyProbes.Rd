\name{matchAffyProbes}
\alias{matchAffyProbes}
\alias{mmProbes}
\title{ Match the probes on an Affymetrix array }
\description{
  Match the individual probes on an Affymetrix array to
  arbitrary targets.
}
\usage{
mmProbes(probes)

matchAffyProbes(probes, targets, chip_type,
                matchmm = TRUE,
                selectMatches = function(x) which(countIndex(x) > 0))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{probes}{ a \code{probetable} object }
  \item{targets}{ a vector of references }
  \item{chip_type}{a name for the chip type.}
  \item{matchmm}{whether to match MM probes or not}
  \item{selectMatches}{a function to select matches}
}
\details{
  
}
\value{
  \code{mmProbes} returns a vector of MM probe sequences.
  
  \code{matchAffyProbes} returns a two-items list:
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

targets <- unlist(lapply(fasta.seq, "[[", "seq"))
names(targets) <- lapply(fasta.seq, "[[", "desc")
names(targets) <- sub("^>.+\\\|(NM[^ \\\|]+|Hs[^ \\\|]+)\\\|? .+$", "\\\1",
                      names(targets))

m <- matchAffyProbes(hgu133aprobe,
                     targets,
                     "HG-U133A")

# see the vignette for building a new environment from the begining
# to the end

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }


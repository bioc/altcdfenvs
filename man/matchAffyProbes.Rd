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
                selectMatches = function(x) which(countIndex(x) > 0),
                ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{probes}{ a \code{probetable} object }
  \item{targets}{ a vector of references }
  \item{chip_type}{a name for the chip type.}
  \item{matchmm}{whether to match MM probes or not}
  \item{selectMatches}{a function to select matches (see Details).}
  \item{...}{further arguments to be passed to \code{matchPDict}.}
}
\details{
  The matching is performed by the function
  \code{\link[Biostrings]{matchPDict}}. The man page for that function
    will indicate what are the options it accepts.
  
  In the case where a large number targets are given, like when
  each target represents a possible mRNA, is it expected to have a
  largely sparse incidence matrix, that is a low number of probes
  matching every target. For that reason, only the index of matching
  probes are associated with each given target, with the function
  \code{selectMatches} giving the definition of what are
  matching probes. The default function just count anything matching,
  but the user can specify a more stringent definition if wanted.
}
\value{
  \code{mmProbes} returns a vector of MM probe sequences.
  
  \code{matchAffyProbes} returns an instance of \code{AffyProbesMatch-class}.
}
\author{ Laurent Gautier }
\note{
%As outlined in Section Details, the structure can be 

}
\seealso{ \code{\link[Biostrings]{matchPDict}} for details on how the
    matching is performed, \code{\link{AffyProbesMatch-class}}
    and \code{\link{buildCdfEnv.biostrings}}}
\examples{

library(hgu133aprobe)

filename <- system.file("exampleData", "sample.fasta",
                        package="altcdfenvs")

fasta.seq <- readFASTA(filename)

targets <- unlist(lapply(fasta.seq, "[[", "seq"))
names(targets) <- lapply(fasta.seq, "[[", "desc")
names(targets) <- sub("^>.+\\\\|(NM[^ \\\\|]+|Hs[^ \\\\|]+)\\\\| ? .+$", "",
                      names(targets))

m <- matchAffyProbes(hgu133aprobe,
                     targets,
                     "HG-U133A")


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }


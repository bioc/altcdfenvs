\name{getxy.probeseq}
\alias{getxy.probeseq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ A function to get the XY coordinates from a probes sequences
  data frame}
\description{
 A function to get the XY coordinates from a probes sequences \code{data.frame}
}
\usage{
getxy.probeseq(ppset.id = NULL, probeseq = NULL, i.row = NULL,
offset.one = TRUE, x.colname = "x", y.colname = "y")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ppset.id}{ The probe sets of interest (a vector of mode \code{character}. }
  \item{probeseq}{ The probe sequence \code{data.frame} (see details). }
  \item{i.row}{ Row indexes in the \code{data.frame} (see details). }
  \item{offset.one}{ Add one to the XY coordinates returned. }
  \item{x.colname, y.colname}{ The probe sequence packages have seen the
    names for the columns in their \code{data.frame}. This parameters
    exists to let us follow these changes.}
}
\details{
  The \code{data.frame} passed as argument \code{probeseq} is expected
  to have (at least) the following columns: \code{Probe.X},
  \code{Probe.Y} and \code{Probe.Set.Name}. When the argument
  \code{ppset.id} is not null, the probe sets
}
\value{
A \code{matrix} of two columns. The first column contains x coordinates,
while the second column contains y coordinates.
}
\author{ Laurent }

\section{Warning }{
  This function should not be confused with the methods \code{index2xy}
  and similar. Here the the XY coordinate come from a \code{data.frame}
  that stores information about an arbitrary number probes on the
  chip. (See the `probe sequence' data packages on Bioconductor, and the
  package \code{matchprobes}).
  
  The methods \code{index2xy} are meant to interact with instances of
  class \code{AffyBatch}.
}

%\seealso{ ~~objects to See Also as \code{\link{~~fun~~}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
}
\keyword{ manip }

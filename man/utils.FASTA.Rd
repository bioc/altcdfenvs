\name{read.FASTA.entry}
\alias{read.FASTA.entry}
\alias{print.FASTA}
\alias{write.FASTA}
\alias{skip.FASTA.entry}
\alias{read.n.FASTA.entries}
\alias{read.n.FASTA.headers}
\alias{read.n.FASTA.sequences}
\alias{read.n.FASTA.entries.split}
\alias{countskip.FASTA.entries}
\alias{read.FASTA.entry}
\alias{grep.FASTA.entry}
\title{Functions to work with FASTA files / connections}
\description{
  Set of function to work with biological sequences stored in FASTA format.
}
\usage{
countskip.FASTA.entries(linebreaks = 3000)
grep.FASTA.entry(pattern, con, ...)
print.FASTA(x, ...)
read.FASTA.entry(con, linebreaks = 3000)
read.n.FASTA.entries(con, n, linebreaks = 3000)
read.n.FASTA.entries.split(con, n, linebreaks = 3000)
read.n.FASTA.headers(con, n, linebreaks = 3000)
read.n.FASTA.sequences(con, n, linebreaks = 3000)
skip.FASTA.entry(con, skip, linebreaks = 3000)
write.FASTA(x, file="data.fasta", append = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{con}{a \link{connection}}
  \item{linebreaks}{(to optimize the parsing, probably safe to leave it
    as it is)}
  \item{n}{number of entries to read}
  \item{pattern}{a pattern (to be passed to the function \code{grep})}
  \item{skip}{number of entries to skip}
  \item{x}{a FASTA sequence object}
  \item{...}{Optional arguments to be forwarded to the function
    \code{print} or to the function \code{grep}}
}
%\details{
%  ~~ If necessary, more details than the __description__  above ~~
%}
\value{
  
}
\author{ Laurent Gautier }
\examples{
}
\keyword{ manip }
\keyword{ IO }
\keyword{ connection }

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
countskip.FASTA.entries(con, linebreaks = 3000)
grep.FASTA.entry(pattern, con, ...)
\method{print}{FASTA}(x, ...)
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
  \item{append}{append to the file (or not)}
  \item{con}{a \link{connection}}
  \item{file}{a file name}
  \item{linebreaks}{(to optimize the parsing, probably safe to leave it
    as it is)}
  \item{n}{number of entries to read}
  \item{pattern}{a pattern (to be passed to the function \code{grep})}
  \item{skip}{number of entries to skip}
  \item{x}{a FASTA sequence object}
  \item{...}{optional arguments to be forwarded to the function
    \code{print} or to the function \code{grep}}
}
\details{
  \code{countskip.FASTA.entries} skips the remaining FASTA entries
  currently remaining in the connection and return the count.
  \code{grep.FASTA.entry} returns the next FASTA entry in the connection
  that matches a given regular expression.
  \code{print.FASTA} prints a FASTA object.
  \code{read.FASTA.entry} reads the next FASTA entry in the connection.
  \code{read.n.FASTA.entries} reads the \code{n} next FASTA entries and
  returns a \code{list} of FASTA objects.
  \code{read.n.FASTA.entries.split} reads the \code{n} next FASTA
  entries and returns a list of two elements: headers and sequences.
  \code{read.n.FASTA.headers} reads the \code{n} next FASTA headers.
  \code{read.n.FASTA.sequences} reads the \code{n} next FASTA sequences.
  \code{skip.FASTA.entry} skips a given number of FASTA entries.
  \code{write.FASTA} write a FASTA object into a connection.
%  ~~ If necessary, more details than the __description__  above ~~
}
\value{
 The value returned depends on the function. See above.
}
\author{ Laurent Gautier }
\examples{
filename <- system.file("exampleData", "sample.fasta",
                        package="altcdfenvs")
con <- file(filename, open="r")

fasta.seq <- grep.FASTA.entry("NM_001544\\\.2", con)
close(con)


print(fasta.seq)
}
\keyword{ manip }
\keyword{ IO }
\keyword{ connection }

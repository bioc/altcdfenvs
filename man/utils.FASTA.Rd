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
\alias{grepFASTA}
\title{Functions to work with FASTA files }
\description{
  Set of function to work with biological sequences stored in FASTA format.
}
\usage{
read.FASTA.entry(con, linebreaks = 3000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{con}{a \link{connection}}
  \item{linebreaks}{(to optimize the parsing, probably safe to leave it
    as it is)}
}
\details{
  ~~ If necessary, more details than the __description__  above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ }

 ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ ~~objects to See Also as \code{\link{~~fun~~}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(con, linebreaks=3000) {
  ## read the next FASTA entry in a connection
  ## (note: the parameters 'linebreaks' should be increased
  ## for very large sequences split in more than 'linebreaks' lines)
  ## return a list of two elements:
  ##  - header: the FASTA header
  ##  - sequence: the sequence
  getnext.FASTA.header <- function(con) {
    line <- readLines(con, n=1)
    
    while(length(line) > 0) {
      if (substr(line, 1, 1) ==  ">") {
        break
      }
      line <- readLines(con, n=1)
    }
    return(line)
  }
  
  bioseq <- vector("list", length=linebreaks)
  i <- 1
  
  header <- getnext.FASTA.header(con)
  
  line <- readLines(con, n=1)
  
  while(length(line) > 0) {
    
    if (substr(line, 1, 1) == ">") {
      pushBack(line, con)
      break
    } else {
      bioseq[[i]] <- line
    }
    i <- i+1
    line <- readLines(con, n=1)
  }
  
  return(list(header=header, sequence=paste(bioseq[1:(i-1)], collapse="")))
  
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

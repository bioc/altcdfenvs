##
## Set of (hopefully convenient) functions to extract
## sequences and headers from FASTA files
##
## Laurent 2003 - under LGPL license


print.FASTA <- function(x, ...) {
  cat("FASTA sequence:\n", ...)
  if (is.null(x$header)) {
    cat("  NULL\n")
  } else {
    cat(paste("  ", substr(x$header, 1, 60), "...\n"), ...)
    cat(paste("  ", substr(x$sequence, 1, 60), "...\n"), ...)
  }
}

write.FASTA <- function(x, file="data.fasta", append = FALSE) {
  cat(x$header, file = file, "\n", sep = "", append = append)
  cat(x$sequence, file = file, "\n", sep = "", append = TRUE)
}

skip.FASTA.entry <- function(con, skip, linebreaks=3000) {
  ## skip FASTA entries in a connection
  for (i in rep(1, skip))
    read.FASTA.entry(con, linebreaks=linebreaks)
}

read.n.FASTA.entries <- function(con, n, linebreaks=3000) {
  ## read n FASTA entries in a connection
  ## return a list of length n
  r.list <- vector("list", length=n)
  for (i in seq(along=r.list))
    r.list[[i]] <- read.FASTA.entry(con)
  return(r.list)
}

read.n.FASTA.headers <- function(con, n, linebreaks=3000) {
  ## read n FASTA headers (skipping the sequences) in a connection
  ## return a vector of mode "character" of length n
  headers <- vector("character", length=n)
  for (i in seq(along=headers))
    headers[i] <- read.FASTA.entry(con)$header
  return(headers)
}

read.n.FASTA.sequences <- function(con, n, linebreaks=3000) {
  ## read n FASTA sequences(skipping the headers) in a connection
  ## return a vector of mode character
  
  seqs <- vector("character", length=n)
  for (i in seq(along=seqs))
    seqs[i] <- read.FASTA.entry(con)$sequence
  return(seqs)
}

read.n.FASTA.entries.split <- function(con, n, linebreaks=3000) {
  ## read n FASTA entries in a connection
  ## return a list of two elements:
  ##   - a vector of headers
  ##   - a vector of sequences
  
  headers <- vector("character", length=n)
  seqs <- vector("character", length=n)
  for (i in seq(along=seqs)) {
    fs <- read.FASTA.entry(con)
    headers[i] <- fs$header
    seqs[i] <- fs$sequence
  }
  r <- list(headers=headers, sequences=seqs)
  return(r)
}

countskip.FASTA.entries <- function(con, linebreaks=3000) {
  ## skip and count the remaining FASTA entries in a connection
  ## (handy to count the entries in a FASTA file)
  ## return an integer
  i <- as.integer(0)
  fs <- read.FASTA.entry(con)
  while(!identical(fs$header, character(0)) && !identical(fs$sequence, NULL)) {
    i <- i+1
    fs <- read.FASTA.entry(con, linebreaks=linebreaks)
  }
  return(i)
}

read.FASTA.entry <- function(con, linebreaks=3000) {
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
  i <- as.integer(1)
  one.integer <- as.integer(1)
  
  header <- getnext.FASTA.header(con)
  
  line <- readLines(con, n=1)
  
  while(length(line) > 0) {
    
    if (substr(line, 1, 1) == ">") {
      pushBack(line, con)
      break
    } else {
      bioseq[[i]] <- line
    }
    i <- i + one.integer
    line <- readLines(con, n=1)
  }
  if (identical(header, character(0))) {
    header <- NULL
    bioseq <- NULL
  } else {
    bioseq <- paste(bioseq[1 : (i-1)], collapse="")
  }
  r <- list(header=header, sequence=bioseq)
  class(r) <- "FASTA"  
  return(r)
}

grep.FASTA.entry <- function(pattern, con, ...) {
  ## grep the first FASTA entry with a header matching the pattern 'pattern'
  ##
  fs <- read.FASTA.entry(con)
  i <- 0
  while (! identical(fs$header, NULL)) {
    i <- i + 1
    if (length(grep(pattern, fs$header, ...)) > 0)
      break
    fs <- read.FASTA.entry(con)
  }
  attr(fs, "i") <- i
  return(fs)
}



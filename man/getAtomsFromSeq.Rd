\name{getAtomsFromSeq}
\alias{getAtomsFromSeq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function computing an atomic composition from amono acid sequence.
}
\description{
Function computing an atomic composition from (naturally occuring) amino acid sequence. 
}
\usage{
getAtomsFromSeq(seq)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seq}{
The character vector of AAString (see Biostrings package) with amino acid sequence. It should contain only letters "A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V" (1-letter symbols of 20 naturally occuring amino acids).
}
}
\details{
The atomic composition is just a summaric atomic composition of 
all amino acids composing the sequence minus (n-1) times the water molecule, where n is a length of given amino acid sequence.
}
\value{
Named list with the following fields with number of correcponding atoms (integer non-negative values):
\itemize{
\item{C},
\item{H},
\item{N},
\item{O},
\item{S}
}
}
\author{
Piotr Dittwald <piotr.dittwald@mimuw.edu.pl>
}

\examples{
  seq1 <-  "AACD"
  aC1 <- getAtomsFromSeq(seq = seq1)
  seq2 <-  AAString("ACCD")
  aC2 <- getAtomsFromSeq(seq = seq2)


}


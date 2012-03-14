\name{calculateIsotopicProbabilities}
\alias{calculateIsotopicProbabilities}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function computing probabilities of aggregated isotopic variants using BRAIN algorithm.
}
\description{
Function computing probabilities of aggregated isotopic variants for chemical components built from carbon, hydrogen, oxygen, nitrogen and sulfur (e.g. peptides). 
}
\usage{calculateIsotopicProbabilities(aC, stopOption = "nrPeaks", 
nrPeaks, coverage, abundantEstim)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{aC}{
List with fields C, H, N, O, S of integer non-negative values (if any field is ommited, then its value is set to 0).
} 
\item{stopOption}{
one of the following strings: "nrPeaks" (default), "coverage", "abundantEstim"
}
  \item{nrPeaks}{
Integer indicating the number of consecutive isotopic variants to be calculated,
starting from the monoisotopic one. This value can always be provided, even if <stopOption>
is not a default setting. In the latter case it is a hard stopping criterion.
}
  \item{coverage}{
Scalar indicating the value of the cumulative aggregated distribution. The
calculations will be stopped after reaching this value.
}
  \item{abundantEstim}{
Integer indicating the number of consecutive isotopic variants to be
calculated, starting from one after the most abundant one. All consecutive isotopic variants before the
most abundant peak are also returned.
}
}
\details{
Remember that the isotopic variants starts from the monoisotopic one.
In case of large chemical molecules, first masses may have very low abundance 
values for the lower mass aggregated values. A sufficient number of peaks should be calculated to reach most abundant isotopic variant.
}
\value{
Probabilities of aggregated isotopic variants (numeric vector)
}
\references{
[Clae] Claesen J., Dittwald P., Burzykowski T. and Valkenborg D. An efficient method to calculate 
the aggregated isotopic distribution and exact center-masses. JASMS, 2012, doi:10.1007/s13361-011-0326-2
}
\author{
Piotr Dittwald <piotr.dittwald@mimuw.edu.pl>
}
\note{
If also masses associated with the aggregated isotopic variants are needed, then the function useBRAIN should be used.
}

\seealso{
\code{\link{useBRAIN}}
 }
\examples{
  nrPeaks = 1000
  aC <-  list(C=23832, H=37816, N=6528, O=7031, S=170)  # Human dynein heavy chain  
  res <- calculateIsotopicProbabilities(aC = aC, stopOption="nrPeaks", 
nrPeaks = nrPeaks)
}


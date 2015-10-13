\name{useBRAIN}
\alias{useBRAIN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function computing probabilities of aggregated isotopic variants and their center-masses using BRAIN algorithm.
}
\description{
Function computing probabilities of isotopic variants and their aggregated masses for chemical components composed of carbon, hydrogen, oxygen, nitrogen 
and sulfur (e.g. peptides). Additionally the function returns also the monoisotopic mass and the average mass of given chemical component.
}
\usage{
useBRAIN(aC, stopOption = "nrPeaks", nrPeaks, coverage, abundantEstim)
}
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
starting from the monoisotopic one. This value can always be provided, even if <stop.option>
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
Function uses recursive formulae based on algebraic Newton-Girard identity (see [Clae]).
}
\value{
Named list with the following fields:
\itemize{
\item{isoDistr}{Probabilities of aggregated isotopic variants (numeric vector)}
\item{masses}{Aggregated masses for isotopic variants (numeric vector)}
\item{monoisotopicMass}{Monoisotopic mass (numeric)}
\item{avgMass}{Average mass - weighted average of the isotopic variants contributing to the most abundant aggregated variant (numeric)}
}
}
\references{
[Clae] Claesen J., Dittwald P., Burzykowski T. and Valkenborg D. An efficient method to calculate 
the aggregated isotopic distribution and exact center-masses. JASMS, 2012, doi:10.1007/s13361-011-0326-2
}
\author{
Piotr Dittwald <piotr.dittwald@mimuw.edu.pl>
}
\note{
Remember that the isotopic variants starts from monoisotopic one. For large chemical molecules, first masses may have very low abundances. So sufficient number of peaks should be calculated to reach most abundant isotopic variant.

If only isotopic probabilities are needed, then the function calculateIsotopicProbabilities should be used.
}

\seealso{
 \code{\link{calculateIsotopicProbabilities}}
}
\examples{
  nrPeaks = 1000  
  aC <-  list(C=23832, H=37816, N=6528, O=7031, S=170)  # Human dynein heavy chain  
  res <-  useBRAIN(aC = aC, stopOption="nrPeaks", nrPeaks = nrPeaks)
}


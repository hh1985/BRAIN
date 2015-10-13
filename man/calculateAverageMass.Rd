\name{calculateAverageMass}
\alias{calculateAverageMass}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function computing theoretical average masses.
}
\description{
Function computing the theoretical average masses for chemical components composed of carbon, hydrogen, oxygen, nitrogen and sulfur (e.g. peptides). 
}
\usage{
calculateAverageMass(aC)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{aC}{
List with fields C, H, N, O, S of integer non-negative values (if any field is ommited, then its value is set to 0).
}
}
\details{
Mass is calculated in Daltons.
}
\value{
Average mass (numeric)
}
\references{
[Clae] Claesen J., Dittwald P., Burzykowski T. and Valkenborg D. An efficient method to calculate 
the aggregated isotopic distribution and exact center-masses. JASMS, 2012, doi:10.1007/s13361-011-0326-2
}
\author{
Piotr Dittwald <piotr.dittwald@mimuw.edu.pl>
}
\seealso{
 \code{\link{useBRAIN}}
}
\examples{
  aC <-  list(C=23832, H=37816, N=6528, O=7031, S=170)  # Human dynein heavy chain
  res <- calculateAverageMass(aC = aC)
}

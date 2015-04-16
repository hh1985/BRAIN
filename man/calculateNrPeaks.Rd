\name{calculateNrPeaks}
\alias{calculateNrPeaks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function computing heuristically the required number of consecutive aggregated isotopic variants.
}
\description{
Function computing heuristically the required number of consecutive aggregated isotopic variants (starting from the monoisotopic mass).

}
\usage{
calculateNrPeaks(aC)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{aC}{
List with fields C, H, N, O, S of integer non-negative values (if any field is ommited, then its value is set to 0).
}
}
\details{
This function uses following rule of thumb: the difference between the theoretical monoisotopic mass and the theoretical average mass is 
computed and multiplied by two. Subsequently, the obtained number is rounded to the nearest integer greater than or equal to 
the multiplied difference. For small molecules, the minimal number of returned variants is five.
}
\value{
Integer number not lower than 5.
}
\references{
[Clae] Claesen J., Dittwald P., Burzykowski T. and Valkenborg D. An efficient method to calculate 
the aggregated isotopic distribution and exact center-masses. JASMS, 2012, doi:10.1007/s13361-011-0326-2
}
\author{
Jurgen Claesen <jurgen.claesen@uhasselt.be>
}
\seealso{
 \code{\link{useBRAIN}}
}
\examples{
  aC <-  list(C=23832, H=37816, N=6528, O=7031, S=170)  # Human dynein heavy chain
  res <- calculateNrPeaks(aC = aC)
}

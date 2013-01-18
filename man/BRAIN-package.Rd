\name{BRAIN-package}
\alias{BRAIN-package}
\alias{BRAIN}
\docType{package}
\title{
Implementation of BRAIN (Baffling Recursive Algorithm for Isotope distributioN calculations)
}
\description{
This package implements BRAIN (Baffling Recursive Algorithm for Isotope distributioN calculations) is described in full details by Claesen et al. The algorithm uses an algebraic approach (Viete's formulas, Newton identities [Macd]) which is 
especially useful for large molecules due to its advantageous scaling properties. 
This version of the package provides functions for calculating the aggregated isotopic distribution and center-masses masses for 
each aggregated isotopic variant for chemical components built from carbon, hydrogen, oxygen, nitrogen and sulfur (e.g. peptides). The natural abundances and 
molecular masses for stable isotopes of C, H, N, O, S are taken from IUPAC 1997 [Rosm]. 
}
\details{
\tabular{ll}{
Package: \tab BRAIN\cr
Type: \tab Package\cr
Version: \tab 1.4.0\cr
Date: \tab 2013-01-18\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
} 
\author{
Piotr Dittwald with contribution of Dirk Valkenborg and Jurgen Claesen

Maintainer: Piotr Dittwald <piotr.dittwald@mimuw.edu.pl>
}
\references{
[Clae] Claesen J., Dittwald P., Burzykowski T. and Valkenborg D. An efficient method to calculate 
the aggregated isotopic distribution and exact center-masses. JASMS, 2012, doi:10.1007/s13361-011-0326-2

[Macd] Macdonald I.G., Symmetric functions and Hall polynomials / by I. G. Macdonald. Clarendon Press; Oxford University Press, Oxford : New York, 1979.

[Rosm] K.J.R. Rosman and P.D.P. Taylor. Isotopic compositions of the elements 1997. Pure and Applied Chemistry, 70(1):217-235, 1998.
}
\keyword{ package }
\examples{
  nrPeaks = 1000
  aC <-  list(C=23832, H=37816, N=6528, O=7031, S=170)  # Human dynein heavy chain  
  res <- useBRAIN(aC = aC, nrPeaks = nrPeaks)
  iso <- res$isoDistr
  masses <- res$masses
  mono <- res$monoisotopicMass
  avgMass <- res$avgMass
}

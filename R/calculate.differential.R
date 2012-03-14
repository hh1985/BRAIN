calculateDifferential <- function(i, aCVec, nrPeaks, isoL, massCoefList){

  isoA <- isoL[[1]]
  isoB <- isoL[[2]]
  isoC <- isoL[[3]]
  isoD <- isoL[[4]]
  isoZ <- isoL[[5]]
  isoPhi <- isoL[[6]]
  isoQ0 <- isoL[[7]]
  isoNrRRoots <- isoL[[8]]
  isoNrCRoots <- isoL[[9]]



  massCoefA <- massCoefList[[1]]
  massCoefB <- massCoefList[[2]]
  massCoefC <- massCoefList[[3]]
  massCoefD <- massCoefList[[4]]
  massCoefZ <- massCoefList[[5]]
  massCoefPhi <- massCoefList[[6]]
  massCoefQ0 <- massCoefList[[7]]
  massCoefNrRRoots <- massCoefList[[8]]
  massCoefNrCRoots <- massCoefList[[9]] 
  
  listAtoms <- massCoefList[[10]]
 
      
  a <- c(isoA, massCoefA[i])

  l <- listAtoms[[i]]

  if (is.null(l$real)){
    b <- isoB
  }else{
    b <- c(isoB, massCoefB[l$real])
  }

  if (is.null(l$comp)){
    c <- isoC
    d <- isoD
    z <- isoZ
    phi <- isoPhi
  }else{
    c <- c(isoC, massCoefC[l$comp])
    d <- c(isoD, massCoefD[l$comp])
    z <- c(isoZ, massCoefZ[l$comp])
    phi <- c(isoPhi, massCoefPhi[l$comp])
  }
  
  q0 <- c(isoQ0, massCoefQ0[i])
  nrRRoots <- c(isoNrRRoots, massCoefNrRRoots[i])
  nrCRoots <- c(isoNrCRoots, massCoefNrCRoots[i])

  len <- length(aCVec) + 1
  
  aCModified <- c(aCVec, 1) - 1 * ((1:len) == i) 
  peaks <- aCVec[i] * peaksFromParameters(aCModified, a, b, c, d, z, phi, q0, nrRRoots, nrCRoots, stopOption = "nrPeaks", nrPeaks)     

}


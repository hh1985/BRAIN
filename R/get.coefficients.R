###library(PolynomF) - needed in function polynom

getCoefficients <- function(coef){


    a <- list(C=NULL, H=NULL, N=NULL, O=NULL, S=NULL)
    b <- list(C=NULL, H=NULL, N=NULL, O=NULL, S=NULL)
    c <- list(C=NULL, H=NULL, N=NULL, O=NULL, S=NULL)
    d <- list(C=NULL, H=NULL, N=NULL, O=NULL, S=NULL)
    q0 <- list(C=NULL, H=NULL, N=NULL, O=NULL, S=NULL)



    ###### atoms coefficients related to their isotopes
    a$C <- (1/coef$C[2])
    b$C <-  - (coef$C[1]/coef$C[2])
    c$C <- NULL
    d$C <- NULL
    q0$C <-coef$C[1]

    a$H <- (1/coef$H[2])
    b$H <-  - (coef$H[1]/coef$H[2])
    c$H <- NULL
    d$H <- NULL
    q0$H <-coef$H[1]

    a$N <- (1/coef$N[2])
    b$N <-  - (coef$N[1]/coef$N[2])
    c$N <- NULL
    d$N <- NULL
    q0$N <-coef$N[1]


    a$O <- (1/coef$O[3])
    b$O <- NULL
    c$O <- -coef$O[2]/(2*coef$O[3])
    d$O <- sqrt(-(coef$O[2])^2+4*coef$O[3]*coef$O[1])/(2*coef$O[3])
    q0$O <- coef$O[1]

    a$S <- (1/coef$S[5])
    b$S <- NULL
    
    rootsS <- sort(solve(polynom(coef$S)))
    c$S <- Re(rootsS[c(1,3)])
    d$S <- Im(rootsS[c(1,3)])
    #c$S and d$S are not exact values (but the error is small)
    # maybe another calculations of root than solving them numerically would be better...
    q0$S <- coef$S[1]

    ##### coefficients together 

    aVec <- c(a$C, a$H, a$N, a$O, a$S)
    bVec <- c(b$C, b$H, b$N, b$O, b$S)
    cVec <- c(c$C, c$H, c$N, c$O, c$S)
    dVec <- c(d$C, d$H, d$N, d$O, d$S)
    zVec<- complex(real = cVec, imaginary = dVec)
    phiVec <- Arg(zVec)
    q0Vec <- c(q0$C, q0$H, q0$N, q0$O, q0$S)

    nrRRoots <- c(length(b$C), length(b$H), length(b$N), length(b$O), length(b$S))
    nrCRoots <- c(length(c$C), length(c$H), length(c$N), length(c$O), length(c$S))

    list(aVec, bVec, cVec, dVec, zVec, phiVec, q0Vec, nrRRoots, nrCRoots)

}


getListAtoms <- function(){
#depents on data from file input.R (i.e. how many real and complex roots are associated with each atom)
 listC <- list(real = 1,comp = NULL)
 listH <- list(real = 2, comp = NULL)
 listN <- list(real = 3, comp = NULL)
 listO <- list(real = NULL, comp = 1)
 listS <- list(real = NULL, comp = c(2:3))

 list(listC, listH, listN, listO, listS)
}

getCoefficientsIso <- function(){

  getCoefficients(listIso) #listIso defined in input.R
}

calculateMassCoefList <- function(){
# listMass and listIso defined in input.R
 massCoef <-  list(C=NULL, H=NULL, N=NULL, O=NULL, S=NULL)

 massCoef$C <- listMass$C * listIso$C
 massCoef$H <- listMass$H * listIso$H
 massCoef$N <- listMass$N * listIso$N
 massCoef$O <- listMass$O * listIso$O
 massCoef$S <- listMass$S * listIso$S  

 massCoefL <- getCoefficients(massCoef)

 ### fragment below assume some properties of polynomial connected with atoms, i.e. how many real and complex roots do they have
 
 massCoefL[[10]] <- getListAtoms()
 massCoefL
 
}


calculateMonoisotopicMass <- function(aC){
  aCVec <- getACVec(aC)
  monomass <- c((listMass$C)[1], (listMass$H)[1], (listMass$N)[1], (listMass$O)[1], (listMass$S)[1])
  sum(aCVec * monomass)
}


calculateAverageMass <- function(aC){

	aCVec <- getACVec(aC)  
	averageC <- sum(listMass$C * listIso$C)
	averageH <- sum(listMass$H * listIso$H)
	averageN <- sum(listMass$N * listIso$N)
	averageO <- sum(listMass$O * listIso$O)
	averageS <- sum(listMass$S * listIso$S)

	averageMass <- averageC * aCVec[1]  + averageH * aCVec[2] + averageN * aCVec[3] + averageO * aCVec[4] + averageS * aCVec[5]
	averageMass
}

getACListValue <- function(val){
  if (is.null(val))
    0
  else
    val 
}


getACVec <- function(aC){
  nrC <- getACListValue(aC$C)
  nrH <- getACListValue(aC$H)
  nrN <- getACListValue(aC$N)
  nrO <- getACListValue(aC$O)
  nrS <- getACListValue(aC$S)
  c(nrC, nrH, nrN, nrO, nrS)
}



getAC <- function(aCVec){  
  list(C=aCVec[1], H=aCVec[2], N=aCVec[3], O=aCVec[4], S=aCVec[5])
}

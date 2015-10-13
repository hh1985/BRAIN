###library(PolynomF) - needed in function polynom

getCoefficients <- function(coef){

    list.null <- list(C=NULL, H=NULL, N=NULL, O=NULL, S=NULL, F = NULL, Br = NULL, P = NULL, Cl=NULL, Na = NULL, I = NULL, K=NULL)	

    a <- list.null
    b <- list.null
    c <- list.null
    d <- list.null
    q0 <- list.null

    ###### atoms coefficients related to their isotopes
    a[["C"]] <- (1/coef[["C"]][2])
    b[["C"]] <-  - (coef[["C"]][1]/coef[["C"]][2])
    c[["C"]] <- NULL
    d[["C"]] <- NULL
    q0[["C"]] <-coef[["C"]][1]

    a[["H"]] <- (1/coef[["H"]][2])
    b[["H"]] <-  - (coef[["H"]][1]/coef[["H"]][2])
    c[["H"]] <- NULL
    d[["H"]] <- NULL
    q0[["H"]] <-coef[["H"]][1]

    a[["N"]] <- (1/coef[["N"]][2])
    b[["N"]] <-  - (coef[["N"]][1]/coef[["N"]][2])
    c[["N"]] <- NULL
    d[["N"]] <- NULL
    q0[["N"]] <-coef[["N"]][1]


    a[["O"]] <- (1/coef[["O"]][3])
    b[["O"]] <- NULL
    c[["O"]] <- -coef[["O"]][2]/(2*coef[["O"]][3])
    d[["O"]] <- sqrt(-(coef[["O"]][2])^2+4*coef[["O"]][3]*coef[["O"]][1])/(2*coef[["O"]][3])
    q0[["O"]] <- coef[["O"]][1]

    a[["S"]] <- (1/coef[["S"]][5])
    b[["S"]] <- NULL
    
    rootsS <- sort(solve(polynom(coef[["S"]])))
    c[["S"]] <- Re(rootsS[c(1,3)])
    d[["S"]] <- Im(rootsS[c(1,3)])
    #c[["S"]] and d[["S"]] are not exact values (but the error is small)
    # maybe another calculations of root than solving them numerically would be better...
    q0[["S"]] <- coef[["S"]][1]




    a[["F"]] <- coef[["F"]][1]
    b[["F"]] <- NULL
    c[["F"]] <- NULL
    d[["F"]] <- NULL
    q0[["F"]] <- coef[["F"]][1]

    a[["Br"]] <- (1/coef[["Br"]][3])
    b[["Br"]] <- NULL    
    c[["Br"]] <- 0
    d[["Br"]] <- sqrt(coef[["Br"]][3]*coef[["Br"]][1])/coef[["Br"]][3]
    q0[["Br"]] <- coef[["Br"]][1]

    a[["P"]] <- coef[["P"]][1]
    b[["P"]] <- NULL
    c[["P"]] <- NULL
    d[["P"]] <- NULL
    q0[["P"]] <- coef[["P"]][1]




    a[["Cl"]] <- (1/coef[["Cl"]][3])
    b[["Cl"]] <- NULL     
    c[["Cl"]] <- 0
    d[["Cl"]] <- sqrt(coef[["Cl"]][3]*coef[["Cl"]][1])/coef[["Cl"]][3]
    q0[["Cl"]] <- coef[["Cl"]][1]

    a[["Na"]] <- coef[["Na"]][1]
    b[["Na"]] <- NULL
    c[["Na"]] <- NULL
    d[["Na"]] <- NULL
    q0[["Na"]] <- coef[["Na"]][1]

    a[["I"]] <- coef[["I"]][1]
    b[["I"]] <- NULL
    c[["I"]] <- NULL
    d[["I"]] <- NULL
    q0[["I"]] <- coef[["I"]][1]


    a[["K"]] <- (1/coef[["K"]][3])
    b[["K"]] <- NULL
    c[["K"]] <- -coef[["K"]][2]/(2*coef[["K"]][3])
    d[["K"]] <- sqrt(-(coef[["K"]][2])^2+4*coef[["K"]][3]*coef[["K"]][1])/(2*coef[["K"]][3])
    q0[["K"]] <- coef[["K"]][1]

  

    ##### coefficients together 

    aVec <- c(a[["C"]], a[["H"]], a[["N"]], a[["O"]], a[["S"]], a[["F"]], a[["Br"]], a[["P"]], a[["Cl"]], a[["Na"]], a[["I"]], a[["K"]])
    bVec <- c(b[["C"]], b[["H"]], b[["N"]], b[["O"]], b[["S"]], b[["F"]], b[["Br"]], b[["P"]], b[["Cl"]], b[["Na"]], b[["I"]], b[["K"]])
    cVec <- c(c[["C"]], c[["H"]], c[["N"]], c[["O"]], c[["S"]], c[["F"]], c[["Br"]], c[["P"]], c[["Cl"]], c[["Na"]], c[["I"]], c[["K"]])
    dVec <- c(d[["C"]], d[["H"]], d[["N"]], d[["O"]], d[["S"]], d[["F"]], d[["Br"]], d[["P"]], d[["Cl"]], d[["Na"]], d[["I"]], d[["K"]])
    zVec<- complex(real = cVec, imaginary = dVec)
    phiVec <- Arg(zVec)
    q0Vec <- c(q0[["C"]], q0[["H"]], q0[["N"]], q0[["O"]], q0[["S"]], q0[["F"]], q0[["Br"]], q0[["P"]], q0[["Cl"]], q0[["Na"]], q0[["I"]], q0[["K"]])

    nrRRoots <- c(length(b[["C"]]), length(b[["H"]]), length(b[["N"]]), length(b[["O"]]), length(b[["S"]]), length(b[["F"]]), length(b[["Br"]]), length(b[["P"]]), length(b[["Cl"]]), length(b[["Na"]]), length(b[["I"]]), length(b[["K"]]))
    nrCRoots <- c(length(c[["C"]]), length(c[["H"]]), length(c[["N"]]), length(c[["O"]]), length(c[["S"]]), length(c[["F"]]), length(c[["Br"]]), length(c[["P"]]), length(c[["Cl"]]), length(c[["Na"]]), length(c[["I"]]), length(c[["K"]]))

    lRes <- list(aVec, bVec, cVec, dVec, zVec, phiVec, q0Vec, nrRRoots, nrCRoots)
    lRes	

}


getListAtoms <- function(){
#depents on data from file input.R (i.e. how many real and complex roots are associated with each atom)
 listC <- list(real = 1,comp = NULL)
 listH <- list(real = 2, comp = NULL)
 listN <- list(real = 3, comp = NULL)
 listO <- list(real = NULL, comp = 1)
 listS <- list(real = NULL, comp = c(2:3))

 listF <- list(real = NULL, comp = NULL)
 listBr <- list(real = NULL, comp = 4)
 listP <- list(real = NULL, comp = NULL)
 listCl <- list(real = NULL, comp = 5)
 listNa <- list(real = NULL, comp = NULL)
 listI <- list(real = NULL, comp = NULL) 

 listK <- list(real = NULL, comp = 6) 


 list(listC, listH, listN, listO, listS, listF, listBr, listP, listCl, listNa, listI, listK)
}

getCoefficientsIso <- function(){

  getCoefficients(listIso) #listIso defined in input.R
}

calculateMassCoefList <- function(){
# listMass and listIso defined in input.R
 massCoef <-  list(C=NULL, H=NULL, N=NULL, O=NULL, S=NULL, F = NULL, Br = NULL, P = NULL, Cl = NULL, Na = NULL, I = NULL, K = NULL)

 massCoef[["C"]] <- listMass[["C"]] * listIso[["C"]]
 massCoef[["H"]] <- listMass[["H"]] * listIso[["H"]]
 massCoef[["N"]] <- listMass[["N"]] * listIso[["N"]]
 massCoef[["O"]] <- listMass[["O"]] * listIso[["O"]]
 massCoef[["S"]] <- listMass[["S"]] * listIso[["S"]]
 massCoef[["F"]] <- listMass[["F"]] * listIso[["F"]]  
 massCoef[["Br"]] <- listMass[["Br"]] * listIso[["Br"]]  
 massCoef[["P"]] <- listMass[["P"]] * listIso[["P"]]  
 massCoef[["Cl"]] <- listMass[["Cl"]] * listIso[["Cl"]]  
 massCoef[["Na"]] <- listMass[["Na"]] * listIso[["Na"]]
 massCoef[["I"]] <- listMass[["I"]] * listIso[["I"]]        
 massCoef[["K"]] <- listMass[["K"]] * listIso[["K"]]        	


 massCoefL <- getCoefficients(massCoef)

 ### fragment below assume some properties of polynomial connected with atoms, i.e. how many real and complex roots do they have
 
 massCoefL[[10]] <- getListAtoms()
 massCoefL
 
}


calculateMonoisotopicMass <- function(aC){
  aCVec <- getACVec(aC)

  
  monomass <- c((listMass[["C"]])[1], (listMass[["H"]])[1], (listMass[["N"]])[1], (listMass[["O"]])[1], (listMass[["S"]])[1], (listMass[["F"]])[1], (listMass[["Br"]])[1], (listMass[["P"]])[1], (listMass[["Cl"]])[1], (listMass[["Na"]])[1], (listMass[["I"]])[1], (listMass[["K"]])[1])
  sum(aCVec * monomass)
}


calculateAverageMass <- function(aC){


###this may be implemented in a loop, also e.g. calculateMassCoefList

	aCVec <- getACVec(aC)  
	averageC <- sum(listMass[["C"]] * listIso[["C"]])
	averageH <- sum(listMass[["H"]] * listIso[["H"]])
	averageN <- sum(listMass[["N"]] * listIso[["N"]])
	averageO <- sum(listMass[["O"]] * listIso[["O"]])
	averageS <- sum(listMass[["S"]] * listIso[["S"]])

	averageF <- sum(listMass[["F"]] * listIso[["F"]])
	averageBr <- sum(listMass[["Br"]] * listIso[["Br"]])
	averageP <- sum(listMass[["P"]] * listIso[["P"]])
	averageCl <- sum(listMass[["Cl"]] * listIso[["Cl"]])
	averageNa <- sum(listMass[["Na"]] * listIso[["Na"]])
	averageI <- sum(listMass[["I"]] * listIso[["I"]])

	averageK <- sum(listMass[["K"]] * listIso[["K"]])


	averageMass <- averageC * aCVec[1]  + averageH * aCVec[2] + averageN * aCVec[3] + averageO * aCVec[4] + averageS * aCVec[5] + averageF * aCVec[6]  + averageBr * aCVec[7] + averageP * aCVec[8] + averageCl * aCVec[9] + averageNa * aCVec[10] + averageI * aCVec[11] + averageK * aCVec[12]
	averageMass
}

getACListValue <- function(val){
  if (is.null(val))
    0
  else
    val 
}


getACVec <- function(aC){
  nrC <- getACListValue(aC[["C"]])
  nrH <- getACListValue(aC[["H"]])
  nrN <- getACListValue(aC[["N"]])
  nrO <- getACListValue(aC[["O"]])
  nrS <- getACListValue(aC[["S"]])

  nrF <- getACListValue(aC[["F"]])
  nrBr <- getACListValue(aC[["Br"]])
  nrP <- getACListValue(aC[["P"]])
  nrCl <- getACListValue(aC[["Cl"]])
  nrNa <- getACListValue(aC[["Na"]])
  nrI <- getACListValue(aC[["I"]])

  nrK <- getACListValue(aC[["K"]])	

  c(nrC, nrH, nrN, nrO, nrS, nrF, nrBr, nrP, nrCl, nrNa, nrI, nrK)
}



getAC <- function(aCVec){  
  list(C=aCVec[1], H=aCVec[2], N=aCVec[3], O=aCVec[4], S=aCVec[5], F = aCVec[6], Br = aCVec[7], P = aCVec[8], Cl = aCVec[9], Na = aCVec[10], I = aCVec[11], K = aCVec[12])
}

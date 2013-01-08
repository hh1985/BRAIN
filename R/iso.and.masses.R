calculateNrPeaks <- function(aC){		
	averageMass <- calculateAverageMass(aC)
	monoMass <-calculateMonoisotopicMass(aC)
	nrPeaks <- max(ceiling((averageMass-monoMass)*2),5)
	nrPeaks
}


checkOption <- function(peaks, i, stopOption, nrPeaks, coverage, abundantEstim){
##checks stop option and stops calculations if needed  
  res <- i
  
  if (stopOption == "coverage") 
    if (sum(peaks[1:i]) >= coverage)
      res <- nrPeaks
    
  
  if (stopOption == "abundantEstim")
     if (which.max(peaks[1:i]) == (i - abundantEstim))
	res <- nrPeaks
      
  res
}




peaksFromParameters <- function(aCVec, a, b, c, d, z, phi, q0, nrRRoots, nrCRoots, stopOption, nrPeaks, coverage, abundantEstim){
##used in calculatePeaks (defined below) and calculateDifferential (calculate.differential.R)

  if (is.null(nrPeaks)){      
      aC <- getAC(aCVec)      
      nrPeaks <- calculateNrPeaks(aC)      
  }
  
  if (!is.na(pmatch(stopOption, "nrPeaks"))) 
        stopOption <- "nrPeaks"
#   print(stopOption)
  STOPS <- c("nrPeaks", "coverage", "abundantEstim")
  stopIdx <- pmatch(stopOption, STOPS)
#   print(stopOption)

  
  if (is.na(stopIdx)) 
      stop("invalid stop option")
  if (stopIdx == -1) 
      stop("ambiguous stop option")

  if (is.na(nrPeaks))
      stop("specify maximal number of calculated peaks")


  if ((stopOption == "coverage"))
      if (is.null(coverage))
	stop("coverage parameter must be specified for this stop option")

  if (stopOption == "abundantEstim")
      if (is.null(abundantEstim))
	stop("abundantEstim parameter must be specified for this stop option")


  v.r <- rep(aCVec, nrRRoots)
  v.c <- rep(aCVec, nrCRoots)
  
  A.r.tmp <- -v.r 
  A.c.tmp <- - 2 * v.c


  calculate.A.r <- function(k){
  l <- length(v.r)
  if (l > 0){
    sum(A.r.tmp)
    }else{
      0
    }  
  }

  calculate.A.c <- function(k){ 
  l <- length(v.c)
  if (l > 0){
    sum(A.c.tmp * cos(k * phi)) 
  }else{
    0
  }  
  }

  calculate.A <- function(k){
    calculate.A.r(k) + calculate.A.c(k)
  }

  ### start of the body of peaks.from.parameters
  ### declarations
  peaks <- numeric(nrPeaks)
  A <- numeric(nrPeaks)

  ### parameters initialization
  peaks[1] <- exp(sum(aCVec*log(q0)))

  ### main loop for calculating peaks
  i <- 2
  length <- 1

  while (i <= nrPeaks){     
#   for (i in 2:nr.peaks){
    A.r.tmp <- A.r.tmp / b
    A.c.tmp <- A.c.tmp / sqrt(c^2 + d^2)
    A[i-1] <- calculate.A(i-1)    
    peaks[i] <- sum(peaks[1:(i-1)] * A[(i-1):1])/(i-1)
    i <- checkOption(peaks, i, stopOption, nrPeaks, coverage, abundantEstim)    
    i <- i + 1
    length <- length + 1   
  } 
  peaks[1:length]
  ### result stored in 'peaks'
}


calculateIsotopicProbabilities <- function(aC, stopOption="nrPeaks", nrPeaks=NULL, coverage=NULL, abundantEstim=NULL){

  calculatePeaks <- function(aCVec, stopOption, nrPeaks, coverage, abundantEstim){
    ### local function

    

    prefixList <- getCoefficientsIso()
    
    prefixPeaks <- peaksFromParameters(aCVec, prefixList[[1]], prefixList[[2]], prefixList[[3]], prefixList[[4]], prefixList[[5]], prefixList[[6]], prefixList[[7]], prefixList[[8]], prefixList[[9]], stopOption=stopOption, nrPeaks=nrPeaks, coverage=coverage, abundantEstim=abundantEstim)
  }
  #### start of calculatePeaks
  aCVec <- getACVec(aC)  
  isoPeaks <- calculatePeaks(aCVec, stopOption, nrPeaks, coverage, abundantEstim) 
}

useBRAIN <- function(aC = aC, stopOption="nrPeaks", nrPeaks=NULL, coverage=NULL, abundantEstim=NULL){
#   print(stopOption)
 iso <- calculateIsotopicProbabilities(aC = aC, stopOption, nrPeaks, coverage, abundantEstim)
  
 isoList <- getCoefficientsIso()
 massCoefList <- calculateMassCoefList()

 nrPeaks <- length(iso) 
 nrConsideredAtoms <- length(isoList[[1]])
 



 aCVec <- getACVec(aC)
 lRes <- lapply((1:nrConsideredAtoms)[aCVec > 0], calculateDifferential, aCVec, nrPeaks, isoList, massCoefList)
 
 for (i in 1:length(lRes)){
      if (i == 1){
	masses <- lRes[[1]]
     }else{
	masses <- masses + lRes[[i]]
    }      
 }
 masses <- masses/iso
 mono = masses[1]
 avgMass <- calculateAverageMass(aC)

 list(isoDistr=iso, masses=masses, monoisotopicMass = mono, avgMass = avgMass)
}


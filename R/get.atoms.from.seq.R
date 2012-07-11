WATER <- c(0,2,0,1,0)

getAtomsFromAminoacid <- function(c){  
  if (c == "I") res <-c(6, 13, 1, 2, 0)
  if (c == "L") res <-c(6, 13, 1, 2, 0)
  if (c == "K") res <-c(6, 14, 2, 2, 0)
  if (c == "M") res <-c(5, 11, 1, 2, 1)
  if (c == "F") res <-c(9, 11, 1, 2, 0)
  if (c == "T") res <-c(4, 9, 1, 3, 0)
  if (c == "W") res <-c(11, 12, 2, 2, 0)
  if (c == "V") res <-c(5, 11, 1, 2, 0)
  if (c == "R") res <-c(6, 14, 4, 2, 0)
  if (c == "H") res <-c(6, 9, 3, 2, 0)
  if (c == "A") res <-c(3, 7, 1, 2, 0)
  if (c == "N") res <-c(4, 8, 2, 3, 0)
  if (c == "D") res <-c(4, 7, 1, 4, 0)
  if (c == "C") res <-c(3, 7, 1, 2, 1)
  if (c == "E") res <-c(5, 9, 1, 4, 0)
  if (c == "Q") res <-c(5, 10, 2, 3, 0)
  if (c == "G") res <-c(2, 5, 1, 2, 0)
  if (c == "P") res <-c(5, 9, 1, 2, 0)
  if (c == "S") res <-c(3, 7, 1, 3, 0)
  if (c == "Y") res <-c(9, 11, 1, 3, 0) 
  res
}


getAtomsFromSeq <- function(seq){
  seqCharacter <- as.character(seq) #if AAString, then converted to character vector
  atoms <- rep(0,5)
  len <- 0
  seq2 <- strsplit(seqCharacter, "")[[1]]
  for (c in seq2){
    atoms <- atoms + getAtomsFromAminoacid(c)
    len <- len + 1
  }
#   subtract water masses
  atoms <- atoms - (len-1) * WATER	
  atomsList <- getAC(atoms)
  atomsList[is.na(atomsList)] <- NULL
  atomsList
}

  
makesubscript <- function(word){
  lettercount <- 0
  letters = strsplit(word, "")
  splitidx <- list()
  
  currentchar = ""
  prevchar = ""
  # For each letter in each formula
  for (j in 1:lengths(letters)){
    lettercount <- lettercount + 1
    
    prevchar = currentchar
    currentchar = letters[[1]][j]
    # print(paste0("letter ", lettercount, ": ", currentchar))
    
    # If current char is number and prev char is a letter
    if (suppressWarnings(isTruthy(as.numeric(currentchar))) &&
        suppressWarnings(!(isTruthy(as.numeric(prevchar))))
    ){
      splitidx <- c(splitidx, lettercount)
      # If current char is a letter followed by number
    } else if (suppressWarnings(!(isTruthy(as.numeric(currentchar)))) &&
               suppressWarnings(isTruthy(as.numeric(prevchar)))){
      splitidx <- c(splitidx, lettercount)
    }
  }
  # print(paste0("Word: ", word))
  # print(paste0("Idx: ", splitidx))
  # cat("\n")
  
  
  finalformula = ""
  numidx = 1
  # For each index put before it either a <sub> or </sub> alternating for subscripts
  for (j in 1:lengths(letters)){
    
    if (is.element(j, splitidx)){
      # If even letter use </sub>
      if (numidx%%2 == 0){
        finalformula = paste0(finalformula, "</sub>", letters[[1]][j])
        
      } else {
        finalformula = paste0(finalformula, "<sub>", letters[[1]][j])
      }
      numidx = numidx + 1
    } else {
      finalformula = paste0(finalformula, letters[[1]][j])
    }
    
  }
  # if (numidx%%2 == 0){
  #   finalformula = paste0(finalformula, "</sub>")
  # }
  
  # print(paste0("Finalformula: ", finalformula))
  return(finalformula)
}
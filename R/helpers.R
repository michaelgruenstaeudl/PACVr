#!/usr/bin/R
#contributors = c("Michael Gruenstaeudl","Nils Jenke")
#email = "m.gruenstaeudl@fu-berlin.de", "nilsj24@zedat.fu-berlin.de"
#version = "2020.01.17.1800"

HistCol <- function(cov, threshold,relative, logScale) {
  
  # Function to generate color vector for histogram data
  # ARGS:
  #       cov:       data.frame of coverage
  #       threshold: numeric value of a specific threshold
  # RETURNS:
  #   color vector
  # Error handling
  if (!is.numeric(threshold) | threshold < 0) {
    warning("User-defined coverage depth threshold must be >=1.")
    stop()
  }
  
  
  if(relative == TRUE & logScale){
    threshold <- mean(cov[,4]) + log(threshold)
    
  }else if(relative == TRUE){
    threshold <- mean(cov[,4]) * threshold
    
  }
  color <- rep("black",nrow(cov))
  ind   <- as.numeric(cov[ ,4]) <= threshold
  color <- replace(color,ind,"red")
  return(color)
}

boolToDeci <- function(boolList) {
  
  out = 0
  boolList <- rev(boolList)
  for (i in 1:length(boolList)){
    out = out + boolList[i]*(2^(i-1))
  }
  return(out)
}
#!/usr/bin/R
#contributors = c("Michael Gruenstaeudl", "Nils Jenke")
#email = "m.gruenstaeudl@fu-berlin.de", "nilsj24@zedat.fu-berlin.de"
#version = "2019.03.15.1800"


In.Interval <- function(end, interval){
  # Check if a given interval is part of another given interval
  # ARGS:
  #   x: Numerical value or vector; begin of interval to check
  #   interval: End Value of an interval
  # RETURNS:
  #   Boolean vector
  # Error handling
  if(!is.numeric(end) == TRUE | !is.numeric(interval) == TRUE) {
    stop("chromEnd and geneEnd must be numerical.\n")
  }
  return(interval >= end) 
} 


Rename_Df <- function(df,str_vec = NULL) {
  # Rename the given df to make it accessible for other functions
  # ARGS:
  #   df: data frame to rename
  #   str_vec: optional colnames
  # RETURNS:
  #   data frame with new colnames
  dfNames <- c("Chromosome","chromStart","chromEnd")
  if(is.null(str_vec) == FALSE) {
    dfNames     <- c(dfNames,str_vec)
  }
  colnames(df)  <- dfNames
  df            <- df[order(df[2]), ]
  row.names(df) <- 1:nrow(df) 
  return(df)
}


HistCol <- function(cov, threshold) {
  # Function to generate color vector for histogram data
  # ARGS:
  #       cov:       data.frame of coverage
  #       threshold: numeric value of a specific threshold
  # RETURNS:
  #   color vector
  # Error handling
  if (!is.numeric(threshold) | threshold < 0) {
    stop("threshold has to be greater than zero")
  }
  color <- rep("black",nrow(cov))
  ind   <- cov[ ,4] < threshold
  color <- replace(color,ind,"red")
  return(color)
}


AssignRegionInfo <- function(toAssign, genomic_regions) {
  # Assign genomic region names to a data.frame depending on chromosome start and chromosome end
  # ARGS:
  #   toAssign: data.frame that contains Chromosome names, chromosome start and chromosome end 
  #   genomic_regions: data.frame that contains genomic regions, chromosome start and chromosome end
  # RETURNS:
  #   data.frame with assigned genomic regions
  if (ncol(genomic_regions) < 3 | ncol(toAssign) < 3) {
    stop("Dataframe needs atleast 3 columns containing chromosome name, chromosome start and chromosome end.\n")
  }
  region.rows   <- colnames(genomic_regions)
  assign.rows   <- colnames(toAssign)
  region.names  <- genomic_regions[ ,1]
  assign.Start  <- toAssign[ ,2]
  genomic.Start <- genomic_regions[ ,2]
  assign.End    <- toAssign[ ,3]
  genomic.End   <- genomic_regions[ ,3]
  levels(toAssign[ ,1]) <- c(region.names,levels(toAssign[ ,1]))
  for (i in length(region.names):1){
    indices <- In.Interval(assign.End,genomic.End[i])
    toAssign[, which(assign.rows == "Chromosome")] <- replace(toAssign[ , which(assign.rows == "Chromosome")], indices, region.names[i])
  }
  return(toAssign)
}


AdjustRegionLocation <- function(toShift, genomic_regions){
  # Shift genome regions so that they fit RCircos validation
  # ARGS:
  #   toShift: data.frame with region names, chromosome begin and chromosome end
  #   genomic_regions: data.frame with region names, chromosome begin and chromosome end
  # RETURNS:
  #   data.frame with shifted regions
  if (ncol(genomic_regions) < 3 | ncol(toShift) < 3) {
    stop("Dataframe needs at least 3 columns containing chromosome name, chromosome start and chromosome end.\n")
  }
  i <- 0
  j <- 0
  if (toShift[1,2] >  0){
    i <- 1
  }
  else if(genomic_regions[1,2] > 0) {
    j <- 1
  }
  toShift[toShift[ ,1] == "LSC",2] = toShift[toShift[ ,1] == "LSC",2] - i - j
  toShift[toShift[ ,1] == "LSC",3] = toShift[toShift[ ,1] == "LSC",3] - i - j
  toShift[toShift[ ,1] == "IRb",2] = toShift[toShift[ ,1] == "IRb",2] - (genomic_regions[1,3] + i + j)
  toShift[toShift[ ,1] == "IRb",3] = toShift[toShift[ ,1] == "IRb",3] - (genomic_regions[1,3] + i + j)
  toShift[toShift[ ,1] == "SSC",2] = toShift[toShift[ ,1] == "SSC",2] - (genomic_regions[2,3] + i + j)
  toShift[toShift[ ,1] == "SSC",3] = toShift[toShift[ ,1] == "SSC",3] - (genomic_regions[2,3] + i + j)
  toShift[toShift[ ,1] == "IRa",2] = toShift[toShift[ ,1] == "IRa",2] - (genomic_regions[3,3] + i + j)
  toShift[toShift[ ,1] == "IRa",3] = toShift[toShift[ ,1] == "IRa",3] - (genomic_regions[3,3] + i + j)
  toShift[which(toShift[ ,2] < 0),2] = 0
  return(toShift)
}

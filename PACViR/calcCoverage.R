#!/usr/bin/R
#contributors = c("Michael Gruenstaeudl","Nils Jenke")
#email = "m.gruenstaeudl@fu-berlin.de", "nilsj24@zedat.fu-berlin.de"
#version = "2019.05.24.1700"

source("helpers.R")
#source("/home/michael_science/git/michaelgruenstaeudl_PACViR/PACViR/helpers.R")

CovCalc <- function(bamFile, windowSize=250, outDir="./PACViR_output/", mosdepthCmd="mosdepth"){
  # Calculates coverage of a given bam file and stores data in data.frame format
  # ARGS:
  #     bamFile: bam file to calculate coverage
  #     windowSize: numeric value to specify the coverage calculation window
  #     outDir: output directory
  #     mosdepthCmd: path to mosdepth
  # RETURNS:
  #     data.frame with region names, chromosome start, chromosome end and coverage calcucation
  if (!is.numeric(windowSize) | windowSize < 0) {
    stop("windowSize has to be greater than zero")
  }
  get_os <- Sys.info()[1]
  
  ### TODO: Windows console
  
  if(get_os == "Windows"){
    system(paste("md", paste(outDir, "\ncoverage-output", sep = "")))
    system(paste(mosdepthCmd, "--by", windowSize, paste(outDir, "\ncoverage-output\noutput", sep = ""), bamFile))
    system(paste("gzip -df", paste(outDir, "\ncoverage-output\noutput.regions.bed.gz", sep = "")))
    cov <- Rename_Df(cov, "coverage")
  } else {
    system(paste("mkdir -p", paste(outDir, "/coverage-output" ,sep = "")))
    system(paste(mosdepthCmd, "--by", windowSize, paste(outDir, "/coverage-output/output", sep = ""), bamFile))
    system(paste("gzip -df", paste(outDir, "/coverage-output/output.regions.bed.gz", sep = "")))
    cov <-read.table(paste(outDir, "coverage-output/output.regions.bed", sep = ""))
    cov <- Rename_Df(cov, "coverage")
  }
  return(cov)
}

SplitCovAtRegionBorders <- function(covData, regionData) {
  # Function to split coverage data that occur in two different regions at
  # the region borders
  # ARGS:
  #   covData: dataframe with gene data
  #   regionData: dataframe with region data
  # RETURNS:
  #   covData dataframe with splitted covDatas
  for (i in 1:nrow(regionData)) {
    for (j in 1:nrow(covData)) {
      if (as.integer(covData[j,2]) >= as.integer(regionData[i,2]) & 
          as.integer(covData[j,3]) >  as.integer(regionData[i,3]) & 
          as.integer(covData[j,2]) <  as.integer(regionData[i,3])){
        covData[nrow(covData)+1,] <- c(as.character(covData[j,1]), regionData[i,3]+1, covData[j,3], covData[j,4])
        covData[j,1] <- regionData[i,1]
        covData[j,3] <- regionData[i,3]
        covData[j,4] <- covData[j,4]
      }
    }
  }
  covData      <- covData[order(as.integer(covData[,2])), ]
  covData[ ,2] <- as.integer(covData[ ,2])
  covData[ ,3] <- as.integer(covData[ ,3])
  return(covData)
}


adjustCoverage <- function(cov, regions) {
  # Shift of coverage regions so that they fit RCircos validation
  # ARGS:
  #   cov: data.frame with region names, chromosome begin, chromosome end and coverage
  #   regions: data.frame with region names, chromosome begin and chromosome end
  # RETURNS:
  #   data.frame with shifted regions
  irb <- as.numeric(cov[cov[ ,1] == 'IRb',2][1])
  ssc <- as.numeric(cov[cov[ ,1] == 'SSC',2][1])
  ira <- as.numeric(cov[cov[ ,1] == 'IRa',2][1])
  cov[cov[ ,1] == 'IRb',2] = as.numeric(cov[cov[ ,1] == 'IRb',2] - irb)
  cov[cov[ ,1] == 'SSC',2] = as.numeric(cov[cov[ ,1] == 'SSC',2] - ssc)
  cov[cov[ ,1] == 'IRa',2] = as.numeric(cov[cov[ ,1] == 'IRa',2] - ira)
  cov[cov[ ,1] == 'IRb',3] = as.numeric(cov[cov[ ,1] == 'IRb',3] - irb)
  cov[cov[ ,1] == 'SSC',3] = as.numeric(cov[cov[ ,1] == 'SSC',3] - ssc)
  cov[cov[ ,1] == 'IRa',3] = as.numeric(cov[cov[ ,1] == 'IRa',3] - ira)
  cov[cov[ ,1] == 'LSC',3][length(cov[cov[ ,1] == 'LSC',3])] = as.numeric(regions[,3][1])
  cov[cov[ ,1] == 'IRb',3][length(cov[cov[ ,1] == 'IRb',3])] = as.numeric(regions[,3][2])
  cov[cov[ ,1] == 'SSC',3][length(cov[cov[ ,1] == 'SSC',3])] = as.numeric(regions[,3][3])
  cov[cov[ ,1] == 'IRa',3][length(cov[cov[ ,1] == 'IRa',3])] = as.numeric(regions[,3][4])
  chromosome <- as.character(cov[,1])
  chromStart <- as.numeric(cov[,2])
  chromEnd <- as.numeric(cov[,3])
  coverage <- as.numeric(cov[,4])
  cov <- data.frame(chromosome,chromStart,chromEnd,coverage)
  return(cov)
}

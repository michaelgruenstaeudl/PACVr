#!/usr/bin/R
#contributors = c("Michael Gruenstaeudl","Nils Jenke")
#email = "m.gruenstaeudl@fu-berlin.de", "nilsj24@zedat.fu-berlin.de"
#version = "2018.10.23.1700"

#source("helpers.R")
source("/home/michael_science/git/michaelgruenstaeudl_PACViR/PACViR/helpers.R")

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
  system(paste("mkdir -p", paste(outDir, "/coverage-output" ,sep = "")))
  system(paste(mosdepthCmd, "--by", windowSize, paste(outDir, "/coverage-output/output", sep = ""), bamFile))
  #system(paste("gzip -df", paste(outDir, "/coverage-output/output.windows.bed.gz", sep = "")))
  system(paste("gzip -df", paste(outDir, "/coverage-output/output.regions.bed.gz", sep = "")))
  #cov <-read.table(paste(outDir, "/coverage-output/output.windows.bed", sep = ""))
  cov <-read.table(paste(outDir, "/coverage-output/output.regions.bed", sep = ""))
  cov <- Rename_Df(cov, "coverage")
  return(cov)
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
  #cov = cov[-nrow(cov),]
  cov[cov[ ,1] == 'IRa',3][length(cov[cov[ ,1] == 'IRa',3])] = as.numeric(regions[,3][4])
  chromosome <- as.character(cov[,1])
  chromStart <- as.numeric(cov[,2])
  chromEnd <- as.numeric(cov[,3])
  coverage <- as.numeric(cov[,4])
  cov <- data.frame(chromosome,chromStart,chromEnd,coverage)
  return(cov)
}

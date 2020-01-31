#!/usr/bin/R
#contributors = c("Michael Gruenstaeudl","Nils Jenke")
#email = "m.gruenstaeudl@fu-berlin.de", "nilsj24@zedat.fu-berlin.de"
#version = "2020.01.17.1800"

CovCalc <- function(bamFile, windowSize=250, tmpDir, logScale, mosdepthCmd="mosdepth"){
  
  # Calculates coverage of a given bam file and stores data in data.frame format
  # ARGS:
  #     bamFile: bam file to calculate coverage
  #     windowSize: numeric value to specify the coverage calculation window
  #     output: name and directory of output file
  #     mosdepthCmd: path to mosdepth
  # RETURNS:
  #     data.frame with region names, chromosome start, chromosome end and coverage calcucation
  if (!is.numeric(windowSize) | windowSize < 0) {
    warning("User-selected window size must be >= 1.")
    stop()
  }
  system2(command=mosdepthCmd, args=c("--by", windowSize,paste(tmpDir, .Platform$file.sep, "coverage", sep = ""), bamFile))
  system2(command="gzip", args=c("-df", paste(tmpDir, .Platform$file.sep, "coverage.regions.bed.gz", sep = "")))
  cov <-read.table(paste(tmpDir, .Platform$file.sep, "coverage.regions.bed", sep = ""))
  colnames(cov) <- c("Chromosome","chromStart","chromEnd","coverage")
  cov$Chromosome <- ""
  
  if(logScale==TRUE){
    cov$coverage <- log(cov$coverage)
  }
  return(cov)
}


DummyCov <- function(chromName, regions, windowSize=250){
  
  # Generates data.frame with dummy coverage values
  # ARGS:
  #     raw_regions
  #     windowSize
  # RETURNS:
  #     data.frame with region names, chromosome start, chromosome end and coverage calcucation
  if (!is.numeric(windowSize) | windowSize < 0) {
    warning("User-selected window size must be >= 1.")
    stop()
  }
  chromLen <- max(as.integer(regions[, "chromEnd"]))
  chromStart <- seq.int(0, chromLen, windowSize)
  chromEnd <- c(seq.int(windowSize, chromLen, windowSize), chromLen)
  Chromosome <- rep("", length(chromStart))
  coverage <- rep(1, length(chromStart))
  cov <- data.frame(Chromosome, chromStart, chromEnd, coverage)
  return(cov)
}

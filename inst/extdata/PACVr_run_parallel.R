#!/usr/bin/env RScript
#contributors=c("Gregory Smith", "Nils Jenke", "Michael Gruenstaeudl")
#email="m_gruenstaeudl@fhsu.edu"
#version="2024.02.14.1130"

library(tcltk)
library(foreach)
library(doParallel)

inDir <- tcltk::tk_choose.dir(default = "~", caption = "Select directory")
inFiles <- list.files(path=inDir, pattern='annotated.gb', full.names=TRUE) 	 # User should modify pattern to target GenBank input files

# Set up parallel backend to use multiple processors
cores_avail <- floor((detectCores()-1)/2)  # only using half of available cores
print(paste("Available cores:", cores_avail))
cl <- makeCluster(cores_avail, type="FORK", outfile=paste0(inDir, "/", "MPI_runs.log"))
registerDoParallel(cl)

run_PACVr <- function(f) {
  inFileDir <- paste0(dirname(f), "/")
  accNum <- gsub("_annotated.gb", "", basename(f))
  print(paste("Processing", accNum))
  
  gbkFile <- paste0(inFileDir, basename(f))
  bamFile <- paste0(inFileDir, accNum, "_mapping_OneMoreLocations.sorted.bam")
  
  tmpName <- paste0(inFileDir, accNum, "_CoverageViz_IRCheckFALSE")
  #sink(paste0(tmpName, ".log"))
  outFile <- paste0(tmpName, ".png")
  PACVr.complete(gbkFile, bamFile, windowSize=250, logScale=FALSE, 
                 threshold=0.5, relative=TRUE, textSize=0.5, 
                 output=outFile)
  #sink()
  
#  tmpName <- paste0(inFileDir, accNum, "_CoverageViz_IRCheck0")
#  sink(paste0(tmpName, ".log"))
#  outFile <- paste0(tmpName, ".png")
#  PACVr.complete(gbkFile, bamFile, windowSize=250, logScale=FALSE, 
#                 threshold=0.5, relative=TRUE, textSize=0.5, 
#                 IRCheck=0, output=outFile)
#  sink()
  
#  tmpName <- paste0(inFileDir, accNum, "_CoverageViz_IRCheck1")
#  sink(paste0(tmpName, ".log"))
#  outFile <- paste0(tmpName, ".png")
#  PACVr.complete(gbkFile, bamFile, windowSize=250, logScale=FALSE, 
#                 threshold=0.5, relative=TRUE, textSize=0.5, 
#                 IRCheck=1, output=outFile)
#  sink()
}

foreach(i=1:length(inFiles), .packages=c("PACVr")) %dopar% {
  run_PACVr(inFiles[i])
}

# Stop parallelization
stopCluster(cl)

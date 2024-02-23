#!/usr/bin/env RScript
#contributors=c("Gregory Smith", "Nils Jenke", "Michael Gruenstaeudl")
#email="m_gruenstaeudl@fhsu.edu"
#version="2024.02.23.0413"

PACVr.read.gb <- function(gbkFile) {
  gbkRaw <- getGbkRaw(gbkFile)
  if (is.null(gbkRaw$char)) {
    return(NULL)
  }
  gbkData <- read.gbWithHandling(gbkRaw)
  return(gbkData)
}

PACVr.verboseInformation <- function(gbkData,
                                     coverageRaw,
                                     analysisSpecs,
                                     plotSpecs) {
  sampleName <- gbkData$sampleName
  verbosePath <- getVerbosePath(sampleName,
                                plotSpecs)
  printCovStats(coverageRaw,
                gbkData$genes,
                gbkData$quadripRegions,
                sampleName,
                analysisSpecs,
                verbosePath)
  if (analysisSpecs$isSyntenyLine) {
    checkIREquality(gbkData$seq,
                    gbkData$quadripRegions,
                    verbosePath,
                    sampleName)
  }
  logger::log_info('Verbose output saved in `{verbosePath}`')
}

PACVr.visualizeWithRCircos <- function(gbkData,
                                       coverage,
                                       analysisSpecs,
                                       plotSpecs) {
  logger::log_info('Generating a visualization of the sequencing coverage')
  isOutput <- plotSpecs$isOutput

  if (isOutput) {
    createVizFile(plotSpecs)
  }

  visualizeWithRCircos(
    gbkData,
    coverage,
    analysisSpecs,
    plotSpecs
  )

  if (isOutput) {
    dev.off()
    logger::log_info('Visualization saved as `{plotSpecs$output}`')
  }
}

#' @title Execute the complete pipeline of \pkg{PACVr}
#' @description This function executes the complete pipeline of \pkg{PACVr} 
#' via a single command.
#'
#' @param gbkFile a character string that specifies the name of, and path to, 
#' the GenBank input file; alternatively, a character string of GenBank data 
#' @param bamFile a character string that specifies the name of, and path to, 
#' the BAM input file
#' @param windowSize a numeric value that specifies window size in which the 
#' coverage is calculated
#' @param logScale a boolean that specifies if the coverage depth is to be 
#' log-transformed before visualizing it
#' @param threshold a numeric value that specifies the threshold for plotting 
#' coverage depth bars in red as opposed to the default black
#' @param IRCheck a numeric value that specifies if tests 
#' for IRs of input genome should be performed, and - if IRs are present - 
#' which line type to be used for visualizing gene synteny between IRs;
#' 0 = IR presence test but no synteny visualization, 
#' 1 = IR presence test and synteny visualization, with ribbon lines between IRs,
#' 2 = IR presence test and synteny visualization, with solid lines between IRs,
#' otherwise = neither IR presence test nor synteny visualization
#' @param relative a boolean that specifies whether the threshold is a relative 
#' value of the average coverage instead of an absolute value
#' @param textSize a numeric value that specifies the relative font size of the 
#' text element in the visualization
#' @param verbose a boolean, that when TRUE, generates additional files with
#' detailed genomic region information
#' @param output a character string that specifies the name of, and path to, 
#' the output file
#' @return A file in pdf format containing a circular visualization of the 
#' input plastid genome and its sequence reads. 
#' As a function, returns 0 in case of visualization success.
#' @export
#' @examples
#' \dontrun{
#' gbkFile <- system.file("extdata", "NC_045072/NC_045072.gb", package="PACVr")
#' bamFile <- system.file("extdata", "NC_045072/NC_045072_subsampled.bam", package="PACVr")
#' outFile <- paste(tempdir(), "/NC_045072__all_reads.pdf", sep="")
#' PACVr.complete(gbkFile=gbkFile, bamFile=bamFile, windowSize=250, logScale=FALSE,
#'                threshold=0.5, IRCheck=1, relative=TRUE, textSize=0.5,
#'                verbose=FALSE, output=outFile)
#' }
#' \dontrun{
#' gbkFile <- system.file("extdata", "MG936619/MG936619.gb", package="PACVr")
#' bamFile <- system.file("extdata", "MG936619/MG936619_subsampled.bam", package="PACVr")
#' outFile <- paste(tempdir(), "/MG936619_CoverageViz.pdf", sep="")
#' PACVr.complete(gbkFile=gbkFile, bamFile=bamFile, windowSize=50, logScale=FALSE,
#'                threshold=0.5, IRCheck=NA, relative=TRUE, textSize=0.5,
#'                verbose=FALSE, output=outFile)
#' }
	
PACVr.complete <- function(gbkFile,
                           bamFile,
                           windowSize=250,
                           logScale=FALSE,
                           threshold=0.5,
                           IRCheck=NA,
                           relative=TRUE,
                           textSize=0.5,
                           verbose=FALSE,
                           output=NA) {
  ######################################################################
  read.gbData <- PACVr.read.gb(gbkFile)
  analysisSpecs <- getAnalysisSpecs(IRCheck,
                                    windowSize)
  gbkData <- PACVr.gbkData(read.gbData,
                           analysisSpecs)
  rm(read.gbData)
  gc()
  if (is.null(gbkData)) {
    logger::log_fatal('Unsuccessful.')
    return(-1)
  }

  ###################################
  plotSpecs <- getPlotSpecs(logScale,
                            threshold,
                            relative,
                            textSize,
                            output)

  ###################################
  coverage <- PACVr.calcCoverage(bamFile,
                                 analysisSpecs$windowSize,
                                 plotSpecs$logScale)

  ###################################
  if (verbose) {
    PACVr.verboseInformation(gbkData,
                             coverage$raw,
                             analysisSpecs,
                             plotSpecs)
  }
  
  ###################################
  PACVr.visualizeWithRCircos(gbkData,
                             coverage$plot,
                             analysisSpecs,
                             plotSpecs)

  ######################################################################
  logger::log_success('Done.')
  ######################################################################
  return(0)
}

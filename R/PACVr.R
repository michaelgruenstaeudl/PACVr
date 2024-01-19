#!/usr/bin/env RScript
#contributors=c("Gregory Smith", "Nils Jenke", "Michael Gruenstaeudl")
#email="m_gruenstaeudl@fhsu.edu"
#version="2024.01.07.2200"

PACVr.read.gb <- function(gbkFile) {
  gbkChar <- getGbkChar(gbkFile)
  if (is.null(gbkChar)) {
    return(NULL)
  }
  gbkData <- read.gbWithHandling(gbkFile, gbkChar)
  return(gbkData)
}

PACVr.parseName <- function (gbkData) {
  return(read.gbSampleName(gbkData))
}

PACVr.parseRegions <- function (gbkData, gbkDataDF) {
  raw_regions <- ExtractAllRegions(gbkDataDF)
  regions <- fillDataFrame(gbkData, raw_regions)
  return(regions)
}

PACVr.parseSource <- function(gbkDataDF) {
  return(parseSource(gbkDataDF))
}

PACVr.parseGenes <- function (gbkDataDF) {
  # This function parses the genes of a GenBank file
  genes <- ExtractAllGenes(gbkDataDF)
  return(genes)
}

PACVr.calcCoverage <-
  function (bamFile, windowSize=250) {
    coverage <- CovCalc(bamFile, windowSize)
    return(coverage)
  }

PACVr.generateIRGeneData <- function(genes, regions,
                                     syntenyLineType) {
  # Parse GenBank file
  if ("IRb" %in% regions[, 4] &&
      "IRa" %in% regions[, 4]) {
    linkData <- GenerateIRSynteny(genes, syntenyLineType)
    return(linkData)
  }
  return(-1)
}

PACVr.verboseInformation <- function(gbkData,
                                     bamFile,
                                     genes,
                                     regions,
                                     output) {
  sampleName <- read.gbSampleName(gbkData)
  # Step 1. Check ...
  if (!is.na(output)) {
    outDir <- dirname(output)
    tmpDir <- file.path(outDir, 
            paste(sampleName["sample_name"],
            ".tmp",
            sep=""))
  } else {
    tmpDir <-
      file.path(".", paste(sampleName["sample_name"],
                   ".tmp",
                   sep=""))
  }
  # Step 2. Check ...
  if (dir.exists(tmpDir) == FALSE) {
    dir.create(tmpDir)
  }
  # Step 3. Write output
  writeTables(regions, bamFile, genes, tmpDir, sampleName)
  checkIREquality(gbkData, regions, tmpDir, sampleName)
}

PACVr.visualizeWithRCircos <- function(gbkData,
                                       genes,
                                       regions,
                                       coverage,
                                       windowSize,
                                       logScale,
                                       threshold,
                                       relative,
                                       linkData,
                                       syntenyLineType,
                                       textSize) {
  # Step 1. Generate plot title
  plotTitle <- read.gbPlotTitle(gbkData)
  # Step 2. Visualize
  visualizeWithRCircos(
    plotTitle,
    genes,
    regions,
    coverage,
    windowSize,
    threshold,
    logScale,
    relative,
    linkData,
    syntenyLineType,
    textSize
  )
}

#' @title Execute the complete pipeline of \pkg{PACVr}
#' @description This function executes the complete pipeline of \pkg{PACVr} via a single command.
#'
#' @param gbkFile a character vector that specifies the name of, and path to, the GenBank input file;
#' alternatively, a string of the GenBank data
#' @param bamFile a character vector that specifies the name of, and path to, the BAM input file
#' @param windowSize a numeric value that specifies window size in which the coverage is calculated
#' @param logScale a boolean that specifies if the coverage depth is to be log-transformed before visualizing it
#' @param threshold a numeric value that specifies the threshold for plotting coverage depth bars in red as opposed to the default black
#' @param syntenyLineType a numeric value of 1 or 2 that specifies the line type for visualizing IR gene synteny; 1 = ribbon lines, 2 = solid lines, otherwise = no line
#' @param relative a boolean that specifies whether the threshold is a relative value of the average coverage instead of an absolute value
#' @param textSize a numeric value that specifies the relative font size of the text element in the visualization
#' @param verbose the decision to generate additional files with detailed genomic region information
#' @param regionsCheck a boolean that specifies if region analysis of genome should be performed; FALSE disables syntenyLineType and verbose
#' @param output a character vector that specifies the name of, and path to, the output file
#' @return A file in pdf format containing a circular visualization of the submitted plastid sample.
#' @export
#' @examples
#'\dontrun{
#' gbkFile <- system.file("extdata",
#'                        "NC_045072/NC_045072.gb",
#'                        package="PACVr")
#' bamFile <- system.file("extdata",
#'                        "NC_045072/NC_045072_PlastomeReadsOnly.sorted.bam",
#'                        package="PACVr")
#' outFile <- paste(tempdir(), "/NC_045072__all_reads.pdf", sep="")
#' PACVr.complete(gbkFile=gbkFile, bamFile=bamFile, windowSize=250, logScale=FALSE,
#'                threshold=0.5, syntenyLineType=1, relative=TRUE, textSize=0.5,
#'                regionsCheck=FALSE, verbose=FALSE, output=outFile
#'                }
PACVr.complete <- function(gbkFile,
                           bamFile,
                           windowSize=250,
                           logScale=FALSE,
                           threshold=0.5,
                           syntenyLineType=1,
                           relative=TRUE,
                           textSize=0.5,
                           regionsCheck=FALSE,
                           verbose=FALSE,
                           output=NA) {
  ######################################################################
  gbkData <- PACVr.read.gb(gbkFile)
  gbkDataDF <- read.gb2DF(gbkData, regionsCheck)
  if (is.null(gbkDataDF)) {
    logger::log_error(paste("No usable data to perform specified analysis"))
    return(NULL)
  }
  
  ###################################
  if (regionsCheck) {
    logger::log_info('Parsing the different genome regions')
    regions <- PACVr.parseRegions(gbkData,
                                  gbkDataDF)
  } else {
    regions <- PACVr.parseSource(gbkDataDF)
  }

  ###################################
  logger::log_info('Parsing the different genes')
  genes <- PACVr.parseGenes(gbkDataDF)

  ###################################
  logger::log_info('Calculating the sequencing coverage')
  coverage <- PACVr.calcCoverage(bamFile,
                                 windowSize)

  ###################################
  linkData <- NULL
  IRCheck <- regionsCheck && isSyntenyLineType(syntenyLineType)
  if (IRCheck) {
    logger::log_info('Inferring the IR regions and the genes within the IRs')
    linkData <- PACVr.generateIRGeneData(genes,
                                         regions,
                                         syntenyLineType)
  }

  ###################################
  if (regionsCheck && verbose) {
      logger::log_info('Generating statistical information on the sequencing coverage')
      PACVr.verboseInformation(gbkData,
                           bamFile,
                           genes,
                           regions,
                           output)
  }
  
  ###################################
  if (!is.na(output)) {
    logger::log_info('Generating a visualization of the sequencing coverage')
    pdf(output, width=10, height=10)
    PACVr.visualizeWithRCircos(
      gbkData,
      genes,
      regions,
      coverage,
      windowSize,
      threshold,
      logScale,
      relative,
      linkData,
      syntenyLineType,
      textSize
    )
    dev.off()
    logger::log_info('Visualization (including coverage) saved as `{output}`')
  } else {
    logger::log_info('No coverage data inferred; generating empty visualization')
    PACVr.visualizeWithRCircos(
      gbkData,
      genes,
      regions,
      coverage,
      windowSize,
      threshold,
      logScale,
      relative,
      linkData,
      syntenyLineType,
      textSize
    )
    dev.off()
    logger::log_info('Visualization (excluding coverage) saved as `{output}`')
  }
  ######################################################################
  logger::log_success('Done.')
  ######################################################################
}

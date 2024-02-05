#!/usr/bin/env RScript
#contributors=c("Gregory Smith", "Nils Jenke", "Michael Gruenstaeudl")
#email="m_gruenstaeudl@fhsu.edu"
#version="2024.02.01.1736"

PACVr.read.gb <- function(gbkFile) {
  gbkRaw <- getGbkRaw(gbkFile)
  if (is.null(gbkRaw$char)) {
    return(NULL)
  }
  gbkData <- read.gbWithHandling(gbkRaw)
  return(gbkData)
}

PACVr.parseName <- function (gbkData) {
  return(read.gbSampleName(gbkData))
}

PACVr.parseQuadripRegions <- function (gbkData, gbkDataDF) {
  raw_quadripRegions <- ParseQuadripartiteStructure(gbkDataDF)
  quadripRegions <- fillDataFrame(gbkData, raw_quadripRegions)
  return(quadripRegions)
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

PACVr.generateIRGeneData <- function(genes, quadripRegions,
                                     syntenyLineType) {
  # Parse GenBank file
  if ("IRb" %in% quadripRegions[, 4] &&
      "IRa" %in% quadripRegions[, 4]) {
    linkData <- GenerateIRSynteny(genes, syntenyLineType)
    return(linkData)
  }
  return(-1)
}

PACVr.verboseInformation <- function(gbkData,
                                     bamFile,
                                     genes,
                                     quadripRegions,
                                     IRCheck,
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
  writeTables(quadripRegions, bamFile, genes, tmpDir, sampleName)
  if (IRCheck) {
    checkIREquality(gbkData, quadripRegions, tmpDir, sampleName)
  }
}

PACVr.visualizeWithRCircos <- function(gbkData,
                                       genes,
                                       quadripRegions,
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
    quadripRegions,
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

PACVr.quadripRegions <- function(gbkData,
                                 gbkDataDF,
                                 isIRCheck) {
  if (isIRCheck) {
    logger::log_info('Parsing the different genome regions')
    quadripRegions <- PACVr.parseQuadripRegions(gbkData,
                                                gbkDataDF)
  } else {
    quadripRegions <- PACVr.parseQuadripRegions(gbkDataDF)
  }
  return(quadripRegions)
}

PACVr.linkData <- function(genes,
                           quadripRegions,
                           IRCheck) {
  linkData <- NULL
  isSyntenyLine <- getIsSyntenyLine(IRCheck)
  if (isSyntenyLine) {
    logger::log_info('Inferring the IR regions and the genes within the IRs')
    linkData <- PACVr.generateIRGeneData(genes,
                                         quadripRegions,
                                         IRCheck)
  }
  return(linkData)
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
#' detailed genomic region information;
#' requires a `IRCheck` value that will perform region analysis
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
  gbkData <- PACVr.read.gb(gbkFile)
  isIRCheck <- getIsIRCheck(IRCheck)
  gbkDataDF <- read.gb2DF(gbkData, isIRCheck)
  if (is.null(gbkDataDF)) {
    logger::log_error(paste("No usable data to perform specified analysis"))
    return(NULL)
  }
  
  ###################################
  quadripRegions <- PACVr.quadripRegions(gbkData,
                                         gbkDataDF,
                                         isIRCheck)

  ###################################
  logger::log_info('Parsing the different genes')
  genes <- PACVr.parseGenes(gbkDataDF)

  ###################################
  logger::log_info('Calculating the sequencing coverage')
  coverage <- PACVr.calcCoverage(bamFile,
                                 windowSize)

  ###################################
  linkData <- PACVr.linkData(genes,
                             quadripRegions,
                             IRCheck)

  ###################################
  if (isIRCheck && verbose) {
      logger::log_info('Generating statistical information on the sequencing coverage')
      PACVr.verboseInformation(gbkData,
                               bamFile,
                               genes,
                               quadripRegions,
                               IRCheck,
                               output)
  } else if (verbose) {
      logger::log_warn(paste0('Verbose output requires `IRCheck` in ',
                              '`', deparse(getIRCheckTypes()), '`'))
  }
  
  ###################################
  if (!is.na(output)) {
    logger::log_info('Generating a visualization of the sequencing coverage')
    pdf(output, width=10, height=10)
    PACVr.visualizeWithRCircos(
      gbkData,
      genes,
      quadripRegions,
      coverage,
      windowSize,
      threshold,
      logScale,
      relative,
      linkData,
      IRCheck,
      textSize
    )
    dev.off()
    logger::log_info('Visualization (including coverage) saved as `{output}`')
  } else {
    logger::log_info('No coverage data inferred; generating empty visualization')
    PACVr.visualizeWithRCircos(
      gbkData,
      genes,
      quadripRegions,
      coverage,
      windowSize,
      threshold,
      logScale,
      relative,
      linkData,
      IRCheck,
      textSize
    )
    dev.off()
    logger::log_info('Visualization (excluding coverage) saved as `{output}`')
  }
  ######################################################################
  logger::log_success('Done.')
  ######################################################################
  return(0)
}

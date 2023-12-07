#!/usr/bin/R
#contributors=c("Gregory Smith", "Nils Jenke", "Michael Gruenstaeudl")
#email="m_gruenstaeudl@fhsu.edu"
#version="2023.11.23.1530"

PACVr.parseName <- function (gbkData) {
  return(read.gbSampleName(gbkData))
}

PACVr.parseRegions <- function (gbkData) {
  raw_regions <- ExtractAllRegions(gbkData)
  regions <- fillDataFrame(gbkData, raw_regions)
  return(regions)
}

PACVr.parseGenes <- function (gbkData) {
  # This function parses the genes of a GenBank file
  genes <- ExtractAllGenes(gbkData)
  return(genes)
}

PACVr.calcCoverage <-
  function (bamFile, regions, windowSize=250) {
    coverage <- CovCalc(bamFile, windowSize)
    return(coverage)
  }

PACVr.generateIRGeneData <- function(gbkData, genes, regions,
                                     syntenyLineType) {
  # Parse GenBank file
  if ("IRb" %in% regions[, 4] &&
      "IRa" %in% regions[, 4] && syntenyLineType < 3) {
    linkData <- GenerateIRSynteny(genes, syntenyLineType)
    return(linkData)
  }
  return(-1)
}

PACVr.verboseInformation <- function(gbkData,
                                     bamFile,
                                     genes,
                                     regions,
                                     output,
                                     sampleName) {
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
#' @param gbkFile a character vector that specifies the name of, and path to, the GenBank input file
#' @param bamFile a character vector that specifies the name of, and path to, the BAM input file
#' @param windowSize a numeric value that specifies window size in which the coverage is calculated
#' @param logScale a boolean that specifies if the coverage depth is to be log-transformed before visualizing it
#' @param threshold a numeric value that specifies the threshold for plotting coverage depth bars in red as opposed to the default black
#' @param syntenyLineType a numeric value of 1, 2 or 3 that specifies the line type for visualizing IR gene synteny; 1 = ribbon lines, 2 = solid lines, 3 = no line
#' @param relative a boolean that specifies whether the threshold is a relative value of the average coverage instead of an absolute value
#' @param textSize a numeric value that specifies the relative font size of the text element in the visualization
#' @param verbose the decision to generate additional files with detailed genomic region information
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
#'                verbose=FALSE, output=outFile
#'                }
PACVr.complete <- function(gbkFile,
                           bamFile,
                           windowSize=250,
                           logScale=FALSE,
                           threshold=0.5,
                           syntenyLineType=1,
                           relative=TRUE,
                           textSize=0.5,
                           verbose=FALSE,
                           output=NA) {
  ######################################################################
  logger::log_info('Reading GenBank flatfile `{gbkFile}`')
  gbkData <- read.gb::read.gb(gbkFile, DNA=TRUE, Type="full", Source="File")
  sampleName <- read.gbSampleName(gbkData)
  
  ###################################
  logger::log_info('Parsing different genome regions and genes')
  regions <- PACVr.parseRegions(gbkData)
  genes <- PACVr.parseGenes(gbkData)

  ###################################
  logger::log_info('Calculating sequencing coverage')
  coverage <- PACVr.calcCoverage(bamFile,
                                 regions,
                                 windowSize)

  ###################################
  logger::log_info('Inferring IR regions and genes within IRs')
  linkData <- PACVr.generateIRGeneData(gbkData,
                                       genes,
                                       regions,
                                       syntenyLineType)
  ###################################
  if (verbose) {
      logger::log_info('Generating statistical information on sequencing coverage')
      PACVr.verboseInformation(gbkData,
                           bamFile,
                           genes,
                           regions,
                           output,
                           sampleName)
  }
  
  ###################################
  if (!is.na(output)) {
    logger::log_info('Generating visualization of sequencing coverage')
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
    logger::log_info('Saved visualization including coverage as `{output}`')
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
    logger::log_info('Saved visualization excluding coverage as `{output}`')
  }
  ######################################################################
  logger::log_success('Done.')
  ######################################################################
}

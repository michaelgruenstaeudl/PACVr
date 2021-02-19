#!/usr/bin/R
#contributors = c("Michael Gruenstaeudl","Nils Jenke")
#email = "m.gruenstaeudl@fu-berlin.de", "nilsj24@zedat.fu-berlin.de"
#version = "2020.07.29.1700"

PACVr.parseName <- function (gbkData) {
  # Parse sample name
  sample_name = c(sample_name=genbankr::accession(gbkData), genome_name=genbankr::seqinfo(gbkData)@genome)
  return(sample_name)
}

PACVr.parseRegions <- function (gbkData) {
  raw_regions <- ExtractAllRegions(gbkData)
  regions <- fillDataFrame(gbkData, raw_regions)
  return(regions)
}

PACVr.parseGenes <- function (gbkData) {
  # Parse GenBank file
  genes <- ExtractAllGenes(gbkData)
  return(genes)
}

PACVr.calcCoverage <- function (bamFile, regions, windowSize = 250) {
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
                                     bam.file,
                                     genes,
                                     regions,
                                     output,
                                     sample_name) {
  if (!is.na(output)) {
    outDir <- dirname(output)
    tmpDir <- file.path(outDir, paste(sample_name["sample_name"], ".tmp", sep = ""))
  } else {
    tmpDir <- file.path(".", paste(sample_name["sample_name"], ".tmp", sep = ""))
  }
  
  
  if (dir.exists(tmpDir) == FALSE) {
    dir.create(tmpDir)
  }
  writeTables(regions, bam.file, genes, tmpDir, sample_name)
  checkIREquality(gbkData, regions, tmpDir, sample_name)
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
  # 1. Generate plot title
  plotTitle <-
    paste(genbankr::sources(gbkData)$organism,
          genbankr::accession(gbkData))
  
  
  # 2. Visualize
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



PACVr.complete <- function(gbk.file,
                           bam.file,
                           windowSize = 250,
                           logScale = FALSE,
                           threshold = 0.5,
                           syntenyLineType = 1,
                           relative = TRUE,
                           textSize = 0.5,
                           verbose = FALSE,
                           output = NA) {
  # 1. Preparatory steps
  gbkData <- genbankr::readGenBank(gbk.file, verbose = FALSE)
  sample_name <- PACVr.parseName(gbkData)
  
  # 2. Conduct operations
  regions <- PACVr.parseRegions(gbkData)
  genes <- PACVr.parseGenes(gbkData)
  coverage <- PACVr.calcCoverage(bam.file, regions, windowSize)
  linkData <- PACVr.generateIRGeneData(gbkData, genes, regions,
                                       syntenyLineType)
  
  # 3. Save files
  if (verbose) {
    PACVr.verboseInformation(gbkData,
                             bam.file,
                             genes,
                             regions,
                             output,
                             sample_name)
  }
  
  # 4. Save plot
  if (!is.na(output)) {
    pdf(output, width = 10, height = 10)
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
  } else {
    # 4. Generate visualization
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
  }
}

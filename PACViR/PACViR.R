#!/usr/bin/R
#contributors = c("Michael Gruenstaeudl", "Nils Jenke")
#email = "m.gruenstaeudl@fu-berlin.de", "nilsj24@zedat.fu-berlin.de"
#version = "2019.03.15.1800"

PACViR.parseRegions <- function (gbkFile) {
    
  # 1. Source custom R functions
    source("parseRegions.R")
    source("helpers.R")
    
  # 2. Parse GenBank file
    gbkData <- genbankr::readGenBank(gbkFile)
    raw_regions <- ExtractAllRegions(gbkData)
    return(raw_regions)
}


PACViR.parseGenes <- function (gbkFile, raw_regions) {
    
  # 1. Source custom R functions
    source("parseGenes.R")
    source("helpers.R")
    
  # 2. Parse GenBank file
    gbkData <- genbankr::readGenBank(gbkFile)
    raw_genes <- ExtractAllGenes(gbkData)
    genes_inclSplitOnes <- SplitGenesAtRegionBorders(raw_genes, raw_regions)
    genes_withRegionInfo <- AssignRegionInfo(genes_inclSplitOnes, raw_regions)
    genes_withUpdRegions <- AdjustRegionLocation(genes_withRegionInfo, raw_regions)
    return(genes_withUpdRegions)
}


PACViR.calcCoverage <- function (bamFile, raw_regions,
                                 windowSize=250,
                                 outDir="./PACViR_output/",
                                 mosdepthCmd="mosdepth") {
    
  # 1. Source custom R functions
    source("calcCoverage.R")
    source("helpers.R")
    
  # 2. Coverage calculation
    raw_coverage <- CovCalc(bamFile, windowSize, outDir, mosdepthCmd)
    cov_withRegionInfo <- AssignRegionInfo(raw_coverage, raw_regions)
    cov_inclSplitOnes <- SplitCovAtRegionBorders(cov_withRegionInfo, raw_regions)
    regions_withUpdRegions <- AdjustRegionLocation(raw_regions, raw_regions)
    cov_withUpdRegions <- adjustCoverage(cov_inclSplitOnes, regions_withUpdRegions)
    return(cov_withUpdRegions)
}


PACViR.generateIRGeneData <- function (genes_withUpdRegions) {
    
  # 1. Source custom R functions
    source("generateIRGeneData.R")
    source("helpers.R")
    
  # 2. Parse GenBank file
    linkData <- GenerateIRGeneData(genes_withUpdRegions)
    return(linkData)
}


PACViR.GenerateHistogramData <- function (cov_withUpdRegions) {
    
  # 1. Source custom R functions
    source("generateHistogramData.R")
    source("helpers.R")
    
  # 2. Parse GenBank file
    lineData <- GenerateHistogramData(cov_withUpdRegions)
    return(lineData)
}


PACViR.visualizeWithRCircos <- function (gbkFile,
                                         genes_withUpdRegions,
                                         regions_withUpdRegions,
                                         cov_withUpdRegions,
                                         threshold=25,
                                         lineData,
                                         linkData) {
    
  # 1. Source custom R functions
    source("visualizeWithRCircos.R")

  # 2. Get gbkData
    gbkData <- genbankr::readGenBank(gbkFile)

  # 3. Calculate average
    avg <- as.integer(unique(lineData[ ,4]))
    
  # 4. Visualize
    visualizeWithRCircos(gbkData, genes_withUpdRegions, regions_withUpdRegions, cov_withUpdRegions, threshold, avg, lineData, linkData)
}

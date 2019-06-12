#!/usr/bin/R
#contributors = c("Michael Gruenstaeudl","Nils Jenke")
#email = "m.gruenstaeudl@fu-berlin.de", "nilsj24@zedat.fu-berlin.de"
#version = "2019.06.12.1530"

PACViR.parseName <- function (gbkFile) {
  
  # Parse sample name
    gbkData <- genbankr::readGenBank(gbkFile)
    sample_name = genbankr::accession(gbkData)
    return(sample_name)
}

PACViR.parseRegions <- function (gbkFile, tmpDir) {
    
  # Parse GenBank file
    gbkData <- genbankr::readGenBank(gbkFile)
    raw_regions <- ExtractAllRegions(gbkData)
    write.table(raw_regions, file = paste(tmpDir, .Platform$file.sep, "regions.PACViR.tmp", sep=""), 
                sep = " ", dec = ".", row.names = TRUE, col.names = TRUE)
    return(raw_regions)
}


PACViR.parseGenes <- function (gbkFile, raw_regions) {
    
  # Parse GenBank file
    gbkData <- genbankr::readGenBank(gbkFile)
    raw_genes <- ExtractAllGenes(gbkData)
    genes_inclSplitOnes <- SplitGenesAtRegionBorders(raw_genes, raw_regions)
    genes_withRegionInfo <- AssignRegionInfo(genes_inclSplitOnes, raw_regions)
    genes_withUpdRegions <- AdjustRegionLocation(genes_withRegionInfo, raw_regions)
    return(genes_withUpdRegions)
}


PACViR.calcCoverage <- function (bamFile, raw_regions,
                                 windowSize=250,
                                 output,
                                 mosdepthCmd="mosdepth") {
    
  # Coverage calculation
    raw_coverage <- CovCalc(bamFile, windowSize, output, mosdepthCmd)
    cov_withRegionInfo <- AssignRegionInfo(raw_coverage, raw_regions)
    cov_inclSplitOnes <- SplitCovAtRegionBorders(cov_withRegionInfo, raw_regions)
    regions_withUpdRegions <- AdjustRegionLocation(raw_regions, raw_regions)
    cov_withUpdRegions <- adjustCoverage(cov_inclSplitOnes, regions_withUpdRegions)
    return(cov_withUpdRegions)
}


PACViR.generateIRGeneData <- function (genes_withUpdRegions) {
    
  # Parse GenBank file
    linkData <- GenerateIRGeneData(genes_withUpdRegions)
    return(linkData)
}


PACViR.GenerateHistogramData <- function (cov_withUpdRegions) {
    
  # Parse GenBank file
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
    
  # 1. Get gbkData
    gbkData <- genbankr::readGenBank(gbkFile)

  # 2. Calculate average
    avg <- as.integer(unique(lineData[ ,4]))
    
  # 3. Visualize
    visualizeWithRCircos(gbkData, genes_withUpdRegions, regions_withUpdRegions, cov_withUpdRegions, threshold, avg, lineData, linkData)
}


PACViR.complete <- function(gbk.file, bam.file, 
                            windowSize = 250, mosdepthCmd = "mosdepth", 
                            threshold = 25, delete = TRUE,
                            output = "./PACViR_output.svg" ) {
  
  # 1. Preparatory steps
  sample_name <- PACViR.parseName(gbk.file)
  outDir <- dirname(output)
  tmpDir <- file.path(outDir, paste(sample_name, ".tmp", sep=""))
  get_os <- Sys.info()[1]
  if(get_os == "Windows"){system(paste("md", tmpDir))} else {system(paste("mkdir -p", tmpDir))}
  
  # 2. Conduct operations
  raw_regions <- PACViR.parseRegions(gbk.file, tmpDir)
  genes_withUpdRegions <- PACViR.parseGenes(gbk.file, raw_regions)
  regions_withUpdRegions <- AdjustRegionLocation(raw_regions, raw_regions)
  cov_withUpdRegions <- PACViR.calcCoverage(bam.file, raw_regions, windowSize, tmpDir, mosdepthCmd)
  linkData <- PACViR.generateIRGeneData(genes_withUpdRegions)
  lineData <- PACViR.GenerateHistogramData(cov_withUpdRegions)
  
  # 3. Save plot
  svg(output)
  PACViR.visualizeWithRCircos(gbk.file, genes_withUpdRegions, regions_withUpdRegions, cov_withUpdRegions, threshold, lineData, linkData)
  dev.off()
  
  # 4. Delete temp files
  if (isTRUE(delete)) {
      system(paste0("rm -r ", tmpDir))
  }
}

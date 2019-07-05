#!/usr/bin/R
#contributors = c("Michael Gruenstaeudl","Nils Jenke")
#email = "m.gruenstaeudl@fu-berlin.de", "nilsj24@zedat.fu-berlin.de"
#version = "2019.07.05.1100"

PACVr.parseName <- function (gbkFile) {
  
  # Parse sample name
    gbkData <- suppressMessages(suppressWarnings(genbankr::readGenBank(gbkFile)))
    sample_name = genbankr::accession(gbkData)
    return(sample_name)
}

PACVr.parseRegions <- function (gbkFile, tmpDir) {
    
  # Parse GenBank file
    gbkData <- suppressMessages(suppressWarnings(genbankr::readGenBank(gbkFile)))
    raw_regions <- ExtractAllRegions(gbkData)
    write.table(raw_regions, file = paste(tmpDir, .Platform$file.sep, "regions.PACVr.tmp", sep=""), 
                sep = " ", dec = ".", row.names = TRUE, col.names = TRUE)
    return(raw_regions)
}


PACVr.parseGenes <- function (gbkFile, raw_regions) {
    
  # Parse GenBank file
    gbkData <- suppressMessages(suppressWarnings(genbankr::readGenBank(gbkFile)))
    raw_genes <- ExtractAllGenes(gbkData)
    genes_inclSplitOnes <- SplitGenesAtRegionBorders(raw_genes, raw_regions)
    genes_withRegionInfo <- AssignRegionInfo(genes_inclSplitOnes, raw_regions)
    genes_withUpdRegions <- AdjustRegionLocation(genes_withRegionInfo, raw_regions)
    return(genes_withUpdRegions)
}


PACVr.calcCoverage <- function (chromName, bamFile, 
                                 raw_regions, windowSize=250,
                                 output,
                                 mosdepthCmd) {
    
  # Coverage calculation
    #mosdepth_present = suppressWarnings(system(paste("command -v", mosdepthCmd), intern=TRUE))
    #if (!(is.null(attr(mosdepth_present, "status")))) {
    mosdepth_present = tryCatch(system2(command="command", args=c("-v", mosdepthCmd), stdout=TRUE), error=function(e) NULL)
    if (is.null(mosdepth_present)) {
      print('The software tool Mosdepth (https://github.com/brentp/mosdepth) was not detected on your system. Please install it, for example via R command: system("conda install -y mosdepth")')
      raw_coverage <- DummyCov(chromName, raw_regions, windowSize)
    } else {
      raw_coverage <- CovCalc(bamFile, windowSize, output, mosdepthCmd)
    }
    cov_withRegionInfo <- AssignRegionInfo(raw_coverage, raw_regions)
    cov_inclSplitOnes <- SplitCovAtRegionBorders(cov_withRegionInfo, raw_regions)
    regions_withUpdRegions <- AdjustRegionLocation(raw_regions, raw_regions)
    cov_withUpdRegions <- adjustCoverage(cov_inclSplitOnes, regions_withUpdRegions)
    return(cov_withUpdRegions)
}


PACVr.generateIRGeneData <- function (genes_withUpdRegions) {
    
  # Parse GenBank file
    linkData <- GenerateIRGeneData(genes_withUpdRegions)
    return(linkData)
}


PACVr.GenerateHistogramData <- function (cov_withUpdRegions) {
    
  # Parse GenBank file
    lineData <- GenerateHistogramData(cov_withUpdRegions)
    return(lineData)
}


PACVr.visualizeWithRCircos <- function (gbkFile,
                                         genes_withUpdRegions,
                                         regions_withUpdRegions,
                                         cov_withUpdRegions,
                                         threshold=25,
                                         lineData,
                                         linkData,
                                         mosdepthCmd) {
    
  # 1. Generate plot title
    gbkData <- suppressMessages(suppressWarnings(genbankr::readGenBank(gbkFile)))
    #mosdepth_present = suppressWarnings(system(paste("command -v", mosdepthCmd), intern=TRUE))
    #if (!(is.null(attr(mosdepth_present, "status")))) {
    mosdepth_present = tryCatch(system2(command="command", args=c("-v", mosdepthCmd), stdout=TRUE), error=function(e) NULL)
    if (is.null(mosdepth_present)) {
      plotTitle <- paste("Dummy data -", genbankr::accession(gbkData), ". Please install mosdepth.")
    } else {
      plotTitle <- genbankr::definition(gbkData)
    }

  # 2. Calculate average
    avg <- as.integer(unique(lineData[ ,4]))
    
  # 3. Visualize
    visualizeWithRCircos(plotTitle, genes_withUpdRegions, regions_withUpdRegions, cov_withUpdRegions, threshold, avg, lineData, linkData)
}


PACVr.complete <- function(gbk.file, bam.file, 
                            windowSize = 250, mosdepthCmd = "mosdepth", 
                            threshold = 25, delete = TRUE,
                            output = "./PACVr_output.pdf" ) {
  
  # 1. Preparatory steps
  sample_name <- PACVr.parseName(gbk.file)
  outDir <- dirname(output)
  tmpDir <- file.path(outDir, paste(sample_name, ".tmp", sep=""))
  get_os <- Sys.info()[1]
  #if(get_os == "Windows"){system(paste("md", tmpDir))} else {system(paste("mkdir -p", tmpDir))}
  if(get_os == "Windows"){system2(command="md", args=tmpDir)} else {system2(command="mkdir", args=c("-p", tmpDir))}
  
  # 2. Conduct operations
  raw_regions <- PACVr.parseRegions(gbk.file, tmpDir)
  genes_withUpdRegions <- PACVr.parseGenes(gbk.file, raw_regions)
  regions_withUpdRegions <- AdjustRegionLocation(raw_regions, raw_regions)
  cov_withUpdRegions <- PACVr.calcCoverage(sample_name, bam.file, raw_regions, windowSize, tmpDir, mosdepthCmd)
  linkData <- PACVr.generateIRGeneData(genes_withUpdRegions)
  lineData <- PACVr.GenerateHistogramData(cov_withUpdRegions)
  
  # 3. Save plot
  pdf(output)
  PACVr.visualizeWithRCircos(gbk.file, genes_withUpdRegions, regions_withUpdRegions, cov_withUpdRegions, threshold, lineData, linkData, mosdepthCmd)
  dev.off()
  
  # 4. Delete temp files
  if (isTRUE(delete)) {
      #system(paste0("rm -r ", tmpDir))
      system2(command="rm", args=c("-r", tmpDir))
  }
}

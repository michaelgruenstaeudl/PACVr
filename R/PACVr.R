#!/usr/bin/R
#contributors = c("Michael Gruenstaeudl","Nils Jenke")
#email = "m.gruenstaeudl@fu-berlin.de", "nilsj24@zedat.fu-berlin.de"
#version = "2020.01.17.1800"


PACVr.parseName <- function (gbkData) {
  # Parse sample name
    sample_name = genbankr::accession(gbkData)
    return(sample_name)
}


PACVr.parseRegions <- function (gbkData, tmpDir) {
    raw_regions <- ExtractAllRegions(gbkData)
    regions <- fillDataFrame(gbkData,raw_regions)
    return(regions)
}


PACVr.parseGenes <- function (gbkData) {
    
  # Parse GenBank file
    genes <- ExtractAllGenes(gbkData)
    return(genes)
}


PACVr.calcCoverage <- function (chromName, bamFile, regions, 
                                windowSize=250, output, logScale=FALSE,
                                mosdepthCmd) {
    
  # Coverage calculation
    mosdepth_present = tryCatch(system2(command="command", args=c("-v", mosdepthCmd), stdout=TRUE), error=function(e) NULL)
    if (is.null(mosdepth_present)) {
      message('The software tool Mosdepth (https://github.com/brentp/mosdepth) was not detected on your system. Please install it, for example via R command: system("conda install -y mosdepth")')
      coverage <- DummyCov(chromName, regions, windowSize)
    } else {
      coverage <- CovCalc(bamFile, windowSize, output, 
                          logScale, mosdepthCmd)
    }
    return(coverage)
}


PACVr.generateIRGeneData <- function(gbkData, genes, regions,
                                     synteny, syntenyLineType) {
  
  # Parse GenBank file
  if("IRb" %in% regions[,4] && "IRa" %in% regions[,4] && synteny == TRUE){
    checkIRSynteny(gbkData, regions)
    linkData <- GenerateIRSynteny(genes, syntenyLineType)
    return(linkData)
  }
  return(-1)
}


PACVr.GenerateHistogramData <- function (regions,coverage) {
    
  # Parse GenBank file
    lineData <- GenerateHistogramData(regions,coverage)
    return(lineData)
}


PACVr.visualizeWithRCircos <- function (gbkData, genes, regions,
                                        coverage, threshold=25, relative,
                                        mosdepthCmd, linkData, lineData, 
                                        syntenyLineType) {
    
  # 1. Generate plot title
    mosdepth_present = tryCatch(system2(command="command", args=c("-v", mosdepthCmd), stdout=TRUE), error=function(e) NULL)
    if (is.null(mosdepth_present)) {
      plotTitle <- paste("Dummy data -", genbankr::sources(gbkData)$organism,genbankr::accession(gbkData), ". Please install mosdepth.")
    } else {
      plotTitle <- paste(genbankr::sources(gbkData)$organism,genbankr::accession(gbkData))
    }

  # 2. Calculate average
    #avg <- as.integer(unique(lineData[ ,4]))
    
  # 3. Visualize
    visualizeWithRCircos(plotTitle, genes, regions, 
                         coverage, threshold, relative, 
                         linkData, lineData, syntenyLineType)
}


PACVr.complete <- function(gbk.file, bam.file, windowSize = 250,
                           mosdepthCmd = "mosdepth", logScale = FALSE, threshold = 0.5,
                           synteny = FALSE, syntenyLineType="1", relative = TRUE, 
                           delete = TRUE, output = NA) {
  
  # 1. Preparatory steps
  #gbkData <- suppressMessages(genbankr::readGenBank(gbk.file,verbose = FALSE))
  gbkData <- genbankr::readGenBank(gbk.file,verbose = FALSE)
  sample_name <- PACVr.parseName(gbkData)
  if (!is.na(output)) {
    outDir <- dirname(output)
    tmpDir <- file.path(outDir, paste(sample_name, ".tmp", sep=""))
  } else {
    tmpDir <- file.path(".", paste(sample_name, ".tmp", sep=""))
  }
  get_os <- Sys.info()[1]
  if(get_os == "Windows"){
      system2(command="md", args=tmpDir)} else {system2(command="mkdir", args=c("-p", tmpDir))
  }
  
  # 2. Conduct operations
  regions <- PACVr.parseRegions(gbkData, tmpDir)
  genes <- PACVr.parseGenes(gbkData)
  coverage <- PACVr.calcCoverage(sample_name, bam.file, regions, 
                                 windowSize, tmpDir, logScale, 
                                 mosdepthCmd)
  linkData <- PACVr.generateIRGeneData(gbkData, genes, regions,
                                       synteny, syntenyLineType)
  lineData <- PACVr.GenerateHistogramData(regions,coverage)
  
  
  # 3. Save plot
  if (!is.na(output)) {
    pdf(output, encoding="ISOLatin9.enc")
    PACVr.visualizeWithRCircos(gbkData, genes, regions, 
                               coverage, threshold, relative, 
                               mosdepthCmd, linkData, lineData, 
                               syntenyLineType)
    dev.off()
  } else {
  # 4. Generate visualization
    PACVr.visualizeWithRCircos(gbkData, genes, regions, 
                               coverage, threshold, relative, 
                               mosdepthCmd, linkData, lineData, 
                               syntenyLineType)
  }
  
  # 5. Delete temp files
  if (isTRUE(delete)) {
    system2(command="rm", args=c("-r", tmpDir))
  }
}

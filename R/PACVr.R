#!/usr/bin/R
#contributors = c("Michael Gruenstaeudl","Nils Jenke")
#email = "m.gruenstaeudl@fu-berlin.de", "nilsj24@zedat.fu-berlin.de"
#version = "2020.07.29.1700"

PACVr.parseName <- function (gbkData) {
  # Parse sample name
    sample_name = genbankr::accession(gbkData)
    return(sample_name)
}

PACVr.parseRegions <- function (gbkData) {
    raw_regions <- ExtractAllRegions(gbkData)
    regions <- fillDataFrame(gbkData,raw_regions)
    return(regions)
}

PACVr.parseGenes <- function (gbkData) {
    
  # Parse GenBank file
    genes <- ExtractAllGenes(gbkData)
    return(genes)
}

PACVr.calcCoverage <- function (bamFile, regions, windowSize=250) {
    
    coverage <- CovCalc(bamFile, windowSize)
    return(coverage)
}


PACVr.generateIRGeneData <- function(gbkData, genes, regions,
                                     syntenyLineType) {
  
  # Parse GenBank file
  if("IRb" %in% regions[,4] && "IRa" %in% regions[,4] && syntenyLineType < 3){
    checkIREquality(gbkData, regions)
    linkData <- GenerateIRSynteny(genes, syntenyLineType)
    return(linkData)
  }
  return(-1)
}


PACVr.visualizeWithRCircos <- function (gbkData, genes, regions,
                                        coverage, windowSize, logScale, 
                                        threshold, relative, linkData, 
                                        syntenyLineType, textSize) {

  # 1. Generate plot title
  plotTitle <- paste(genbankr::sources(gbkData)$organism,genbankr::accession(gbkData))
  
  
  # 2. Visualize
    visualizeWithRCircos(plotTitle, genes, regions, 
                         coverage, windowSize, threshold,
                         logScale, relative, linkData, 
                         syntenyLineType, textSize)
}



PACVr.complete <- function(gbk.file, bam.file, windowSize = 250,
                           logScale = FALSE, threshold = 0.5, syntenyLineType = 1, 
                           relative = TRUE, textSize=0.5, delete = TRUE, 
                           output = NA) {
  
  # 1. Preparatory steps
  gbkData <- genbankr::readGenBank(gbk.file,verbose = FALSE)
  sample_name <- PACVr.parseName(gbkData)
  
  # 2. Conduct operations
  regions <- PACVr.parseRegions(gbkData)
  genes <- PACVr.parseGenes(gbkData)
  coverage <- PACVr.calcCoverage(bam.file, regions, windowSize)
  linkData <- PACVr.generateIRGeneData(gbkData, genes, regions,
                                       syntenyLineType)
  
  # 3. Save files
  if (!is.na(output)) {
    outDir <- dirname(output)
    tmpDir <- file.path(outDir, paste(sample_name, ".tmp", sep=""))
  } else {
    tmpDir <- file.path(".", paste(sample_name, ".tmp", sep=""))
    }
  
  if (dir.exists(tmpDir) == FALSE) {
    dir.create(tmpDir)
  }
  
  if (relative == TRUE) {
    coverage$lowCoverage <- coverage$coverage < mean(coverage$coverage) * threshold
  } else { 
    coverage$lowCoverage <- coverage$coverage < mean(coverage$coverage) 
    }
  coverage$lowCoverage[coverage$lowCoverage == TRUE] <- "*"
  coverage$lowCoverage[coverage$lowCoverage == FALSE] <- ""
  
  write.csv(coverage[,c("chromStart","chromEnd","coverage","lowCoverage")],
            paste(tmpDir, .Platform$file.sep, tools::file_path_sans_ext(basename(bam.file)),"_coverage.regions.bed", sep=""), 
            row.names = FALSE, quote = FALSE)
  write.csv(coverage[coverage$lowCoverage == "*",c("chromStart","chromEnd","coverage","lowCoverage")], 
            paste(tmpDir, .Platform$file.sep, sample_name, "_low_coverage.csv", sep=""), 
            row.names = FALSE, quote = FALSE)

  
  # 4. Save plot
  if (!is.na(output)) {
    pdf(output, width=10, height = 10)
    PACVr.visualizeWithRCircos(gbkData, genes, regions, 
                               coverage, windowSize, threshold,
                               logScale, relative, linkData, 
                               syntenyLineType, textSize)
    dev.off()
  } else {
  # 4. Generate visualization
    PACVr.visualizeWithRCircos(gbkData, genes, regions, 
                               coverage, windowSize, threshold,
                               logScale, relative, linkData,
                               syntenyLineType, textSize)
    dev.off()
  }

  
  # 5. Delete temp files
  if (isTRUE(delete)) {
    unlink(tmpDir,recursive = TRUE)
  }
}

#!/usr/bin/env RScript
#contributors=c("Gregory Smith", "Nils Jenke", "Michael Gruenstaeudl")
#email="m_gruenstaeudl@fhsu.edu"
#version="2024.05.16.1704"

writeCovTables <- function(covData, sample_name, dir) {
  writeStatsTable(covData$ir_genes, sample_name, dir, "coverage.genes")
  writeStatsTable(covData$ir_regions, sample_name, dir, "coverage.regions")
  writeStatsTable(covData$ir_noncoding, sample_name, dir, "coverage.noncoding")
}

writeStatsTable <- function(df, sample_name, dir, fileName) {
  write.table(
    df,
    paste(dir, .Platform$file.sep, sample_name["sample_name"], "_", fileName, ".tsv", sep = ""),
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )
}

writeSumTables <- function(summary, sample_name, dir) {
  sumField2File <- getSumField2File()
  for (sumField in names(summary)) {
    sumFile <- sumField2File[[sumField]]
    writeStatsTable(summary[[sumField]], sample_name, dir, sumFile)
  }
}

getSumField2File <- function() {
  sumField2File <- list(
    "regions_summary" = "summary.regions",
    "genes_summary" = "coverage.summary.genes",
    "noncoding_summary" = "coverage.summary.noncoding"
  )
  return(sumField2File)
}

printCovStats <- function(gbkData,
                          coverageRaw,
                          analysisSpecs,
                          outputSpecs) {
  sampleName <- gbkData$sampleName
  seqnames <- unname(sampleName[sampleName %in% names(coverageRaw)])
  if (length(seqnames) == 0) {
    logger::log_error("Neither `ACCESSION` nor `VERSION` matches BAM sample name")
    return(NULL)
  }

  logger::log_info('Generating statistical information on the sequencing coverage')
  quadripRegions <- gbkData$quadripRegions
  genes <- gbkData$genes
  covData <- getCovData(quadripRegions, genes)
  covData <- filter_IR_genes(quadripRegions, coverageRaw, seqnames, covData, analysisSpecs)
  covData <- filter_IR_noncoding(quadripRegions, coverageRaw, seqnames, covData, analysisSpecs)
  covData <- filter_IR_regions(coverageRaw, seqnames, covData, analysisSpecs)
  covData <- setLowCoverages(covData, analysisSpecs)

  # Writing values to output table
  statsFilePath <- outputSpecs$statsFilePath
  writeCovTables(covData,
                 sampleName,
                 statsFilePath)

  return(covData)
}

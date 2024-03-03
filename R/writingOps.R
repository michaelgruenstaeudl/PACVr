#!/usr/bin/env RScript
#contributors=c("Gregory Smith", "Nils Jenke", "Michael Gruenstaeudl")
#email="m_gruenstaeudl@fhsu.edu"
#version="2024.03.03.0204"

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

writeCovSumTables <- function(covSummaries, sample_name, dir) {
  writeStatsTable(covSummaries$genes_summary, sample_name, dir, "coverage.summary.genes")
  writeStatsTable(covSummaries$regions_summary, sample_name, dir, "coverage.summary.regions")
  writeStatsTable(covSummaries$noncoding_summary, sample_name, dir, "coverage.summary.noncoding")
}

printCovStats <- function(coverageRaw,
                          genes,
                          quadripRegions,
                          sampleName,
                          analysisSpecs,
                          dir) {
  seqnames <- unname(sampleName[sampleName %in% names(coverageRaw)])
  if (length(seqnames) == 0) {
    logger::log_error("Neither `ACCESSION` nor `VERSION` matches BAM sample name")
    return(-1)
  }

  logger::log_info('Generating statistical information on the sequencing coverage')
  covData <- getCovData(quadripRegions, genes, analysisSpecs)
  covData <- filter_IR_genes(quadripRegions, coverageRaw, seqnames, covData)
  covData <- filter_IR_noncoding(quadripRegions, coverageRaw, seqnames, covData)
  covData <- filter_IR_regions(coverageRaw, seqnames, covData)
  covData <- setLowCoverage(covData)

  # Writing values to output table
  writeCovTables(covData, sampleName, dir)

  # Getting and writing summarized coverage data, grouped by `quadripRegions`
  covSummaries <- getCovSummaries(covData,
                                  analysisSpecs)
  writeCovSumTables(covSummaries, sampleName, dir)
}

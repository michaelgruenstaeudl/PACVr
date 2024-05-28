#!/usr/bin/env RScript
#contributors=c("Gregory Smith", "Nils Jenke", "Michael Gruenstaeudl")
#email="m_gruenstaeudl@fhsu.edu"
#version="2024.05.25.0158"

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

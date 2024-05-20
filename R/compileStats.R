#!/usr/bin/env RScript
#contributors=c("Gregory Smith", "Nils Jenke", "Michael Gruenstaeudl")
#email="m_gruenstaeudl@fhsu.edu"
#version="2024.05.16.1704"

PACVr.compileCovStats <- function(gbkData,
                                  coverageRaw,
                                  analysisSpecs,
                                  outputSpecs) {
  outputSpecs$makeStatsFolder()

  covData <- printCovStats(gbkData,
                           coverageRaw,
                           analysisSpecs,
                           outputSpecs)

  # Count of N nucleotides for source
  regionsSummary <- getAmbigCounts(gbkData,
                                   analysisSpecs)

  # Add count of mismatches between IRs
  regions_name <- analysisSpecs$regions_name
  if (analysisSpecs$isSyntenyLine) {
    IR_mismatches <- checkIREquality(gbkData,
                                     analysisSpecs)
    regionsSummary <- dplyr::full_join(regionsSummary,
                                       IR_mismatches,
                                       regions_name)
  }

  # Summarized coverage data, grouped by `quadripRegions`;
  # add previous results to this
  if (!is.null(covData)) {
    summary <- getCovSummaries(covData,
                               analysisSpecs)
    summary$regions_summary <- dplyr::full_join(summary$regions_summary,
                                                regionsSummary,
                                                regions_name)
  } else {
    summary <- list(
      regions_summary = regionsSummary
    )
  }

  writeSumTables(summary,
                 gbkData$sampleName,
                 outputSpecs$statsFilePath)

  logger::log_info('Coverage statistics saved to `{outputSpecs$statsFilePath}`')
}

getAmbigCounts <- function(gbkData, analysisSpecs) {
  ambigCounts <- data.frame(
    N_count = unname(Biostrings::alphabetFrequency(gbkData$sequences)[, "N"])
  )
  ambigCounts[analysisSpecs$regions_name] <- "Complete_genome"
  return(ambigCounts)
}

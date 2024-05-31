#!/usr/bin/env RScript
#contributors=c("Gregory Smith", "Nils Jenke", "Michael Gruenstaeudl")
#email="m_gruenstaeudl@fhsu.edu"
#version="2024.05.31.0123"

PACVr.compileCovStats <- function(gbkData,
                                  coverage) {
  coverage$makeStatsFolder()
  coverage$printCovStats(gbkData)
  coverage$printCovSumStats(gbkData)
  logger::log_info('Coverage statistics saved to `{coverage$outputSpecs$statsFilePath}`')
}

getAmbigCounts <- function(gbkData, analysisSpecs) {
  ambigCounts <- data.frame(
    N_count = unname(Biostrings::alphabetFrequency(gbkData$sequences)[, "N"])
  )
  ambigCounts[analysisSpecs$regions_name] <- "Complete_genome"
  return(ambigCounts)
}

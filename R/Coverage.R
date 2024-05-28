#!/usr/bin/env RScript
#contributors=c("Gregory Smith", "Nils Jenke", "Michael Gruenstaeudl")
#email="m_gruenstaeudl@fhsu.edu"
#version="2024.05.28.2100"

Coverage <- R6Class("Coverage",
  public = list(
    # fields
    analysisSpecs = NULL,
    outputSpecs = NULL,
    coverageRaw = NULL,
    coveragePlot = NULL,
    seqNames = NULL,
    covData = NULL,
    sumData = NULL,

    # constructor
    initialize = function(bamFile,
                          analysisSpecs = NULL,
                          outputSpecs = NULL) {
      private$setAnalysisSpecs(analysisSpecs)
      private$setOutputSpecs(outputSpecs)
      private$setBaseCoverages(bamFile)
    },

    # public methods

    # precondition: `outputSpecs` is set
    makeStatsFolder = function() {
      self$outputSpecs$makeStatsFolder()
    },

    printCovStats = function(gbkData) {
      self$setSeqNames(gbkData)
      if (length(self$seqNames) == 0) {
        logger::log_error("Neither `ACCESSION` nor `VERSION` matches BAM sample name")
        return()
      }

      logger::log_info('Generating statistical information on the sequencing coverage')
      self$setCovData(gbkData)
      self$writeCovTables(gbkData)
    },

    # precondition: `coverageRaw` is set
    setSeqNames = function(gbkData) {
      sampleName <- gbkData$sampleName
      self$seqNames <- unname(sampleName[sampleName %in% names(self$coverageRaw)])
    },

    # precondition: `analysisSpecs`, `coverageRaw`, and `seqNames` are set
    setCovData = function(gbkData) {
      quadripRegions <- gbkData$quadripRegions
      genes <- gbkData$genes
      seqnames <- self$seqNames
      analysisSpecs <- self$analysisSpecs
      coverageRaw <- self$coverageRaw

      covData <- getCovData(quadripRegions, genes)
      covData <- filter_IR_genes(quadripRegions, coverageRaw, seqnames, covData, analysisSpecs)
      covData <- filter_IR_noncoding(quadripRegions, coverageRaw, seqnames, covData, analysisSpecs)
      covData <- filter_IR_regions(coverageRaw, seqnames, covData, analysisSpecs)
      covData <- setLowCoverages(covData, analysisSpecs)
      self$covData <- covData
    },

    # precondition: `covData` and `outputSpecs` are set
    writeCovTables = function(gbkData) {
      writeCovTables(self$covData,
                     gbkData$sampleName,
                     self$outputSpecs$statsFilePath)
    },

    printCovSumStats = function(gbkData) {
      self$setSumData(gbkData)
      self$writeSumTables(gbkData)
    },

    setSumData = function(gbkData) {
      # Count of N nucleotides for source
      self$setAmbigCounts(gbkData)

      # Add count of mismatches between IRs
      self$setIRMismatches(gbkData)

      # Summarized coverage data, grouped by `quadripRegions`;
      # add previous results to this
      self$setCovSummary()
    },

    # precondition: `analysisSpecs` is set
    setAmbigCounts = function(gbkData) {
      self$sumData <- getAmbigCounts(gbkData,
                                     self$analysisSpecs)
    },

    # precondition: `analysisSpecs` and `sumData` are set
    setIRMismatches = function(gbkData) {
      analysisSpecs <- self$analysisSpecs
      if (!analysisSpecs$isSyntenyLine) {
        return()
      }

      IR_mismatches <- checkIREquality(gbkData,
                                       analysisSpecs)
      self$sumData <- dplyr::full_join(self$sumData,
                                       IR_mismatches,
                                       analysisSpecs$regions_name)
    },

    # precondition: `analysisSpecs`, `covData`, and `sumData` are set
    setCovSummary = function() {
      covData <- self$covData
      analysisSpecs <- self$analysisSpecs

      if (!is.null(covData)) {
        summary <- getCovSummaries(covData,
                                   analysisSpecs)
        summary$regions_summary <- dplyr::full_join(summary$regions_summary,
                                                    self$sumData,
                                                    analysisSpecs$regions_name)
      } else {
        summary <- list(
          regions_summary = self$sumData
        )
      }
      self$sumData <- summary
    },

    # precondition: `outputSpecs` and `sumData` are set
    writeSumTables = function(gbkData) {
      writeSumTables(self$sumData,
                     gbkData$sampleName,
                     self$outputSpecs$statsFilePath)
    }
  ),

  # private setters for constructor
  private = list(
    setAnalysisSpecs = function(analysisSpecs) {
      if(is.null(analysisSpecs)) {
        self$analysisSpecs <- AnalysisSpecs$new()
      } else {
        self$analysisSpecs <- analysisSpecs
      }
    },

    setOutputSpecs = function(outputSpecs) {
      if(is.null(outputSpecs)) {
        self$outputSpecs <- OutputSpecs$new()
      } else {
        self$outputSpecs <- outputSpecs
      }
    },

    # precondition: `analysisSpecs` is set
    setBaseCoverages = function(bamFile) {
      coverage <- PACVr.calcCoverage(bamFile,
                                     self$analysisSpecs$windowSize,
                                     self$outputSpecs$logScale)
      self$coverageRaw <- coverage$raw
      self$coveragePlot <- coverage$plot
    }
  )
)

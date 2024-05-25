#!/usr/bin/env RScript
#contributors=c("Gregory Smith", "Nils Jenke", "Michael Gruenstaeudl")
#email="m_gruenstaeudl@fhsu.edu"
#version="2024.05.24.2053"

GBKData <- R6Class("GBKData",
  public = list(
    # fields
    analysisSpecs = NULL,
    genes = NULL,
    sequences = NULL,
    lengths = NULL,
    sampleName = NULL,
    plotTitle = NULL,
    quadripRegions = NULL,
    sourceRegion = NULL,
    linkData = NULL,
    
    # constructor
    initialize = function(gbkFile,
                          analysisSpecs = NULL) {
      read.gbData <- PACVr.read.gb(gbkFile)
      private$setAnalysisSpecs(analysisSpecs)

      # main derivative of `read.gb` data
      gbkSeqFeatures <- read.gbSeqFeaturesAdapt(read.gbData,
                                                self$analysisSpecs)
      if (is.null(gbkSeqFeatures)) {
        logger::log_fatal('Parsing of any sequence features unsuccessful.')
        return(NULL)
      }

      # other derivatives of `read.gb` data
      private$setSequences(read.gbData)
      private$setLengths()
      private$setSampleName(read.gbData)
      private$setPlotTitle(read.gbData)

      # `read.gbData` no longer needed
      rm(read.gbData)
      gc()

      # `gbkSeqFeatures` derivatives
      private$setSourceRegion(gbkSeqFeatures)
      self$setQuadripRegions(gbkSeqFeatures)
      private$setGenes(gbkSeqFeatures)
      private$setIRCheckFields()
      private$setLinkData()

      # `gbkSeqFeatures` no longer needed
      rm(gbkSeqFeatures)
      gc()
    },

    # public setter
    # precondition: `analysisSpecs` is set
    setQuadripRegions = function(gbkSeqFeatures) {
      if (self$analysisSpecs$isIRCheck) {
        logger::log_info('Parsing the quadripartite genome structure')
        self$quadripRegions <- PACVr.parseQuadripRegions(self$lengths,
                                                         gbkSeqFeatures,
                                                         self$analysisSpecs)
      } else {
        self$quadripRegions <- self$sourceRegion
      }
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

    setSequences = function(read.gbData) {
      self$sequences <- read.gbSequence(read.gbData)
    },

    # precondition: `sequences` is set
    setLengths = function() {
      self$lengths <- read.gbLengths(self$sequences)
    },

    setSampleName = function(read.gbData) {
      self$sampleName <- read.gbSampleName(read.gbData)
    },

    setPlotTitle = function(read.gbData) {
      self$plotTitle <- read.gbPlotTitle(read.gbData)
    },

    setSourceRegion = function(gbkSeqFeatures) {
      self$sourceRegion <- PACVr.parseSource(gbkSeqFeatures)
    },

    setGenes = function(gbkSeqFeatures) {
      self$genes <- PACVr.parseGenes(gbkSeqFeatures)
    },

    # precondition: `quadripRegions` is set
    setIRCheckFields = function() {
      IRBandRequirements <- getIRBandRequirements()
      missingBands <- IRBandRequirements[!(IRBandRequirements %in% self$quadripRegions$Band)]
      isSyntenyLine <- ifelse(length(missingBands) == 0, TRUE, FALSE)
      if (self$analysisSpecs$isSyntenyLine && !isSyntenyLine) {
        logger::log_warn(paste0("Unable to test synteny; missing IR band(s): ",
                                "`",
                                paste0(missingBands, collapse = "`, `"),
                                "`"))
        self$analysisSpecs$setIRCheckFields(0)
      }
    },

    # precondition: `genes` and `analysisSpecs` are set
    setLinkData = function() {
      # Parse GenBank file
      if (self$analysisSpecs$isSyntenyLine) {
        self$linkData  <- GenerateIRSynteny(self$genes,
                                            self$analysisSpecs)
      }
    }
  )
)

# "static" "fields" and "methods" used for the class
getIRBandRequirements <- function() {
  IRBandRequirements <- c("IRa", "IRb")
  return(IRBandRequirements)
}

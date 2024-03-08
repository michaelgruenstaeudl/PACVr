#!/usr/bin/env RScript
#contributors=c("Gregory Smith", "Nils Jenke", "Michael Gruenstaeudl")
#email="m_gruenstaeudl@fhsu.edu"
#version="2024.03.03.0509"

AnalysisSpecs <- R6Class("AnalysisSpecs",
  public = list(
    # fields
    syntenyLineType = NULL,
    isIRCheck = FALSE,
    windowSize = 250,
    isSyntenyLine = FALSE,
    regions_name = "Source",
    regions_start = "srcStart",
    regions_end = "srcEnd",
    
    # constructor
    initialize = function(IRCheck = NA,
                          windowSize = 250) {
      self$setIRCheckFields(IRCheck)
      private$setWindowSize(windowSize)
    },

    # public setter
    setIRCheckFields = function(IRCheck) {
      private$setSyntenyLineType(IRCheck)
      private$setIsIRCheck(IRCheck)
      private$setIsSyntenyLine()
      private$setCovFields()
    }
  ),

  # private setters for constructor
  private = list(
    setSyntenyLineType = function(IRCheck) {
      syntenyLineTypes <- getSyntenyLineTypes()
      if (is.numeric(IRCheck) &&
          IRCheck %in% syntenyLineTypes) {
        self$syntenyLineType <- IRCheck
      } else {
        self$syntenyLineType <- NULL
      }
    },

    setIsIRCheck = function(IRCheck) {
      IRCheckTypes <- getIRCheckTypes()
      self$isIRCheck <- (IRCheck %in% IRCheckTypes)
    },

    setWindowSize = function(windowSize) {
      windowSize <- filterPosNumeric(windowSize)
      if (is.null(windowSize)) {
        logger::log_warn("Using default value for `windowSize`: 250")
      } else {
        self$windowSize <- windowSize
      }
    },

    # precondition: `syntenyLineType` is set
    setIsSyntenyLine = function() {
      self$isSyntenyLine <- !is.null(self$syntenyLineType)
    },

    # precondition: `isIRCheck` is set
    setCovFields = function() {
      if (self$isIRCheck) {
        self$regions_name <- "Chromosome"
        self$regions_start <- "chromStart"
        self$regions_end <- "chromEnd"
      }
    }
  )
)

# "static" "fields" and "methods" used for the class
getSyntenyLineTypes <- function() {
  syntenyLineTypes <- c(1, 2)
  return(syntenyLineTypes)
}

getIRCheckTypes <- function() {
  IRCheckTypes <- c(0, 1, 2)
  return(IRCheckTypes)
}

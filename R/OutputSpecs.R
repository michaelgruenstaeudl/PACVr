#!/usr/bin/env RScript
#contributors=c("Gregory Smith", "Nils Jenke", "Michael Gruenstaeudl")
#email="m_gruenstaeudl@fhsu.edu"
#version="2024.03.23.2208"

OutputSpecs <- R6Class("OutputSpecs",
  public = list(
    # fields
    logScale = FALSE,
    threshold = 0.5,
    relative = TRUE,
    textSize = 0.5,
    output = NULL,
    outputType = NULL,
    isOutput = FALSE,
    statsFilePath = NULL,
    
    # constructor
    initialize = function(logScale = FALSE,
                          threshold = 0.5,
                          relative = TRUE,
                          textSize = 0.5,
                          output = NULL,
                          sampleName = NULL) {
      private$setLogScale(logScale)
      private$setThreshold(threshold)
      private$setRelative(relative)
      private$setTextSize(textSize)
      private$setOutputFields(output)
      private$setStatsFilePath(sampleName)
    },

    makeStatsFolder = function() {
      if (!dir.exists(self$statsFilePath)) {
        dir.create(self$statsFilePath)
      }
    }
  ),

  # private setters for constructor
  private = list(
    setLogScale = function(logScale) {
      logScale <- filterLogical(logScale)
      if (is.null(logScale)) {
        logger::log_warn("Using default value for `logScale`: FALSE")
      } else {
        self$logScale <- logScale
      }
    },

    setThreshold = function(threshold) {
      threshold <- filterPosNumeric(threshold)
      if (is.null(threshold)) {
        logger::log_warn("Using default value for `threshold`: 0.5")
      } else {
        self$threshold <- threshold
      }
    },

    setRelative = function(relative) {
      relative <- filterLogical(relative)
      if (is.null(relative)) {
        logger::log_warn("Using default value for `relative`: TRUE")
      } else {
        self$relative <- relative
      }
    },

    setTextSize = function(textSize) {
      textSize <- filterPosNumeric(textSize)
      if (is.null(textSize)) {
        logger::log_warn("Using default value for `textSize`: 0.5")
      } else {
        self$textSize <- textSize
      }
    },

    setOutputFields = function(output) {
      outputTypes <- paste(getOutputTypes(), collapse = "|")
      outputPattern <- sprintf("^(?:.+\\.)(%s)$", outputTypes)
      outputMatch <- regexec(outputPattern, output, ignore.case = TRUE)
      outputVec <- regmatches(output, outputMatch)
    
      # non-char `output` or non-match for char `output`
      if ((length(outputVec) == 0) || (length(outputVec[[1]]) == 0)) {
        logger::log_info("No `output` file detected")
      } else {
        self$output <- outputVec[[1]][1]
        self$outputType <- tolower(outputVec[[1]][2])
        self$isOutput <- TRUE
      }
    },

    # precondition: `isOutput` and `output` are set
    setStatsFilePath = function(sampleName) {
      # Step 1. Check ...
      if (self$isOutput) {
        outDir <- file.path(dirname(self$output),
                            paste(sampleName["sample_name"],
                                  ".tmp",
                                  sep=""))
      } else {
        outDir <-
          file.path(".", paste(sampleName["sample_name"],
                               ".tmp",
                               sep=""))
      }
      self$statsFilePath <- outDir
    }
  )
)

# "static" "fields" and "methods" used for the class
getOutputTypes <- function() {
  outputTypes <- c("pdf",
                   "png")
  return(outputTypes)
}

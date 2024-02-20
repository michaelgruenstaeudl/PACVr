#!/usr/bin/env RScript
#contributors=c("Gregory Smith", "Nils Jenke", "Michael Gruenstaeudl")
#email="m_gruenstaeudl@fhsu.edu"
#version="2024.02.01.1736"

HistCol <- function(cov, threshold, relative, logScale) {
  # Function to generate color vector for histogram data
  # ARGS:
  #       cov:       data.frame of coverage
  #       threshold: numeric value of a specific threshold
  # RETURNS:
  #   color vector
  # Error handling
  if (!is.numeric(threshold) | threshold < 0) {
    warning("User-defined coverage depth threshold must be >=1.")
    stop()
  }
  if (relative == TRUE & logScale) {
    threshold <- mean(cov[, 4]) + log(threshold)
  } else if (relative == TRUE) {
    threshold <- mean(cov[, 4]) * threshold
  }
  color <- rep("black", nrow(cov))
  ind <- as.numeric(cov[, 4]) <= threshold
  color <- replace(color, ind, "red")
  return(color)
}
    
boolToDeci <- function(boolList) {
  out = 0
  boolList <- rev(boolList)
  for (i in 1:length(boolList)) {
    out = out + boolList[i] * (2 ^ (i - 1))
  }
  return(out)
}

validateColors <- function(colorsToValidate) {
  colorNames <- colors()
  unsupportedColors <- colorsToValidate[!(colorsToValidate %in% colorNames)]
  if (length(unsupportedColors) > 0) {
    stop("Unsupported R plot color defined.")
  }
}

getAnalysisSpecs <- function(IRCheck,
                             windowSize) {
  analysisSpecs <- list(
    syntenyLineType = getSyntenyLineType(IRCheck),
    isIRCheck = getIsIRCheck(IRCheck),
    windowSize = filterPosNumeric(windowSize)
  )
  analysisSpecs$isSyntenyLine <- !is.null(analysisSpecs$syntenyLineType)
  return(analysisSpecs)
}

getPlotSpecs <- function(logScale,
                         threshold,
                         relative,
                         textSize,
                         output) {
  plotSpecs <- list(
    logScale = filterLogical(logScale),
    threshold = filterPosNumeric(threshold),
    relative = filterLogical(relative),
    textSize = filterPosNumeric(textSize),
    output = getOutput(output)
  )
  plotSpecs$isOutput <- !is.null(plotSpecs$output)
  return(plotSpecs)
}

filterByType <- function(x, typeFun) {
  if (typeFun(x)) {
    return(x)
  } else {
    return(NULL)
  }
}

filterPosNumeric <- function(x) {
  return (
    filterByType(x, is.pos.numeric)
  )
}

filterLogical <- function(x) {
  return (
    filterByType(x, is.logical)
  )
}

is.pos.numeric <- function(x) {
  return (
    is.numeric(x) && (x > 0)
  )
}

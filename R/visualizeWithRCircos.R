#!/usr/bin/env RScript
#contributors=c("Gregory Smith", "Nils Jenke", "Michael Gruenstaeudl")
#email="m_gruenstaeudl@fhsu.edu"
#version="2024.02.01.1736"

visualizeWithRCircos <- function(gbkData,
                                 coverage,
                                 analysisSpecs,
                                 plotSpecs) {
  regions <- gbkData$quadripRegions
  genes <- gbkData$genes

  # STEP 1. RCIRCOS INITIALIZATION
  RCircosInit(regions)
  
  # STEP 2. SET PARAMETER FOR IDEOGRAM
  setPlotParams(genes,
                regions,
                coverage,
                plotSpecs)

  # STEP 3. GRAPHIC DEVICE INITIALIZATION
  RCircos::RCircos.Set.Plot.Area()
  RCircos::RCircos.Chromosome.Ideogram.Plot()

  # STEP 4. GENERATE PLOT
  logger::log_info('  Generating RCircos plot')
  positions <- plotMain(genes,
                        coverage)
  
  # STEP 5. OPTIONAL PLOTS
  averageLines <- NULL
  if (analysisSpecs$isIRCheck) {
    plotRegionNames(regions)
    averageLines <- plotAverageLines(regions,
                                     coverage,
                                     analysisSpecs$windowSize,
                                     positions)
  }

  if (analysisSpecs$isSyntenyLine) {
    plotIRLinks(gbkData$linkData,
                analysisSpecs$syntenyLineType)
  }

  # STEP 6. GENERATE TITLE AND LEGEND
  logger::log_info('  Generating title and legend for visualization')
  graphics::title(paste(gbkData$plotTitle), line = -4.5, cex.main = 0.8)
  addLegend(coverage,
            averageLines,
            plotSpecs)
  
}

RCircosInit <- function(regions) {
  suppressMessages(
    RCircos::RCircos.Set.Core.Components(
      cyto.info      =  regions,
      chr.exclude    =  NULL,
      tracks.inside  =  0,
      tracks.outside =  0
    )
  )
}

setPlotParams <- function(genes,
                          regions,
                          coverage,
                          plotSpecs) {
  textSize <- plotSpecs$textSize

  RCircosEnvironment.params <- RCircos::RCircos.Get.Plot.Parameters()
  RCircosEnvironment.params$base.per.unit <- 1
  RCircosEnvironment.params$chrom.paddings <- 1
  RCircosEnvironment.params$track.height <- 0.07
  RCircosEnvironment.params$text.size <- textSize
  RCircosEnvironment.params$track.background <- "gray71"
  RCircosEnvironment.params$sub.tracks <- 4
  RCircosEnvironment.params$char.width <-
    6000000 * (max(regions$chromEnd) / (52669 + 310 * (nrow(genes)))) / textSize

  RCircosEnvironment.params$hist.colors <- HistCol(coverage,
                                                   plotSpecs$threshold,
                                                   plotSpecs$relative,
                                                   plotSpecs$logScale)
  RCircosEnvironment.params$line.color <- "yellow3"
  RCircosEnvironment.params$chrom.width <- 0.05
  RCircosEnvironment.params$track.in.start <- 1.08
  RCircosEnvironment.params$track.out.start <- 1.5
  RCircosEnvironment.params$radius.len <- 3

  PACVr.Reset.Plot.Parameters(RCircosEnvironment.params)
  RCircosEnvironment.cyto <- RCircos::RCircos.Get.Plot.Ideogram()
  RCircosEnvironment.cyto$ChrColor <- "black"
  RCircos::RCircos.Reset.Plot.Ideogram(RCircosEnvironment.cyto)
}

plotMain <- function(genes, coverage) {
  PACVr.Ideogram.Tick.Plot(
    tick.num = 10,
    track.for.ticks = 2,
    add.text.size = 0.1
  )
  
  PACVr.Gene.Connector.Plot(
    genomic.data = genes,
    track.num = 1,
    side = "in"
    #inside.pos = inside.pos, outside.pos = outside.pos)
  )
  
  PACVr.Gene.Name.Plot(
    gene.data = genes,
    name.col = 4,
    track.num = 2,
    side = "in"
  )
  
  positions <- list(
    outside.pos = RCircos::RCircos.Get.Plot.Boundary(track.num = 5, "in")[1],
    inside.pos = RCircos::RCircos.Get.Plot.Boundary(track.num = 6, "in")[2]
  )
  
  PACVr.Histogram.Plot(
    hist.data = coverage,
    data.col = 4,
    track.num = 5,
    side = "in",
    outside.pos = positions$outside.pos,
    inside.pos = positions$inside.pos
  )
  return(positions)
}

addLegend <- function(coverage,
                      averageLines,
                      plotSpecs) {
  legendParams <- getLegendParams(coverage,
                                  averageLines,
                                  plotSpecs)
  graphics::legend(
    x = legendParams$x,
    y = legendParams$y,
    legend = legendParams$legend,
    pch = c(15, 15, NA, legendParams$partSpace),
    lty = c(NA, NA, legendParams$avgLine, legendParams$partSpace),
    lwd = 2,
    col = c("black", "red", legendParams$avgColor, legendParams$partSpace),
    cex = 0.5,
    bty = "n"
  )
}

getLegendParams <- function(coverage,
                            averageLines,
                            plotSpecs) {
  meanCoverage <- mean(coverage[, 4])
  threshold <- plotSpecs$threshold
  if (plotSpecs$relative == TRUE) {
    absolute <- trunc(meanCoverage * threshold)
    perc <- threshold * 100
    
    legendParams <- list(
      x = -1.6,
      y = -1.2,
      vals = c(absolute, 
               threshold * 100, 
               absolute,
               perc)
    )
  } else {
    absolute <- round(threshold / trunc(meanCoverage) * 100)
    
    legendParams <- list(
      x = "bottomleft",
      y = NULL,
      vals = c(threshold, 
               absolute, 
               threshold,
               absolute)
    )
  }
  
  legendParams$legend <- c(
    paste(
      "Coverage > ",
      legendParams$vals[[1]],
      "X ",
      "(=",
      legendParams$vals[[2]],
      "% of avg. cov.)",
      sep = ""
    ),
    as.expression(bquote(
      "Coverage" <= .(
      paste(legendParams$vals[[3]], "X (=", legendParams$vals[[4]], "% of avg. cov.)", sep = "")
    )))
  )
  
  if (is.vector(averageLines)) {
    legendParams$legend <- c(legendParams$legend,
                             "Average Coverage:",
                             averageLines)
    legendParams$avgLine <- 1
    legendParams$avgColor <- "yellow3"
    legendParams$partSpace <- rep(NA, length(averageLines))
  }
  
  return(legendParams)
}

getOutputFields <- function(output) {
  outputTypes <- paste(getOutputTypes(), collapse = "|")
  outputPattern <- sprintf("^(?:.+\\.)(%s)$", outputTypes)
  outputMatch <- regexec(outputPattern, output, ignore.case = TRUE)
  outputVec <- regmatches(output, outputMatch)

  # non-char `output` or non-match for char `output`
  if ((length(outputVec) == 0) || (length(outputVec[[1]]) == 0)) {
    output <- NULL
    outputType <- NULL
    isOutput <- FALSE
  } else {
    output <- outputVec[[1]][1]
    outputType <- tolower(outputVec[[1]][2])
    isOutput <- TRUE
  }

  outputFields <- list(
    output = output,
    outputType = outputType,
    isOutput = isOutput
  )
  return(outputFields)
}

createVizFile <- function(plotSpecs) {
  output <- plotSpecs$output
  outputType <- plotSpecs$outputType

  if (outputType == "pdf") {
    pdf(output,
        width=10,
        height=10)
  } else if (outputType == "png") {
    png(output,
        width=480,
        height=480)
  }
}

getOutputTypes <- function() {
  outputTypes <- c("pdf",
                   "png")
  return(outputTypes)
}

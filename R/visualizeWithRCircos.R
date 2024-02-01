#!/usr/bin/env RScript
#contributors=c("Gregory Smith", "Nils Jenke", "Michael Gruenstaeudl")
#email="m_gruenstaeudl@fhsu.edu"
#version="2024.02.01.0322"

visualizeWithRCircos <- function(plotTitle,
                                 genes,
                                 regions,
                                 coverage,
                                 windowSize,
                                 logScale,
                                 threshold,
                                 relative,
                                 linkData,
                                 syntenyLineType = NA,
                                 textSize) {

  if (logScale == TRUE) {
    coverage$coverage <- log(cov$coverage)
    #coverage$coverage <- log(coverage$coverage)
  }
  coverage$Chromosome <- ""
  
  # STEP 1. RCIRCOS INITIALIZATION
  RCircosInit(regions)
  
  # STEP 2. SET PARAMETER FOR IDEOGRAM
  setPlotParams(genes, regions, coverage, logScale, 
                threshold, relative, textSize)

  # STEP 3. GRAPHIC DEVICE INITIALIZATION
  RCircos::RCircos.Set.Plot.Area()
  RCircos::RCircos.Chromosome.Ideogram.Plot()

  # STEP 4. GENERATE PLOT
  logger::log_info('  Generating RCircos plot')
  positions <- plotMain(genes, coverage)
  
  # STEP 5. OPTIONAL PLOTS
  averageLines <- NULL
  if (isRealRegions(regions)) {
    plotRegionNames(regions)
    averageLines <- plotAverageLines(regions, coverage, windowSize, positions)
  }

  if (is.data.frame(linkData)) {
    plotIRLinks(linkData, syntenyLineType)
  }
  
  # STEP 6. GENERATE TITLE AND LEGEND
  logger::log_info('  Generating title and legend for visualization')
  graphics::title(paste(plotTitle), line = -4.5, cex.main = 0.8)
  addLegend(relative, coverage, threshold, averageLines)
  
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
                          logScale,
                          threshold,
                          relative,
                          textSize) {
  RCircosEnvironment.params <- RCircos::RCircos.Get.Plot.Parameters()
  RCircosEnvironment.params$base.per.unit <- 1
  RCircosEnvironment.params$chrom.paddings <- 1
  RCircosEnvironment.params$track.height <- 0.07
  RCircosEnvironment.params$text.size <- textSize
  RCircosEnvironment.params$track.background <- "gray71"
  RCircosEnvironment.params$sub.tracks <- 4
  RCircosEnvironment.params$char.width <-
    6000000 * (max(regions$chromEnd) / (52669 + 310 * (nrow(genes)))) / textSize

  RCircosEnvironment.params$hist.colors <- HistCol(coverage, threshold, relative, logScale)
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

addLegend <- function(relative, coverage, threshold, averageLines) {
  legendParams <- getLegendParams(relative, coverage, threshold, averageLines)
  
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

getLegendParams <- function(relative, coverage, threshold, averageLines) {
  meanCoverage <- mean(coverage[, 4])
  if (relative == TRUE) {
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

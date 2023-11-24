#!/usr/bin/R
#contributors=c("Gregory Smith", "Nils Jenke", "Michael Gruenstaeudl")
#email="m_gruenstaeudl@fhsu.edu"
#version="2023.11.23.1530"

 
#' Title
#'
#' @param plotTitle
#' @param genes
#' @param regions
#' @param coverage
#' @param threshold
#'
#' @return
#' @export
#'
#' @examples
visualizeWithRCircos <- function(plotTitle,
                                 genes,
                                 regions,
                                 coverage,
                                 windowSize,
                                 logScale,
                                 threshold,
                                 relative,
                                 linkData,
                                 syntenyLineType = 3,
                                 textSize) {

  # Generates the visualization of genome data and their tracks
  # ARGS:
  #   plotTitle: character string
  #   genes_withUpdRegions: data frame that contains the genomic region, gene start, gene end and gene names
  #   regions_withUpdRegions: data frame that contains the genomic region, region start, region end and two dummy columns
  #   cov_withUpdRegions: data frame that contains the genomic region, coverage start, coverage end and coverage value
  #   threshold: numeric value that indicate how many bases are covered at a given threshold
  #   avg:  numeric value of the average coverage value
  #   lineData: data frame that contains IRb region information, gene start (IRb), gene end (IRb), gene names, IRa region
  #             information, gene start (IRa), gene end (IRb)and gene names.
  #   linkData: data frame that contains genomic region, coverage start, coverage end and coverage value
  # RETURNS:
  #   ---

  if (logScale == TRUE) {
    coverage$coverage <- log(cov$coverage)
    #coverage$coverage <- log(coverage$coverage)
  }
  coverage$Chromosome <- ""
  
  # STEP 1. RCIRCOS INITIALIZATION
  
  #library(RCircos)  # NOT ALLOWED BY CRAN
  #RCircosEnvironment <- get("RCircos.Env", envir = globalenv())  # NOT ALLOWED BY CRAN
  RCircosEnvironment <- RCircos::RCircos.Env
 
  suppressMessages(
    RCircos::RCircos.Set.Core.Components(
      cyto.info      = regions,
      chr.exclude    =  NULL,
      tracks.inside  =  0,
      tracks.outside =  0
    )
  )
  
  # STEP 2. SET PARAMETER FOR IDEOGRAM
  RCircosEnvironment.params <- RCircos::RCircos.Get.Plot.Parameters()
  RCircosEnvironment.params$base.per.unit <- 1
  RCircosEnvironment.params$chrom.paddings <- 1
  RCircosEnvironment.params$track.height <- 0.07
  RCircosEnvironment.params$text.size <- textSize
  RCircosEnvironment.params$track.background <- "gray71"
  RCircosEnvironment.params$sub.tracks <- 4
  RCircosEnvironment.params$char.width <-
    6000000 * (max(regions$chromEnd) / (52669 + 310 * (nrow(genes)))) / textSize

  # TO DO - Please check why the below lines produce the warning message:
  suppressMessages({
    suppressWarnings({
  RCircosEnvironment.params$hist.color <- HistCol(coverage, threshold, relative, logScale)
  RCircosEnvironment.params$line.color <- "yellow3"
  RCircosEnvironment.params$chrom.width <- 0.05
  RCircosEnvironment.params$track.in.start <- 1.08
  RCircosEnvironment.params$track.out.start <- 1.5
  RCircosEnvironment.params$radius.len <- 3
  PACVr.Reset.Plot.Parameters(RCircosEnvironment.params)
  RCircosEnvironment.cyto <- RCircos::RCircos.Get.Plot.Ideogram()  
  # The above lines causes message:
  #Warning message:
  #  In !parameters$text.color %in% colorNames || !parameters$hist.color %in%  :
  #  'length(x) = 644 > 1' in coercion to 'logical(1)'
    })
  })
  
  RCircosEnvironment.cyto$ChrColor <- "black"
  RCircos::RCircos.Reset.Plot.Ideogram(RCircosEnvironment.cyto)
  
  # STEP 3. GRAPHIC DEVICE INITIALIZATION
  suppressMessages(RCircos::RCircos.Set.Plot.Area())
  suppressMessages(RCircos::RCircos.Chromosome.Ideogram.Plot())
  
  # STEP 4. GENERATE PLOT
  logger::log_info('  Generating RCircos plot')
  PACVr.Ideogram.Tick.Plot(
    tick.num = 10,
    track.for.ticks = 2,
    add.text.size = 0.1
  )
  
  suppressMessages(PACVr.Gene.Connector.Plot(
    genomic.data = genes,
    track.num = 1,
    side = "in"
    #inside.pos = inside.pos, outside.pos = outside.pos)
  ))
  
  #  suppressMessages(
  PACVr.Gene.Name.Plot(
    gene.data = genes,
    name.col = 4,
    track.num = 2,
    side = "in"
  )
  #  )
  
  PACVr.Gene.Name.Plot(
    gene.data = regions,
    name.col = 4,
    track.num = 1,
    side = "out",
    rotate = 90,
    correction = 0.2,
    add.text.size = 0.2
  )
  
  outside.pos <- RCircos::RCircos.Get.Plot.Boundary(track.num = 5, "in")[1]
  inside.pos <- RCircos::RCircos.Get.Plot.Boundary(track.num = 6, "in")[2]
  
  suppressMessages(
    PACVr.Histogram.Plot(
      hist.data = coverage,
      data.col = 4,
      track.num = 5,
      side = "in",
      outside.pos = outside.pos,
      inside.pos = inside.pos
    )
  )
  
  averageLines <- c()
  for (i in 1:nrow(regions)) {
    lineData <-
      GenerateHistogramData(regions[i,], coverage, windowSize, (i == nrow(regions)))
    averageLines <-
      c(averageLines, paste(regions[i, 4], ": ", trunc(lineData[1, 4]), "X", sep = ""))
    suppressMessages(
      PACVr.Line.Plot(
        line.data = lineData,
        data.col = 4,
        track.num = 5,
        side = "in",
        min.value = min(coverage[4]),
        max.value = max(coverage[4]),
        outside.pos = outside.pos,
        inside.pos = inside.pos,
        genomic.columns = 3,
        is.sorted = TRUE
      )
    )
  }
  
  if (is.data.frame(linkData) == TRUE) {
    if (syntenyLineType == 1) {
      suppressMessages(
        RCircos::RCircos.Ribbon.Plot(
          ribbon.data = linkData,
          track.num = 7,
          by.chromosome = FALSE,
          genomic.columns = 3,
          twist = TRUE
        )
      )
    }
    
    else if (syntenyLineType == 2) {
      suppressMessages(
        RCircos::RCircos.Link.Plot(
          link.data = linkData,
          track.num = 7,
          by.chromosome = FALSE,
          genomic.columns = 3,
          lineWidth = rep(0.5, nrow(linkData))
        )
      )
    }
  }
  
  # STEP 5. GENERATE TITLE AND LEGEND
  logger::log_info('  Generating title and legend for visualization')
  graphics::title(paste(plotTitle), line = -4.5, cex.main = 0.8)
  if (relative == TRUE) {
    absolute <- trunc(mean(coverage[, 4]) * threshold)
    perc <- threshold * 100
    graphics::legend(
      x = -1.6,
      y = -1.2,
      legend = c(
        paste(
          "Coverage > ",
          trunc(mean(coverage[, 4]) * threshold),
          "X ",
          "(=",
          threshold * 100,
          "% of avg. cov.)",
          sep = ""
        ),
        as.expression(bquote("Coverage" <= .(
          paste(" ", absolute, "X (=", perc, "% of avg. cov.)", sep = "")
        ))),
        "Average Coverage:",
        averageLines
      ),
      pch = c(15, 15, NA, rep(NA, length(averageLines))),
      lty = c(NA, NA, 1, rep(NA, length(averageLines))),
      lwd = 2,
      col = c("black", "red", "yellow3", rep(NA, length(averageLines))),
      cex = 0.5,
      bty = "n"
    )
  } else {
    absolute <- round(threshold / trunc(mean(coverage[, 4])) * 100)
    graphics::legend(
      "bottomleft",
      legend = c(
        paste(
          "Coverage > ",
          threshold,
          "X ",
          "(=",
          round(threshold / trunc(mean(coverage[, 4])) * 100),
          "% of avg. cov.)",
          sep = ""
        ),
        as.expression(bquote("Coverage" <= .(
          paste(threshold, "X (=", absolute, "% of avg. cov.)", sep = "")
        ))),
        "Average Coverage:",
        averageLines
      ),
      pch = c(15, 15, NA, rep(NA, length(averageLines))),
      lty = c(NA, NA, 1, rep(NA, length(averageLines))),
      lwd = 2,
      col = c("black", "red", "yellow3", rep(NA, length(averageLines))),
      cex = 0.5,
      bty = "n"
    )
  }
}

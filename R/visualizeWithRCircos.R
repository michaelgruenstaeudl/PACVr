#!/usr/bin/R
#contributors = c("Michael Gruenstaeudl","Nils Jenke")
#email = "m.gruenstaeudl@fu-berlin.de", "nilsj24@zedat.fu-berlin.de"
#version = "2020.01.17.1800"

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
visualizeWithRCircos <- function(plotTitle, genes, regions, 
                                 coverage, threshold=25,relative,
                                 linkData, syntenyLineType=3) {
  
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

  # 1. RCIRCOS INITIALIZATION

  # See for explanation: https://stackoverflow.com/questions/56875962/r-package-transferring-environment-from-imported-package/56894153#56894153
  RCircos.Env <- RCircos::RCircos.Env

  suppressMessages(
    RCircos.Set.Core.Components(cyto.info      = regions, 
                                chr.exclude    =  NULL,
                                tracks.inside  =  9, 
                                tracks.outside =  0)
  )

  # 2. SET PARAMETER FOR IDEOGRAM

  rcircos.params <- RCircos.Get.Plot.Parameters()
  rcircos.params$base.per.unit <- 1
  rcircos.params$chrom.paddings <- 1
  rcircos.params$track.height <- 0.07#0.07
  rcircos.params$text.size <- 0.1
  #rcircos.params$plot.radius <- 2
  #rcircos.params$radius.len <- 2
  rcircos.params$track.background <- "gray71"
  rcircos.params$sub.tracks <- 4
  rcircos.params$char.width <- 16666667*(max(regions$chromEnd)/50000)
  rcircos.params$hist.color <- HistCol(coverage, threshold,relative)
  rcircos.params$line.color <- "yellow3"
  RCircos.Reset.Plot.Parameters(rcircos.params)
  
  rcircos.cyto <- RCircos.Get.Plot.Ideogram()
  rcircos.cyto$ChrColor <- "black"
  RCircos.Reset.Plot.Ideogram(rcircos.cyto)
  
  
  # 3. GRAPHIC DEVICE INITIALIZATION
  suppressMessages(
    RCircos.Set.Plot.Area()
  )
  suppressMessages(
    RCircos.Chromosome.Ideogram.Plot()
  )
  
  # 4. GENERATE PLOT
  PACVr.Ideogram.Tick.Plot(tick.num=10, track.for.ticks=2, text.size=0.4)
  
  suppressMessages(
    PACVr.Gene.Connector.Plot(genomic.data=genes, track.num=1, side="in")
  )
  
  suppressMessages(
    PACVr.Gene.Name.Plot(gene.data=genes, name.col=4, track.num=2,
                         side="in", text.size=0.3
                        )
  )
  
  PACVr.Gene.Name.Plot(gene.data=regions, name.col=4, track.num=1, 
                       side="out", rotate=90, correction=0.2,
                       text.size=0.5)
  
  outside.pos <- RCircos.Get.Plot.Boundary(track.num = 6, "in")[1]
  inside.pos <- RCircos.Get.Plot.Boundary(track.num = 9, "in")[2]
  
  suppressMessages(
    PACVr.Histogram.Plot(hist.data=coverage, data.col= 4, track.num= 9,
                         side= "in", outside.pos=outside.pos,inside.pos=inside.pos
                        )
 )

  averageLines <- c()
  for (i in 1:nrow(regions)){
    lineData <- GenerateHistogramData(regions[i,],coverage, windowSize, (i == nrow(regions)))
    averageLines <- c(averageLines, paste(regions[i,4],": ",trunc(lineData[1,4]),"X", sep = ""))
    suppressMessages(
      PACVr.Line.Plot(line.data = lineData, data.col = 4, track.num = 9,
                      side = "in", min.value = min(coverage[4]), max.value = max(coverage[4]),
                      outside.pos=outside.pos,inside.pos=inside.pos,genomic.columns = 3, 
                      is.sorted = TRUE
      )
    )
  }
  
  if(is.data.frame(linkData) == TRUE){
    if(syntenyLineType == 1){
      suppressMessages(
        RCircos.Ribbon.Plot(ribbon.data=linkData, track.num=10, by.chromosome=FALSE,
                            genomic.columns=3, twist=TRUE
                           )
      )
    }
  
    else if(syntenyLineType == 2){
      suppressMessages(
        RCircos.Link.Plot(link.data=linkData, track.num=10, by.chromosome=FALSE,
                          genomic.columns=3, lineWidth=rep(0.5, nrow(linkData))
                         )
      )
    }
  }
  
  # 5. GENERATE TITLE AND LEGEND
  title(paste(plotTitle),line = -1.75, cex.main = 0.8)
  if (relative == TRUE) {
    absolute <- trunc(mean(coverage[,4]) * threshold)
    perc <- threshold*100
    legend("bottomleft", legend = c(paste("Coverage > ", trunc(mean(coverage[,4]) * threshold),"X ", "(=",threshold*100,"% of avg. cov.)", sep = ""), 
                                    as.expression(bquote("Coverage"<=.(paste(" ",absolute,"X (=",perc,"% of avg. cov.)",sep="")))),
                                    "Average Coverage:",
                                    averageLines),
           pch = c(15, 15, NA, rep(NA,length(averageLines))), lty = c(NA, NA, 1, rep(NA,length(averageLines))), lwd = 2,
           col = c("black", "red", "yellow3",rep(NA,length(averageLines))), cex = 0.5, bty = "n"
    )
  } else {
    absolute <- round(threshold/trunc(mean(coverage[,4]))*100)
    legend("bottomleft", legend = c(paste("Coverage > ", threshold,"X ", "(=",round(threshold/trunc(mean(coverage[,4]))*100),"% of avg. cov.)", sep = ""), 
                                    as.expression(bquote("Coverage"<=.(paste(threshold, "X (=",absolute,"% of avg. cov.)",sep="")))),
                                    "Average Coverage:",
                                    averageLines),
           pch = c(15, 15, NA, rep(NA,length(averageLines))), lty = c(NA, NA, 1, rep(NA,length(averageLines))), lwd = 2,
           col = c("black", "red", "yellow3",rep(NA,length(averageLines))), cex = 0.5, bty = "n"
    )
  }
}
#!/usr/bin/R
#contributors = c("Michael Gruenstaeudl","Nils Jenke")
#email = "m.gruenstaeudl@fu-berlin.de", "nilsj24@zedat.fu-berlin.de"
#version = "2018.10.23.1700"

#source("helpers.R")
source("/home/michael_science/git/michaelgruenstaeudl_PACViR/PACViR/helpers.R")

visualizeWithRCircos <- function(gbkData, genes_withUpdRegions, regions_withUpdRegions, cov_withUpdRegions, threshold=25, avg, lineData, linkData) {
  # foo bar baz
  # ARGS:
  #   gbkData: foo bar baz
  #   genes_withUpdRegions: foo bar baz
  #   regions_withUpdRegions: foo bar baz
  #   cov_withUpdRegions: foo bar baz
  #   threshold: foo bar baz
  #   avg: foo bar baz
  #   lineData: foo bar baz
  #   linkData: foo bar baz
  # RETURNS:
  #   ---

  # 1. RCIRCOS INITIALIZATION
  RCircos::RCircos.Set.Core.Components(cyto.info      =  regions_withUpdRegions, 
                                       chr.exclude    =  NULL,
                                       tracks.inside  =  8, 
                                       tracks.outside =  0)

  # 2. SET PARAMETER FOR IDEOGRAM
  params <- RCircos::RCircos.Get.Plot.Parameters()
  params$base.per.unit <- 225
  params$track.background <- NULL
  params$sub.tracks <- 1
  params$hist.width <- 1
  params$hist.color <- HistCol(cov_withUpdRegions, threshold)
  params$line.color <- "yellow3"
  RCircos::RCircos.Reset.Plot.Parameters(params)
  
  # 3. GRAPHIC DEVICE INITIALIZATION
  RCircos::RCircos.Set.Plot.Area()
  RCircos::RCircos.Chromosome.Ideogram.Plot()
  
  # 4. GENERATE TITLE AND LEGEND
  title(paste(genbankr::definition(gbkData)),line = -1)
  legend(x= 1.5, y= 2,
         legend = c(paste("Coverage >", threshold), 
                    paste("Coverage <", threshold),
                    paste("Average Coverage:"),
                    paste("-LSC:", avg[1]),
                    paste("-IRb:", avg[2]),
                    paste("-SSC:", avg[3]),
                    paste("-IRa:", avg[4])),
         pch    = c(15, 15, NA, NA, NA, NA, NA),
         lty    = c(NA, NA, 1, NA, NA, NA, NA),
         lwd    = 2,
         col    = c("black", "red", "yellow3"),
         cex    = 0.6 
        )
  
  # 5. GENERATE PLOT
  RCircos::RCircos.Gene.Connector.Plot(genomic.data = genes_withUpdRegions,
                                       track.num    = 1,
                                       side         = "in"
                                      )
  
  RCircos::RCircos.Gene.Name.Plot(gene.data  = genes_withUpdRegions,
                                  name.col   = 4,
                                  track.num  = 2,
                                  side       = "in"
                                  )
  
  RCircos::RCircos.Histogram.Plot(hist.data     = cov_withUpdRegions,
                                  data.col      = 4,
                                  track.num     = 5, 
                                  side          = "in",
                                  outside.pos   = RCircos::RCircos.Get.Plot.Boundary(track.num = 5, "in")[1],
                                  inside.pos    = RCircos::RCircos.Get.Plot.Boundary(track.num = 5, "in")[1]-0.3
                                 )
  
  RCircos::RCircos.Line.Plot(line.data       = lineData,
                             data.col        = 4,
                             track.num       = 5,
                             side            = "in",
                             min.value       = min(cov_withUpdRegions[4]),
                             max.value       = max(cov_withUpdRegions[4]),
                             outside.pos     = RCircos::RCircos.Get.Plot.Boundary(track.num = 5, "in")[1],
                             inside.pos      = RCircos::RCircos.Get.Plot.Boundary(track.num = 5, "in")[1]-0.3,
                             genomic.columns = 3,
                             is.sorted       = TRUE
                            )

  RCircos::RCircos.Link.Plot(link.data     =  linkData, 
                             track.num     = 8,
                             by.chromosome = TRUE
                            )

}

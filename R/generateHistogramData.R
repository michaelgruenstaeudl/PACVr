#!/usr/bin/R
#contributors = c("Michael Gruenstaeudl","Nils Jenke")
#email = "m.gruenstaeudl@fu-berlin.de", "nilsj24@zedat.fu-berlin.de"
#version = "2020.01.17.1800"

GenerateHistogramData <- function(regions,coverage) {
  
  # Function to generate line data for RCircos.Line.Plot
  # ARGS:
  #   coverage: data.frame of coverage
  # RETURNS:
  #   data.frame with region means to plot over histogram data
  # Error handling
  avgLine <- data.frame(chromosome="",start=min(regions[,2]),stop=max(regions[,3]),seg.mean=mean(coverage[,4]))

  return(avgLine)
}

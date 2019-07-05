#!/usr/bin/R
#contributors = c("Michael Gruenstaeudl","Nils Jenke")
#email = "m.gruenstaeudl@fu-berlin.de", "nilsj24@zedat.fu-berlin.de"
#version = "2019.07.05.1100"

GenerateHistogramData <- function(lineData) {
  # Function to generate line data for RCircos.Line.Plot
  # ARGS:
  #   lineData: data.frame of coverage
  # RETURNS:
  #   data.frame with region means to plot over histogram data
  # Error handling
  lsc <- lineData[lineData[,1] == "LSC",4]
  irb <- lineData[lineData[,1] == "IRb",4]
  ssc <- lineData[lineData[,1] == "SSC",4]
  ira <- lineData[lineData[,1] == "IRa",4]
  lineData[lineData[,1] == "LSC",4] <- mean(lsc)
  lineData[lineData[,1] == "IRb",4] <- mean(irb)
  lineData[lineData[,1] == "SSC",4] <- mean(ssc)
  lineData[lineData[,1] == "IRa",4] <- mean(ira)
  return(lineData)
}

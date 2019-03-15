#!/usr/bin/R
#contributors = c("Michael Gruenstaeudl", "Nils Jenke")
#email = "m.gruenstaeudl@fu-berlin.de", "nilsj24@zedat.fu-berlin.de"
#version = "2019.03.15.1800"

source("helpers.R")
#source("/home/michael_science/git/michaelgruenstaeudl_PACViR/PACViR/helpers.R")

ExtractAllRegions <- function(gbkData) {
  # Function to extract specific regions from Genbank input data
  # ARGS:
  #   gbkData: Genbank input data parsed by genbankr
  # RETURNS:
  #   regions in data frame format
  region <- BiocGenerics::as.data.frame(genbankr::otherFeatures(gbkData))
  region <- subset(region, grepl("IRb|IRa|repeat|inverted", region[,"note"], ignore.case = FALSE))
  if (nrow(region) > 2) {
    region <- region[order(region[,4],decreasing = TRUE),]
    region <- region[1:2,]
  }
  if (nrow(region) != 2) {
    stop("Something wrong with IRb/IRa\n")
  }
  region <- region[,1:3]
  region <- region[order(region[,3],decreasing = FALSE),]
  region[1] <- c("IRb","IRa")
  LSC = c(1,region[1,2]-1)
  SSC = c(region[1,3]+1,region[2,2]-1)
  region[3, ] <- c("LSC", LSC)
  region[4, ] <- c("SSC", SSC)
  region[ ,2] <- as.integer(region[ ,2])
  region[ ,3] <- as.integer(region[ ,3])
  region <- region[order(region[ ,3], decreasing = FALSE), ]
  if (nrow(region) != 4) {
    stop("Could not find every region; Genbank file must at least contain notes on IRb/IRa.\n")
  }
  region <- cbind(region, Band  = c("","","",""), Stain = c("","","",""))
  region <- Rename_Df(region, c("Band","Stain"))
  return(region)
}

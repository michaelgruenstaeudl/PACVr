#!/usr/bin/R
#contributors = c("Michael Gruenstaeudl","Nils Jenke")
#email = "m.gruenstaeudl@fu-berlin.de", "nilsj24@zedat.fu-berlin.de"
#version = "2019.06.11.1530"

source("helpers.R")
#source("/home/michael_science/git/michaelgruenstaeudl_PACViR/PACViR/helpers.R")

filter <- function(allRegions_L, where) {
  out = subset(allRegions_L, grepl("IRb|IRa|repeat|inverted", allRegions_L[,where], ignore.case = FALSE))
  if (nrow(out) < 1) {warning(paste("Inverted repeat info not found under qualifier ", where))}
  return(out)
}

ExtractAllRegions <- function(gbkData) {
  # Function to extract specific regions from Genbank input data
  # ARGS:
  #   gbkData: Genbank input data parsed by genbankr
  # RETURNS:
  #   regions in data frame format
  allRegions_L <- BiocGenerics::as.data.frame(genbankr::otherFeatures(gbkData))
  region_L <- tryCatch(filter(allRegions_L,"note"), warning = function(w) filter(allRegions_L,"standard_name"))
  if (nrow(region_L) > 2) {
    region_L <- region_L[order(region_L[,4],decreasing = TRUE),]
    region_L <- region_L[1:2,]
  }
  if (nrow(region_L) != 2) {
    stop("ERROR: The number of detected inverted repeat regions is not two\n")
  }
  region_L <- region_L[,1:3]
  region_L <- region_L[order(region_L[,3],decreasing = FALSE),]
  region_L[1] <- c("IRb","IRa")
  LSC <- c(1,region_L[1,2]-1)
  SSC <- c(region_L[1,3]+1,region_L[2,2]-1)
  region_L[3, ] <- c("LSC", LSC)
  region_L[4, ] <- c("SSC", SSC)
  region_L[ ,2] <- as.integer(region_L[ ,2])
  region_L[ ,3] <- as.integer(region_L[ ,3])
  region_L <- region_L[order(region_L[ ,3], decreasing = FALSE), ]
  if (nrow(region_L) != 4) {
    stop("ERROR: Could not find all necessary regions; Genbank file must at least contain note-qualifiers on IRb and IRa.\n")
  }
  region_L <- cbind(region_L, Band  = c("","","",""), Stain = c("","","",""))
  region_L <- Rename_Df(region_L, c("Band","Stain"))
  return(region_L)
}

#!/usr/bin/R
#contributors = c("Michael Gruenstaeudl","Nils Jenke")
#email = "m.gruenstaeudl@fu-berlin.de", "nilsj24@zedat.fu-berlin.de"
#version = "2018.10.23.1700"

library("optparse")

CmdLineArgs <- function() {
  # Function to generate command line arguments
  # ARGS:---
  # RETURNS:
  #   user inputs 
  option_list <- list(
              # -g/i do not work due to R version bug
                  make_option(opt_str = c("-k","--gbkFile"), 
                              type    = "character", 
                              default = NULL, 
                              dest    = "gbkFile",
                              help    = "GenBank filename", 
                              metavar = "character"),
                  make_option(opt_str = c("-b","--bamFile"), 
                              type    = "character", 
                              default = NULL, 
                              dest    = "bamFile",
                              help    = "BAM filename", 
                              metavar = "character"),
                  make_option(opt_str = c("-r","--windowSize"), 
                              type    = "numeric", 
                              default = 250, 
                              dest    = "windowSize",
                              help    = "window size to calculate the coverage [default = %default]", 
                              metavar = "integer"),
                  make_option(opt_str = c("-m","--mosdepthCmd"), 
                              type    = "character", 
                              default = "mosdepth", 
                              dest    = "mosdepthCmd",
                              help    = "command to execute mosdepth on system [default = %default]", 
                              metavar = "character"),
                  make_option(opt_str = c("-t","--threshold"), 
                              type    = "numeric", 
                              default = 25, 
                              dest    = "threshold",
                              help    = "threshold for plotting coverage at different color [default = %default]", 
                              metavar = "integer"),
                  make_option(opt_str = c("-o","--outDir"), 
                              type    = "character", 
                              default = "./PACViR_output/", 
                              dest    = "outDir",
                              help    = "name of output directory [default= %default]", 
                              metavar = "character"))
  opt_parse <- optparse::OptionParser(option_list=option_list)
  opt <- optparse::parse_args(opt_parse)
  if (is.null(opt$gbkFile)){
    print_help(opt_parse)
    stop("No .gb file supplied", call.=FALSE)
  }
  if (is.null(opt$bamFile)) {
    print_help(opt_parse)
    stop("No .bam file supplied", call.=FALSE)
  }
  return(opt)
}

opt <- CmdLineArgs()

########################################################################

#source("PACViR.R")
#source("helpers.R")
source("/home/michael_science/git/michaelgruenstaeudl_PACViR/PACViR/PACViR.R")
source("/home/michael_science/git/michaelgruenstaeudl_PACViR/PACViR/helpers.R")

raw_regions <- PACViR.parseRegions(opt$gbkFile)

genes_withUpdRegions <- PACViR.parseGenes(opt$gbkFile, raw_regions)

regions_withUpdRegions <- AdjustRegionLocation(raw_regions, raw_regions)

cov_withUpdRegions <- PACViR.calcCoverage(opt$bamFile, raw_regions, opt$windowSize, opt$outDir, opt$mosdepthCmd)

linkData <- PACViR.generateIRGeneData(genes_withUpdRegions)

lineData <- PACViR.GenerateHistogramData(cov_withUpdRegions)

svg(paste(opt$outDir, "test.svg", sep=""))
PACViR.visualizeWithRCircos(opt$gbkFile, genes_withUpdRegions, regions_withUpdRegions, cov_withUpdRegions, opt$threshold, lineData, linkData)
dev.off()

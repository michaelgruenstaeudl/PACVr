#!/usr/bin/env RScript
#contributors=c("Gregory Smith", "Nils Jenke", "Michael Gruenstaeudl")
#email="m_gruenstaeudl@fhsu.edu"
#version="2024.02.05.2100"

library("optparse")

########################################################################

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
                              help    = "a text string that specifies the name of, and path to, the GenBank input file", 
                              metavar = "character"),
                  make_option(opt_str = c("-b","--bamFile"), 
                              type    = "character", 
                              default = NULL, 
                              dest    = "bamFile",
                              help    = "a text string that specifies the name of, and path to, the BAM input file", 
                              metavar = "character"),
                  make_option(opt_str = c("-w","--windowSize"), 
                              type    = "numeric", 
                              default = 250, 
                              dest    = "windowSize",
                              help    = "a numeric value that specifies window size in which the coverage is calculated [default = %default]", 
                              metavar = "integer"),
                  make_option(opt_str = c("-l","--logScale"), 
                              type    = "logical", 
                              default = FALSE, 
                              dest    = "logScale",
                              help    = "a boolean that specifies if the coverage depth is to be log-transformed before visualizing it [default = %default]", 
                              metavar = "logical"),
                  make_option(opt_str = c("-t","--threshold"), 
                              type    = "numeric", 
                              default = 0.5, 
                              dest    = "threshold",
                              help    = "a numeric value that specifies the threshold for plotting coverage depth bars in red as opposed to the default black [default = %default]", 
                              metavar = "integer"),
                  make_option(opt_str = c("-i","--IRCheck"),
                              type    = "numeric", 
                              default = 1,
                              dest    = "IRCheck",
                              help    = paste("a numeric value that specifies if tests for IRs of input genome",
                                              "should be performed, and - if IRs are present - which line type",
                                              "to be used  for visualizing gene synteny between IRs;",
                                              "0 = IR presence test but no synteny visualization,",
                                              "1 = IR presence test and synteny visualization, with ribbon lines between IRs,",
                                              "2 = IR presence test and synteny visualization, with solid lines between IRs,",
                                              "otherwise = neither IR presence test nor synteny visualization",
                                              "[default = %default]"),
                              metavar = "integer"),
                  make_option(opt_str = c("-r","--relative"), 
                              type    = "logical", 
                              default = TRUE, 
                              dest    = "relative",
                              help    = "a boolean that specifies whether the threshold is a relative value of the average coverage instead of an absolute value [default = %default]", 
                              metavar = "logical"),
                  make_option(opt_str = c("-x","--textSize"), 
                              type    = "numeric", 
                              default = 0.5, 
                              dest    = "textSize",
                              help    = "a numeric value that specifies the relative font size of the text element in the visualization [default = %default]", 
                              metavar = "integer"),
                  make_option(opt_str = c("-c","--tabularCovStats"), 
                              type    = "logical", 
                              default = FALSE, 
                              dest    = "tabularCovStats",
                              help    = paste("a boolean, that when TRUE, generates additional files with",
                                              "detailed genomic region information"),
                              metavar = "logical"),
                  make_option(opt_str = c("-o","--output"), 
                              type    = "character", 
                              default = "./PACVr_output.pdf", 
                              dest    = "output",
                              help    = "a text string that specifies the name of, and path to, the output file [default= %default]", 
                              metavar = "character"))
  
  opt_parse <- optparse::OptionParser(option_list=option_list)
  opt <- optparse::parse_args(opt_parse)
  if (is.null(opt$gbkFile)){
    print_help(opt_parse)
    warning("No .gb file supplied", call.=FALSE)
    stop("Error")
  }
  if (is.null(opt$bamFile)) {
    print_help(opt_parse)
    warning("No .bam file supplied", call.=FALSE)
    stop("Error")
  }
  return(opt)
}

opt <- CmdLineArgs()

########################################################################

require("PACVr")
PACVr.complete(gbkFile = opt$gbkFile,
               bamFile = opt$bamFile,
               windowSize = opt$windowSize,
               logScale = opt$logScale,
               threshold = opt$threshold,
               IRCheck = opt$IRCheck,
               relative = opt$relative,
               textSize = opt$textSize,
               tabularCovStats = opt$tabularCovStats,
               output = opt$output)

########################################################################


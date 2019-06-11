#!/usr/bin/R
#contributors = c("Michael Gruenstaeudl","Nils Jenke")
#email = "m.gruenstaeudl@fu-berlin.de", "nilsj24@zedat.fu-berlin.de"
#version = "2019.06.05.1730"

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
                  make_option(opt_str = c("-d","--delete"), 
                              type    = "logical", 
                              default = TRUE, 
                              dest    = "delete",
                              help    = "decision to delete temporary files upon program execution [default = %default]", 
                              metavar = "logical"),
                  make_option(opt_str = c("-o","--output"), 
                              type    = "character", 
                              default = "./PACViR_output.svg", 
                              dest    = "output",
                              help    = "name of output file [default= %default]", 
                              metavar = "character"))
  
  opt_parse <- optparse::OptionParser(option_list=option_list)
  opt <- optparse::parse_args(opt_parse)
  if (is.null(opt$gbkFile)){
    print_help(opt_parse)
    stop("ERROR: No .gb file supplied", call.=FALSE)
  }
  if (is.null(opt$bamFile)) {
    print_help(opt_parse)
    stop("ERROR: No .bam file supplied", call.=FALSE)
  }
  return(opt)
}

opt <- CmdLineArgs()

########################################################################

source("PACViR.R")
PACViR.complete(gbk.file = opt$gbkFile, bam.file = opt$bamFile, 
                windowSize = opt$windowSize, mosdepthCmd = opt$mosdepthCmd, 
                threshold = opt$threshold, delete = opt$delete, 
                output = opt$output)


###########
# TESTING #
###########

# MH800324
source("PACViR.R")
PACViR.complete(
gbk.file="/home/michael_science/research/01_manuscripts/01_First_Author/01_PACViR/05_ANALYSES/02_improve_visualizations/MH800324.gb", 
bam.file="/home/michael_science/research/01_manuscripts/01_First_Author/01_PACViR/05_ANALYSES/02_improve_visualizations/MH800324_PlastomeReadsOnly.sorted.bam", 
windowSize=250, 
threshold=15, 
delete=FALSE,
output="/home/michael_science/research/01_manuscripts/01_First_Author/01_PACViR/05_ANALYSES/02_improve_visualizations/MH800324_AssemblyCoverage_viz.svg"
)

# MH899017
source("PACViR.R")
PACViR.complete(
gbk.file="/home/michael_science/research/01_manuscripts/01_First_Author/01_PACViR/05_ANALYSES/02_improve_visualizations/MH899017.gb", 
bam.file="/home/michael_science/research/01_manuscripts/01_First_Author/01_PACViR/05_ANALYSES/02_improve_visualizations/MH899017_PlastomeReadsOnly.sorted.bam", 
windowSize=250, 
threshold=15, 
delete=FALSE,
output="/home/michael_science/research/01_manuscripts/01_First_Author/01_PACViR/05_ANALYSES/02_improve_visualizations/MH899017_AssemblyCoverage_viz.svg"
)

# MH161174
source("PACViR.R")
PACViR.complete(
gbk.file="/home/michael_science/research/01_manuscripts/01_First_Author/01_PACViR/05_ANALYSES/02_improve_visualizations/MH161174.gb", 
bam.file="/home/michael_science/research/01_manuscripts/01_First_Author/01_PACViR/05_ANALYSES/02_improve_visualizations/MH161174_PlastomeReadsOnly.sorted.bam", 
windowSize=250, 
threshold=15, 
delete=FALSE,
output="/home/michael_science/research/01_manuscripts/01_First_Author/01_PACViR/05_ANALYSES/02_improve_visualizations/MH161174_AssemblyCoverage_viz.svg"
)

# MH161175
source("PACViR.R")
PACViR.complete(
gbk.file="/home/michael_science/research/01_manuscripts/01_First_Author/01_PACViR/05_ANALYSES/02_improve_visualizations/MH161175.gb", 
bam.file="/home/michael_science/research/01_manuscripts/01_First_Author/01_PACViR/05_ANALYSES/02_improve_visualizations/MH161175_PlastomeReadsOnly.sorted.bam", 
windowSize=250, 
threshold=15, 
delete=FALSE,
output="/home/michael_science/research/01_manuscripts/01_First_Author/01_PACViR/05_ANALYSES/02_improve_visualizations/MH161175_AssemblyCoverage_viz.svg"
)

q("no")

###########
###########

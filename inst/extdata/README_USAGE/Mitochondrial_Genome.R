# In R:
library(PACVr)
gbkFile <- system.file("extdata", "MG936619/MG936619.gb", package="PACVr")
bamFile <- system.file("extdata", "MG936619/MG936619_subsampled.bam",
                       package="PACVr")

outFile <- paste(tempdir(), "/MG936619_CoverageViz.pdf", sep="")
#outFile <- "../MG936619_CoverageViz.pdf"  # on R-Studio for Windows
#outFile <- "~/MG936619_CoverageViz.pdf"   # on R-Studio for Linux

## ONLY COVERAGE VALUES, NO REGION INDICATORS ##
exitStatus <- PACVr.complete(gbkFile, bamFile, windowSize=50, logScale=FALSE,
                             threshold=0.5, relative=TRUE, textSize=0.5,
                             output=outFile)

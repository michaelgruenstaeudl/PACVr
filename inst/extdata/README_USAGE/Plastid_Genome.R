# In R:
library(PACVr)
gbkFile <- system.file("extdata", "NC_045072/NC_045072.gb", package="PACVr")
bamFile <- system.file("extdata", "NC_045072/NC_045072_subsampled.bam",
                       package="PACVr")

outFile <- paste(tempdir(), "/NC_045072_CoverageViz.pdf", sep="")
#outFile <- "../NC_045072_CoverageViz.pdf"  # on R-Studio for Windows
#outFile <- "~/NC_045072_CoverageViz.pdf"   # on R-Studio for Linux

## ONLY COVERAGE VALUES, NO REGION INDICATORS ##
PACVr.complete(gbkFile, bamFile, windowSize=250, logScale=FALSE,
               threshold=0.5, relative=TRUE, textSize=0.5,
               output=outFile)

## COVERAGE VALUES PLUS REGION INDICATORS ##
PACVr.complete(gbkFile, bamFile, windowSize=250, logScale=FALSE,
               threshold=0.5, relative=TRUE, textSize=0.5,
               regionsCheck=0, output=outFile)

## COVERAGE VALUES PLUS REGION INDICATORS PLUS IR SYNTENY LINES ##
PACVr.complete(gbkFile, bamFile, windowSize=250, logScale=FALSE,
               threshold=0.5, relative=TRUE, textSize=0.5,
               regionsCheck=1, output=outFile)

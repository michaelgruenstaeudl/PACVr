# In R:
library(PACVr)
gbkFile <- system.file("extdata", "NC_045072/NC_045072.gb", package="PACVr")
bamFile <- system.file("extdata", "NC_045072/NC_045072_subsampled.bam",
                       package="PACVr")

outFile <- paste(tempdir(), "/NC_045072_CoverageViz.pdf", sep="")
#outFile <- "../NC_045072_CoverageViz.pdf"  # on R-Studio for Windows
#outFile <- "~/NC_045072_CoverageViz.pdf"   # on R-Studio for Linux
exitStatusVec <- c()  # useful for large automated jobs

## ONLY COVERAGE VALUES, NO REGION INDICATORS ##
exitStatus <- PACVr.complete(gbkFile, bamFile, windowSize=250, logScale=FALSE,
                             threshold=0.5, relative=TRUE, textSize=0.5,
                             output=outFile)
exitStatusVec <- c(exitStatusVec, exitStatus)

## COVERAGE VALUES PLUS REGION INDICATORS ##
exitStatus <- PACVr.complete(gbkFile, bamFile, windowSize=250, logScale=FALSE,
                             threshold=0.5, relative=TRUE, textSize=0.5,
                             regionsCheck=0, output=outFile)
exitStatusVec <- c(exitStatusVec, exitStatus)

## COVERAGE VALUES PLUS REGION INDICATORS PLUS IR SYNTENY LINES ##
exitStatus <- PACVr.complete(gbkFile, bamFile, windowSize=250, logScale=FALSE,
                             threshold=0.5, relative=TRUE, textSize=0.5,
                             regionsCheck=1, output=outFile)
exitStatusVec <- c(exitStatusVec, exitStatus)

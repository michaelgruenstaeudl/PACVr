\name{PACViR-package}
\alias{PACViR-package}
\alias{PACViR}
\docType{package}

\title{
  Plastome Assembly Coverage Visualization in R
}

\description{
  This function executes the complete \pkg{PACViR} pipeline.
}

\usage{
PACViR.complete(gbk.file, bam.file, 
                windowSize = 250, mosdepthCmd = "mosdepth", 
                threshold = 25, outDir = "./PACViR_output/")
}
\arguments{
\item{gbk.file}{a character vector that specifies the path to a .gbk file}

\item{bam.file}{a character vector that specifies the path to a .bam file}

\item{windowSize}{numeric value that is used to calculate the average coverage in a certain window size}

\item{mosdepthCmd}{a character vector that includes the path to execute mosdepth}

\item{threshold}{numeric value that indicate how many bases are covered at a given threshold}

\item{outDir}{a character vector that specifies the output path}

}

\details{ This is the main function of \pkg{PACViR} to execute the full pipeline. It takes a genbank file and a BAM file that are used to visualize plastome data in a circular way using the \pkg{RCircos} package. Those two input arguments are mandatory. The visualization contains gene names, gene locations and coverage data displayed as histrogram. The pipeline will produce a pdf file format containing the circular visualization of the plastome data. Additionally to this it will give the mosdepth output. 

}

\value{ invisible NULL.
}

\references{

}

\note{
This pipeline only works for plastome data.
}


\author{
  Michael Gruenstaeudl, Nils Jenke

  Maintainer: Michael Gruenstaeudl \email{mi.gruenstaeudl@gmail.com}
}



\examples{
#path to .gbk file
gbk.file <- 'gbk/DAS01.gb'

#path to .bam file
bam.file <- 'bam/DAS01.bam'

windowSize <-100 

mosdepthCmd <- 'mosdepth'

threshold <- 55

outDir <- 'output/'

PACViR.complete(gbk.file, bam.file, windowSize, mosdepthCmd, threshold, outDir)

}

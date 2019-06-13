\name{PACViR.complete}
\alias{PACViR.complete}

\title{
Execute the complete PACViR pipeline via a single command
}


\description{
This function executes the complete \pkg{PACViR} pipeline.
}

\usage{
PACViR.complete(gbk.file,
                bam.file,
                windowSize=250,
                mosdepthCmd='mosdepth',
                threshold=25,
                delete=TRUE,
                output='./PACViR_output.svg')
}


\arguments{
\item{gbk.file}{a character vector that specifies the name of, and path to, the GenBank input file}
\item{bam.file}{a character vector that specifies the name of, and path to, the BAM input file}
\item{windowSize}{a numeric value that specifies window size in which the coverage is calculated}
\item{mosdepthCmd}{a character vector that specifies the command to execute mosdepth on system}
\item{threshold}{a numeric value that specifies the threshold for plotting coverage at different color}
\item{delete}{the decision to delete temporary files upon program execution}
\item{output}{a character vector that specifies the name of, and path to, the output file}
}

\details{
This is the main function of \pkg{PACViR} to execute the full pipeline. It takes a Genbank file and a BAM file that are used to visualize plastome data in a circular way using the \pkg{RCircos} package. Those two input arguments are mandatory. The visualization contains gene names, gene locations and coverage data displayed as histrogram. The pipeline will produce a svg file format containing the circular visualization of the plastome data. Additionally, it will give the mosdepth output.
}

\value{ invisible NULL.
}

\references{
Gruenstaeudl M., Jenke N. (2019) PACViR: Plastome Assembly Coverage Visualization in R. R package version 0.4.
}


\author{
Michael Gruenstaeudl, Nils Jenke
Maintainer: Michael Gruenstaeudl \email{m.gruenstaeudl@fu-berlin.de}
}


\examples{

# MH161174
gbkFile <- system.file("extdata", "MH161174/MH161174.gb", package="PACViR")
bamFile <- system.file("extdata", "MH161174/MH161174_PlastomeReadsOnly.sorted.bam", 
                       package="PACViR")
outFile <- paste(getwd(), "/MH161174_AssemblyCoverage_viz.svg", sep="")
PACViR.complete(gbk.file=gbkFile, bam.file=bamFile, windowSize=250, 
                mosdepthCmd='mosdepth', threshold=15, delete=FALSE, output=outFile)

# MH899017
gbkFile <- system.file("extdata", "MH899017/MH899017.gb", package="PACViR")
bamFile <- system.file("extdata", "MH899017/MH899017_PlastomeReadsOnly.sorted.bam", 
                       package="PACViR")
outFile <- paste(getwd(), "/MH899017_AssemblyCoverage_viz.svg", sep="")
PACViR.complete(gbk.file=gbkFile, bam.file=bamFile, windowSize=250, 
                mosdepthCmd='mosdepth', threshold=15, delete=FALSE, output=outFile)

}

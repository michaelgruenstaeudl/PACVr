\name{PACViR-package}
\alias{PACViR-package}
\alias{PACViR}
\docType{package}

\title{
  Plastome Assembly Coverage Visualization in R
}

\description{
  \pkg{PACViR} provides functions to read GenBank and BAM file formats to visualize plastome data in a circular way using     RCircos-package. This is meant to be a further quality check for NGS downstream analysis. Especially for assembly and       annotation steps.
}

\note{

  \bold{Software Requirements}
  
  Mandatory requirements for \pkg{PACViR}:

  External installation of Mosdepth (>= 0.2.2) following \url{https://github.com/brentp/mosdepth}.
  
  
  
  \bold{Dataset Requirements}
  
  The user-supplied data set must contain Plastome data only with information about the exact postion of the repeat           regions. Additionally to that those regions have to be named as 'IRa', 'IRb' or 'Inverted Repeat a', 'Inverted Repeat b'.   The sequence length of the genome should be within 100kb - 200kb, otherwise it may lead to an erroneous visualization.


  \bold{Input File Requirements}

  In order to execute \pkg{PACViR}, a user must provide paths to two different types of input files: (a) a file in            GenBank format and (b) a file in BAM format. The GenBank Flat File Format needs to be exactly as described in the sample    GenBank record of NCBI (\url{https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html}). The BAM File Format has to fit the specifications described in the Sequence Alignment/Map Format Specification documentation (\url{https://samtools.github.io/hts-specs/SAMv1.pdf}). Note that the accompanying BAM index file is mandatory for the use of \pkg{PACViR} and must, thus, be present. 
}

\author{
  Michael Gruenstaeudl, Nils Jenke

  Maintainer: Michael Gruenstaeudl \email{mi.gruenstaeudl@gmail.com}, Nils Jenke \email{nilsj24@zedat.fu-berlin.de}
}

\references{

}
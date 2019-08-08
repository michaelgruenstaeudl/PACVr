\name{PACVr-package}
\alias{PACVr-package}
\alias{PACVr}
\docType{package}

\title{
  Plastome Assembly Coverage Visualization in R
}

\description{
  \pkg{PACVr} visualizes the coverage depth of a complete plastid genome as well as the equality of its inverted repeat regions in relation to the circular, quadripartite genome structure and the location of individual genes.
}

\note{
  \bold{Software Dependencies}
  
  For full functionality, \pkg{PACVr} requires the installation of mosdepth (\url{https://github.com/brentp/mosdepth}) on the system.
  
  \bold{Input Requirements}
  
  The input to \pkg{PACVr} consists of two input files that contain information on genome sequence and structure as well as coverage depth. Specifically, users must provide (a) a file in GenBank flat file format that complies with the GenBank record specifications (\url{https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html}), and (b) a file in BAM format that complies with the specifications described in the Sequence Alignment/Map Format documentation (\url{https://samtools.github.io/hts-specs/SAMv1.pdf}) and is accompanied by an ancillary index file.
  
  \bold{Data Requirements}

  The user-supplied GenBank flat file must contain a sequence record of a complete plastid genome that contains a sequence features for each of the inverted repeat regions. Specifically, the inverted repeat regions must be named 'IRa' and 'IRb', or 'Inverted Repeat a' and 'Inverted Repeat b'. Moreover, the total sequence length of the sequence record must be between 100kb and 200kb.
  
}

\author{
  Michael Gruenstaeudl, Nils Jenke

  Maintainer: Michael Gruenstaeudl \email{m.gruenstaeudl@fu-berlin.de}
}

\references{
  Gruenstaeudl, M. and Jenke, N. (2019). PACVr: Plastome Assembly Coverage Visualization in R. bioRxiv 697821; doi: https://doi.org/10.1101/697821
}

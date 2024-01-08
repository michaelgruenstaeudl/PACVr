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
  \bold{Input Requirements}
  The input to \pkg{PACVr} consists of two input files that contain information on (a) genome sequence and structure, and (b) on coverage depth. Specifically, users must provide (a) a file in GenBank flat file format that complies with the GenBank record specifications (\url{https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html}), and (b) a file in BAM format that complies with the specifications described in the Sequence Alignment/Map Format documentation (\url{https://samtools.github.io/hts-specs/SAMv1.pdf}) and is accompanied by an ancillary index file.
  
  \bold{Data Requirements}
  The user-supplied GenBank flat file must contain a sequence record of a complete plastid genome. For effective visualizations, the sequence record should contain a sequence feature for each of the inverted repeat regions, whereby the regions should be named 'IRa' and 'IRb' or 'Inverted Repeat a' and 'Inverted Repeat b', respectively. For effective visualizations, the total sequence length of the complete genome should be between 50 kb and 250 kb.
}

\author{
  Gregory Smith, Nils Jenke, Michael Gruenstaeudl
  Maintainer: Michael Gruenstaeudl \email{m_gruenstaeudl@fhsu.edu}
}

\references{
  Gruenstaeudl, M. and Jenke, N. (2020). PACVr: plastome assembly coverage visualization in R. BMC Bioinformatics 21: 207. doi: https://doi.org/10.1186/s12859-020-3475-0
}

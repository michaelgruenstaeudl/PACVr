\name{PACVr-package}
\alias{PACVr-package}
\alias{PACVr}
\docType{package}

\title{
  Plastome Assembly Coverage Visualization in R
}

\description{
  \pkg{PACVr} is a user-friendly software tool to visualize the coverage depth of a complete plastid genome as well as the equality of its IR regions while simultaneously accounting for the circular, quadripartite structure of the genome, thus providing an aid to optimize the process of plastid genome assembly.
}

\note{
  \bold{Software Requirements}
  
  Optional dependency for \pkg{PACVr}: Mosdepth (installation via \url{https://github.com/brentp/mosdepth})
  
  \bold{Dataset Requirements}
  
  The user-supplied data set must contain plastome data with information about the exact postion of the repeat regions. Those regions have to be named as 'IRa', 'IRb' or 'Inverted Repeat a', 'Inverted Repeat b'. The sequence length of the genome should be within 100kb-200kb; otherwise it may lead to an erroneous visualization.
  
  
  \bold{Input File Requirements}
  
  In order to execute \pkg{PACVr}, a user must provide paths to two different types of input files: (a) a file in GenBank format and (b) a file in BAM format. The GenBank flat file format needs to comply with the GenBank record specifications of NCBI (\url{https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html}). The BAM file format needs to comply with the specifications described in the Sequence Alignment/Map Format documentation (\url{https://samtools.github.io/hts-specs/SAMv1.pdf}). Please note that the accompanying BAM index file is mandatory for the use of \pkg{PACVr} and must, thus, be present.
}

\author{
  Michael Gruenstaeudl, Nils Jenke

  Maintainer: Michael Gruenstaeudl \email{m.gruenstaeudl@fu-berlin.de}
}

\references{
  Gruenstaeudl M., Jenke N. (2019) PACVr: Plastome Assembly Coverage Visualization in R. R package version 0.7.5. Manuscript in prep.
}

#!/usr/bin/R
#contributors=c("Gregory Smith", "Nils Jenke", "Michael Gruenstaeudl")
#email="m_gruenstaeudl@fhsu.edu"
#version="2023.11.23.1530"

CovCalc <- function(bamFile, windowSize = 250) {
  # Calculates coverage of a given bam file and stores data in data.frame format
  # ARGS:
  #     bamFile: bam file to calculate coverage
  #     windowSize: numeric value to specify the coverage calculation window
  #     output: name and directory of output file
  # RETURNS:
  #     data.frame with region names, chromosome start, chromosome end and coverage calcucation
  logger::log_info('  Coverage calculation for `{bamFile}`')
  if (!is.numeric(windowSize) | windowSize < 0) {
    warning("User-selected window size must be >= 1.")
    stop()
  }
  cov <- GenomicAlignments::coverage(bamFile)
  bins <-
    GenomicRanges::tileGenome(
      sum(IRanges::runLength(cov)),
      tilewidth = windowSize,
      cut.last.tile.in.chrom = TRUE
    )
  cov <- GenomicRanges::binnedAverage(bins, cov, "coverage")
  cov <-
    as.data.frame(cov)[c("seqnames", "start", "end", "coverage")]
  colnames(cov) <-
    c("Chromosome", "chromStart", "chromEnd", "coverage")
  cov$coverage <- ceiling(as.numeric(cov$coverage))
  return(cov)
}

GenerateIRSynteny <- function(genes, syntenyLineType) {
  logger::log_info('  Testing gene synteny in IRs')
  n_occur <- data.frame(table(genes[, 4]), stringsAsFactors = FALSE)
  n_occur <- n_occur[n_occur$Freq == 2,]
  ir_synteny <- c()
  if (syntenyLineType == "1") {
    for (gene in n_occur$Var1) {
      duplicateGene <- genes[which(gene == genes$gene), 1:3]
      ir_synteny <-
        rbind(
          ir_synteny,
          cbind(duplicateGene[1,], duplicateGene[2,], stringsAsFactors = FALSE),
          stringsAsFactors = FALSE
        )
    }
  } else if (syntenyLineType == "2") {
    for (gene in n_occur$Var1) {
      duplicateGene <- genes[which(gene == genes$gene), 1:3]
      duplicateGene[1, 2] <-
        mean(as.numeric(duplicateGene[1, 2:3]))
      duplicateGene[1, 3] <- duplicateGene[1, 2]
      duplicateGene[2, 2] <-
        mean(as.numeric(duplicateGene[2, 2:3]))
      duplicateGene[2, 3] <- duplicateGene[2, 2]
      ir_synteny <-
        rbind(
          ir_synteny,
          cbind(duplicateGene[1,], duplicateGene[2,], stringsAsFactors = FALSE),
          stringsAsFactors = FALSE
        )
    }
  }
  ir_synteny$PlotColor <- "dodgerblue4"
  return(ir_synteny)
}

GenerateHistogramData <-
  function(region, coverage, windowSize, lastOne) {
    # Function to generate line data for RCircos.Line.Plot
    # ARGS:
    #   coverage: data.frame of coverage
    # RETURNS:
    #   data.frame with region means to plot over histogram data
    # Error handling
    logger::log_info('  Generating histogram data for region `{region[4]}`')
    if (lastOne) {
      coverage <-
        coverage[(floor(region[1, 2] / windowSize) + 1):ceiling(region[1, 3] / windowSize),]
    } else {
      coverage <-
        coverage[(floor(region[1, 2] / windowSize) + 1):floor(region[1, 3] / windowSize) +
                   1,]
    }
    coverage[, 4] <- mean(coverage[, 4])
    return(coverage)
  }

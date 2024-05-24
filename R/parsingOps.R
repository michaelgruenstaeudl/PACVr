#!/usr/bin/env RScript
#contributors=c("Gregory Smith", "Nils Jenke", "Michael Gruenstaeudl")
#email="m_gruenstaeudl@fhsu.edu"
#version="2024.05.24.2053"

PACVr.parseGenes <- function (gbkSeqFeatures) {
  # Function to extract gene information from Genbank flatfile data
  # ARGS:
  #   gbkData (resulting data frame from parsing read.gb object)
  # RETURNS:
  #   genes in data frame format
  logger::log_info('  Extracting information on genes')
  gene_L <- read.gbGenes(gbkSeqFeatures)
  gene_L <- gene_L[, c(1:3, which(colnames(gene_L) == "gene"))]
  colnames(gene_L) <- c("Chromosome", "chromStart", "chromEnd", "gene")
  gene_L$Chromosome <- ""
  gene_L <- gene_L[order(gene_L$chromStart),]
  row.names(gene_L) <- 1:nrow(gene_L)
  return(gene_L)
}

PACVr.parseQuadripRegions <- function (gbkLengths,
                                       gbkSeqFeatures,
                                       analysisSpecs) {
  raw_quadripRegions <- ParseQuadripartiteStructure(gbkSeqFeatures)
  if (is.null(raw_quadripRegions)) {
    logger::log_warn("No regions identified; using full genome instead.")
    analysisSpecs$setIRCheckFields(NA)
    quadripRegions <- PACVr.parseSource(gbkSeqFeatures)
  } else {
    quadripRegions <- fillDataFrame(gbkLengths, raw_quadripRegions)
  }
  return(quadripRegions)
}

PACVr.calcCoverage <- function (bamFile,
                                windowSize=250,
                                logScale) {
  logger::log_info('Calculating the sequencing coverage for `{bamFile}`')
  coverageRaw <- GenomicAlignments::coverage(bamFile)
  coveragePlot <- CovCalc(coverageRaw,
                          windowSize,
                          logScale)
  coverage <- list(
    raw = coverageRaw,
    plot = coveragePlot
  )
  return(coverage)
}

PACVr.parseSource <- function(gbkSeqFeatures) {
  type <-
    Band <-
    Stain <-
    Chromosome <-
    start <-
    end <-
    seqnames <-
    NULL
  source <- gbkSeqFeatures %>%
    dplyr::filter(type=="source") %>%
    dplyr::rename(chromStart = start,
                  chromEnd = end,
                  Band = seqnames) %>%
    dplyr::mutate(Chromosome = "",
                  Stain = "gpos75",
                  Band = sub(" ", "_", Band)) %>%
    dplyr::select(dplyr::all_of(c("Chromosome", "chromStart",
                                  "chromEnd", "Band", "Stain")))
  return(source)
}

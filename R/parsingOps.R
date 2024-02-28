#!/usr/bin/env RScript
#contributors=c("Gregory Smith", "Nils Jenke", "Michael Gruenstaeudl")
#email="m_gruenstaeudl@fhsu.edu"
#version="2024.02.28.0231"

ExtractAllGenes <- function(gbkSeqFeatures) {
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

PACVr.gbkData <- function(read.gbData, analysisSpecs) {
  gbkSeqFeatures <- read.gbSeqFeatures(read.gbData, analysisSpecs)
  if (is.null(gbkSeqFeatures)) {
    logger::log_error(paste("No usable sequence features detected to perform specified analysis"))
    return(NULL)
  }

  # derived from gbkData
  gbkSeq <- read.gbSequence(read.gbData)
  lengths <- read.gbLengths(gbkSeq)
  sampleName <- read.gbSampleName(read.gbData)
  plotTitle <- read.gbPlotTitle(read.gbData)
  rm(read.gbData)
  gc()

  # primarily derived from gbkSeqFeatures
  quadripRegions <- PACVr.quadripRegions(lengths,
                                         gbkSeqFeatures,
                                         analysisSpecs$isIRCheck)
  genes <- PACVr.parseGenes(gbkSeqFeatures)
  linkData <- PACVr.linkData(genes,
                             quadripRegions,
                             analysisSpecs$syntenyLineType)
  rm(gbkSeqFeatures)
  gc()

  gbkData <- list(
    genes = genes,
    seq = gbkSeq,
    lengths = lengths,
    sampleName = sampleName,
    plotTitle = plotTitle,
    quadripRegions = quadripRegions,
    linkData = linkData
  )
  return(gbkData)
}

PACVr.parseGenes <- function (gbkSeqFeatures) {
  # This function parses the genes of a GenBank file
  logger::log_info('Parsing the different genes')
  genes <- ExtractAllGenes(gbkSeqFeatures)
  return(genes)
}

PACVr.quadripRegions <- function(gbkLengths,
                                 gbkSeqFeatures,
                                 isIRCheck) {
  if (isIRCheck) {
    logger::log_info('Parsing the quadripartite genome structure')
    quadripRegions <- PACVr.parseQuadripRegions(gbkLengths,
                                                gbkSeqFeatures)
  } else {
    quadripRegions <- PACVr.parseSource(gbkSeqFeatures)
  }
  return(quadripRegions)
}

PACVr.parseQuadripRegions <- function (gbkLengths, gbkSeqFeatures) {
  raw_quadripRegions <- ParseQuadripartiteStructure(gbkSeqFeatures)
  if (is.null(raw_quadripRegions)) {
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

PACVr.generateIRGeneData <- function(genes, quadripRegions,
                                     syntenyLineType) {
  # Parse GenBank file
  if ("IRb" %in% quadripRegions[, 4] &&
    "IRa" %in% quadripRegions[, 4]) {
    linkData <- GenerateIRSynteny(genes, syntenyLineType)
    return(linkData)
  }
  return(-1)
}

PACVr.linkData <- function(genes,
                           quadripRegions,
                           syntenyLineType) {
  linkData <- NULL
  if (!is.null(syntenyLineType)) {
    logger::log_info('Inferring the IR regions and the genes within the IRs')
    linkData <- PACVr.generateIRGeneData(genes,
                                         quadripRegions,
                                         syntenyLineType)
  }
  return(linkData)
}

PACVr.parseSource <- function(gbkSeqFeatures) {
  logger::log_info("No regions identified; using full genome instead.")

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

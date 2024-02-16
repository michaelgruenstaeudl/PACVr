#!/usr/bin/env RScript
#contributors=c("Gregory Smith", "Nils Jenke", "Michael Gruenstaeudl")
#email="m_gruenstaeudl@fhsu.edu"
#version="2024.02.01.1736"

ExtractAllGenes <- function(gbkDataDF) {
  # Function to extract gene information from Genbank flatfile data
  # ARGS:
  #   gbkData (resulting data frame from parsing read.gb object)
  # RETURNS:
  #   genes in data frame format
  logger::log_info('  Extracting information on genes')
  gene_L <- read.gbGenes(gbkDataDF)
  gene_L <- gene_L[, c(1:3, which(colnames(gene_L) == "gene"))]
  colnames(gene_L) <- c("Chromosome", "chromStart", "chromEnd", "gene")
  gene_L$Chromosome <- ""
  gene_L <- gene_L[order(gene_L$chromStart),]
  row.names(gene_L) <- 1:nrow(gene_L)
  return(gene_L)
}

parseSource <- function(gbkDataDF) {
PACVr.parseGenes <- function (gbkDataDF) {
  # This function parses the genes of a GenBank file
  logger::log_info('Parsing the different genes')
  genes <- ExtractAllGenes(gbkDataDF)
  return(genes)
}

PACVr.quadripRegions <- function(gbkLengths,
                                 gbkDataDF,
                                 isIRCheck) {
  if (isIRCheck) {
    logger::log_info('Parsing the quadripartite genome structure')
    quadripRegions <- PACVr.parseQuadripRegions(gbkLengths,
                                                gbkDataDF)
  } else {
    quadripRegions <- PACVr.parseSource(gbkDataDF)
  }
  return(quadripRegions)
}

PACVr.parseQuadripRegions <- function (gbkLengths, gbkDataDF) {
  raw_quadripRegions <- ParseQuadripartiteStructure(gbkDataDF)
  quadripRegions <- fillDataFrame(gbkLengths, raw_quadripRegions)
  return(quadripRegions)
}

PACVr.calcCoverage <- function (coverageRaw,
                                windowSize=250,
                                logScale) {
  logger::log_info('Calculating the sequencing coverage')
  coverage <- CovCalc(coverageRaw, windowSize)
  if (logScale == TRUE) {
    coverage$coverage <- log(cov$coverage)
    #coverage$coverage <- log(coverage$coverage)
  }
  coverage$Chromosome <- ""
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

PACVr.parseSource <- function(gbkDataDF) {
  type <-
    Band <-
    Stain <-
    Chromosome <-
    start <-
    end <-
    seqnames <-
    NULL
  source <- gbkDataDF %>%
    dplyr::filter(type=="source") %>%
    dplyr::rename(chromStart = start,
                  chromEnd = end,
                  Band = seqnames) %>%
    dplyr::mutate(Chromosome = "",
                  Stain = "gpos75") %>%
    dplyr::select(dplyr::all_of(c("Chromosome", "chromStart",
                                  "chromEnd", "Band", "Stain")))
  return(source)
}

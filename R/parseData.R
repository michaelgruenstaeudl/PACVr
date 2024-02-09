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

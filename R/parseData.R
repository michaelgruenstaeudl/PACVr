#!/usr/bin/R
#contributors=c("Gregory Smith", "Nils Jenke", "Michael Gruenstaeudl")
#email="m_gruenstaeudl@fhsu.edu"
#version="2023.11.23.1530"

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

#!/usr/bin/R
#contributors = c("Michael Gruenstaeudl","Nils Jenke")
#email = "m.gruenstaeudl@fu-berlin.de", "nilsj24@zedat.fu-berlin.de"
#version = "2020.01.17.1800"

ExtractAllGenes <- function(gbkData) {
  
  # Function to extract genes from genbank file
  # ARGS:
  #   gbkData: genbank file parsed by genbankr
  # RETURNS:
  #   genes in data frame format
  gene_L <- BiocGenerics::as.data.frame(genbankr::genes(gbkData))
  gene_L <- gene_L[ ,c(1:3,which(colnames(gene_L) == "gene"))]
  colnames(gene_L) <- c("Chromosome","chromStart","ChromEnd","gene")
  gene_L$Chromosome <- ""
  return(gene_L)
}

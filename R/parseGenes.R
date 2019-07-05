#!/usr/bin/R
#contributors = c("Michael Gruenstaeudl","Nils Jenke")
#email = "m.gruenstaeudl@fu-berlin.de", "nilsj24@zedat.fu-berlin.de"
#version = "2019.07.05.1100"

ExtractAllGenes <- function(gbkData) {
  # Function to extract genes from genbank file
  # ARGS:
  #   gbkData: genbank file parsed by genbankr
  # RETURNS:
  #   genes in data frame format
  gene_L <- BiocGenerics::as.data.frame(genbankr::genes(gbkData))
  gene_L <- gene_L[ ,c(1:3,which(colnames(gene_L) == "gene"))]
  gene_L <- Rename_Df(gene_L,
                    "gene")
  return(gene_L)
}


SplitGenesAtRegionBorders <- function(geneData, regionData) {
  # Function to split those genes that occur in two different regions at
  # the region borders
  # ARGS:
  #   geneData: dataframe with gene data
  #   regionData: dataframe with region data
  # RETURNS:
  #   geneData dataframe with splitted geneDatas
  for (i in 1:nrow(regionData)) {
    for (j in 1:nrow(geneData)) {
      if (as.integer(geneData[j,2]) >= as.integer(regionData[i,2]) & 
          as.integer(geneData[j,3]) >  as.integer(regionData[i,3]) & 
          as.integer(geneData[j,2]) <  as.integer(regionData[i,3])){
        geneData[nrow(geneData)+1,] <- c(as.character(geneData[j,1]), regionData[i,3]+1, geneData[j,3], paste(geneData[j,4],"_",regionData[i+1,1],sep = ""))
        geneData[j,3] <- regionData[i,3]
        geneData[j,4] <- paste(geneData[j,4],"_",regionData[i,1],sep = "")
      }
    }
  }
  geneData      <- geneData[order(as.integer(geneData[,2])), ]
  geneData[ ,2] <- as.integer(geneData[ ,2])
  geneData[ ,3] <- as.integer(geneData[ ,3])
  return(geneData)
}

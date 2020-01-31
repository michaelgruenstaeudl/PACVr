#!/usr/bin/R
#contributors = c("Michael Gruenstaeudl","Nils Jenke")
#email = "m.gruenstaeudl@fu-berlin.de", "nilsj24@zedat.fu-berlin.de"
#version = "2020.01.17.1800"

checkIRSynteny <- function(gbkData, regions){
  
  gbkSeq <- genbankr::getSeq(gbkData)
  repeatB <- as.numeric(regions[which(regions[,4] == "IRb"),2:3])
  repeatA <- as.numeric(regions[which(regions[,4] == "IRa"),2:3])
  if(repeatB[2] - repeatB[1] != repeatA[2] - repeatA[1]){
    message("Inverted repeats differ in sequence length")
  }
  if(gbkSeq[[1]][repeatB[1]:repeatB[2]] != Biostrings::reverseComplement(gbkSeq[[1]][repeatA[1]:repeatA[2]])){
    message("Inverted repeats differ in sequence")
  }
}

GenerateIRSynteny <- function(genes, syntenyLineType) {
  
  n_occur <- data.frame(table(genes[,4]),stringsAsFactors = FALSE)
  n_occur <- n_occur[n_occur$Freq == 2,]
  ir_synteny <- c()
  
  if (syntenyLineType == "1"){
    for (gene in n_occur$Var1) {
      duplicateGene <- genes[which(gene == genes$gene),1:3]
      ir_synteny <- rbind(ir_synteny,cbind(duplicateGene[1,],duplicateGene[2,],stringsAsFactors=FALSE),stringsAsFactors=FALSE) 
    }
  } else if (syntenyLineType == "2"){
      for (gene in n_occur$Var1) {
        duplicateGene <- genes[which(gene == genes$gene),1:3]
        duplicateGene[1,2] <- mean(as.numeric(duplicateGene[1,2:3]))
        duplicateGene[1,3] <- duplicateGene[1,2]
        duplicateGene[2,2] <- mean(as.numeric(duplicateGene[2,2:3]))
        duplicateGene[2,3] <- duplicateGene[2,2]
        ir_synteny <- rbind(ir_synteny,cbind(duplicateGene[1,],duplicateGene[2,],stringsAsFactors=FALSE),stringsAsFactors=FALSE) 
      }
    }
  ir_synteny$PlotColor <- "dodgerblue4"
  return(ir_synteny)
}

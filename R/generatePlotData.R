#!/usr/bin/R
#contributors = c("Michael Gruenstaeudl","Nils Jenke")
#email = "m.gruenstaeudl@fu-berlin.de", "nilsj24@zedat.fu-berlin.de"
#version = "2020.07.29.1700"

CovCalc <- function(bamFile, windowSize=250){
  
  # Calculates coverage of a given bam file and stores data in data.frame format
  # ARGS:
  #     bamFile: bam file to calculate coverage
  #     windowSize: numeric value to specify the coverage calculation window
  #     output: name and directory of output file
  # RETURNS:
  #     data.frame with region names, chromosome start, chromosome end and coverage calcucation
  if (!is.numeric(windowSize) | windowSize < 0) {
    warning("User-selected window size must be >= 1.")
    stop()
  }
  
  cov <- GenomicAlignments::coverage(bamFile)
  bins <- GenomicRanges::tileGenome(sum(IRanges::runLength(cov)), tilewidth=windowSize, cut.last.tile.in.chrom=TRUE)
  cov <- GenomicRanges::binnedAverage(bins, cov, "coverage")
  cov <- as.data.frame(cov)[c("seqnames","start","end","coverage")]
  colnames(cov) <- c("Chromosome","chromStart","chromEnd","coverage")
  cov$coverage <- ceiling(as.numeric(cov$coverage))
  
  return(cov)
}

checkIREquality <- function(gbkData, regions){
  
  gbkSeq <- genbankr::getSeq(gbkData)

  repeatB <- as.numeric(regions[which(regions[,4] == "IRb"),2:3])
  repeatA <- as.numeric(regions[which(regions[,4] == "IRa"),2:3])
  if(repeatB[2] - repeatB[1] != repeatA[2] - repeatA[1]){
    message("WARNING: Inverted repeats differ in sequence length")
    message(paste("The IRb has a total lengths of: ", repeatB[2]-repeatB[1], " bp", sep=""))
    message(paste("The IRa has a total lengths of: ", repeatA[2]-repeatA[1], " bp", sep=""))
  }
  if(gbkSeq[[1]][repeatB[1]:repeatB[2]] != Biostrings::reverseComplement(gbkSeq[[1]][repeatA[1]:repeatA[2]])){
    IRa_seq <- Biostrings::DNAString(gbkSeq[[1]][repeatB[1]:repeatB[2]])
    IRa_seq <- split(IRa_seq, ceiling(seq_along(IRa_seq)/10000))
    IRb_seq <- Biostrings::DNAString(Biostrings::reverseComplement(gbkSeq[[1]][repeatA[1]:repeatA[2]]))
    IRb_seq <- split(IRb_seq, ceiling(seq_along(IRb_seq)/10000))
    
    IR_diff_SNPS <- c()
    IR_diff_gaps <- c()
    for(i in  1:min(length(IRa_seq), length(IRb_seq))){
      subst_mat <- Biostrings::nucleotideSubstitutionMatrix(match=1, mismatch=-3, baseOnly=TRUE)
      globalAlign <- tryCatch(
        {Biostrings::pairwiseAlignment(IRa_seq[[i]], IRb_seq[[i]], substitutionMatrix=subst_mat, gapOpening=5, gapExtension=2)},
        error = function(e) {return(NULL)})
      if(is.null(globalAlign)) break
      IR_diff_SNPS <- c(IR_diff_SNPS, which(strsplit(Biostrings::compareStrings(globalAlign), "")[[1]]=="?"))
      IR_diff_gaps <- c(IR_diff_gaps, which(strsplit(Biostrings::compareStrings(globalAlign), "")[[1]]=="-"))
    }
    message("WARNING: Inverted repeats differ in sequence")
    if(length(IR_diff_SNPS)>0){
      message(paste("When aligned, the IRs differ through a total of ", length(IR_diff_SNPS)," SNPS. These SNPS are located at the following nucleotide positions: ", paste(unlist(IR_diff_SNPS), collapse=" "), sep=""))
    }
    if(length(IR_diff_gaps)>0){
      message(paste("When aligned, the IRs differ through a total of ", length(IR_diff_gaps)," gaps. These gaps are located at the following nucleotide positions: ", paste(unlist(IR_diff_gaps), collapse=" "), sep=""))
    }
  }
  message("Proceeding with coverage depth visualization, but without quadripartite genome structure ...")
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

GenerateHistogramData <- function(region,coverage, windowSize, lastOne) {
  
  # Function to generate line data for RCircos.Line.Plot
  # ARGS:
  #   coverage: data.frame of coverage
  # RETURNS:
  #   data.frame with region means to plot over histogram data
  # Error handling
  if (lastOne){
    coverage <- coverage[(floor(region[1,2] / windowSize)+1):ceiling(region[1,3] / windowSize),]
  } else {
    coverage <- coverage[(floor(region[1,2] / windowSize)+1):floor(region[1,3] / windowSize)+1,]
  }
  coverage[,4] <- mean(coverage[,4])
  return(coverage)
}

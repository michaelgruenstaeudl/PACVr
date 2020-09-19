#!/usr/bin/R
#contributors = c("Michael Gruenstaeudl","Nils Jenke")
#email = "m.gruenstaeudl@fu-berlin.de", "nilsj24@zedat.fu-berlin.de"
#version = "2020.07.29.1700"

HistCol <- function(cov, threshold, relative, logScale) {
  
  # Function to generate color vector for histogram data
  # ARGS:
  #       cov:       data.frame of coverage
  #       threshold: numeric value of a specific threshold
  # RETURNS:
  #   color vector
  # Error handling
  if (!is.numeric(threshold) | threshold < 0) {
    warning("User-defined coverage depth threshold must be >=1.")
    stop()
  }
  
  if(relative == TRUE & logScale){
    threshold <- mean(cov[,4]) + log(threshold)
    
  }else if(relative == TRUE){
    threshold <- mean(cov[,4]) * threshold
    
  }
  color <- rep("black",nrow(cov))
  ind   <- as.numeric(cov[ ,4]) <= threshold
  color <- replace(color,ind,"red")
  return(color)
}

boolToDeci <- function(boolList) {
  
  out = 0
  boolList <- rev(boolList)
  for (i in 1:length(boolList)){
    out = out + boolList[i]*(2^(i-1))
  }
  return(out)
}

writeTables <- function(regions, genes, cov, relative, threshold, dir, sample) {
  ir_regions <- IRanges::IRanges(start=regions$chromStart, end=regions$chromEnd, names = regions$Band) 
  ir_genes <- IRanges::IRanges(start=genes$chromStart, end=genes$chromEnd, names = genes$gene)
  ir_noncoding <- IRanges::gaps(ir_genes, start = min(regions$chromStart), end = max(regions$chromEnd))
  noncoding <- as.data.frame(ir_noncoding)[,c("start","end")]
  ir_cov <- IRanges::IRanges(start=cov$chromStart, end=cov$chromEnd)
  
  overlaps_genes <- IRanges::findOverlaps(query = ir_genes, subject = ir_regions, type="any", select = "all")
  genes$Chromosome <- sapply(overlaps_genes, function(x) regions$Band[x])
  genes$Chromosome <- gsub("^c\\(\"|\"|\"|\"\\)$", "",as.character(genes$Chromosome))
  
  overlaps_cov <- IRanges::findOverlaps(query = ir_cov, subject = ir_regions, type="any", select = "all")
  cov$Chromosome <- sapply(overlaps_cov, function(x) regions$Band[x])
  cov$Chromosome <- gsub("^c\\(\"|\"|\"|\"\\)$", "",as.character(cov$Chromosome))
  
  overlaps_genes_cov <- IRanges::findOverlaps(query = ir_cov, subject = ir_genes, type = "any")
  genes["coverage"] <- ceiling(tapply(cov$coverage[overlaps_genes_cov@from], overlaps_genes_cov@to, mean))
  
  overlaps_noncoding_cov <- IRanges::findOverlaps(query = ir_cov, subject = ir_noncoding, type = "any")
  noncoding["coverage"] <- ceiling(tapply(cov$coverage[overlaps_noncoding_cov@from], overlaps_noncoding_cov@to, mean))
  
  overlaps_noncoding <- IRanges::findOverlaps(query = ir_noncoding, subject = ir_regions, type="any", select = "all")
  noncoding$Chromosome <- sapply(overlaps_noncoding, function(x) regions$Band[x])
  noncoding$Chromosome <- gsub("^c\\(\"|\"|\"|\"\\)$", "",as.character(noncoding$Chromosome))
  
  noncoding <- noncoding[,c("Chromosome", "start", "end","coverage")]
  colnames(noncoding) <- c("Chromosome", "chrStart", "chrEnd", "coverage")
  
  if (relative == TRUE) {
    genes$lowCoverage <- genes$coverage < mean(genes$coverage) * threshold
    cov$lowCoverage <- cov$coverage < mean(cov$coverage) * threshold
    noncoding$lowCoverage <- noncoding$coverage < mean(noncoding$coverage) * threshold
  } else {
    genes$lowCoverage <- genes$coverage < mean(genes$coverage)
    cov$lowCoverage <- cov$coverage < mean(cov$coverage)
    noncoding$lowCoverage <- noncoding$coverage < mean(noncoding$coverage)
    
  }
  genes$lowCoverage[genes$lowCoverage == TRUE] <- "*"
  genes$lowCoverage[genes$lowCoverage == FALSE] <- ""
  cov$lowCoverage[cov$lowCoverage == TRUE] <- "*"
  cov$lowCoverage[cov$lowCoverage == FALSE] <- ""
  noncoding$lowCoverage[noncoding$lowCoverage == TRUE] <- "*"
  noncoding$lowCoverage[noncoding$lowCoverage == FALSE] <- ""
  
  write.csv(genes, paste(dir, .Platform$file.sep, sample,"_coverage.genes.bed", sep=""), 
            row.names = FALSE, quote = FALSE)
  write.csv(cov, paste(dir, .Platform$file.sep, sample,"_coverage.regions.bed", sep=""), 
            row.names = FALSE, quote = FALSE)
  write.csv(noncoding, paste(dir, .Platform$file.sep, sample,"_coverage.noncoding.bed", sep=""), 
            row.names = FALSE, quote = FALSE)
}

checkIREquality <- function(gbkData, regions, dir, sample){
  gbkSeq <- genbankr::getSeq(gbkData)
  if("IRb" %in% regions[,4] && "IRa" %in% regions[,4]){
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
    write.csv(data.frame(Number_N = Biostrings::alphabetFrequency(gbkSeq)["N"], Mismatches = length(IR_diff_SNPS)+length(IR_diff_gaps)), paste(dir, .Platform$file.sep, sample,"_IR_quality.csv", sep=""), 
              row.names = FALSE, quote = FALSE)
  } else {
    write.csv(data.frame(Number_N = Biostrings::alphabetFrequency(gbkSeq)["N"], Mismatches = NA), paste(dir, .Platform$file.sep, sample,"_IR_quality.csv", sep=""), 
              row.names = FALSE, quote = FALSE)
  }
}
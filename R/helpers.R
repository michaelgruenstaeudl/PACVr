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
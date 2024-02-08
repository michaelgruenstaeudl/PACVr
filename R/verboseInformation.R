#!/usr/bin/env RScript
#contributors=c("Gregory Smith", "Nils Jenke", "Michael Gruenstaeudl")
#email="m_gruenstaeudl@fhsu.edu"
#version="2024.02.08.2300"

printCovStats <- function(gbkData,
                          bamFile,
                          genes,
                          quadripRegions,
                          analysisSpecs,
                          output) {
  sampleName <- read.gbSampleName(gbkData)
  # Step 1. Check ...
  if (!is.na(output)) {
    outDir <- dirname(output)
    tmpDir <- file.path(outDir,
                        paste(sampleName["sample_name"],
                              ".tmp",
                              sep=""))
  } else {
    tmpDir <-
      file.path(".", paste(sampleName["sample_name"],
                           ".tmp",
                           sep=""))
  }
  # Step 2. Check ...
  if (dir.exists(tmpDir) == FALSE) {
    dir.create(tmpDir)
  }
  # Step 3. Write output
  printCovValsAsTable(quadripRegions, bamFile, genes, tmpDir, sampleName)
  if (!is.null(analysisSpecs$syntenyLineType)) {
    checkIREquality(gbkData, quadripRegions, tmpDir, sampleName)
  }
}

printCovValsAsTable <- function(regions, bamFile, genes, dir, sample_name) {
  covData <- getCovData(regions, genes)
  covData <- update.ir_genes(regions, bamFile, sample_name, covData)
  covData <- update.ir_noncoding(regions, bamFile, sample_name, covData)
  covData <- update.ir_regions(bamFile, sample_name, covData)
  covData <- setLowCoverage(covData)

  # Writing values to output table
  writeCovTables(covData, sample_name, dir)
}

getCovData <- function(regions, genes) {
  ir_regions <- IRanges::IRanges(
    start = regions$chromStart,
    end = regions$chromEnd,
    names = regions$Band
  )
  ir_genes <- unlist(IRanges::slidingWindows(
    IRanges::reduce(
      IRanges::IRanges(
        start = genes$chromStart,
        end = genes$chromEnd,
        names = genes$gene
      )
    ),
    width = 250L,
    step = 250L
  ))
  ir_noncoding <- unlist(IRanges::slidingWindows(
    IRanges::reduce(IRanges::gaps(
      ir_genes,
      start = min(regions$chromStart),
      end = max(regions$chromEnd)
    )),
    width = 250L,
    step = 250L
  ))
  overlaps_genes <- IRanges::findOverlaps(
    query = ir_genes,
    subject = ir_regions,
    type = "any",
    select = "all"
  )
  overlaps_noncoding <- IRanges::findOverlaps(
    query = ir_noncoding,
    subject = ir_regions,
    type = "any",
    select = "all"
  )

  covData <- list(
    ir_regions = ir_regions,
    ir_genes = ir_genes,
    ir_noncoding = ir_noncoding,
    overlaps_genes = overlaps_genes,
    overlaps_noncoding = overlaps_noncoding
  )
  return(covData)
}

update.ir_genes <- function(regions, bamFile, sample_name, covData) {
  ir_genes <- GenomicRanges::GRanges(seqnames = sample_name["genome_name"], covData$ir_genes)
  ir_genes <- GenomicRanges::binnedAverage(ir_genes, GenomicAlignments::coverage(bamFile), "coverage")
  ir_genes <- as.data.frame(ir_genes)[c("seqnames", "start", "end", "coverage")]
  colnames(ir_genes) <- c("Chromosome", "chromStart", "chromEnd", "coverage")
  ir_genes$coverage <- ceiling(as.numeric(ir_genes$coverage))
  ir_genes$Chromosome <- sapply(covData$overlaps_genes, function(x) regions$Band[x])
  ir_genes$Chromosome <- gsub("^c\\(\"|\"|\"|\"\\)$", "", as.character(ir_genes$Chromosome))

  covData$ir_genes <- ir_genes
  return(covData)
}

update.ir_noncoding <- function(regions, bamFile, sample_name, covData) {
  ir_noncoding <- GenomicRanges::GRanges(seqnames = sample_name["genome_name"], covData$ir_noncoding)
  ir_noncoding <- GenomicRanges::binnedAverage(ir_noncoding, GenomicAlignments::coverage(bamFile), "coverage")
  ir_noncoding <- as.data.frame(ir_noncoding)[c("seqnames", "start", "end", "coverage")]
  colnames(ir_noncoding) <- c("Chromosome", "chromStart", "chromEnd", "coverage")
  ir_noncoding$coverage <- ceiling(as.numeric(ir_noncoding$coverage))
  ir_noncoding$Chromosome <- sapply(covData$overlaps_noncoding, function(x) regions$Band[x])
  ir_noncoding$Chromosome <- gsub("^c\\(\"|\"|\"|\"\\)$", "", as.character(ir_noncoding$Chromosome))

  covData$ir_noncoding <- ir_noncoding
  return(covData)
}

update.ir_regions <- function(bamFile, sample_name, covData) {
  ir_regions <- unlist(IRanges::slidingWindows(covData$ir_regions, width = 250L, step = 250L))
  ir_regions <- GenomicRanges::GRanges(seqnames = sample_name["genome_name"], ir_regions)
  ir_regions <- GenomicRanges::binnedAverage(ir_regions, GenomicAlignments::coverage(bamFile), "coverage")
  chr <- ir_regions@ranges@NAMES
  ir_regions <- as.data.frame(ir_regions, row.names = NULL)[c("seqnames", "start", "end", "coverage")]
  ir_regions["seqnames"] <- chr
  colnames(ir_regions) <- c("Chromosome", "chromStart", "chromEnd", "coverage")
  ir_regions$coverage <- ceiling(as.numeric(ir_regions$coverage))

  covData$ir_regions <- ir_regions
  return(covData)
}

setLowCoverage <- function(covData) {
  # ir_regions
  ir_regions <- covData$ir_regions
  cov_regions <-
    aggregate(
      coverage ~ Chromosome,
      data = ir_regions,
      FUN = function(x)
        ceiling(mean(x) - sd(x))
    )
  ir_regions$lowCoverage <- with(ir_regions, ir_regions$coverage < cov_regions$coverage[match(Chromosome, cov_regions$Chromosome)])
  ir_regions$lowCoverage[ir_regions$lowCoverage == TRUE] <- "*"
  ir_regions$lowCoverage[ir_regions$lowCoverage == FALSE] <- ""
  covData$ir_regions <- ir_regions

  # ir_genes
  ir_genes <- covData$ir_genes
  ir_genes$lowCoverage <- ir_genes$coverage < mean(ir_genes$coverage) - sd(ir_genes$coverage)
  ir_genes$lowCoverage[ir_genes$lowCoverage == TRUE] <- "*"
  ir_genes$lowCoverage[ir_genes$lowCoverage == FALSE] <- ""
  covData$ir_genes <- ir_genes

  # ir_noncoding
  ir_noncoding <- covData$ir_noncoding
  ir_noncoding$lowCoverage <- ir_noncoding$coverage < mean(ir_noncoding$coverage) - sd(ir_noncoding$coverage)
  ir_noncoding$lowCoverage[ir_noncoding$lowCoverage == TRUE] <- "*"
  ir_noncoding$lowCoverage[ir_noncoding$lowCoverage == FALSE] <- ""
  covData$ir_noncoding <- ir_noncoding

  return(covData)
}

writeCovTables <- function(covData, sample_name, dir) {
  write.table(
    covData$ir_genes,
    paste(dir, .Platform$file.sep, sample_name["sample_name"], "_coverage.genes.tsv", sep = ""),
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )
  write.table(
    covData$ir_regions,
    paste(dir, .Platform$file.sep, sample_name["sample_name"], "_coverage.regions.tsv", sep = ""),
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )
  write.table(
    covData$ir_noncoding,
    paste(dir, .Platform$file.sep, sample_name["sample_name"], "_coverage.noncoding.tsv", sep = ""),
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )
}

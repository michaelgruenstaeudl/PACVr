#!/usr/bin/env RScript
#contributors=c("Gregory Smith", "Nils Jenke", "Michael Gruenstaeudl")
#email="m_gruenstaeudl@fhsu.edu"
#version="2024.02.08.2300"

getVerbosePath <- function(sampleName, output) {
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
  return(tmpDir)
}

printCovStats <- function(bamFile,
                          genes,
                          quadripRegions,
                          sampleName,
                          analysisSpecs,
                          dir) {
  logger::log_info('Generating statistical information on the sequencing coverage')
  covData <- getCovData(quadripRegions, genes, analysisSpecs)
  covData <- update.ir_genes(quadripRegions, bamFile, sampleName, covData)
  covData <- update.ir_noncoding(quadripRegions, bamFile, sampleName, covData)
  covData <- update.ir_regions(bamFile, sampleName, covData)
  covData <- setLowCoverage(covData)

  # Writing values to output table
  writeCovTables(covData, sampleName, dir)

  # Getting and writing summarized coverage data, grouped by `quadripRegions`
  covSummaries <- getCovSummaries(covData)
  writeCovSumTables(covSummaries, sampleName, dir)
}

getCovData <- function(regions, genes, analysisSpecs) {
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

  if (analysisSpecs$isIRCheck) {
    regions_name <- "Chromosome"
    regions_start <- "chromStart"
    regions_end <- "chromEnd"
  } else {
    regions_name <- "Source"
    regions_start <- "srcStart"
    regions_end <- "srcEnd"
  }

  covData <- list(
    ir_regions = ir_regions,
    ir_genes = ir_genes,
    ir_noncoding = ir_noncoding,
    overlaps_genes = overlaps_genes,
    overlaps_noncoding = overlaps_noncoding,
    regions_name = regions_name,
    regions_start = regions_start,
    regions_end = regions_end
  )
  return(covData)
}

update.ir_genes <- function(regions, bamFile, sample_name, covData) {
  ir_genes <- GenomicRanges::GRanges(seqnames = sample_name["genome_name"], covData$ir_genes)
  ir_genes <- GenomicRanges::binnedAverage(ir_genes, GenomicAlignments::coverage(bamFile), "coverage")
  ir_genes <- as.data.frame(ir_genes)[c("seqnames", "start", "end", "coverage")]
  regions_name <- covData$regions_name
  colnames(ir_genes) <- c(regions_name, covData$regions_start, covData$regions_end, "coverage")
  ir_genes$coverage <- ceiling(as.numeric(ir_genes$coverage))
  ir_genes[[regions_name]] <- sapply(covData$overlaps_genes, function(x) regions$Band[x])
  ir_genes[[regions_name]] <- gsub("^c\\(\"|\"|\"|\"\\)$", "", as.character(ir_genes[[regions_name]]))

  covData$ir_genes <- ir_genes
  return(covData)
}

update.ir_noncoding <- function(regions, bamFile, sample_name, covData) {
  ir_noncoding <- GenomicRanges::GRanges(seqnames = sample_name["genome_name"], covData$ir_noncoding)
  ir_noncoding <- GenomicRanges::binnedAverage(ir_noncoding, GenomicAlignments::coverage(bamFile), "coverage")
  ir_noncoding <- as.data.frame(ir_noncoding)[c("seqnames", "start", "end", "coverage")]
  regions_name <- covData$regions_name
  colnames(ir_noncoding) <- c(regions_name, covData$regions_start, covData$regions_end, "coverage")
  ir_noncoding$coverage <- ceiling(as.numeric(ir_noncoding$coverage))
  ir_noncoding[[regions_name]] <- sapply(covData$overlaps_noncoding, function(x) regions$Band[x])
  ir_noncoding[[regions_name]] <- gsub("^c\\(\"|\"|\"|\"\\)$", "", as.character(ir_noncoding[[regions_name]]))

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
  colnames(ir_regions) <- c(covData$regions_name, covData$regions_start, covData$regions_end, "coverage")
  ir_regions$coverage <- ceiling(as.numeric(ir_regions$coverage))

  covData$ir_regions <- ir_regions
  return(covData)
}

setLowCoverage <- function(covData) {
  # ir_regions
  ir_regions <- covData$ir_regions
  regions_name <- covData$regions_name
  aggFormula <- as.formula(paste("coverage ~", regions_name))
  cov_regions <-
    aggregate(
      aggFormula,
      data = ir_regions,
      FUN = function(x)
        ceiling(mean(x) - sd(x))
    )
  ir_regions$lowCoverage <- ir_regions$coverage < cov_regions$coverage[match(ir_regions[[regions_name]],
                                                                             cov_regions[[regions_name]])]
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
  writeVerboseTable(covData$ir_genes, sample_name, dir, "coverage.genes")
  writeVerboseTable(covData$ir_regions, sample_name, dir, "coverage.regions")
  writeVerboseTable(covData$ir_noncoding, sample_name, dir, "coverage.noncoding")
}

writeVerboseTable <- function(df, sample_name, dir, fileName) {
  write.table(
    df,
    paste(dir, .Platform$file.sep, sample_name["sample_name"], "_", fileName, ".tsv", sep = ""),
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )
}

# adapted from `nilsj9/PlastidSequenceCoverage`
getCovSummaries <- function(covData) {
  regions_name <- covData$regions_name

  covSummaries <- getCovDepths(covData, regions_name)
  regions_evenness <- getCovEvenness(covData$ir_regions, regions_name)
  covSummaries$regions_summary <- dplyr::full_join(covSummaries$regions_summary,
                                                   regions_evenness,
                                                   regions_name)
  return(covSummaries)
}

getCovDepths <- function(covData, regions_name) {
  regions_depth <- getCovDepth(covData$ir_regions, regions_name)
  genes_depth <- getCovDepth(covData$ir_genes, regions_name)
  noncoding_depth <- getCovDepth(covData$ir_noncoding, regions_name)
  
  covDepths <- list(
    regions_summary = regions_depth,
    genes_summary = genes_depth,
    noncoding_summary = noncoding_depth
  )
  return(covDepths)
}

getCovDepth <- function(covDataField, regions_name) {
  covDepth <- covDataField %>%
    groupByRegionsName(regions_name) %>%
    dplyr::summarise(
      lowCoverage = sum(lowCoverage == "*", na.rm = TRUE),
      .groups = "drop"
    )
  return(covDepth)
}

getCovEvenness <- function(covDataField, regions_name) {
  covEvenness <- covDataField %>%
    groupByRegionsName(regions_name) %>%
    dplyr::summarise(
      evenness = evennessScore(coverage),
      .groups = "drop"
    )
  return(covEvenness)
}

groupByRegionsName <- function(df, regions_name) {
  return(
    dplyr::group_by(df, dplyr::across(dplyr::all_of(regions_name)))
  )
}

evennessScore <- function(coverage) {
  coverage_mean <- round(mean(coverage))
  D2 <- coverage[coverage <= coverage_mean]
  E <- 1 - (length(D2) - sum(D2) / coverage_mean) / length(coverage)
  return(E)
}

writeCovSumTables <- function(covSummaries, sample_name, dir) {
  writeVerboseTable(covSummaries$genes_summary, sample_name, dir, "coverage.summary.genes")
  writeVerboseTable(covSummaries$regions_summary, sample_name, dir, "coverage.summary.regions")
  writeVerboseTable(covSummaries$noncoding_summary, sample_name, dir, "coverage.summary.noncoding")
}

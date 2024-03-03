#!/usr/bin/env RScript
#contributors=c("Gregory Smith", "Nils Jenke", "Michael Gruenstaeudl")
#email="m_gruenstaeudl@fhsu.edu"
#version="2024.03.03.0433"

CovCalc <- function(coverageRaw,
                    windowSize = 250,
                    logScale) {
  # Calculates coverage of a given bam file and stores data in data.frame format
  # ARGS:
  #     coverageRaw: coverage data from `GenomicAlignments::coverage()` on bam file
  #     windowSize: numeric value to specify the coverage calculation window
  # RETURNS:
  #     data.frame with region names, chromosome start, chromosome end and coverage calcucation
  if (!is.numeric(windowSize) | windowSize < 0) {
    logger::log_error("User-selected window size must be >= 1.")
    stop() # Should 'stop()' be replaced with 'return(NULL)' ?
  }
  bins <- GenomicRanges::tileGenome(
      sum(IRanges::runLength(coverageRaw)),
      tilewidth = windowSize,
      cut.last.tile.in.chrom = TRUE
    )
  cov <- GenomicRanges::binnedAverage(bins, coverageRaw, "coverage")
  cov <- as.data.frame(cov)[c("seqnames", "start", "end", "coverage")]
  colnames(cov) <- c("Chromosome", "chromStart", "chromEnd", "coverage")
  cov$coverage <- ceiling(as.numeric(cov$coverage))

  if (logScale) {
    cov$coverage <- log(cov$coverage)
  }
  cov$Chromosome <- ""

  return(cov)
}

getStatsFilePath <- function(sampleName,
                           plotSpecs) {
  # Step 1. Check ...
  if (plotSpecs$isOutput) {
    outDir <- dirname(plotSpecs$output)
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

filter_IR_genes <- function(regions, coverageRaw, seqnames, covData) {
  ir_genes <- GenomicRanges::GRanges(seqnames = seqnames, covData$ir_genes)
  ir_genes <- GenomicRanges::binnedAverage(ir_genes, coverageRaw, "coverage")
  ir_genes <- as.data.frame(ir_genes)[c("seqnames", "start", "end", "coverage")]
  regions_name <- covData$regions_name
  colnames(ir_genes) <- c(regions_name, covData$regions_start, covData$regions_end, "coverage")
  ir_genes$coverage <- ceiling(as.numeric(ir_genes$coverage))
  ir_genes[[regions_name]] <- sapply(covData$overlaps_genes, function(x) regions$Band[x])
  ir_genes[[regions_name]] <- gsub("^c\\(\"|\"|\"|\"\\)$", "", as.character(ir_genes[[regions_name]]))

  covData$ir_genes <- ir_genes
  return(covData)
}

filter_IR_noncoding <- function(regions, coverageRaw, seqnames, covData) {
  ir_noncoding <- GenomicRanges::GRanges(seqnames = seqnames, covData$ir_noncoding)
  ir_noncoding <- GenomicRanges::binnedAverage(ir_noncoding, coverageRaw, "coverage")
  ir_noncoding <- as.data.frame(ir_noncoding)[c("seqnames", "start", "end", "coverage")]
  regions_name <- covData$regions_name
  colnames(ir_noncoding) <- c(regions_name, covData$regions_start, covData$regions_end, "coverage")
  ir_noncoding$coverage <- ceiling(as.numeric(ir_noncoding$coverage))
  ir_noncoding[[regions_name]] <- sapply(covData$overlaps_noncoding, function(x) regions$Band[x])
  ir_noncoding[[regions_name]] <- gsub("^c\\(\"|\"|\"|\"\\)$", "", as.character(ir_noncoding[[regions_name]]))

  covData$ir_noncoding <- ir_noncoding
  return(covData)
}

filter_IR_regions <- function(coverageRaw, seqnames, covData) {
  ir_regions <- unlist(IRanges::slidingWindows(covData$ir_regions, width = 250L, step = 250L))
  ir_regions <- GenomicRanges::GRanges(seqnames = seqnames, ir_regions)
  ir_regions <- GenomicRanges::binnedAverage(ir_regions, coverageRaw, "coverage")
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
  aggFormula <- stats::as.formula(paste("coverage ~", regions_name))
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

# adapted from `nilsj9/PlastidSequenceCoverage`
getCovSummaries <- function(covData,
                            analysisSpecs) {
  regions_name <- covData$regions_name

  covData <- filterCovData(covData,
                           analysisSpecs)
  covSummaries <- getCovDepths(covData,
                               regions_name)
  covSummaries <- updateRegionsSummary(covSummaries,
                                       covData$ir_regions,
                                       regions_name)
  return(covSummaries)
}

updateRegionsSummary <- function(covSummaries,
                                 covDataRegions,
                                 regions_name) {
  covSumRegions <- covSummaries$regions_summary
  regions_evenness <- getCovEvenness(covDataRegions,
                                     regions_name)
  covSumRegions <- dplyr::full_join(covSumRegions,
                                    regions_evenness,
                                    regions_name)

  if (regions_name != "Source") {
    genome_summary <- getGenomeSummary(covDataRegions,
                                       regions_name)
    covSumRegions <- dplyr::bind_rows(covSumRegions,
                                      genome_summary)
  }

  covSummaries$regions_summary <- covSumRegions
  return(covSummaries)
}

filterCovData <- function(covData,
                          analysisSpecs) {
  regions_name <- covData$regions_name
  regions_start <- covData$regions_start
  regions_end <- covData$regions_end
  windowSize <- analysisSpecs$windowSize

  covData$ir_regions <- filterCovDataField(covData$ir_regions,
                                           regions_name,
                                           regions_start,
                                           regions_end,
                                           windowSize)
  covData$ir_genes <- filterCovDataField(covData$ir_genes,
                                         regions_name,
                                         regions_start,
                                         regions_end,
                                         windowSize)
  covData$ir_noncoding <- filterCovDataField(covData$ir_noncoding,
                                             regions_name,
                                             regions_start,
                                             regions_end,
                                             windowSize)
  return(covData)
}

filterCovDataField <- function(covDataField,
                               regions_name,
                               regions_start,
                               regions_end,
                               windowSize) {
  sizeThreshold <- windowSize - 1
  covDataField["length"] <- covDataField[regions_end] - covDataField[regions_start]
  covDataField <- covDataField[covDataField$length >= sizeThreshold, ]

  # update junction naming
  covDataField[[regions_name]] <- gsub("^(\\w+),\\s+", "Junction_\\1_", covDataField[[regions_name]])
  covDataField[[regions_name]] <- gsub("(\\w+),\\s+", "\\1_", covDataField[[regions_name]])
  return(covDataField)
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

getCovDepth <- function(covDataField, regions_name = NULL) {
  if (!is.null(regions_name)) {
    covDataField <- covDataField %>%
      groupByRegionsName(regions_name)
  }
  covDepth <- covDataField %>%
    calcCovDepth()
  return(covDepth)
}

calcCovDepth <- function(df) {
  lowCoverage <-
    lowCovWin_abs <-
    length <-
    regionLen <-
    NULL

  return(
    df %>%
    dplyr::summarise(
      lowCovWin_abs = sum(lowCoverage == "*", na.rm = TRUE),
      regionLen = sum(length, na.rm = TRUE),
      .groups = "drop") %>%
    addLowCovWin_relToRegionLen()
  )
}

addLowCovWin_relToRegionLen <- function(df) {
  lowCovWin_abs <-
    lowCovWin_relToRegionLen <-
    regionLen <-
    NULL

  return (
    df %>%
    dplyr::mutate(lowCovWin_relToRegionLen = lowCovWin_abs / regionLen)
  )
}

getCovEvenness <- function(covDataField, regions_name = NULL) {
  if (!is.null(regions_name)) {
    covDataField <- covDataField %>%
      groupByRegionsName(regions_name)
  }
  covEvenness <- covDataField %>%
    calcCovEvenness()
  return(covEvenness)
}

calcCovEvenness <- function(df) {
  coverage <- NULL

  return(
    df %>%
    dplyr::summarise(
      evenness = evennessScore(coverage),
      .groups = "drop"
    )
  )
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

getGenomeSummary <- function(covDataField, regions_name) {
  genome_depth <- getCovDepth(covDataField)
  genome_evenness <- getCovEvenness(covDataField)

  genome_summary <- genome_depth %>%
    dplyr::bind_cols(genome_evenness)
  genome_summary[regions_name] <- "Complete_genome"
  return(genome_summary)
}


#!/usr/bin/env RScript
#contributors=c("Gregory Smith", "Nils Jenke", "Michael Gruenstaeudl")
#email="m_gruenstaeudl@fhsu.edu"
#version="2024.01.07.2200"

read.gb2DF <- function(gbkData, regionsCheck) {
  fileDF <- data.frame()
  for (sample in gbkData) {
    sampleDF <- parseFeatures(sample$FEATURES, regionsCheck)
    if (!is.null(sampleDF)) {
      fileDF <- dplyr::bind_rows(fileDF, sampleDF)
    }
  }
  # signal no usable data
  if (ncol(fileDF) == 0 && nrow(fileDF) == 0) {
    return(NULL)
  }
  return(fileDF)
}

parseFeatures <- function(features, regionsCheck) {
  sampleDF <- data.frame()
  for (feature in features) {
    feature <- parseFeature(feature)
    if (!is.null(feature)) {
      sampleDF <- dplyr::bind_rows(sampleDF, feature)
    }
  }
  subsetCols <- checkFeatureQualifiers(sampleDF, regionsCheck)
  if (is.null(subsetCols)) {
    return(NULL)
  }
  type <- NULL
  source <- sampleDF %>%
              dplyr::filter(type=="source")
  sampleDF <- sampleDF %>% 
                dplyr::mutate(seqnames = as.factor(source[, "organism"])) %>%
                dplyr::select(dplyr::all_of(subsetCols))
  return(sampleDF)
}

parseFeature <- function(feature) {
  # first check if feature is desired
  locAndTypeIndex <- 1
  type <- feature[locAndTypeIndex,1]
  if (isIgnoredFeature(type)) {
    return(NULL)
  }
  
  # transform source data frame making sure final result
  # is still data frame
  feature <- as.data.frame(t(feature))
  colNames <- feature[1,]
  feature <- as.data.frame(feature[-1,])
  colnames(feature) <- colNames
  rownames(feature) <- NULL
  
  # fix sequence location(s) and feature type
  locationsStr <- feature[1,locAndTypeIndex]
  feature <- feature %>% 
              dplyr::rename_with(~ "locations", 
                                 .cols = dplyr::all_of(locAndTypeIndex)) %>%
              dplyr::mutate(type = type)
  
  # create derivative "start" and "end" variables from "locations"
  # and create feature copies for multiple sequences
  locationsList <- 
    regmatches(locationsStr, gregexpr("\\d+", locationsStr))[[1]] %>%
    as.integer()
  feature <- feature %>%
              dplyr::mutate(start = NA, end = NA)
  for (pairIndex in 1:(length(locationsList) / 2)) {
    if (pairIndex > 1) {
      feature <- dplyr::bind_rows(feature, feature[pairIndex-1,])
    }
    startIndex = pairIndex * 2 - 1
    feature[pairIndex, "start"] = locationsList[startIndex]
    feature[pairIndex, "end"]   = locationsList[startIndex+1]
  }
  return(feature)
}

read.gbSeq <- function(gbkData) {
  sampleSequences <- c()
  for (sample in gbkData) {
    sampleSequences <- c(sampleSequences, sample$ORIGIN)
  }
  return(Biostrings::DNAStringSet(sampleSequences))
}

read.gbGenes <- function(gbkDataDF) {
  type <- NULL
  subsetCols <- c("seqnames", "start", "end", "gene")
  gene_L <- gbkDataDF %>%
              dplyr::filter(type=="gene") %>% 
              dplyr::select(dplyr::all_of(subsetCols))
  rownames(gene_L) <- NULL
  return(gene_L)
}

read.gbOther <- function(gbkDataDF) {
  type <- NULL
  subsetCols <- c("seqnames", "start", "end", 
                  "gene", "note", "standard_name")
  regions <- gbkDataDF %>%
              dplyr::filter(!type %in% c("gene", "exon", "transcript",
                                  "CDS", "variant")) %>% 
              dplyr::select(dplyr::all_of(subsetCols))
  rownames(regions) <- NULL
  return(regions)
}

read.gbLengths <- function(gbkData) {
  sampleLengths <- c()
  for (sample in gbkData) {
    sampleLengths <- c(sampleLengths, nchar(sample$ORIGIN))
  }
  return(sampleLengths)
}

read.gbSampleName <- function(gbkData) {
  sampleNames <- c()
  for (sample in gbkData) {
    sampleNames <- c(sampleNames, 
                   c(sample_name = sample$ACCESSION, 
                     genome_name = sample$VERSION))
  }
  return(sampleNames)
}

read.gbPlotTitle <- function(gbkData) {
  plotTitles <- c()
  for (sample in gbkData) {
    plotTitles <- c(plotTitles, paste(sample$DEFINITION, sample$ACCESSION))
  }
  return(plotTitles)
}

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
  if (relative == TRUE & logScale) {
    threshold <- mean(cov[, 4]) + log(threshold)
  } else if (relative == TRUE) {
    threshold <- mean(cov[, 4]) * threshold
  }
  color <- rep("black", nrow(cov))
  ind <- as.numeric(cov[, 4]) <= threshold
  color <- replace(color, ind, "red")
  return(color)
}
    
boolToDeci <- function(boolList) {
  out = 0
  boolList <- rev(boolList)
  for (i in 1:length(boolList)) {
    out = out + boolList[i] * (2 ^ (i - 1))
  }
  return(out)
}

writeTables <- function(regions, bamFile, genes, dir, sample_name) {
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

	# Operations on genes
    ir_genes <- GenomicRanges::GRanges(seqnames = sample_name["genome_name"], ir_genes)
    ir_genes <- GenomicRanges::binnedAverage(ir_genes, GenomicAlignments::coverage(bamFile), "coverage")
    ir_genes <- as.data.frame(ir_genes)[c("seqnames", "start", "end", "coverage")]
    colnames(ir_genes) <- c("Chromosome", "chromStart", "chromEnd", "coverage")
    ir_genes$coverage <- ceiling(as.numeric(ir_genes$coverage))
    ir_genes$Chromosome <- sapply(overlaps_genes, function(x) regions$Band[x])
    ir_genes$Chromosome <- gsub("^c\\(\"|\"|\"|\"\\)$", "", as.character(ir_genes$Chromosome))
   
    # Operations on non-coding regions
    ir_noncoding <- GenomicRanges::GRanges(seqnames = sample_name["genome_name"], ir_noncoding)
    ir_noncoding <- GenomicRanges::binnedAverage(ir_noncoding, GenomicAlignments::coverage(bamFile), "coverage")
    ir_noncoding <- as.data.frame(ir_noncoding)[c("seqnames", "start", "end", "coverage")]
    colnames(ir_noncoding) <- c("Chromosome", "chromStart", "chromEnd", "coverage")
    ir_noncoding$coverage <- ceiling(as.numeric(ir_noncoding$coverage))
    ir_noncoding$Chromosome <- sapply(overlaps_noncoding, function(x) regions$Band[x])
    ir_noncoding$Chromosome <- gsub("^c\\(\"|\"|\"|\"\\)$", "", as.character(ir_noncoding$Chromosome))

    # Operations on IRs
    ir_regions <- unlist(IRanges::slidingWindows(ir_regions, width = 250L, step = 250L))
    ir_regions <- GenomicRanges::GRanges(seqnames = sample_name["genome_name"], ir_regions)
    ir_regions <- GenomicRanges::binnedAverage(ir_regions, GenomicAlignments::coverage(bamFile), "coverage")
    chr <- ir_regions@ranges@NAMES
    ir_regions <- as.data.frame(ir_regions, row.names = NULL)[c("seqnames", "start", "end", "coverage")]
    ir_regions["seqnames"] <- chr
    colnames(ir_regions) <- c("Chromosome", "chromStart", "chromEnd", "coverage")
    ir_regions$coverage <- ceiling(as.numeric(ir_regions$coverage))
    
    cov_regions <-
      aggregate(
        coverage ~ Chromosome,
        data = ir_regions,
        FUN = function(x)
          ceiling(mean(x) - sd(x))
      )
    ir_regions$lowCoverage <- with(ir_regions, ir_regions$coverage < cov_regions$coverage[match(Chromosome, cov_regions$Chromosome)])
    ir_genes$lowCoverage <- ir_genes$coverage < mean(ir_genes$coverage) - sd(ir_genes$coverage)
    ir_noncoding$lowCoverage <- ir_noncoding$coverage < mean(ir_noncoding$coverage) - sd(ir_noncoding$coverage)
    
    ir_genes$lowCoverage[ir_genes$lowCoverage == TRUE] <- "*"
    ir_genes$lowCoverage[ir_genes$lowCoverage == FALSE] <- ""
    ir_regions$lowCoverage[ir_regions$lowCoverage == TRUE] <- "*"
    ir_regions$lowCoverage[ir_regions$lowCoverage == FALSE] <- ""
    ir_noncoding$lowCoverage[ir_noncoding$lowCoverage == TRUE] <- "*"
    ir_noncoding$lowCoverage[ir_noncoding$lowCoverage == FALSE] <- ""
    
    # Writing values to output table
    write.table(
      ir_genes,
      paste(dir, .Platform$file.sep, sample_name["sample_name"], "_coverage.genes.bed", sep = ""),
      row.names = FALSE,
      quote = FALSE,
      sep = "\t"
    )
    write.table(
      ir_regions,
      paste(dir, .Platform$file.sep, sample_name["sample_name"], "_coverage.regions.bed", sep = ""),
      row.names = FALSE,
      quote = FALSE,
      sep = "\t"
    )
    write.table(
      ir_noncoding,
      paste(dir, .Platform$file.sep, sample_name["sample_name"], "_coverage.noncoding.bed", sep = ""),
      row.names = FALSE,
      quote = FALSE,
      sep = "\t"
    )
 }

validateColors <- function(colorsToValidate) {
  colorNames <- colors()
  unsupportedColors <- colorsToValidate[!(colorsToValidate %in% colorNames)]
  if (length(unsupportedColors) > 0) {
    stop("Unsupported R plot color defined.")
  }
}

checkFeatureQualifiers <- function(sampleDF, regionsCheck) {
  subsetCols <- c("start", "end", "gene", 
                  "note", "type")
  if (regionsCheck) {
    subsetCols <- c(subsetCols, "standard_name")
  }
  missingCols <- subsetCols[!(subsetCols %in% colnames(sampleDF))]
  if (length(missingCols) > 0) {
    logger::log_warn(paste0("Unable to analyze sample as specified; ", 
                            "missing feature qualifiers: ",
                            "'",
                            paste(missingCols, collapse = "', '"),
                            "'"))
    return(NULL)
  }
  subsetCols <- c(subsetCols, "seqnames")
  return(subsetCols)
}

isIgnoredFeature <- function(featureName) {
  ignoredFeatures <- c("D-loop")
  return(featureName %in% ignoredFeatures)
}

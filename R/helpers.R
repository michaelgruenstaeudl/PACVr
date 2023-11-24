#!/usr/bin/R
#contributors=c("Gregory Smith", "Nils Jenke", "Michael Gruenstaeudl")
#email="m_gruenstaeudl@fhsu.edu"
#version="2023.11.21.2100"

`%>%` <- dplyr::`%>%

read.gb2DF <- function(gbkData) {
  fileDF <- data.frame()
  for (sample in gbkData) {
    sampleDF <- parseFeatures(sample$FEATURES)
    fileDF <- dplyr::bind_rows(fileDF, sampleDF)
  }
  return(fileDF)
}

parseFeatures <- function(features) {
  # first feature is sample information
  source <- parseFeature(features[[1]])
  features <- features[-1]

  sampleDF <- data.frame()
  for (feature in features) {
    feature <- parseFeature(feature)
    sampleDF <- dplyr::bind_rows(sampleDF, feature)
  }
  sampleDF <- sampleDF %>% 
                dplyr::mutate(seqnames = as.factor(source[, "organism"]))
  return(sampleDF)
}

parseFeature <- function(feature) {
  # transform source data frame
  feature <- as.data.frame(t(feature))
  colnames(feature) <- feature[1,]
  feature <- feature[-1,]
  rownames(feature) <- NULL
  
  # fix sequence location(s) and feature type
  locAndTypeIndex <- 1
  locationsStr <- feature[1,locAndTypeIndex]
  type <- names(feature)[locAndTypeIndex]
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

read.gbGenes <- function(gbkData) {
  bgkDataDF <- read.gb2DF(gbkData)
  gene_L <- bgkDataDF %>% 
              dplyr::filter(type=="gene") %>% 
              dplyr::select(all_of(c("seqnames", "start", "end", "gene")))
  rownames(gene_L) <- NULL
  return(gene_L)
}

read.gbOther <- function(gbkData) {
  bgkDataDF <- read.gb2DF(gbkData)
  regions <- bgkDataDF %>% 
              dplyr::filter(!type %in% c("gene", "exon", "transcript",
                                  "CDS", "variant")) %>% 
              dplyr::select(all_of(c("seqnames", "start", "end", 
                              "gene", "note", "standard_name")))
  rownames(regions) <- NULL
  return(regions)
}

read.gbLengths <- function(gbkData) {
  sampleLengths <- c()
  for (sample in gbkData) {
    sampleLengths <- c(sampleLengths, 
                       nchar(sample$ORIGIN))
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
    plotTitles <- c(plotTitles, 
                    paste(sample$DEFINITION, 
                          sample$ACCESSION))
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
  ind   <- as.numeric(cov[, 4]) <= threshold
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

writeTables <-
  function(regions,
           bamFile,
           genes,
           dir,
           sample_name) {
    ir_regions <-
      IRanges::IRanges(
        start = regions$chromStart,
        end = regions$chromEnd,
        names = regions$Band
      )
    ir_genes <-
      unlist(IRanges::slidingWindows(
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
    ir_noncoding <-
      unlist(IRanges::slidingWindows(
        IRanges::reduce(IRanges::gaps(
          ir_genes,
          start = min(regions$chromStart),
          end = max(regions$chromEnd)
        )),
        width = 250L,
        step = 250L
      ))
    overlaps_genes <-
      IRanges::findOverlaps(
        query = ir_genes,
        subject = ir_regions,
        type = "any",
        select = "all"
      )

    ir_genes <-
      GenomicRanges::GRanges(seqnames = sample_name["genome_name"], ir_genes)
    ir_genes <-
      GenomicRanges::binnedAverage(ir_genes,
                                   GenomicAlignments::coverage(bamFile),
                                   "coverage")
    ir_genes <-
      as.data.frame(ir_genes)[c("seqnames", "start", "end", "coverage")]
    colnames(ir_genes) <-
      c("Chromosome", "chromStart", "chromEnd", "coverage")
    ir_genes$coverage <- ceiling(as.numeric(ir_genes$coverage))
    ir_genes$Chromosome <-
      sapply(overlaps_genes, function(x)
        regions$Band[x])
    ir_genes$Chromosome <-
      gsub("^c\\(\"|\"|\"|\"\\)$",
           "",
           as.character(ir_genes$Chromosome))
    
    overlaps_noncoding <-
      IRanges::findOverlaps(
        query = ir_noncoding,
        subject = ir_regions,
        type = "any",
        select = "all"
      )
    
    ir_noncoding <-
      GenomicRanges::GRanges(seqnames = sample_name["genome_name"], ir_noncoding)
    ir_noncoding <-
      GenomicRanges::binnedAverage(ir_noncoding,
                                   GenomicAlignments::coverage(bamFile),
                                   "coverage")
    ir_noncoding <-
      as.data.frame(ir_noncoding)[c("seqnames", "start", "end", "coverage")]
    colnames(ir_noncoding) <-
      c("Chromosome", "chromStart", "chromEnd", "coverage")
    ir_noncoding$coverage <-
      ceiling(as.numeric(ir_noncoding$coverage))
    ir_noncoding$Chromosome <-
      sapply(overlaps_noncoding, function(x)
        regions$Band[x])
    ir_noncoding$Chromosome <-
      gsub("^c\\(\"|\"|\"|\"\\)$",
           "",
           as.character(ir_noncoding$Chromosome))
    
    ir_regions <-
      unlist(IRanges::slidingWindows(ir_regions, width = 250L, step = 250L))
    ir_regions <-
      GenomicRanges::GRanges(seqnames = sample_name["genome_name"], ir_regions)
    ir_regions <-
      GenomicRanges::binnedAverage(ir_regions,
                                   GenomicAlignments::coverage(bamFile),
                                   "coverage")
    chr <- ir_regions@ranges@NAMES
    ir_regions <-
      as.data.frame(ir_regions, row.names = NULL)[c("seqnames", "start", "end", "coverage")]
    ir_regions["seqnames"] <- chr
    colnames(ir_regions) <-
      c("Chromosome", "chromStart", "chromEnd", "coverage")
    ir_regions$coverage <- ceiling(as.numeric(ir_regions$coverage))
    
    cov_regions <-
      aggregate(
        coverage ~ Chromosome,
        data = ir_regions,
        FUN = function(x)
          ceiling(mean(x) - sd(x))
      )
    ir_regions$lowCoverage <-
      with(ir_regions, ir_regions$coverage < cov_regions$coverage[match(Chromosome, cov_regions$Chromosome)])
    
    ir_genes$lowCoverage <-
      ir_genes$coverage < mean(ir_genes$coverage) - sd(ir_genes$coverage)
    ir_noncoding$lowCoverage <-
      ir_noncoding$coverage < mean(ir_noncoding$coverage) - sd(ir_noncoding$coverage)
    
    ir_genes$lowCoverage[ir_genes$lowCoverage == TRUE] <- "*"
    ir_genes$lowCoverage[ir_genes$lowCoverage == FALSE] <- ""
    ir_regions$lowCoverage[ir_regions$lowCoverage == TRUE] <- "*"
    ir_regions$lowCoverage[ir_regions$lowCoverage == FALSE] <- ""
    ir_noncoding$lowCoverage[ir_noncoding$lowCoverage == TRUE] <- "*"
    ir_noncoding$lowCoverage[ir_noncoding$lowCoverage == FALSE] <- ""
    
    write.table(
      ir_genes,
      paste(
        dir,
        .Platform$file.sep,
        sample_name["sample_name"],
        "_coverage.genes.bed",
        sep = ""
      ),
      row.names = FALSE,
      quote = FALSE,
      sep = "\t"
    )
    write.table(
      ir_regions,
      paste(
        dir,
        .Platform$file.sep,
        sample_name["sample_name"],
        "_coverage.regions.bed",
        sep = ""
      ),
      row.names = FALSE,
      quote = FALSE,
      sep = "\t"
    )
    write.table(
      ir_noncoding,
      paste(
        dir,
        .Platform$file.sep,
        sample_name["sample_name"],
        "_coverage.noncoding.bed",
        sep = ""
      ),
      row.names = FALSE,
      quote = FALSE,
      sep = "\t"
    )
  }

checkIREquality <- function(gbkData, regions, dir, sample_name) {
  gbkSeq <- read.gbSeq(gbkData)
  if ("IRb" %in% regions[, 4] && "IRa" %in% regions[, 4]) {
    repeatB <- as.numeric(regions[which(regions[, 4] == "IRb"), 2:3])
    repeatA <-
      as.numeric(regions[which(regions[, 4] == "IRa"), 2:3])
    IR_diff_SNPS <- c()
    IR_diff_gaps <- c()
    if (repeatB[2] - repeatB[1] != repeatA[2] - repeatA[1]) {
      message("WARNING: Inverted repeats differ in sequence length")
      message(paste(
        "The IRb has a total lengths of: ",
        repeatB[2] - repeatB[1],
        " bp",
        sep = ""
      ))
      message(paste(
        "The IRa has a total lengths of: ",
        repeatA[2] - repeatA[1],
        " bp",
        sep = ""
      ))
    }
    if (gbkSeq[[1]][repeatB[1]:repeatB[2]] != Biostrings::reverseComplement(gbkSeq[[1]][repeatA[1]:repeatA[2]])) {
      IRa_seq <- Biostrings::DNAString(gbkSeq[[1]][repeatB[1]:repeatB[2]])
      IRa_seq <- split(IRa_seq, ceiling(seq_along(IRa_seq) / 10000))
      IRb_seq <- Biostrings::DNAString(Biostrings::reverseComplement(gbkSeq[[1]][repeatA[1]:repeatA[2]]))
      IRb_seq <- split(IRb_seq, ceiling(seq_along(IRb_seq) / 10000))
      
      for (i in  1:min(length(IRa_seq), length(IRb_seq))) {
        subst_mat <-
          Biostrings::nucleotideSubstitutionMatrix(match = 1,
                                                   mismatch = -3,
                                                   baseOnly = TRUE)
        globalAlign <- tryCatch({
          Biostrings::pairwiseAlignment(
            IRa_seq[[i]],
            IRb_seq[[i]],
            substitutionMatrix = subst_mat,
            gapOpening = 5,
            gapExtension = 2
          )
        },
        error = function(e) {
          return(NULL)
        })
        if (is.null(globalAlign))
          break
        IR_diff_SNPS <-
          c(IR_diff_SNPS, which(strsplit(
            Biostrings::compareStrings(globalAlign), ""
          )[[1]] == "?"))
        IR_diff_gaps <-
          c(IR_diff_gaps, which(strsplit(
            Biostrings::compareStrings(globalAlign), ""
          )[[1]] == "-"))
      }
      message("WARNING: Inverted repeats differ in sequence")
      if (length(IR_diff_SNPS) > 0) {
        message(
          paste(
            "When aligned, the IRs differ through a total of ",
            length(IR_diff_SNPS),
            " SNPS. These SNPS are located at the following nucleotide positions: ",
            paste(unlist(IR_diff_SNPS), collapse = " "),
            sep = ""
          )
        )
      }
      if (length(IR_diff_gaps) > 0) {
        message(
          paste(
            "When aligned, the IRs differ through a total of ",
            length(IR_diff_gaps),
            " gaps. These gaps are located at the following nucleotide positions: ",
            paste(unlist(IR_diff_gaps), collapse = " "),
            sep = ""
          )
        )
      }
      message(
        "Proceeding with coverage depth visualization, but without quadripartite genome structure ..."
      )
    }
    write.csv(
      data.frame(
        Number_N = unname(Biostrings::alphabetFrequency(gbkSeq)[, "N"]),
        Mismatches = length(IR_diff_SNPS) + length(IR_diff_gaps)
      ),
      paste0(dir, .Platform$file.sep, sample_name["sample_name"], "_IR_quality.csv"),
      row.names = FALSE,
      quote = FALSE
    )
  } else {
    write.csv(
      data.frame(
        Number_N = unname(Biostrings::alphabetFrequency(gbkSeq)[, "N"]),
        Mismatches = NA
      ),
      paste0(dir, .Platform$file.sep, sample_name["sample_name"], "_IR_quality.csv"),
      row.names = FALSE,
      quote = FALSE
    )
  }
}

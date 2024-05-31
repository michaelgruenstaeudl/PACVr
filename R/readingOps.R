#!/usr/bin/env RScript
#contributors=c("Gregory Smith", "Nils Jenke", "Michael Gruenstaeudl")
#email="m_gruenstaeudl@fhsu.edu"
#version="2024.05.31.0457"

read.gbSeqFeaturesAdapt <- function(gbkData, analysisSpecs) {
  gbkSeqFeatures <- read.gbSeqFeatures(gbkData,
                                       analysisSpecs)
  if (is.null(gbkSeqFeatures)) {
    logger::log_warn(paste("No samples match specifications;",
                           "retrying with basic coverage requirements."))
    analysisSpecs$setIRCheckFields(NA)
    gbkSeqFeatures <- read.gbSeqFeatures(gbkData,
                                         analysisSpecs)
  }
  return(gbkSeqFeatures)
}

read.gbSeqFeatures <- function(gbkData, analysisSpecs) {
  fileDF <- data.frame()
  for (sample in gbkData) {
    sampleDF <- parseFeatures(sample$FEATURES, analysisSpecs)
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

parseFeatures <- function(features, analysisSpecs) {
  sampleDF <- data.frame()
  for (feature in features) {
    feature <- parseFeature(feature)
    if (!is.null(feature)) {
      sampleDF <- dplyr::bind_rows(sampleDF, feature)
    }
  }
  # check if can we can use the sample
  subsetCols <- checkFeatureQualifiers(sampleDF, analysisSpecs)
  if (is.null(subsetCols)) {
    return(NULL)
  }
  rownames(sampleDF) <- NULL

  # create derivative "start" and "end" variables from "locations",
  # create "seqname", and subset according to analysis needs
  . <- NULL
  sampleDF <- sampleDF %>%
    addStartEnd() %>%
    addSeqname() %>%
    subset4Analysis(., subsetCols)
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
  colNames <- sub("^/", "", feature[1,])
  feature <- as.data.frame(feature[-1,])
  colnames(feature) <- colNames
  colNames <- names(feature)

  # fix sequence location(s) and feature type
  . <- NULL
  feature <- feature %>%
    stats::setNames(make.unique(names(.))) %>%
    dplyr::rename(locations = {{ type }}) %>%
    dplyr::mutate(type = type)

  # combine duplicate qualifiers of same name
  feature <- combineDupQuals(feature, colNames)

  return(feature)
}

combineDupQuals <- function(feature, colNames) {
  dupQualifiers <- colNames[duplicated(colNames)]
  if (length(dupQualifiers) == 0) {
    return(feature)
  }

  # combine into single string the columns that had the same name before
  # setNames() was applied
  for (qualifier in dupQualifiers) {
    uniteMatch <- paste0("^", qualifier, "$", 
                         "|",
                         "^", qualifier, "\\.\\d+$")
    uniteCols <- grep(uniteMatch, names(feature), value = TRUE)
    feature <- feature %>%
      tidyr::unite("uniteCol", dplyr::all_of(uniteCols), sep = "; ") %>%
      dplyr::select(-dplyr::matches(uniteMatch)) %>%
      dplyr::rename_with(~qualifier, dplyr::matches("uniteCol"))
  }
  return(feature)
}

checkFeatureQualifiers <- function(sampleDF, analysisSpecs) {
  subsetData <- getSubsetData(sampleDF, analysisSpecs)
  if (length(subsetData$missingCols) > 0) {
    logger::log_warn(paste0("Unable to analyze sample as specified; ",
                            "missing feature qualifiers: ",
                            "'",
                            paste(subsetData$missingCols, collapse = "', '"),
                            "'"))
    return(NULL)
  }
  # add future generated columns that are needed
  subsetCols <- c(subsetData$subsetCols, "start", "end", "seqnames")
  return(subsetCols)
}

addStartEnd <- function(sampleDF) {
  start <-
    end <-
      feature_index <-
        . <-
          NULL

  # keep track of feature ordering in case of multiple sequences
  # for features
  sampleDF <- sampleDF %>%
    dplyr::mutate(feature_index = as.integer(rownames(.)))

  locationsList <- getLocationsList(sampleDF)
  for (index in 1:length(locationsList)) {
    featureLocations <- locationsList[[index]]
    locationsCount <- length(featureLocations)
    if (locationsCount == 0) {
      next
    }

    sampleDF[index, c("start", "end")] <- c(featureLocations[[1]][[1]],
                                            featureLocations[[1]][[2]])
    # create feature copies for multiple sequences
    if (locationsCount > 1) {
      sampleDF <- addMultiSeqFeats(sampleDF, featureLocations, index)
    }
  }

  sampleDF <- sampleDF %>%
    dplyr::mutate(start = as.integer(start),
                  end = as.integer(end)) %>%
    dplyr::arrange(feature_index) %>%
    dplyr::select(-feature_index)
  return(sampleDF)
}

addMultiSeqFeats <- function(sampleDF, featureLocations, index) {
  . <- NULL

  featureLocations <- featureLocations[-1]
  featureCopies <- sampleDF %>%
    dplyr::slice(index) %>%
    replicate(length(featureLocations),
              ., simplify=FALSE) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(start = featureLocations[[1]][[1]],
                  end = featureLocations[[1]][[2]])
  sampleDF <- sampleDF %>%
    dplyr::bind_rows(., featureCopies)
  return(sampleDF)
}

addSeqname <- function(sampleDF) {
  type <- NULL
  source <- sampleDF %>%
    dplyr::filter(type=="source")
  sampleDF <- sampleDF %>%
    dplyr::mutate(seqnames = as.factor(source[, "organism"]))
  return(sampleDF)
}

subset4Analysis <- function(sampleDF, subsetCols) {
  sampleDF <- sampleDF %>%
    dplyr::select(dplyr::all_of(subsetCols))
  return(sampleDF)
}

getLocationsList <- function(sampleDF, removeRemotes=TRUE) {
  locationsList <- sampleDF %>%
    dplyr::pull("locations") %>%
    lapply(normalizeLocation) %>%
    lapply(function(x) unlist(strsplit(x, ",")))
  if (removeRemotes) {
    locationsList <- lapply(locationsList,
                            function(x) x[!grepl(":", x)])
  }
  locationsList <- lapply(locationsList,
                          function(x) strsplit(x, ":|\\.{2,}"))
  return(locationsList)
}

# based on standard listed by INSDC
# https://www.insdc.org/submitting-standards/feature-table/
normalizeLocation <- function(locationsStr) {
  # standardize delimiter
  cleanStr <- gsub("(?![A-Za-z]+)(\\d+)\\.(\\d+)(?!:)",
                   "\\1..\\2", locationsStr, perl=TRUE)
  cleanStr <- gsub("(?![A-Za-z]+)(\\d+)\\^(\\d+)(?!:)",
                   "\\1..\\2", cleanStr, perl=TRUE)
  # remove unknown boundary identifiers
  cleanStr <- gsub("[<>](\\d+)",
                   "\\1", cleanStr)
  # remove operators
  cleanStr <- gsub("complement|join|order|[\\(\\)]",
                   "", cleanStr)
  # standardize single bases
  cleanStr <- gsub("(^|,)(.+:)?(\\d+)($|,)",
                   "\\1\\2\\3..\\3\\4", cleanStr, perl=TRUE)
  return(cleanStr)
}

getSubsetCols <- function(analysisSpecs) {
  subsetCols <- c("gene",
                  "type")
  if (analysisSpecs$isIRCheck) {
    subsetCols <- c(subsetCols,
                    "note",
                    "standard_name")
  }
  return(subsetCols)
}

getSubsetData <- function(sampleDF, analysisSpecs) {
  subsetCols <- getSubsetCols(analysisSpecs)
  missingCols <- subsetCols[!(subsetCols %in% colnames(sampleDF))]
  if (analysisSpecs$isIRCheck && ("standard_name" %in% missingCols)) {
    logger::log_info("Using `note` for IR name qualifier")
    subsetCols <- subsetCols[subsetCols != "standard_name"]
    missingCols <- missingCols[missingCols != "standard_name"]
  }
  subsetData <- list(
    subsetCols = subsetCols,
    missingCols = missingCols
  )
  return(subsetData)
}

isIgnoredFeature <- function(featureName) {
  ignoredFeatures <- c("D-loop")
  return(featureName %in% ignoredFeatures)
}

read.gbSequence <- function(gbkData) {
  sampleSequences <- c()
  for (sample in gbkData) {
    sampleSequences <- c(sampleSequences, sample$ORIGIN)
  }
  return(Biostrings::DNAStringSet(sampleSequences))
}

read.gbGenes <- function(gbkSeqFeatures) {
  type <- NULL
  subsetCols <- c("seqnames", "start", "end", "gene")
  gene_L <- gbkSeqFeatures %>%
    dplyr::filter(type=="gene") %>%
    dplyr::select(dplyr::all_of(subsetCols))
  rownames(gene_L) <- NULL
  return(gene_L)
}

read.gbOther <- function(gbkSeqFeatures) {
  type <- NULL
  regions <- gbkSeqFeatures %>%
    dplyr::filter(!type %in% c("gene", "exon", "transcript",
                               "CDS", "variant"))
  rownames(regions) <- NULL
  return(regions)
}

read.gbLengths <- function(gbSeqs) {
  sampleLengths <- Biostrings::width(gbSeqs)
  return(sampleLengths)
}

read.gbSampleName <- function(gbkData) {
  sampleNames <- c()
  for (sample in gbkData) {
    sampleNames <- c(sampleNames,
                     c(sample_name = sample$VERSION,
                       genome_name = sample$ACCESSION))
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

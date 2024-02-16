#!/usr/bin/env RScript
#contributors=c("Gregory Smith", "Nils Jenke", "Michael Gruenstaeudl")
#email="m_gruenstaeudl@fhsu.edu"
#version="2024.02.01.1736"

FilterByKeywords <- function(allRegions, where) {
  # Function to filter list based on genomic keywords
  # ARGS:
  #   ...
  # RETURNS:
  #   ...
  out = subset(
    allRegions,
    grepl(
      "^IR|repeat|invert|^LSC|^SSC|[large,long]\\ssingle\\scopy|[short,small]\\ssingle\\scopy",
      allRegions[, where],
      ignore.case=TRUE
    )
  )
  if (nrow(out) < 1) {
    warning(paste("Inverted repeat info not found for qualifier ", where, ".", sep = ""))
  }
  return(out)
}

ParseQuadripartiteStructure <- function(gbkDataDF) {
  # Function to extract the quadripartite region information from 
  # Genbank flatfile data
  # ARGS:
  #   gbkDataDF (resulting data frame from parsing read.gb object)
  # RETURNS:
  #   regions in data frame format
  logger::log_info('  Extracting information on genomic regions')
  allRegions <- read.gbOther(gbkDataDF)
  colNames <- colnames(allRegions)
  if ("standard_name" %in% colNames) {
    filterWhere <- "standard_name"
  } else if ("note" %in% colNames) {
    filterWhere <- "note"
  }
  quadripRegions <- FilterByKeywords(allRegions, filterWhere)
  quadripRegions <- quadripRegions[, c("start", "end", "note")]
  colnames(quadripRegions) <- c("chromStart", "chromEnd", "Band")

  quadripRegions$Chromosome <- ""
  quadripRegions$Stain <- "gpos100"
  quadripRegions <- quadripRegions[c("Chromosome", "chromStart", "chromEnd", "Band", "Stain")]
  quadripRegions <- quadripRegions[order(quadripRegions[, 3], decreasing=FALSE),]
  quadripRegions$Band[which(grepl("LSC|large|long", quadripRegions$Band, ignore.case=TRUE) == TRUE)] <- "LSC"
  quadripRegions$Band[which(grepl("SSC|small|short", quadripRegions$Band, ignore.case=TRUE) == TRUE)] <- "SSC"
  quadripRegions$Band[which(grepl("IRa|\\sa|IR1|\\s1", quadripRegions$Band, ignore.case=TRUE) == TRUE)] <- "###A"
  quadripRegions$Band[which(grepl("IRb|\\sb|IR2|\\s2", quadripRegions$Band, ignore.case=TRUE) == TRUE)] <- "###B"
  quadripRegions$Band[which(grepl("IR|invert|repeat", quadripRegions$Band, ignore.case=TRUE) == TRUE)] <- "IR"
  quadripRegions$Band[which(grepl("###A", quadripRegions$Band, ignore.case=TRUE) == TRUE)] <- "IRa"
  quadripRegions$Band[which(grepl("###B", quadripRegions$Band, ignore.case=TRUE) == TRUE)] <- "IRb"
  row.names(quadripRegions) <- 1:nrow(quadripRegions)
  quadripRegions <- quadripRegions[order(quadripRegions[, 3], decreasing=FALSE),]
  return(quadripRegions)
}

fillDataFrame <- function(gbkLengths, quadripRegions) {
  # Function to annotate plastid genome with quadripartite regions 
  # based on their position within the genome
  # ARGS:
  #   gbkData (i.e., GenBank flatfile data as parsed by read.gb())
  # RETURNS:
  #   ...
  logger::log_info('  Annotating plastid genome with quadripartite regions')
  seqLength <- gbkLengths
  if ((nrow(quadripRegions) == 0) || (quadripRegions[1, 2] == -1)) {
    quadripRegions[1,] <- c("", as.numeric(1), as.numeric(seqLength), "NA", "gpos100")
    quadripRegions[, 2] <- as.numeric(quadripRegions[, 2])
    quadripRegions[, 3] <- as.numeric(quadripRegions[, 3])
    return(quadripRegions)
  } else {
    start <- 1
    for (i in 1:nrow(quadripRegions)) {
      if (quadripRegions[i, 2] > start) {
        quadripRegions[nrow(quadripRegions) + 1,] <- c("", start, as.numeric(quadripRegions[i, 2]) - 1, "NA", "gpos100")
      }
      start <- as.numeric(quadripRegions[i, 3]) + 1
    }
    if (start - 1 < seqLength) {
      quadripRegions[nrow(quadripRegions) + 1,] <- c("", start, seqLength, "NA", "gpos100")
    }
    quadripRegions <- quadripRegions[order(as.numeric(quadripRegions[, 2]), decreasing=FALSE),]
    row.names(quadripRegions) <- 1:nrow(quadripRegions)
    quadripRegions[, 2] <- as.numeric(quadripRegions[, 2])
    quadripRegions[, 3] <- as.numeric(quadripRegions[, 3])
    
    regionAvail <- boolToDeci(c("LSC", "IRb", "SSC", "IRa") %in% quadripRegions[, 4])
    quadripRegions[, 6] <- quadripRegions[, 3] - quadripRegions[, 2]
    
    if (regionAvail == 5) {
      # only IRa and IRb
      quadripRegions[which(quadripRegions[, 4] != "NA"), 6] <- 0
      quadripRegions[which(quadripRegions[, 6] == max(quadripRegions[, 6])), 4] <- "LSC"
      quadripRegions[which(quadripRegions[, 4] != "NA"), 6] <- 0
      quadripRegions[which(quadripRegions[, 6] == max(quadripRegions[, 6])), 4] <- "SSC"
      message("Annotation for LSC and SSC were automatically added")
    } else if (regionAvail == 7) {
      # only IRa, SSC and IRb
      quadripRegions[which(quadripRegions[, 4] != "NA"), 6] <- 0
      quadripRegions[which(quadripRegions[, 6] == max(quadripRegions[, 6])), 4] <- "LSC"
      message("Annotation for LSC was automatically added")
    } else if (regionAvail == 10) {
      # only LSC and SSC
      IRs <- data.frame(table(quadripRegions[which(quadripRegions[, 4] == "NA"), 6]), stringsAsFactors=FALSE)
      IRs <- IRs[IRs$Freq == 2, 1]
      if (length(IRs) >= 1) {
        IRs <- max(as.numeric(as.character(IRs)))
        quadripRegions[which(quadripRegions[, 6] == IRs), 4] <- c("IRb", "IRa")
        message("Annotation for IRb and IRa were automatically added")
      }
    } else if (regionAvail == 11) {
      # only LSC, SSC and IRa
      quadripRegions[which(quadripRegions[, 4] == "NA" & quadripRegions[, 6] == quadripRegions[which(quadripRegions[, 4] == "IRa"), 6]), 4] <- "IRb"
      message("Annotation for IRb was automatically added")
    } else if (regionAvail == 13) {
      # only LSC, IRb and IRa
      quadripRegions[which(quadripRegions[, 4] != "NA"), 6] <- 0
      quadripRegions[which(quadripRegions[, 6] == max(quadripRegions[, 6])), 4] <- "SSC"
      message("Annotation for SSC was automatically added")
    } else if (regionAvail == 14) {
      # only LSC, IRb and SSC
      quadripRegions[which(quadripRegions[, 4] == "NA" & quadripRegions[, 6] == quadripRegions[which(quadripRegions[, 4] == "IRb"), 6]), 4] <- "IRa"
      message("Annotation for IRa was automatically added")
    }
    quadripRegions <- quadripRegions[-6]
    quadripRegions$Stain[which(quadripRegions$Band == "LSC")] <- "gpos75"
    quadripRegions$Stain[which(quadripRegions$Band == "SSC")] <- "gpos50"
    quadripRegions$Stain[which(quadripRegions$Band == "IRa")] <- "gpos25"
    quadripRegions$Stain[which(quadripRegions$Band == "IRb")] <- "gpos25"
    return(quadripRegions)
  }
}

plotAverageLines <- function(quadripRegions, coverage, windowSize, positions) {
  averageLines <- c()
  for (i in 1:nrow(quadripRegions)) {
    lineData <- GenerateHistogramData(quadripRegions[i,], coverage, windowSize, (i == nrow(quadripRegions)))
    averageLines <- c(averageLines, paste(quadripRegions[i, 4], ": ", trunc(lineData[1, 4]), "X", sep = ""))
    PACVr.Line.Plot(
      line.data = lineData,
      data.col = 4,
      track.num = 5,
      side = "in",
      min.value = min(coverage[4]),
      max.value = max(coverage[4]),
      outside.pos = positions$outside.pos,
      inside.pos = positions$inside.pos,
      genomic.columns = 3,
      is.sorted = TRUE
    )
  }
  return(averageLines)
}

plotRegionNames <- function(quadripRegions) {
  PACVr.Gene.Name.Plot(
    gene.data = quadripRegions,
    name.col = 4,
    track.num = 1,
    side = "out",
    rotate = 90,
    correction = 0.2,
    add.text.size = 0.2
  )
}

isRealRegions <- function(quadripRegions) {
  return(nrow(quadripRegions) > 1)
}

getIsIRCheck <- function(IRCheck) {
  IRCheckTypes <- getIRCheckTypes()
  return(IRCheck %in% IRCheckTypes)
}

getIRCheckTypes <- function() {
  IRCheckTypes <- c(0, 1, 2)
  return(IRCheckTypes)
}

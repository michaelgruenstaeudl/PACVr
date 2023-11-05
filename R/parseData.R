#!/usr/bin/R
#contributors=c("Michael Gruenstaeudl", "Nils Jenke")
#email="m_gruenstaeudl@fhsu.edu"
#version="2023.11.04.2200"

filter <- function(allRegions, where) {
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

ExtractAllGenes <- function(gbkData) {
  # Function to extract genes from genbank file
  # ARGS:
  #   gbkData (i.e., GenBank flatfile parsed by genbankr)
  # RETURNS:
  #   genes in data frame format
  gene_L <- BiocGenerics::as.data.frame(genbankr::genes(gbkData))        # Use of genbankr
  gene_L <- gene_L[, c(1:3, which(colnames(gene_L) == "gene"))]
  colnames(gene_L) <-
    c("Chromosome", "chromStart", "chromEnd", "gene")
  gene_L$Chromosome <- ""
  gene_L <- gene_L[order(gene_L$chromStart),]
  row.names(gene_L) <- 1:nrow(gene_L)
  return(gene_L)
}

ExtractAllRegions <- function(gbkData) {
  # Function to extract specific regions from Genbank input data
  # ARGS:
  #   gbkData (i.e., GenBank flatfile parsed by genbankr)
  # RETURNS:
  #   regions in data frame format
  allRegions <- BiocGenerics::as.data.frame(genbankr::otherFeatures(gbkData))        # Use of genbankr
  regions <-
    tryCatch(
      tryCatch(
        filter(allRegions, "note"),
        warning = function(w)
          filter(allRegions, "standard_name")
      ),
      error = function(e)
        return(
          data.frame(
            start = c(-1),
            end = c(-1),
            note = c("empty"),
            stringsAsFactors=FALSE
          )
        )
    )
  regions <- regions[, c("start", "end", "note")]
  colnames(regions) <- c("chromStart", "chromEnd", "Band")
  regions$Chromosome <- ""
  regions$Stain <- "gpos100"
  regions <- regions[c("Chromosome", "chromStart", "chromEnd", "Band", "Stain")]
  regions <- regions[order(regions[, 3], decreasing=FALSE),]
  regions$Band[which(grepl("LSC|large|long", regions$Band, ignore.case=TRUE) == TRUE)] <- "LSC"
  regions$Band[which(grepl("SSC|small|short", regions$Band, ignore.case=TRUE) == TRUE)] <- "SSC"
  regions$Band[which(grepl("IRa|\\sa|IR1|\\s1", regions$Band, ignore.case=TRUE) == TRUE)] <- "###A"
  regions$Band[which(grepl("IRb|\\sb|IR2|\\s2", regions$Band, ignore.case=TRUE) == TRUE)] <- "###B"
  regions$Band[which(grepl("IR|invert|repeat", regions$Band, ignore.case=TRUE) == TRUE)] <- "IR"
  regions$Band[which(grepl("###A", regions$Band, ignore.case=TRUE) == TRUE)] <- "IRa"
  regions$Band[which(grepl("###B", regions$Band, ignore.case=TRUE) == TRUE)] <- "IRb"
  row.names(regions) <- 1:nrow(regions)
  regions <- regions[order(regions[, 3], decreasing=FALSE),]
  return(regions)
}


fillDataFrame <- function(gbkData, regions) {
  seqlength <- genbankr::seqinfo(gbkData)@seqlengths                    # Use of genbankr
  if ((nrow(regions) == 0) || (regions[1, 2] == -1)) {
    regions[1,] <-
      c("", as.numeric(1), as.numeric(seqlength), "NA", "gpos100")
    regions[, 2] <- as.numeric(regions[, 2])
    regions[, 3] <- as.numeric(regions[, 3])
    return(regions)
  } else {
    start <- 1
    for (i in 1:nrow(regions)) {
      if (regions[i, 2] > start) {
        regions[nrow(regions) + 1,] <- c("", start, as.numeric(regions[i, 2]) - 1, "NA", "gpos100")
      }
      start <- as.numeric(regions[i, 3]) + 1
    }
    if (start - 1 < seqlength) {
      regions[nrow(regions) + 1,] <- c("", start, seqlength, "NA", "gpos100")
    }
    regions <-
      regions[order(as.numeric(regions[, 2]), decreasing=FALSE),]
    row.names(regions) <- 1:nrow(regions)
    regions[, 2] <- as.numeric(regions[, 2])
    regions[, 3] <- as.numeric(regions[, 3])
    
    regionAvail <- boolToDeci(c("LSC", "IRb", "SSC", "IRa") %in% regions[, 4])
    regions[, 6] <- regions[, 3] - regions[, 2]
    
    if (regionAvail == 5) {
      # only IRa and IRb
      regions[which(regions[, 4] != "NA"), 6] <- 0
      regions[which(regions[, 6] == max(regions[, 6])), 4] <- "LSC"
      regions[which(regions[, 4] != "NA"), 6] <- 0
      regions[which(regions[, 6] == max(regions[, 6])), 4] <- "SSC"
      message("Annotation for LSC and SSC were automatically added")
    } else if (regionAvail == 7) {
      # only IRa, SSC and IRb
      regions[which(regions[, 4] != "NA"), 6] <- 0
      regions[which(regions[, 6] == max(regions[, 6])), 4] <- "LSC"
      message("Annotation for LSC was automatically added")
    } else if (regionAvail == 10) {
      # only LSC and SSC
      IRs <- data.frame(table(regions[which(regions[, 4] == "NA"), 6]), stringsAsFactors=FALSE)
      IRs <- IRs[IRs$Freq == 2, 1]
      if (length(IRs) >= 1) {
        IRs <- max(as.numeric(as.character(IRs)))
        regions[which(regions[, 6] == IRs), 4] <- c("IRb", "IRa")
        message("Annotation for IRb and IRa were automatically added")
      }
    } else if (regionAvail == 11) {
      # only LSC, SSC and IRa
      regions[which(regions[, 4] == "NA" &
                    regions[, 6] == regions[which(regions[, 4] == "IRa"), 6]), 4] <- "IRb"
      message("Annotation for IRb was automatically added")
    } else if (regionAvail == 13) {
      # only LSC, IRb and IRa
      regions[which(regions[, 4] != "NA"), 6] <- 0
      regions[which(regions[, 6] == max(regions[, 6])), 4] <- "SSC"
      message("Annotation for SSC was automatically added")
    } else if (regionAvail == 14) {
      # only LSC, IRb and SSC
      regions[which(regions[, 4] == "NA" &
                    regions[, 6] == regions[which(regions[, 4] == "IRb"), 6]), 4] <- "IRa"
      message("Annotation for IRa was automatically added")
    }
    regions <- regions[-6]
    regions$Stain[which(regions$Band == "LSC")] <- "gpos75"
    regions$Stain[which(regions$Band == "SSC")] <- "gpos50"
    regions$Stain[which(regions$Band == "IRa")] <- "gpos25"
    regions$Stain[which(regions$Band == "IRb")] <- "gpos25"
    return(regions)
  }
}

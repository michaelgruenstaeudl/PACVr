#!/usr/bin/env RScript
#contributors=c("Gregory Smith", "Nils Jenke", "Michael Gruenstaeudl")
#email="m_gruenstaeudl@fhsu.edu"
#version="2024.05.20.1821"

read.gbWithHandling <- function(gbkRaw, count=0) {
  gbkData <- tryCatch({
    suppressMessages(
      read.gb::read.gb(gbkRaw$char, DNA=TRUE, Type="full", Source="Char")
    )
  },
    error = function(e) {
      if (conditionMessage(e) == "dim(X) must have a positive length") {
        logger::log_warn(
          paste("Encountered an unqualified feature in read.gb;",
                "permanent fix w/",
                "devtools::install_github(\"alephnull7/read.gb\")")
        )
        return(read.gbOnFix(gbkRaw, count))
      } else {
        logger::log_error(
          logger::skip_formatter(
            paste("read.gb encountered an unexpected error:", e)
          )
        )
        return(NULL)
      }
    })
  return(gbkData)
}

read.gbOnFix <- function(gbkRaw, count) {
  count <- count+1
  if (count > 1) {
    logger::log_error("Unable to fix file for read.gb")
    return(NULL)
  }
  
  gbkRaw <- fixGbkChar(gbkRaw)
  logger::log_info("Retrying read.gb with fixed data")
  gbkData <- read.gbWithHandling(gbkRaw, count)
  return(gbkData)
}

# adapted from read.gb::read.gb
fixGbkChar <- function(gbkRaw) {
  gbkFile <- gbkRaw$file
  gbkChar <- gbkRaw$char
  
  # intended directory of fixed file if file was provided
  if (!is.null(gbkFile)) {
    gbkFileFixed <- gsub("\\.gb$", "-PACVr.gb",  gbkFile)
  } else {
    gbkFileFixed <- NULL
  }

  # used fixed file if already created
  if (!is.null(gbkFileFixed) && file.exists(gbkFileFixed)) {
    return(getGbkRaw(gbkFileFixed))
  }

  ## Separation of reports :
  gbkCharFixed <- ""
  SampleS <- gregexpr("LOCUS {2,}", gbkChar)[[1]]
  SampleE <- gregexpr("(?<!:)//", gbkChar, perl = T)[[1]]
  
  ## Treatment of reports
  for(k in 1:length(SampleS)){
    startIndex <- SampleS[k]
    endIndex <- SampleE[k]+1
    Temp <- substr(gbkChar, startIndex, endIndex)
    y <- Reorganize.report(Temp)
    y <- fixReferences(y)
    Temp <- Collapse.report(y)
    gbkCharFixed <- paste0(gbkCharFixed, Temp)
  }
  
  if (!is.null(gbkFileFixed)) {
    writeLines(gbkCharFixed, gbkFileFixed)
    logger::log_info('Fixed GenBank file saved as `{gbkFileFixed}`')
  }
  
  gbkRaw$char <- gbkCharFixed
  return(gbkRaw)
}

# adapted from read.gb::Reorganize.report
Reorganize.report <- function(Temp) {
  Fb <- gregexpr("FEATURES", Temp)[[1]][1]
  Fe <- gregexpr("ORIGIN", Temp, ignore.case = F)[[1]][1]
  if(Fe == -1){Fe <- gregexpr("CONTIG", Temp, ignore.case = F)[[1]][1]}
  
  y <- list()
  y$BEFORE_FEATURES <- substr(Temp, 1, Fb-1)
  y$FEATURES <- substr(Temp, Fb, Fe-1)
  y$AFTER_FEATURES <- substring(Temp, Fe, nchar(Temp))
  return(y)
}

Collapse.report <- function(y) {
  Temp <- paste0(y, collapse = "")
  return(Temp)
}

# adapted from read.gb::Reference.sep
fixReferences <- function(y) {
  TempF <- y$FEATURES
  Tag <- getFeatureTags()

  ## Separation of Feature items :
  TagChr <- c()
  for(i in 1:length(Tag)){
    if(grepl(Tag[i], TempF, fixed = F)[[1]] == TRUE){
      P <- as.numeric(gregexpr(Tag[i], TempF, fixed = F)[[1]])
      for(j in 1:length(P)){
        if(P[j] != -1){
          TagChr <- c(TagChr, P[j])
        }
      }
    }
  }
  TagChr <- c(TagChr, nchar(TempF))
  TagChr <- sort(TagChr)
  
  spaceLength <- getSpaceLength(TempF)
  Feature <- c(substr(TempF, 1, TagChr[1]-1))
  for(k in 2:length(TagChr)){
    startIndex <- TagChr[k-1]
    endIndex <- TagChr[k]-1
    Feat <- substr(TempF, startIndex, endIndex)
    if (!grepl("\\/", Feat)) {
      Feat <- paste0(Feat, 
                     paste(rep(" ", spaceLength), collapse = ""),
                     "/note=\"PACVr-placeholder\"\n")
    }
    Feature <- c(Feature, Feat)
  }
  Feature <- c(Feature, 
               substr(TempF, TagChr[length(TagChr)], nchar(TempF)))
  y$FEATURES <- paste0(Feature, collapse = "")
  return(y)
}

getSpaceLength <- function(TempF) {
  spaceMatch <- regexpr("\\n +\\/", TempF)
  if (spaceMatch == -1) {
    return(21)
  }
  spaceLength <- attr(spaceMatch, "match.length")
  spaceLength <- spaceLength - 2
  return(spaceLength)
}

# taken from read.gb::Reference.sep
getFeatureTags <- function() {
  Tag <- c(" +assembly_gap {2,}", " +C_region {2,}", 
           " +CDS {2,}", " +centromere {2,}", 
           " +D-loop {2,}", " +D_segment {2,}", 
           " +exon {2,}", " +gap {2,}", 
           " +gene {2,}", " +iDNA {2,}", 
           " +intron {2,}", " +J_segment {2,}", 
           " +mat_peptide {2,}", " +misc_binding {2,}", 
           " +misc_difference {2,}", " +misc_feature {2,}", 
           " +misc_recomb {2,}", " +misc_RNA {2,}", 
           " +misc_structure {2,}", " +mobile_element {2,}", 
           " +modified_base {2,}", " +mRNA {2,}", 
           " +ncRNA {2,}", " +N_region {2,}", 
           " +old_sequence {2,}", " +operon {2,}", 
           " +oriT {2,}", " +polyA_site {2,}", 
           " +precursor_RNA {2,}", " +prim_transcript {2,}", 
           " +primer_bind {2,}", " +propeptide {2,}", 
           " +protein_bind {2,}", " +regulatory {2,}", 
           " +repeat_region {2,}", " +rep_origin {2,}", 
           " +rRNA {2,}", " +S_region {2,}", 
           " +sig_peptide {2,}", " +source {2,}", 
           " +stem_loop {2,}", " +STS {2,}", 
           " +telomere {2,}", " +tmRNA {2,}", 
           " +transit_peptide {2,}", " +tRNA {2,}", 
           " +unsure {2,}", " +V_region {2,}", 
           " +V_segment {2,}", " +variation {2,}", 
           " +3[[:punct:]]UTR +"," +5[[:punct:]]UTR +")
  return(Tag)
}

getGbkRaw <- function(gbkFile) {
  if (!is.character(gbkFile)) {
    logger::log_error("Parameter `gbkFile` is not a character string")
    gbkChar <- NULL
  } else if (grepl("\\.gb$", gbkFile) && file.exists(gbkFile)) {
    logger::log_info('Reading GenBank flatfile `{gbkFile}`')
    gbkChar <- tryCatch({
      suppressWarnings(readChar(gbkFile, 
                                file.info(gbkFile)$size, 
                                nchars = 99999999)
                       )
    },
    error = function(e) {
      logger::log_error(paste0('Error when reading GenBank flatfile',
                              '`{gbkFile}`: "{e}"')
                        )
      return(NULL)
    })
  } else if (grepl("\\.gb$", gbkFile)) {
    logger::log_error("GenBank flatfile does not exist `{gbkFile}`")
    gbkChar <- NULL
  } else if (grepl("\\.\\w+$", gbkFile)) {
    logger::log_error("Non-GenBank file provided `{gbkFile}`")
    gbkChar <- NULL
  } else {
    logger::log_info(paste0("Reading GenBank string \"",
                            substr(gbkFile, 1, 20),
                            "...\"")
                    )
    gbkChar <- gbkFile
    gbkFile <- NULL
  }
  gbkRaw <- list(file = gbkFile,
                 char = gbkChar)
  return(gbkRaw)
}

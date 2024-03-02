#!/usr/bin/env RScript
#contributors=c("Gregory Smith", "Nils Jenke", "Michael Gruenstaeudl")
#email="m_gruenstaeudl@fhsu.edu"
#version="2024.03.02.0100"

GBKData <- R6::R6Class("GBKData",
  public = list(
    # fields
    analysisSpecs = NULL,
    read.gb = NULL,
    features = NULL,
    genes = NULL,
    sequences = NULL,
    lengths = NULL,
    sampleName = NULL,
    plotTitle = NULL,
    quadripRegions = NULL,
    linkData = NULL,
    
    # constructor
    initialize = function(gbkFile,
                          analysisSpecs) {
      self$read.gb <- PACVr.read.gb(gbkFile)
      self$analysisSpecs <- analysisSpecs

      # main derivative of `read.gb` data
      private$setFeatures()
      if (is.null(self$features)) {
        logger::log_fatal('Parsing of any sequence features unsuccessful.')
        return(NULL)
      }

      # other derivatives of `read.gb` data
      private$setSequences()
      private$setLengths()
      private$setSampleName()
      private$setPlotTitle()

      # `read.gb` no longer needed
      self$read.gb <- NULL
      gc()

      # `features` derivatives
      private$setQuadripRegions()
      private$setGenes()
      private$setLinkData()

      # `features` no longer needed
      self$features <- NULL
      gc()
    }
  ),

  # private setters for constructor
  private = list(
    # precondition: `read.gb` and `analysisSpecs` are set
    setFeatures = function() {
      fileDF <- data.frame()
      for (sample in self$read.gb) {
        sampleDF <- parseFeatures(sample$FEATURES, self$analysisSpecs)
        if (!is.null(sampleDF)) {
          fileDF <- dplyr::bind_rows(fileDF, sampleDF)
        }
      }
      if (ncol(fileDF) > 0 || nrow(fileDF) > 0) {
        self$features <- fileDF
      }
    },

    # precondition: `read.gb` is set
    setSequences = function() {
      sampleSequences <- NULL
      for (sample in self$read.gb) {
        sampleSequences <- c(sampleSequences, sample$ORIGIN)
      }
      self$sequences <- Biostrings::DNAStringSet(sampleSequences)
    },

    # precondition: `seq` is set
    setLengths = function() {
      self$lengths <- Biostrings::width(self$sequences)
    },

    # precondition: `read.gb` is set
    setSampleName = function() {
      sampleNames <- NULL
      for (sample in self$read.gb) {
        sampleNames <- c(sampleNames,
                         c(sample_name = sample$VERSION,
                           genome_name = sample$ACCESSION))
      }
      self$sampleName <- sampleNames
    },

    # precondition: `read.gb` is set
    setPlotTitle = function() {
      plotTitles <- NULL
      for (sample in self$read.gb) {
        plotTitles <- c(plotTitles, paste(sample$DEFINITION, sample$ACCESSION))
      }
      self$plotTitle <- plotTitles
    },

    # precondition: `features` and `analysisSpecs` are set
    setQuadripRegions = function() {
      if (self$analysisSpecs$isIRCheck) {
        logger::log_info('Parsing the quadripartite genome structure')
        quadripRegions <- PACVr.parseQuadripRegions(self$lengths,
                                                    self$features)
      } else {
        quadripRegions <- PACVr.parseSource(self$features)
      }
      self$quadripRegions <- quadripRegions
    },

    # precondition: `features` is set
    setGenes = function() {
      logger::log_info('  Extracting information on genes')
      gene_L <- read.gbGenes(self$features)
      gene_L <- gene_L[, c(1:3, which(colnames(gene_L) == "gene"))]
      colnames(gene_L) <- c("Chromosome", "chromStart", "chromEnd", "gene")
      gene_L$Chromosome <- ""
      gene_L <- gene_L[order(gene_L$chromStart),]
      row.names(gene_L) <- 1:nrow(gene_L)
      self$genes <- gene_L
    },

    # precondition: `genes`, `quadripRegions`, and `analysisSpecs` are set
    setLinkData = function() {
      if (self$analysisSpecs$isSyntenyLine) {
        logger::log_info('Inferring the IR regions and the genes within the IRs')
        self$linkData <- PACVr.generateIRGeneData(self$genes,
                                                  self$quadripRegions,
                                                  self$analysisSpecs$syntenyLineType)
      }
    }
  )
)

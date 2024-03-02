#!/usr/bin/env RScript
#contributors=c("Gregory Smith", "Nils Jenke", "Michael Gruenstaeudl")
#email="m_gruenstaeudl@fhsu.edu"
#version="2024.03.02.0100"

GBKData <- R6::R6Class("GBKData",
  public = list(
    # fields
    analysisSpecs = NULL,
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
      read.gbData <- PACVr.read.gb(gbkFile)
      self$analysisSpecs <- analysisSpecs

      # main derivative of `read.gb` data
      gbkSeqFeatures <- read.gbSeqFeaturesAdapt(read.gbData, analysisSpecs)
      if (is.null(gbkSeqFeatures)) {
        logger::log_fatal('Parsing of any sequence features unsuccessful.')
        return(NULL)
      }

      # other derivatives of `read.gb` data
      private$setSequences(read.gbData)
      private$setLengths()
      private$setSampleName(read.gbData)
      private$setPlotTitle(read.gbData)

      # `read.gbData` no longer needed
      rm(read.gbData)
      gc()

      # `gbkSeqFeatures` derivatives
      private$setQuadripRegions(gbkSeqFeatures)
      private$setGenes(gbkSeqFeatures)
      private$setIRCheckFields()
      private$setLinkData()

      # `gbkSeqFeatures` no longer needed
      rm(gbkSeqFeatures)
      gc()
    }
  ),

  # private setters for constructor
  private = list(
    setSequences = function(read.gbData) {
      sampleSequences <- NULL
      for (sample in read.gbData) {
        sampleSequences <- c(sampleSequences, sample$ORIGIN)
      }
      self$sequences <- Biostrings::DNAStringSet(sampleSequences)
    },

    # precondition: `seq` is set
    setLengths = function() {
      self$lengths <- Biostrings::width(self$sequences)
    },

    setSampleName = function(read.gbData) {
      sampleNames <- NULL
      for (sample in read.gbData) {
        sampleNames <- c(sampleNames,
                         c(sample_name = sample$VERSION,
                           genome_name = sample$ACCESSION))
      }
      self$sampleName <- sampleNames
    },

    setPlotTitle = function(read.gbData) {
      plotTitles <- NULL
      for (sample in read.gbData) {
        plotTitles <- c(plotTitles, paste(sample$DEFINITION, sample$ACCESSION))
      }
      self$plotTitle <- plotTitles
    },

    # precondition: `analysisSpecs` is set
    setQuadripRegions = function(gbkSeqFeatures) {
      if (self$analysisSpecs$isIRCheck) {
        logger::log_info('Parsing the quadripartite genome structure')
        quadripRegions <- PACVr.parseQuadripRegions(self$lengths,
                                                    gbkSeqFeatures,
                                                    self$analysisSpecs)
      } else {
        quadripRegions <- PACVr.parseSource(gbkSeqFeatures)
      }
      self$quadripRegions <- quadripRegions
    },

    setGenes = function(gbkSeqFeatures) {
      logger::log_info('  Extracting information on genes')
      gene_L <- read.gbGenes(gbkSeqFeatures)
      gene_L <- gene_L[, c(1:3, which(colnames(gene_L) == "gene"))]
      colnames(gene_L) <- c("Chromosome", "chromStart", "chromEnd", "gene")
      gene_L$Chromosome <- ""
      gene_L <- gene_L[order(gene_L$chromStart),]
      row.names(gene_L) <- 1:nrow(gene_L)
      self$genes <- gene_L
    },

    # precondition: `quadripRegions` is set
    setIRCheckFields = function() {
      if (self$analysisSpecs$isSyntenyLine &&
          (!("IRb" %in% self$quadripRegions[, 4]) ||
           !("IRa" %in% self$quadripRegions[, 4]) )) {
        logger::log_warn("Unable to find synteny: missing `IRa` or `IRb`")
        self$analysisSpecs$setIRCheckFields(0)
      }
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

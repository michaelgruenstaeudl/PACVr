#!/usr/bin/env RScript
#contributors=c("Gregory Smith", "Nils Jenke", "Michael Gruenstaeudl")
#email="m_gruenstaeudl@fhsu.edu"
#version="2024.05.19.1300"

# pacman::p_load loads packages if they have been installed 
# and installs if they are missing
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, vroom)

### INIT NEEDED DIRECTORIES ###
set_script_dir <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (interactive()) {
    script_dir <- getwd()
  } else if (length(match) > 0) {
    script_dir <- normalizePath(dirname(sub(needle, "", cmdArgs[match])))
  } else {
    script_dir <- normalizePath(dirname(sys.frames()[[1]]$ofile))
  }
  setwd(script_dir)
}

set_script_dir()


### METADATA ###
create_metadata_df <- function() {
  # collect files
  metadata_files <- list.files(
    path = inDir_data,
    full.names = TRUE,
    recursive = TRUE,
    pattern = "metadata.csv"
  )

  # create df
  meta_data_list <- purrr::map(metadata_files, ~ {
    data <- vroom::vroom(
      .x,
      delim = ",",
      col_types = c(
        avgLength = "i",
        Model = "c",
        `Assembly Method` = "c",
        `Sequencing Technology` = "c"
      )
    )
    data$Accession <- gsub("_metadata.csv", "", basename(.x))
    data
  })
  meta_data <- bind_rows(meta_data_list)

  return(meta_data)
}

normalize_metadata_avgLength <- function(meta_data) {
  # remove NAs
  meta_data <- meta_data[, colSums(is.na(meta_data)) != nrow(meta_data)] %>%
    arrange(Accession)

  # calculate new avgLength for samples with multiple SRA metadata
  meta_data <- meta_data %>%
    group_by(Accession) %>%
    summarize(newAvgLength = round(mean(avgLength, na.rm = TRUE))) %>%
    inner_join(meta_data, ., by = "Accession") %>%
    distinct(Accession, .keep_all = TRUE) %>%
    select(-avgLength) %>%
    rename(avgLength = newAvgLength)

  return(meta_data)
}

join_metadata_with_samplelist <- function(meta_data) {
  sample_list <- read.csv(inFn_sampleList, header = TRUE)
  meta_data <- inner_join(meta_data, sample_list, by = "Accession")
  return(meta_data)
}

normalize_metadata_names <- function(meta_data) {
  # resolve spelling mistake
  if ("Eragrostis_tef" %in% meta_data$Samples) {
	  meta_data[which(meta_data$Samples == "Eragrostis_tef"), ]$`Assembly Method` <- "Geneious"
  }

  # unify assembly methods
  nameList <- c(
    "Velvet",
    "Consed",
    "Newbler",
    "Spades",
    "CLC",
    "SOAP",
    "Mira",
    "Ray",
    "Geneious",
    "Obitools",
    "Flash",
    "PriceTI",
    "Spades",
    "Yasra",
    "GS",
    "abyss",
    "Allpath"
  )
  assignList <- c(
    "Velvet",
    "Consed",
    "Newbler Assembler",
    "SPAdes",
    "CLC Assembly",
    "SOAPdenovo",
    "MIRA",
    "Ray",
    "Geneious",
    "OBITools",
    "FLASH",
    "PriceTI",
    "SPAdes",
    "YASRA",
    "GS De novo assembler",
    "ABySS",
    "ALLPATHS-LG"
  )

  for (i in 1:length(nameList)) {
	if (nameList[i] %in% meta_data$`Assembly Method`) {
	  meta_data$`Assembly Method`[grepl(nameList[i], meta_data$`Assembly Method`, ignore.case = TRUE)] <- assignList[i]
    }
  }

  # standardize sampling methods
  meta_data$Model <- gsub("Illumina ", "", meta_data$Model)
  meta_data$Model <- gsub("Genome Analyzer", "GA", meta_data$Model)

  # annotate avgLength
  meta_data$avgLength <-
    ifelse(
      meta_data$avgLength <= 300,
      paste0(meta_data$avgLength, " (s)"),
      paste0(meta_data$avgLength, " (m)")
    )

  return(meta_data)
}

get_metadata_df <- function() {
  meta_data <- create_metadata_df() %>%
    normalize_metadata_avgLength() %>%
    join_metadata_with_samplelist() %>%
    normalize_metadata_names()

  return(meta_data)
}

### REGIONS COVERAGE STATS ###

create_cov_sum_df <- function(file_pattern) {
  # collect files
  cov_sum_files <- list.files(
    path = inDir_data,
    full.names = TRUE,
    recursive = TRUE,
    pattern = file_pattern
  )

  # create df
  cov_sum_data_list <- purrr::map(cov_sum_files, ~ {
    data <- vroom::vroom(
      .x,
      delim = "\t",
      col_types = c(
        lowCovWin_abs = "i",
        regionLen = "i",
        lowCovWin_relToRegionLen = "n",
        E_score = "n",
        N_count = "i"
      )
    )
    data$Accession <- gsub(file_pattern, "", basename(.x))
    data
  })
  cov_sum_data <- bind_rows(cov_sum_data_list)

  return(cov_sum_data)
}

transform_regions_sum <- function(regions_sum) {
  regions_wide <- regions_sum %>%
    pivot_wider(
      id_cols = Accession,
      names_from = Chromosome,
      values_from = lowCovWin_relToRegionLen
    )

  complete_sum <- regions_sum %>%
    filter(Chromosome == "Complete_genome") %>%
    select(E_score, N_count, IR_mismatches, Accession)

  regions_sum <- inner_join(regions_wide, complete_sum, by = "Accession")
  return(regions_sum)
}

transform_cov_sum <- function(cov_sum, col_name) {
  col_name_sym <- ensym(col_name)
  
  cov_sum <- cov_sum %>%
    group_by(Accession) %>%
    summarise(!!col_name_sym := sum(lowCovWin_abs) / sum(regionLen))

  return(cov_sum)
}

get_regions_sum_df <- function() {
  regions_sum_df <- create_cov_sum_df(".[0-9]_summary.regions.tsv$") %>%
    transform_regions_sum()
  return(regions_sum_df)
}

get_coding_sum_df <- function() {
  coding_sum_df <- create_cov_sum_df(".[0-9]_coverage.summary.genes.tsv$") %>%
    transform_cov_sum("coding")
  return(coding_sum_df)
}

get_noncoding_sum_df <- function() {
  noncoding_sum_df <- create_cov_sum_df(".[0-9]_coverage.summary.noncoding.tsv$") %>%
    transform_cov_sum("noncoding")
  return(noncoding_sum_df)
}

### COMBINE METADATA AND COVERAGE ###

get_cov_df <- function() {
  meta_data <- get_metadata_df()
  regions_sum <- get_regions_sum_df()
  coding_sum <- get_coding_sum_df()
  noncoding_sum <- get_noncoding_sum_df()

  df_list <- list(meta_data, regions_sum, coding_sum, noncoding_sum)
  cov_data <- reduce(df_list, inner_join, by = "Accession")
  return(cov_data)
}

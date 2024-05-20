#!/usr/bin/env RScript
#contributors=c("Gregory Smith", "Nils Jenke", "Michael Gruenstaeudl")
#email="m_gruenstaeudl@fhsu.edu"
#version="2024.05.19.1300"

source("coverage_data_assembly.R")

input_path <- "input/samples_list.csv"
metadata_extraction_all <- function() {
  set_script_dir()
  samples_list <- read.csv(input_path, stringsAsFactors = FALSE)

  for (index in 1:nrow(samples_list)) {
    sample <- samples_list[index,]
    sra <- sample["SRA"]
    accession <- sample["Accession"]
    
    pacvr_cmd <- c("./metadata_extraction.sh", accession, sra)
    system2(command = "bash", args = pacvr_cmd)
  }
}

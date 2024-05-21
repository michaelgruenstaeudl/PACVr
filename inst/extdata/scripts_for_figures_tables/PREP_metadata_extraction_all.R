#!/usr/bin/env RScript
#contributors=c("Gregory Smith", "Nils Jenke", "Michael Gruenstaeudl")
#email="m_gruenstaeudl@fhsu.edu"
#version="2024.05.19.1300"

source("PREP_coverage_data_assembly.R")

metadata_extraction_all <- function() {
  set_script_dir()
  JenkeEtAl2024_samples_list <- read.csv(inFn_sampleList, stringsAsFactors = FALSE)

  for (index in 1:nrow(JenkeEtAl2024_samples_list)) {
    sample <- JenkeEtAl2024_samples_list[index,]
    sra <- sample["SRA"]
    accession <- sample["Accession"]
    
    pacvr_cmd <- c("./PREP_metadata_extraction.sh", inDir_data, accession, sra)
    system2(command = "bash", args = pacvr_cmd)
  }
}

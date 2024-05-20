#!/usr/bin/env RScript
#contributors=c("Gregory Smith", "Nils Jenke", "Michael Gruenstaeudl")
#email="m_gruenstaeudl@fhsu.edu"
#version="2024.05.19.1300"

source("metadata_extraction_all.R")
source("coverage_tables.R")
source("coverage_figures.R")

main <- function() {
  metadata_extraction_all()
  dir.create("images")
  create_tables()
  create_figures()
}

main()

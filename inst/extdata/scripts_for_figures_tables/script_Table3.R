#!/usr/bin/env RScript
#contributors=c("Gregory Smith", "Nils Jenke", "Michael Gruenstaeudl")
#email="m_gruenstaeudl@fhsu.edu"
#version="2024.05.19.1300"

library(tidyverse)
library(xtable)
library(rstatix)

source("metadata_extraction_all.R")
source("coverage_data_assembly.R")
source("coverage_data_preparation.R")

# TABLE 3 - summary of Kruskal-Wallis tests
create_table_3 <- function(figure_data) {
  kruskal_results <- get_kruskal_results(figure_data)
  xtab_3 <- xtable(kruskal_results, digits = 3)
  print(
    xtab_3,
    file = "./images/JenkeEtAl2024_Table_3.tex",
    compress = FALSE,
    include.rownames = FALSE
  )
}

# Execute commands
metadata_extraction_all()
cov_data <- get_cov_df()
figure_data <- get_figure_data(cov_data)
dir.create("images")
create_table_3(figure_data)

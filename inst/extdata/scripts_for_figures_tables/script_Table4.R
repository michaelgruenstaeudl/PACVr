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

# TABLE 4 - summary of Wilcoxon tests
create_table_4 <- function(figure_data) {
  table_4 <- get_wilcox_results(figure_data) %>%
    filter(variable == "low_coverage") %>%
    select(-c(variable,statistic)) %>%
    rename("State 1" = group1,
           "State 2" = group2,
           d_c = effsize)
  xtab_4 <- xtable(table_4, digits = 3)
  print(
    xtab_4,
    file = "./images/JenkeEtAl2024_Table_4.tex",
    compress = FALSE,
    include.rownames = FALSE
  )
}

# Execute commands
metadata_extraction_all()
cov_data <- get_cov_df()
figure_data <- get_figure_data(cov_data)
dir.create("images")
create_table_4(figure_data)

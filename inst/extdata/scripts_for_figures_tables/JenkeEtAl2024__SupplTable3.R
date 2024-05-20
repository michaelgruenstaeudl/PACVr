#!/usr/bin/env RScript
#contributors=c("Gregory Smith", "Nils Jenke", "Michael Gruenstaeudl")
#email="m_gruenstaeudl@fhsu.edu"
#version="2024.05.19.1300"

library(tcltk) # For dialog boxes
library(tidyverse)
library(xtable)
library(rstatix)

w_dir = tk_choose.dir(caption='Select directory with the .gb and .bam files')

source("PREP_metadata_extraction_all.R")
source("PREP_coverage_data_assembly.R")
source("PREP_coverage_data_preparation.R")

# SUPPLEMENTARY TABLE 3 - pairwise Wilcoxon tests for evenness
create_supp_table_3 <- function(figure_data) {
  supp_table_3 <- get_wilcox_results(figure_data) %>%
    filter(variable == "E_score") %>%
    filter(group1 != "medium") %>%
    filter(group1 != "n.s.") %>%
    filter(group2 != "n.s.") %>%
    select(-c(variable,p,statistic)) %>%
    rename("State 1" = group1,
           "State 2" = group2,
           d_c = effsize)
  supp_xtab_3 <- xtable(supp_table_3, digits = 3)
  print(
    supp_xtab_3,
    file = "./images/JenkeEtAl2024_SupplTable_3.tex",
    compress = FALSE,
    include.rownames = FALSE
  )
}

# Execute commands
metadata_extraction_all()
cov_data <- get_cov_df()
figure_data <- get_figure_data(cov_data)
dir.create("images")
create_supp_table_3(figure_data)

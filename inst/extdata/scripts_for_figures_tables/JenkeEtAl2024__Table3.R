#!/usr/bin/env RScript
#contributors=c("Gregory Smith", "Nils Jenke", "Michael Gruenstaeudl")
#email="m_gruenstaeudl@fhsu.edu"
#version="2024.05.26.1500"

# pacman::p_load loads packages if they have been installed 
# and installs if they are missing
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tcltk, tidyverse, xtable, rstatix)

# Select input directory and input files
inDir_data = tk_choose.dir(default="~", caption='Select directory with the .gb and .bam files')
inFn_sampleList = tk_choose.files(caption='Select samples list file in csv-format (e.g., input/JenkeEtAl2024_samples_list.csv)')

source("PREP_metadata_extraction_all.R")
source("PREP_coverage_data_assembly.R")
source("PREP_coverage_data_preparation.R")

# TABLE 3 - summary of Kruskal-Wallis tests
create_table_3 <- function(figure_data) {
  kruskal_results <- get_kruskal_results(figure_data) %>% 
    mutate(p = pval_asterisk(format(round(p, digits=3), nsmall=3))) %>% 
    mutate(d_c = effectsize_symbol(format(round(d_c, digits=3), nsmall=3)))

  xtab_3 <- xtable(kruskal_results, digits=3)
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

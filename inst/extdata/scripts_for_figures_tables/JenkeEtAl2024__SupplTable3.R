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
           d_c = effsize) %>% 
    mutate(p.adj = pval_asterisk(format(round(p.adj, digits=3), nsmall=3))) %>% 
    mutate(d_c = effectsize_symbol(format(round(d_c, digits=3), nsmall=3)))

  supp_xtab_3 <- xtable(supp_table_3, digits=3)
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

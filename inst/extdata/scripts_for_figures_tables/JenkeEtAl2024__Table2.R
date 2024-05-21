#!/usr/bin/env RScript
#contributors=c("Gregory Smith", "Nils Jenke", "Michael Gruenstaeudl")
#email="m_gruenstaeudl@fhsu.edu"
#version="2024.05.19.1300"

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

# TABLE 2 - descriptive statistics for the partitions
create_table_2 <- function(cov_data) {
  table_2 <-
    cov_data %>% select(E_score, LSC, IRb, SSC, IRa, 
                        coding, noncoding, N_count, IR_mismatches) %>%
    rename("E-score" = E_score,
           Ns = N_count,
           mismatches = IR_mismatches) %>%
    rstatix::get_summary_stats(show = c("n", "min", "max", "q1", "q3",
                                        "median", "mean", "sd")) %>%
    mutate("NA" = nrow(cov_data) - .$n) %>%
    select(variable, n, "NA", min, max, q1, q3, median, mean, sd)
  
  xtab_2 <- xtable(table_2, digits = 3)
  print(
    xtab_2,
    file = "./images/JenkeEtAl2024_Table_2.tex",
    compress = FALSE,
    include.rownames = FALSE
  )
}

# Execute commands
metadata_extraction_all()
cov_data <- get_cov_df()
dir.create("images")
create_table_2(cov_data)

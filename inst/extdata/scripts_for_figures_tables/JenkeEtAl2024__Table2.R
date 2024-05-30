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

# TABLE 2 - descriptive statistics for the partitions
create_table_2 <- function(cov_data) {
  # prepare source dataframe for rstatix
  table_2 <-
    cov_data %>% select(LSC, IRb, SSC, IRa,
                        coding, noncoding,
                        E_score,
                        N_count, IR_mismatches)
  wrsd_cols <- c(
    "LSC", "IRb", "SSC", "IRa",
    "coding", "noncoding"
  )
  for (col_name in wrsd_cols) {
    table_2 <- handle_outliers(table_2, col_name, 3, 1)
  }

  # create summary table
  table_2 <-
    table_2 %>%
    rename("E-score" = E_score,
           Ns = N_count,
           mismatches = IR_mismatches) %>%
    rstatix::get_summary_stats(show = c("n", "min", "max", "q1", "q3",
                                        "median", "mean", "sd")) %>%
    mutate("NA" = nrow(cov_data) - .$n) %>%
    select(variable, n, "NA", min, max, q1, q3, median, mean, sd)
  
  xtab_2 <- xtable(table_2, digits=c(0,0,0,0,2,2,2,2,2,2,2))
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

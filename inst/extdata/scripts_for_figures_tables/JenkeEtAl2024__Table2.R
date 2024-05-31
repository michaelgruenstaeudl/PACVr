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
  # long WRSD table with labeled outliers
  outlier_label_df <- cov_data %>%
    select(Accession, LSC, IRb, SSC, IRa, coding, noncoding) %>%
    pivot_longer(
      cols = -Accession,
      names_to = "regions", 
      values_to = "low_coverage",
      names_prefix = "original_"
    ) %>%
    handle_cov_outliers("regions", "low_coverage", 2)

  # summarized outliers
  outlier_df <- outlier_label_df %>%
    filter(outlier == TRUE) %>%
    group_by(regions) %>%
    summarize(outliers = sum(n())) %>%
    rename(variable = regions)

  # non-WRSD metrics
  asym_qual_df <- cov_data %>%
    select(Accession, N_count, IR_mismatches, E_score) %>%
    rename(
      "E-score" = E_score,
      Ns = N_count,
      mismatches = IR_mismatches
    )

  # create summary stats table
  table_2 <- outlier_label_df %>%
    filter(outlier == FALSE) %>%
    pivot_wider(
      id_cols = Accession,
      names_from = regions,
      values_from = low_coverage
    ) %>%
    full_join(., asym_qual_df) %>%
    select(
      LSC, IRb, SSC, IRa,
      coding, noncoding,
      "E-score",
      Ns, mismatches
    ) %>%
    rstatix::get_summary_stats(
      show = c("n", "min", "max", "q1", "q3",
               "median", "mean", "sd")
    ) %>%
    full_join(., outlier_df) %>%
    mutate("NA" = nrow(cov_data) - .$n - ifelse(is.na(outliers), 0, outliers)) %>%
    select(variable, n, "NA", outliers, min, max, q1, q3, median, mean, sd)
  
  xtab_2 <- xtable(table_2, digits=c(0,0,0,0,0,2,2,2,2,2,2,2))
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

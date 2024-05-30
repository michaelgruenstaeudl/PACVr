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

# TABLE 1 - overview of sequencing platform and assembly software evenness
get_table_1_subtable <- function(cov_data, grouping) {
  subtable <- cov_data %>%
    rename(State = !!ensym(grouping)) %>%
    handle_cov_outliers("State", "E_score", 2) %>%
    group_by(State) %>%
    mutate(count = sum(outlier == FALSE)) %>%
    ungroup() %>%
    mutate(State = replace_na(State, "missing")) %>%
    mutate(State := ifelse(count < 5, "missing", State)) %>%
    select(-count) %>%
    group_by(State) %>%
    summarize(count = ifelse(
                first(State) == "missing",
                sum(n()),
                sum(outlier == FALSE)
              ), 
              evenness = ifelse (
                first(State) == "missing",
                mean(E_score, na.rm = TRUE),
                mean(E_score[outlier == FALSE], na.rm = TRUE)
              ), 
              std = ifelse (
                first(State) == "missing",
                sd(E_score, na.rm = TRUE),
                sd(E_score[outlier == FALSE], na.rm = TRUE)
              ),
              outliers = ifelse(
                first(State) == "missing",
                NA,
                sum(outlier == TRUE))
              ) %>%
    mutate(percentage = (count / sum(count)) * 100) %>%
    mutate(cumulative_percentage = cumsum(percentage))
}

create_table_1 <- function(cov_data) {
  subtable_platform <- get_table_1_subtable(cov_data, "SequencingMethod")
  subtable_assembly <- get_table_1_subtable(cov_data, "AssemblyMethod")
  table_1 <- bind_rows(subtable_platform, subtable_assembly) %>%
    select(State, count, 
           percentage, cumulative_percentage, 
           outliers, 
           evenness, std)
  xtab_1 <- xtable(table_1, digits=c(0,0,0,1,1,0,2,2))
  print(
    xtab_1,
    file = "./images/JenkeEtAl2024_Table_1.tex",
    compress = FALSE,
    include.rownames = FALSE
  )
}

# Execute commands
metadata_extraction_all()
cov_data <- get_cov_df()
dir.create("images")
create_table_1(cov_data)

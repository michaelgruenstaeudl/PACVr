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

# TABLE 1 - overview of sequencing platform and assembly software evenness
get_table_1_subtable <- function(cov_data, grouping) {
  subtable <- cov_data %>%
    rename(State = !!ensym(grouping)) %>%
    group_by(State) %>%
    mutate(count = n()) %>%
    ungroup() %>%
    mutate(State := ifelse(count < 5, "missing", State)) %>%
    select(-count) %>%
    group_by(State) %>%
    summarize(count = n(), evenness = mean(E_score), std = sd(E_score)) %>%
    mutate(percentage = (count / sum(count)) * 100) %>%
    mutate(cumulative_percentage = cumsum(percentage))
}

create_table_1 <- function(cov_data) {
  subtable_platform <- get_table_1_subtable(cov_data, "Model")
  subtable_assembly <- get_table_1_subtable(cov_data, "Assembly Method")
  table_1 <- bind_rows(subtable_platform, subtable_assembly) %>%
    select(State, count, percentage, cumulative_percentage, evenness, std)
  xtab_1 <- xtable(table_1)
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

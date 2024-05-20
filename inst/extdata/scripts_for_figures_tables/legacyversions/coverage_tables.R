#!/usr/bin/env RScript
#contributors=c("Gregory Smith", "Nils Jenke", "Michael Gruenstaeudl")
#email="m_gruenstaeudl@fhsu.edu"
#version="2024.05.19.1300"

library(tidyverse)
library(xtable)
library(rstatix)
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
    file = "./images/table_1.tex",
    compress = FALSE,
    include.rownames = FALSE
  )
}

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
    file = "./images/table_2.tex",
    compress = FALSE,
    include.rownames = FALSE
  )
}

# TABLE 3 - summary of Kruskal-Wallis tests
create_table_3 <- function(figure_data) {
  kruskal_results <- get_kruskal_results(figure_data)
  xtab_3 <- xtable(kruskal_results, digits = 3)
  print(
    xtab_3,
    file = "./images/table_3.tex",
    compress = FALSE,
    include.rownames = FALSE
  )
}

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
    file = "./images/table_4.tex",
    compress = FALSE,
    include.rownames = FALSE
  )
}

# SUPPLEMENTARY TABLE 2 - sample metadata
to_string_without_scientific <- function(x) {
  if (is.numeric(x)) {
    return(format(x, scientific = FALSE))
  } else {
    return(x)
  }
}

create_supp_table_2 <- function(cov_data) {
  supp_table_2 <- cov_data %>%
    select(Accession, Model, "Assembly Method", E_score, 
           LSC, IRb, SSC, IRa, coding, noncoding, 
           N_count, IR_mismatches) %>%
    rename(Sample = Accession,
           "Seq. platform" = Model,
           "Asm. software" = "Assembly Method",
           "E-score" = E_score,
           Ns = N_count,
           Mism. = IR_mismatches) %>%
    arrange(Sample) %>%
    mutate_if(is.numeric, round, digits = 4) %>%
    mutate_if(is.numeric, to_string_without_scientific) %>%
    mutate_all(~ifelse(is.na(.), "-", .)) %>%
    mutate_all(~str_replace_all(., "^\\s*NA$", "-")) %>%
    mutate_all(~str_replace_all(., "(?i)ass\\w*$", "Asm."))
  
  supp_xtab_2 <- xtable(supp_table_2)
  print(
    supp_xtab_2,
    file = "./images/supp_table_2.tex",
    compress = FALSE,
    include.rownames = FALSE
  )
}

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
    file = "./images/supp_table_3.tex",
    compress = FALSE,
    include.rownames = FALSE
  )
}

### CREATE ALL TABLES ###
create_tables <- function() {
  cov_data <- get_cov_df()
  figure_data <- get_figure_data(cov_data)

  create_table_1(cov_data)
  create_table_2(cov_data)
  create_table_3(figure_data)
  create_table_4(figure_data)
  create_supp_table_2(cov_data)
  create_supp_table_3(figure_data)
}

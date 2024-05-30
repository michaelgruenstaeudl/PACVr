#!/usr/bin/env RScript
#contributors=c("Gregory Smith", "Nils Jenke", "Michael Gruenstaeudl")
#email="m_gruenstaeudl@fhsu.edu"
#version="2024.05.19.1300"

# pacman::p_load loads packages if they have been installed 
# and installs if they are missing
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, rstatix)

# Is the number of regions with significantly reduced coverage significantly
# different in coding versus non-coding regions?
wilcox_coding_data <- function(cov_data) {
  wilcox_df <- cov_data %>%
    select(coding, noncoding) %>%
    pivot_longer(
      cols = everything(), 
      names_to = "sequences", 
      values_to = "low_coverage") %>%
    arrange(sequences) %>%
    handle_cov_outliers("sequences", "low_coverage")
}

# Is evenness of coverage significantly different across different read lengths?
wilcox_length_data <- function(cov_data) {
  wilcox_df <- cov_data %>%
    select(avgLength, E_score) %>%
    filter(!is.na(avgLength)) %>%
    mutate(Readlength = ifelse(avgLength >= 100 & avgLength <= 300,
                               "short",
                               "medium"))
}

# Is the number of regions with significantly reduced coverage significantly
# different across the four partitions of the plastid genome (i.e., LSC, IR, SSC)?
kruskal_regions_data <- function(cov_data) {
  kruskal_df <- cov_data %>%
    select(LSC, IRb, SSC, IRa) %>%
    pivot_longer(
      cols = everything(), 
      names_to = "regions", 
      values_to = "low_coverage") %>%
    arrange(regions) %>%
    handle_cov_outliers("regions", "low_coverage")
}

# Is evenness of coverage significantly different across different assembly software tools?
kruskal_AssemblyMethod_data <- function(cov_data) {
  kruskal_df <- cov_data %>%
    select("AssemblyMethod", E_score) %>%
    filter(!is.na(AssemblyMethod)) %>%
    handle_cov_outliers("AssemblyMethod", "E_score") %>%
    group_by(AssemblyMethod) %>%
    filter(n() >= 5) %>%
    ungroup()
}

# Is evenness of coverage significantly different across different sequencing forms?
kruskal_SequencingMethod_data <- function(cov_data) {
  kruskal_df <- cov_data %>%
    select(SequencingMethod, E_score) %>%
    filter(!is.na(SequencingMethod)) %>%
    handle_cov_outliers("SequencingMethod", "E_score") %>%
    group_by(SequencingMethod) %>%
    filter(n() >= 5) %>% 
    ungroup()
}

# Collect needed data for figures
get_figure_data <- function(cov_data) {
  figure_data <- list(
    wilcox_coding = wilcox_coding_data(cov_data),
    wilcox_length = wilcox_length_data(cov_data),
    kruskal_regions = kruskal_regions_data(cov_data),
    kruskal_AssemblyMethod = kruskal_AssemblyMethod_data(cov_data),
    kruskal_SequencingMethod = kruskal_SequencingMethod_data(cov_data)
  )
}

# Kruskal-Wallis test results
get_kruskal_results <- function(figure_data) {
  kruskal_df <- figure_data$kruskal_regions
  kruskal_df2 <- figure_data$kruskal_AssemblyMethod
  kruskal_df3 <- figure_data$kruskal_SequencingMethod
  
  kruskal_results <-
    rstatix::kruskal_test(kruskal_df[!is.na(kruskal_df$low_coverage), ], low_coverage ~ regions) %>%
    bind_cols(., select(
      rstatix::kruskal_effsize(kruskal_df, low_coverage ~ regions),
      effsize
    )) %>%
    bind_rows(., bind_cols(
      rstatix::kruskal_test(kruskal_df3[!is.na(kruskal_df3$E_score), ], E_score ~ SequencingMethod),
      select(rstatix::kruskal_effsize(kruskal_df3, E_score ~ SequencingMethod), effsize)
    )) %>%
    bind_rows(., bind_cols(
      rstatix::kruskal_test(kruskal_df2[!is.na(kruskal_df2$E_score), ], E_score ~ AssemblyMethod),
      select(
        rstatix::kruskal_effsize(kruskal_df2, E_score ~ AssemblyMethod),
        effsize
      )
    )) %>%
    mutate(variable = c("regions", "AssemblyMethod", "SequencingMethod")) %>%
    select(variable, n, p, effsize) %>%
    rename(d_c = effsize)
}

# Wilcoxon rank-sum test results
effSize_conv <- function(x) {
  return(2 * x / (sqrt(1 - x ^ 2)))
}

get_wilcox_results <- function(figure_data) {
  wilcox_df <- figure_data$wilcox_coding
  wilcox_df2 <- figure_data$wilcox_length
  
  kruskal_df <- figure_data$kruskal_regions
  kruskal_df2 <- figure_data$kruskal_AssemblyMethod
  kruskal_df3 <- figure_data$kruskal_SequencingMethod
  
  wilcox_results <-
    rstatix::pairwise_wilcox_test(kruskal_df[!is.na(kruskal_df$low_coverage), ], low_coverage ~ regions, p.adjust.method = "BH") %>%
    bind_cols(., select(
      rstatix::wilcox_effsize(kruskal_df, low_coverage ~ regions),
      effsize
    )) %>%
    bind_rows(., bind_cols(
      rstatix::wilcox_test(wilcox_df[!is.na(wilcox_df$low_coverage), ], low_coverage ~ sequences),
      (select(
        rstatix::wilcox_effsize(wilcox_df, low_coverage ~ sequences),
        effsize
      ))
    )) %>%
    bind_rows(., bind_cols(
      rstatix::pairwise_wilcox_test(kruskal_df3[!is.na(kruskal_df3$E_score), ], E_score ~ SequencingMethod, p.adjust.method = "BH"),
      (select(
        rstatix::wilcox_effsize(kruskal_df3, E_score ~ SequencingMethod),
        effsize
      ))
    )) %>%
    bind_rows(., bind_cols(
      rstatix::pairwise_wilcox_test(kruskal_df2[!is.na(kruskal_df2$E_score), ], E_score ~ AssemblyMethod, p.adjust.method = "BH"),
      (select(
        rstatix::wilcox_effsize(kruskal_df2, E_score ~ AssemblyMethod),
        effsize
      ))
    )) %>%
    bind_rows(., bind_cols(
      rstatix::wilcox_test(wilcox_df2[!is.na(wilcox_df2$E_score), ], E_score ~ Readlength),
      (select(
        rstatix::wilcox_effsize(wilcox_df2, E_score ~ Readlength),
        effsize
      ))
    )) %>%
    mutate(effsize = round(effSize_conv(effsize), 3)) %>%
    rename(variable = .y.) %>%
    select(-p.adj.signif)
}

# Outlier filtering functions
handle_cov_outliers <- function(cov_sum, names, values, remove = TRUE) {
  cov_sum <- cov_sum %>%
    group_by(!!sym(names)) %>%
    group_modify(~ handle_outliers(.x, values, 3, remove)) %>%
    ungroup()
  return(cov_sum)
}

handle_outliers <- function(df, col_name, k, remove) {
  bounds <- find_tukey_fences(df, col_name, k)

  col_name_sym <- sym(col_name)
  outlier_expr <- expr(!!col_name_sym >= bounds[1] & !!col_name_sym <= bounds[2])
  if (remove) {
    df <- df %>%
      filter(!!outlier_expr)
  } else {
    df <- df %>%
      mutate(!!col_name_sym := ifelse(!!outlier_expr, !!col_name_sym, NA))
  }
  return(df)
}

find_tukey_fences <- function(df, col_name, k = 1.5) {
  quartile <- quantile(df[[col_name]], probs = c(0.25, 0.75), na.rm = TRUE)
  iqr <- diff(quartile)
  lower_bound <- quartile[1] - k * iqr
  upper_bound <- quartile[2] + k * iqr

  bounds <- c(lower_bound, upper_bound)
  return(bounds)
}

# Adding asterisks to p-values
pval_asterisk <- function(x){
  symbols   = c(" ***"," **"," *","")
  cutoffs = c(0, 0.001, 0.01, 0.05, 1)
  index <- findInterval(x, cutoffs, left.open=T, rightmost.closed=T)
  paste(x, symbols[index], sep="")
}

# Adding letters to effect sizes
effectsize_symbol <- function(x){
  symbols   = c(""," s"," m"," l")
  cutoffs = c(0, 0.2, 0.5, 0.8, 10)
  index <- findInterval(x, cutoffs, left.open=T, rightmost.closed=T)
  paste(x, symbols[index], sep="")
}

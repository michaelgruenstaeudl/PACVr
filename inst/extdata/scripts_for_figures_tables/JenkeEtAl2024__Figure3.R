#!/usr/bin/env RScript
#contributors=c("Gregory Smith", "Nils Jenke", "Michael Gruenstaeudl")
#email="m_gruenstaeudl@fhsu.edu"
#version="2024.05.26.1500"

# pacman::p_load loads packages if they have been installed 
# and installs if they are missing
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tcltk, tidyverse, ggpubr)

# Select input directory and input files
inDir_data = tk_choose.dir(default="~", caption='Select directory with the .gb and .bam files')
inFn_sampleList = tk_choose.files(caption='Select samples list file in csv-format (e.g., input/JenkeEtAl2024_samples_list.csv)')

source("PREP_metadata_extraction_all.R")
source("PREP_coverage_data_assembly.R")
source("PREP_coverage_data_preparation.R")

# FIGURE 3 - regression plot for evenness explained by assembly quality metrics
get_asm_qual_figure <- function(df, x_name, x_lab) {
  x_name_sym <- sym(x_name)
  df_sub <- df %>%
    select(E_score, !!x_name_sym) %>%
    filter(!is.na(!!x_name_sym)) %>%
    mutate(!!x_name_sym := 1 / (!!x_name_sym + 1))
  
  asm_qual_figure <- ggscatter(
    df_sub,
    x = x_name,
    y = "E_score",
    add = "reg.line",
    add.params = list(color = "blue", fill = "lightgray"),
    conf.int = TRUE,
    cor.coef = TRUE,
    cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n"),
    xlab = x_lab,
    ylab = "E-score",
    xlim = c(0, 1.01)
  ) +
    theme(plot.margin = unit(c(0.25, 0.5, 0, 0), "cm")) +
    font("axis.title", size = 10) +
    font("axis.text", size = 10) +
    #theme(axis.text.y = element_blank()) +
    stat_cor(
      method = "spearman",
      p.accuracy = 0.001,
      r.accuracy = 0.01,
      label.y = 1.1,
      label.x = 0.2
    )
}

create_figure_3 <- function(cov_data) {
  n_count_figure <- get_asm_qual_figure(cov_data, 
                                        "N_count", 
                                        "Assembly quality (#Ns)")
  IR_mismatches_figure <- get_asm_qual_figure(cov_data, 
                                              "IR_mismatches", 
                                              "Assembly quality (#Mismatches)")
  ggarrange(
    n_count_figure,
    IR_mismatches_figure,
    nrow = 1,
    ncol = 2,
    labels = c("A", "B"),
    label.y = 1.02,
    label.x = -0.02
  )
  ggsave(
    filename = "JenkeEtAl2024_Figure_3.pdf",
    path = "./images",
    device = 'pdf',
    dpi = 700,
    width = 15.5,
    height = 7.66,
    units = ("cm")
  )
}

# Execute commands
metadata_extraction_all()
cov_data <- get_cov_df()
dir.create("images")
create_figure_3(cov_data)

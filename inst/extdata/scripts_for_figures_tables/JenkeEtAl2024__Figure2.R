#!/usr/bin/env RScript
#contributors=c("Gregory Smith", "Nils Jenke", "Michael Gruenstaeudl")
#email="m_gruenstaeudl@fhsu.edu"
#version="2024.05.19.1300"

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

# FIGURE 2 - E-score and WRSD distributions by groupings
get_wrsd_figure <- function(df, grouping) {
  wrsd_label <- "WRSD per kb"
  y_axis_max_lim <- 0.005
  
  wrsd_figure <- ggpubr::ggboxplot(
    df,
    x = grouping,
    y = "low_coverage",
    color = grouping,
    add = "jitter",
    bxp.errorbar = TRUE,
    legend = "none",
    ylab = wrsd_label
  ) +
    rremove("xlab") +
    theme(plot.margin = unit(c(0.25, 0.5, 0.5, 0), "cm"))
  ggpar(wrsd_figure, ylim = c(0, y_axis_max_lim))
  wrsd_figure <- wrsd_figure +
    font("axis.text", size = 10) +
    font("ylab", size = 10)
}

get_escore_figure <- function(df, grouping) {
  escore_figure <- ggpubr::ggboxplot(
    df,
    x = grouping,
    y = "E_score",
    color = grouping,
    add = "jitter",
    legend = "none",
    bxp.errorbar = TRUE,
    ylab = "E-score"
  ) +
    rotate_x_text(25) +
    rremove("xlab") +
    theme(plot.margin = unit(c(0.25, 0.5, 0, 0), "cm"))
  ggpar(escore_figure, ylim = c(0, 1.0))
  escore_figure <- escore_figure +
    font("axis.text", size = 10) +
    font("ylab", size = 10)
}

create_figure_2 <- function(figure_data) {
  coding_figure <- get_wrsd_figure(figure_data$wilcox_coding, "sequences")
  partition_figure <- get_wrsd_figure(figure_data$kruskal_regions, "regions")
  assembly_figure <- get_escore_figure(figure_data$kruskal_assembly, "Assembly")
  model_figure <- get_escore_figure(figure_data$kruskal_model, "SequencingMethod")
  
  ggarrange(
    partition_figure,
    coding_figure,
    model_figure,
    assembly_figure,
    nrow = 2,
    ncol = 2,
    labels = c("A", "B", "C", "D"),
    label.y = 1.015,
    label.x = -0.025,
    heights = c(1, 1.12, 0.99)
  )
  ggsave(
    filename = "JenkeEtAl2024_Figure_2.pdf",
    path = "./images",
    device = 'pdf',
    dpi = 700,
    width = 15.50,
    height = 20.75,
    units = ("cm")
  )
}

# Execute commands
metadata_extraction_all()
cov_data <- get_cov_df()
figure_data <- get_figure_data(cov_data)
dir.create("images")
create_figure_2(figure_data)

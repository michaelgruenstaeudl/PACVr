#!/usr/bin/env RScript
#contributors=c("Gregory Smith", "Nils Jenke", "Michael Gruenstaeudl")
#email="m_gruenstaeudl@fhsu.edu"
#version="2024.05.19.1300"

# pacman::p_load loads packages if they have been installed 
# and installs if they are missing
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tcltk, tidyverse, ggpubr, optparse)

# Select input directory and input files
inDir_data = tk_choose.dir(default="~", caption='Select directory with the .gb and .bam files')
inFn_sampleList = tk_choose.files(caption='Select samples list file in csv-format (e.g., input/JenkeEtAl2024_samples_list.csv)')

source("PREP_metadata_extraction_all.R")
source("PREP_coverage_data_assembly.R")
source("PREP_coverage_data_preparation.R")

# FIGURE 1A - evenness outliers
create_figure_1a <- function(cov_data) {
  outlier <- ggpubr::ggboxplot(
    cov_data,
    x = NULL,
    y = "E_score",
    add = "jitter",
    bxp.errorbar = TRUE,
    label = "Accession",
    font.label = list(size = 10, color = "red"),
    label.select = list(top.down = 13),
    repel = TRUE,
    xlab = "Samples",
    ylab = "E-score"
  ) +
    rremove("x.text") +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    font("y.text", size = 10) +
    font("xylab", size = 10)
  outlier$layers[[3]]$aes_params$size <- 0.5
  ggsave(
    filename = "JenkeEtAl2024_Figure_1A.pdf",
    path = "./images",
    device = 'pdf',
    dpi = 700,
    width = 7.748895,
    height = 7.700765,
    units = ("cm")
  )
}

# Execute commands
metadata_extraction_all()
cov_data <- get_cov_df()
dir.create("images")
create_figure_1a(cov_data)

#!/usr/bin/env RScript
#contributors=c("Gregory Smith", "Nils Jenke", "Michael Gruenstaeudl")
#email="m_gruenstaeudl@fhsu.edu"
#version="2024.05.19.1300"

# pacman::p_load loads packages if they have been installed 
# and installs if they are missing
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tcltk, tidyverse, ggplot2, ggrepel, ggpubr, optparse, png, grid)

# Select input directory and input files
inDir_data = tk_choose.dir(default="~", caption='Select directory with the .gb and .bam files')
inFn_sampleList = tk_choose.files(caption='Select samples list file in csv-format (e.g., input/JenkeEtAl2024_samples_list.csv)')
inFn_figure1b = tk_choose.files(caption='Select image file used for figure 1b (e.g., input/NC_026562_CoverageViz_crop.png)')

source("PREP_metadata_extraction_all.R")
source("PREP_coverage_data_assembly.R")
source("PREP_coverage_data_preparation.R")

# FIGURE 1A - evenness outliers
create_figure_1a <- function(cov_data) {
  # manually jitter values
  set.seed(2)
  cov_data <- cov_data %>%
    mutate(jittered_x = jitter(rep(1, n()), amount = 0.2))
  top_outliers <- cov_data %>% 
    slice_min(E_score, n = 13)
  
  # create ggplot
  outlier <- ggplot(cov_data, aes(x = 1, y = E_score)) +
    geom_boxplot(outlier.shape = NA, width = 0.5, height) +
    geom_point(aes(x = jittered_x), alpha = 0.5, size = 0.75) +
    geom_text_repel(
      data = top_outliers,
      aes(x = jittered_x, y = E_score, label = Accession),
      size = 2,
      color = "red",
      nudge_x = 0.1,
      nudge_y = -0.1
    ) +
    theme_minimal() +
    theme(
      plot.margin = unit(c(1, 1, 1, 1), "cm"),
      text = element_text(size = 10)
    ) +
    labs(x = "Samples", y = "E-score") +
    rremove("x.text")
}

# FIGURE 1B - previously created viz of outlier
create_figure_1b <- function() {
  img <- readPNG(inFn_figure1b)
  img_grob <- rasterGrob(img, width = unit(1, "npc"), height = unit(1, "npc"))
  png_plot <- ggplot() +
    annotation_custom(img_grob) +
    theme_void()
}

# FIGURE 1 - arranged display of evenness outliers and an outlier example
create_figure_1 <- function(cov_data) {
  figure_1a <- create_figure_1a(cov_data)
  figure_1b <- create_figure_1b()

  ggarrange(
    figure_1a,
    figure_1b,
    nrow = 1,
    ncol = 2,
    labels = c("A", "B"),
    label.y = 1.02,
    label.x = -0.02
  )
  ggsave(
    filename = "JenkeEtAl2024_Figure_1.pdf",
    path = "./images",
    device = 'pdf',
    dpi = 700,
    width = 18,
    height = 7.700765,
    units = ("cm")
  )
}

# Execute commands
metadata_extraction_all()
cov_data <- get_cov_df()
dir.create("images")
create_figure_1(cov_data)

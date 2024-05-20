#!/usr/bin/env RScript
#contributors=c("Gregory Smith", "Nils Jenke", "Michael Gruenstaeudl")
#email="m_gruenstaeudl@fhsu.edu"
#version="2024.05.19.1300"

library(tidyverse)
library(ggpubr)
source("coverage_data_assembly.R")
source("coverage_data_preparation.R")

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
    xlab = "Samples"
  ) +
    rremove("x.text") +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    font("y.text", size = 10) +
    font("xylab", size = 10)
  outlier$layers[[3]]$aes_params$size <- 0.5
  ggsave(
    filename = "figure_1a.pdf",
    path = "./images",
    device = 'pdf',
    dpi = 700,
    width = 7.748895,
    height = 7.700765,
    units = ("cm")
  )
}

# FIGURE 2 - E-score and WRSD distributions by groupings
get_wrsd_figure <- function(df, grouping) {
  wrsd_label <- "WRSD-Length Ratio"
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
    theme(plot.margin = unit(c(0, 0, 0.5, 0), "cm"))
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
    theme(plot.margin = unit(c(0, 0, 0.5, 0), "cm"))
  ggpar(escore_figure, ylim = c(0, 1.0))
  escore_figure <- escore_figure +
    font("axis.text", size = 10) +
    font("ylab", size = 10)
}

create_figure_2 <- function(figure_data) {
  coding_figure <- get_wrsd_figure(figure_data$wilcox_coding, "sequences")
  partition_figure <- get_wrsd_figure(figure_data$kruskal_regions, "regions")
  assembly_figure <- get_escore_figure(figure_data$kruskal_assembly, "Assembly")
  model_figure <- get_escore_figure(figure_data$kruskal_model, "Model")
  
  ggarrange(
    partition_figure,
    coding_figure,
    model_figure,
    assembly_figure,
    nrow = 2,
    ncol = 2,
    labels = c("A", "B", "C", "D"),
    label.y = 1.02,
    label.x = -0.02,
    heights = c(1, 1.1156625968355375019018071372623, 0.99)
  )
  ggsave(
    filename = "figure_2.pdf",
    path = "./images",
    device = 'pdf',
    dpi = 700,
    width = 15.49779,
    height = 20.75,
    units = ("cm")
  )
}

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
    theme(plot.margin = unit(c(0.05, 0, 0, 0), "cm")) +
    font("axis.title", size = 10) +
    font("axis.text", size = 10) +
    theme(axis.text.y = element_blank()) +
    stat_cor(
      method = "spearman",
      p.accuracy = 0.001,
      r.accuracy = 0.01,
      label.y = 1.025,
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
    filename = "figure_3.pdf",
    path = "./images",
    device = 'pdf',
    dpi = 700,
    width = 15.49779,
    height = 7.654888,
    units = ("cm")
  )
}

### CREATE ALL FIGURES ###
create_figures <- function() {
  cov_data <- get_cov_df()
  figure_data <- get_figure_data(cov_data)

  create_figure_1a(cov_data)
  create_figure_2(figure_data)
  create_figure_3(cov_data)
}

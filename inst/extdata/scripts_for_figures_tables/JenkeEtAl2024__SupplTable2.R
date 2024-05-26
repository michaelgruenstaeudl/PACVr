#!/usr/bin/env RScript
#contributors=c("Gregory Smith", "Nils Jenke", "Michael Gruenstaeudl")
#email="m_gruenstaeudl@fhsu.edu"
#version="2024.05.19.1300"

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
    select(Accession, SequencingMethod, AssemblyMethod, E_score, 
           LSC, IRb, SSC, IRa, coding, noncoding, 
           N_count, IR_mismatches) %>%
    rename(Sample = Accession,
           "Seq. platform" = SequencingMethod,
           "Asm. software" = AssemblyMethod,
           "E-score" = E_score,
           Ns = N_count,
           Mism. = IR_mismatches) %>%
    arrange(Sample) %>%
    mutate_if(is.numeric, round, digits = 2) %>%
    mutate_if(is.numeric, to_string_without_scientific) %>%
    mutate_all(~ifelse(is.na(.), "-", .)) %>%
    mutate_all(~str_replace_all(., "^\\s*NA$", "-")) %>%
    mutate_all(~str_replace_all(., "(?i)ass\\w*$", "Asm."))
  
  supp_xtab_2 <- xtable(supp_table_2)
  print(
    supp_xtab_2,
    file = "./images/JenkeEtAl2024_SupplTable_2.tex",
    compress = FALSE,
    include.rownames = FALSE
  )
}

# Execute commands
metadata_extraction_all()
cov_data <- get_cov_df()
dir.create("images")
create_supp_table_2(cov_data)

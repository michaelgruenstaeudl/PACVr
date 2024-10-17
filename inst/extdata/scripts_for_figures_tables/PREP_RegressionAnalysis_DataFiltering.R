#!/usr/bin/env RScript
#contributors=c("Gregory Smith", "Michael Gruenstaeudl")
#email="g_smith10@mail.fhsu.edu"
#version="2024.10.16.1800"

########################################################################
#               DEPENDENCIES, DEFINING INPUT/OUTPUT                    #
########################################################################

if (!require("pacman")) install.packages("pacman")
pacman::p_load(tcltk)

# Select input directory and input files
inDir_data <- tk_choose.dir(default="~", caption='Select the directory that contains the .gb and .bam files of this analysis')
inFn_sampleList <- tk_choose.files(caption='Select the sample list in csv-format (e.g., input/JenkeEtAl2024_samples_list.csv)')
outDir_fn <- paste(getwd(), '/data_for_regression_analyses_', format(Sys.time(), "%Y_%m_%d"),'.rds', sep='')

# Find and load script dependency
assembly_name <- "PREP_coverage_data_assembly.R"
assembly_path <- list.files(path = getwd(), pattern = assembly_name, recursive = TRUE, full.names = TRUE)
preparation_name <- "PREP_coverage_data_preparation.R"
preparation_path <- list.files(path = getwd(), pattern = preparation_name, recursive = TRUE, full.names = TRUE)

# tryCatch necessary if there are multiple R files with the same name on computer
tryCatch(
    expr = {source(assembly_path)},
	error = function(e){source(assembly_path[[grep("git", assembly_path)]])}
)
tryCatch(
    expr = {source(preparation_path)},
	error = function(e){source(preparation_path[[grep("git", preparation_path)]])}
)

########################################################################
#                             FUNCTIONS                                #
########################################################################

get_long_cov_df <- function() {
  metadata_df <- get_metadata_df(FALSE)

  noncoding_pattern <- ".[0-9]_coverage.summary.noncoding.tsv$"
  noncoding_df <- create_cov_sum_df(noncoding_pattern) %>%
    mutate(RegionSubset = "noncoding")

  coding_pattern <- ".[0-9]_coverage.summary.genes.tsv$"
  coding_df <- create_cov_sum_df(coding_pattern) %>%
    mutate(RegionSubset = "coding")

  regions_pattern <- ".[0-9]_summary.regions.tsv$"
  regions_df <- create_cov_sum_df(regions_pattern) %>%
    mutate(RegionSubset = "all") %>%
    mutate(Chromosome = replace(Chromosome, Chromosome == "Complete_genome", "Unpartitioned"))

  coverage_df <- bind_rows(list(noncoding_df, coding_df, regions_df)) %>%
    inner_join(metadata_df, by = "Accession") %>%
    mutate(Chromosome = relevel(as.factor(Chromosome), ref = "Unpartitioned")) %>%
    mutate(Accession = as.factor(Accession)) %>%
    mutate(RegionSubset = relevel(as.factor(RegionSubset), ref = "all")) %>%
    mutate(AssemblyMethod = as.factor(AssemblyMethod)) %>%
    mutate(SequencingMethod = as.factor(SequencingMethod)) %>%
    mutate(IR_mismatches = as.integer(IR_mismatches)) %>%
    rename(lowCovWin = lowCovWin_abs) %>%
    mutate(N_count = 1 / (N_count + 1)) %>%
    mutate(IR_mismatches = 1/ (IR_mismatches + 1))
}

########################################################################
#                                MAIN                                  #
########################################################################

## Get all coverage data including available metadata
cov_df <- get_long_cov_df()

## Save to R object
saveRDS(cov_df, file=outDir_fn)

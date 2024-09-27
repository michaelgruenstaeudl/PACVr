if (!require("pacman")) install.packages("pacman")
pacman::p_load(tcltk, tidyverse, tidymodels, vip, rpart.plot)

# Select input directory and input files
inDir_data <- tk_choose.dir(default="~", caption='Select directory with the .gb and .bam files')
inFn_sampleList <- tk_choose.files(caption='Select samples list file in csv-format (e.g., input/JenkeEtAl2024_samples_list.csv)')

# Find and load script dependency
assembly_name <- "PREP_coverage_data_assembly.R"
assembly_path <- list.files(path = getwd(), pattern = assembly_name, recursive = TRUE, full.names = TRUE)
preparation_name <- "PREP_coverage_data_preparation.R"
preparation_path <- list.files(path = getwd(), pattern = preparation_name, recursive = TRUE, full.names = TRUE)
source(assembly_path)
source(preparation_path)


### FUNCTIONS ###
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

filter_most_metadata <- function(df) {
  df %>% select(-TITLE, -SRA, -Samples, -Accession)
}

filter_for_coverage <- function(df) {
  df %>% filter_most_metadata() %>%
    select(-N_count, -IR_mismatches, -E_score, -AssemblyMethod, -SequencingMethod, -lowCovWin_perKilobase)
}

filter_region_x_subset <- function(df) {
  df %>%
    filter_for_coverage() %>%
    filter(!is.na(Chromosome)) %>%
    filter(!str_detect(Chromosome, "Junction")) %>%
    filter(Chromosome != "Unpartitioned") %>%
    filter(RegionSubset != "all")
}

filter_complete_subset <- function(df) {
  df %>%
    filter_for_coverage() %>%
    filter(Chromosome == "Unpartitioned") %>%
    select(-Chromosome) %>%
    filter(RegionSubset != "all")
}

filter_complete_region <- function(df) {
  df %>%
    filter_for_coverage() %>%
    filter(RegionSubset == "all") %>%
    select(-RegionSubset) %>%
    filter(Chromosome != "Unpartitioned")
}

filter_complete_genome <- function(df) {
  df %>%
    filter_most_metadata() %>%
    select(-lowCovWin_perKilobase) %>%
    filter(Chromosome == "Unpartitioned") %>%
    filter(RegionSubset == "all") %>%
    select(-Chromosome, -RegionSubset)
}

filter_small_method_classes <- function(df) {
  df %>%
    group_by(AssemblyMethod) %>%
    filter(n() >= 5) %>%
    ungroup() %>%
    group_by(SequencingMethod) %>%
    filter(n() >= 5) %>%
    ungroup()
}

split_train_test <- function(df, prop = 0.8) {
  set.seed(42)
  initial_split(
    df %>%
      drop_na(),
    prop = prop
  )
}

get_tree_fit <- function(df, response_name, grid_levels = 5) {
  # The below code is broadly borrowed from
  # https://www.tidymodels.org/start/tuning/

  tree_model <- decision_tree(
    cost_complexity = tune(),
    min_n = tune(),
    tree_depth = tune()
  ) %>%
    set_engine("rpart") %>%
    set_mode("regression")

  tree_grid <- grid_regular(
    cost_complexity(),
    min_n(),
    tree_depth(),
    levels = grid_levels
  )

  split <- df %>%
    split_train_test()
  set.seed(321)
  folds <- split %>%
    training() %>%
    vfold_cv(v = 10)

  formula <- as.formula(paste(response_name, "~ ."))
  tree_wf <- workflow() %>%
    add_model(tree_model) %>%
    add_formula(formula)

  set.seed(123)
  tree_res <- tree_wf %>%
    tune_grid(resamples = folds, grid = tree_grid)

  best_tree <- tree_res %>%
    select_best(metric = "rmse")

  final_wf <- tree_wf %>%
    finalize_workflow(best_tree)

  final_fit <- final_wf %>%
    last_fit(split)
}


## Get all coverage data including available metadata
cov_df <- get_long_cov_df()


### COVERAGE ANALYSIS ###

# The data that looks at only region or subset (coding-noncoding) coverage
# is available by `filter_complete_region()` and `filter_complete_subset()`,
# but conceptually provides less detail as data by region-subset,
# and from preliminary models does not appear to result in statistically significant differences.
# Instead we will directly proceed to region-subset analysis.

## Mutually disjoint region-subset coverage data
reg_sub_df <- cov_df %>%
  filter_region_x_subset()

## Preliminary linear regression model
reg_sub_lm_fit <- linear_reg() %>%
  fit(lowCovWin ~ . + Chromosome * RegionSubset, data = reg_sub_df)
tidy(reg_sub_lm_fit)
summary(reg_sub_lm_fit$fit)

# When controlling for both region and coding-noncoding subset,
# including an interaction term for the two, we generally see a statistically significant differences
# in these two dimensions.
# Specifically, compared to the reference ChromosomeIRa:RegionSubsetcoding,
# the other region-subsets are statistically different, except for IRb,
# which is not significally different in either coding or noncoding.
# We could perform additional diagnostics tests to iterate on the model to be more statistically valid,
# or regularize the model with something like LASSO or ridge regression, but for now let's move onto
# other model types as the motivation for the hypothesis tests seems supported.

## Decision tree tuned to minimize RMSE on a training subset of `reg_sub_df`
reg_sub_tree_fit <- get_tree_fit(reg_sub_df, "lowCovWin")

# Performance on testing subset
reg_sub_tree_fit %>%
  collect_metrics()

# Structure of tree and importance of predictors
reg_sub_tree <- extract_workflow(reg_sub_tree_fit)
reg_sub_tree
reg_sub_tree %>%
  extract_fit_engine() %>%
  rpart.plot(roundint = FALSE)
reg_sub_tree %>%
  extract_fit_parsnip() %>%
  vip()

# We see that the resulting model has an incredibly similar explanatory power (R^2) to the linear regression,
# and a slightly lower RMSE. The importance of the variables in the decision tree are, in order from greatest
# to least: region length, region, region subset, and average sample length. Inspecting the final tree created,
# we notice that while both coding-noncoding are important for prediction, only the SSC region was utilized.
# This leads up to believe that although the effect of region on coverage is not as great as subset type,
# both are important in reducing variance.

## Decision tree for just region and region subset
rs_2_sub_df <- reg_sub_df %>%
  handle_outliers("lowCovWin", 3, 0) %>%
  select(lowCovWin, Chromosome, RegionSubset)

rs_2_tree_fit <- get_tree_fit(rs_2_sub_df, "lowCovWin")

rs_2_tree_fit %>%
  collect_metrics()
rs_2_tree <- extract_workflow(rs_2_tree_fit)
rs_2_tree
rs_2_tree %>%
  extract_fit_engine() %>%
  rpart.plot(roundint = FALSE)
rs_2_tree %>%
  extract_fit_parsnip() %>%
  vip()

# Focusing on region and region subset, in addition to removing outliers, supports the results from above.
# That is, that the coding-noncoding split and the SSC region are most important for predicting
# coverage, but we get a broader view of how region is used beyond just SSC.


## EVENNESS ANALYSIS ##

# Linear regression models were attempted on the evenness data, with no statistical
# significance for any of the variables of interest. We instead continue with non-linear
# models that more likely represent the relationship between these variables and evenness,
# if such a relationship exists.

## Examine complete genomes for evenness
genome_df <- cov_df %>%
  filter_complete_genome()

genome_filter_df <- genome_df %>%
  filter_small_method_classes() %>%
  handle_outliers("E_score", 3, 0) %>%
  select(-lowCovWin, -AssemblyMethod)

genome_tree_fit <- get_tree_fit(genome_filter_df, "E_score")
genome_tree_fit %>%
  collect_metrics()
genome_tree <- extract_workflow(genome_tree_fit)
genome_tree
genome_tree %>%
  extract_fit_engine() %>%
  rpart.plot(roundint = FALSE)
genome_tree %>%
  extract_fit_parsnip() %>%
  vip()

# To prevent removing too many observations we remove the assembly method from this
# model, and to prevent over-reliance on the highly correlated low coverage count this is also removed.
# Of the tree variables of interest that are present here, only sequencing method appears to have importance.
# We will continue will single variable decision tree models for all four variables of interest.

## Assembly method decision tree
assembly_df <- genome_df %>%
  filter_small_method_classes() %>%
  handle_cov_outliers("AssemblyMethod", "E_score") %>%
  select(E_score, AssemblyMethod)

assembly_tree_fit <- get_tree_fit(assembly_df, "E_score")
assembly_tree_fit %>% collect_metrics()
assembly_tree <- extract_workflow(assembly_tree_fit)
assembly_tree
assembly_tree %>%
  extract_fit_engine() %>%
  rpart.plot(roundint = FALSE)
assembly_tree %>%
  extract_fit_parsnip() %>%
  vip()

# When considering only assembly methods with 5 or more observations and
# removing in-class outliers, there does appear a very minor importance
# on assembly method, especially given the low R^2. If we elect to not
# remove outliers or remove outliers based on the entire data set,
# the assembly method is not determined to have any importance.
# If we keep in the assembly methods with low counts, this artificially increases
# the importance of the method and has a major impact on the explanatory value (R^2).

## Sequencing method decision tree
sequencing_df <- genome_df %>%
  filter_small_method_classes() %>%
  handle_cov_outliers("SequencingMethod", "E_score") %>%
  select(E_score, SequencingMethod)

sequencing_tree_fit <- get_tree_fit(sequencing_df, "E_score")
sequencing_tree_fit %>% collect_metrics()
sequencing_tree <- extract_workflow(sequencing_tree_fit)
sequencing_tree
sequencing_tree %>%
  extract_fit_engine() %>%
  rpart.plot(roundint = FALSE)
sequencing_tree %>%
  extract_fit_parsnip() %>%
  vip()

# When considering only sequencing methods with 5 or more observations and
# removing in-class outliers, there does appear a moderate importance
# on assembly method, with the moderate R^2 supporting this.

## Ambiguous count

# A linear regression model is already being used for the analysis.
# Applying a decision tree to this metric (as a count or normalized)
# finds no importance of this variable.

# IR mismatches

# A linear regression model is already being used for the analysis.
# Applying a decision tree to this metric (as a count or normalized)
# finds no importance of this variable.

## BONUS - Read length decision tree
length_df <- genome_df %>%
  filter_small_method_classes() %>%
  handle_outliers("E_score", 3, 0) %>%
  select(E_score, avgLength)

length_tree_fit <- get_tree_fit(length_df, "E_score")
length_tree_fit %>% collect_metrics()
length_tree <- extract_workflow(length_tree_fit)
length_tree
length_tree %>%
  extract_fit_engine() %>%
  rpart.plot(roundint = FALSE)
length_tree %>%
  extract_fit_parsnip() %>%
  vip()

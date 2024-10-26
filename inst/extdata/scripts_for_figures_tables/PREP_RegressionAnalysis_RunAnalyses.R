########################################################################
#                DEPENDENCIES, DEFINING INPUT/OUTPUT                   #
########################################################################

if (!require("pacman")) install.packages("pacman")
pacman::p_load(svglite, tcltk, tidyverse, tidymodels, vip, rpart.plot)

inData_fn <- tk_choose.files(default="~", caption='Select the RDS file generated via script `PREP_RegressionAnalysis_DataFiltering`', multi=FALSE)
cov_df <- readRDS(file=inData_fn)

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
#                            FUNCTIONS                                 #
########################################################################

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

plot_regression_tree <- function(tree, descr) {
    descr = gsub(" ", "_", descr)
    svglite(paste0("DecisionTree__Effects_of_", descr,".svg"), width=8, height=4)
    tree %>%
        extract_fit_engine() %>%
        rpart.plot(roundint = FALSE)
    dev.off()
}

plot_variable_importance <- function(tree, descr) {
    descr = gsub(" ", "_", descr)
    my_plot <- tree %>%
                extract_fit_parsnip() %>%
                vip(aesthetics = list(alpha=0.5))
    my_plot <- my_plot + theme_bw()  + 
        theme(
          axis.text.x=element_text(size=6),
          axis.text.y=element_text(size=6),
          axis.title.x=element_text(size=8),
          axis.title.y=element_text(size=8)
          )
    ggsave(filename=paste0("VariableImportance__Effects_of_", descr,".svg"), 
	       plot=my_plot, device='svg', dpi=300, width=8, height=4, units=("cm"))
}

########################################################################
#                               MAIN                                   #
########################################################################

#########################
### COVERAGE ANALYSIS ###
#########################

# The coverage data subdivided by only quadripartite genome region or coding/noncoding division 
# is available by `filter_complete_region()` and `filter_complete_subset()`,
# but conceptually provides less detail as data by quadripartite region-coding/noncoding division,
# and from preliminary models does not appear to result in statistically significant differences.
# Instead we will directly proceed to region-subset analysis.

## Mutually disjoint quadripartite region-coding/noncoding division coverage data
reg_sub_df <- cov_df %>%
  filter_region_x_subset()

## Preliminary LINEAR REGRESSION MODEL
reg_sub_lm_fit <- linear_reg() %>%
  fit(lowCovWin ~ . + Chromosome * RegionSubset, data = reg_sub_df)
tidy(reg_sub_lm_fit)
summary(reg_sub_lm_fit$fit)

# When controlling for both region and coding/noncoding division,
# including an interaction term for the two, we generally see statistically significant differences
# in these two dimensions.
# Specifically, compared to the reference ChromosomeIRa:RegionSubsetcoding,
# the other region-subsets are statistically different, except for IRb,
# which is not significally different in either coding or noncoding.
# We could perform additional diagnostics tests to iterate on the model to be more statistically valid,
# or regularize the model with something like LASSO or ridge regression, but for now let's move onto
# other model types as the motivation for the hypothesis tests seems supported.

########################################################################
## Effects of the QUADRIPARTITE REGION and CODING-NONCODING DIVISION on SEQ COVERAGE
########################################################################
descr = "QUADRIPARTITE REGION and CODING-NONCODING DIVISION on SEQ COVERAGE"

## DECISION TREE tuned to minimize RMSE on a training subset of `reg_sub_df`
reg_sub_tree_fit <- get_tree_fit(reg_sub_df, "lowCovWin")

# Performance on testing subset
reg_sub_tree_fit %>% collect_metrics()

# Structure of tree and importance of predictors
reg_sub_tree <- extract_workflow(reg_sub_tree_fit)
reg_sub_tree

# Decision-Tree plot
plot_regression_tree(reg_sub_tree, descr)

# Variable-Importance plot: plotting importance scores for the model predictors
plot_variable_importance(reg_sub_tree, descr)

# RESULTS:
# We see that the resulting model has an incredibly similar explanatory power (R^2) to the linear regression,
# and a slightly lower RMSE. The importance of the variables in the decision tree are, in order from greatest
# to least: region length, region, region subset, and average sample length. Inspecting the final tree created,
# we notice that while both coding-noncoding are important for prediction, only the SSC region was utilized.
# This leads up to believe that although the effect of region on coverage is not as great as subset type,
# both are important in reducing variance.


########################################################################
## Effects of only the QUADRIPARTITE REGION on SEQ COVERAGE
########################################################################
descr = "only QUADRIPARTITE REGION on SEQ COVERAGE"

rs_2_sub_df <- reg_sub_df %>%
  handle_outliers("lowCovWin", 3, 0) %>%
  select(lowCovWin, Chromosome, RegionSubset)

rs_2_tree_fit <- get_tree_fit(rs_2_sub_df, "lowCovWin")

rs_2_tree_fit %>% collect_metrics()
rs_2_tree <- extract_workflow(rs_2_tree_fit)
rs_2_tree

# Decision-Tree plot
plot_regression_tree(rs_2_tree, descr)

# Variable-Importance plot: plotting importance scores for the model predictors
plot_variable_importance(rs_2_tree, descr)

# RESULTS:
# Focusing on region and region subset, in addition to removing outliers, supports the results from above.
# That is, that the coding-noncoding split and the SSC region are most important for predicting
# coverage, but we get a broader view of how region is used beyond just SSC.

########################################################################


#######################
## EVENNESS ANALYSIS ##
#######################

# Linear regression models were attempted on the evenness data, with no statistical
# significance for any of the variables of interest. 

# We instead continue with NON-LINEAR MODELS that more likely represent 
# the relationship between these variables and evenness, if such a 
# relationship exists.

########################################################################
## Effects of FOUR VARIABLES (AMBIGUOUS NUCL, IR MISMATCHES, 
## SEQ PLATFORM, and ASSEMBLY METHOD) on SEQ EVENNESS
########################################################################
descr = "FOUR VARIABLES on SEQ EVENNESS"

genome_df <- cov_df %>%
  filter_complete_genome()

genome_filter_df <- genome_df %>%
  filter_small_method_classes() %>%
  handle_outliers("E_score", 3, 0) %>%
  select(-lowCovWin, -AssemblyMethod)

genome_tree_fit <- get_tree_fit(genome_filter_df, "E_score")
genome_tree_fit %>% collect_metrics()
genome_tree <- extract_workflow(genome_tree_fit)
genome_tree

# Decision-Tree plot
plot_regression_tree(genome_tree, descr)

# Variable-Importance plot: plotting importance scores for the model predictors
plot_variable_importance(genome_tree, descr)

# RESULTS:
# To prevent removing too many observations we remove the assembly method from this
# model, and to prevent over-reliance on the highly correlated low coverage count this is also removed.
# Of the tree variables of interest that are present here, only sequencing method appears to have importance.
# We will continue will single variable decision tree models for all four variables of interest.


########################################################################
## Effects of only ASSEMBLY METHOD on SEQ EVENNESS
########################################################################
descr = "only ASSEMBLY METHOD on SEQ EVENNESS"

asmMethod_df <- genome_df %>%
  filter_small_method_classes() %>%
  handle_cov_outliers("AssemblyMethod", "E_score") %>%
  select(E_score, AssemblyMethod)

asmMethod_tree_fit <- get_tree_fit(asmMethod_df, "E_score")
asmMethod_tree_fit %>% collect_metrics()
asmMethod_tree <- extract_workflow(asmMethod_tree_fit)
asmMethod_tree

# Decision-Tree plot
plot_regression_tree(asmMethod_tree, descr)

# Variable-Importance plot: plotting importance scores for the model predictors
plot_variable_importance(asmMethod_tree, descr)

# RESULTS:
# When considering only assembly methods with 5 or more observations and
# removing in-class outliers, there does appear a very minor importance
# on assembly method, especially given the low R^2. If we elect to not
# remove outliers or remove outliers based on the entire data set,
# the assembly method is not determined to have any importance.
# If we keep in the assembly methods with low counts, this artificially increases
# the importance of the method and has a major impact on the explanatory value (R^2).


########################################################################
## Effects of only SEQ PLATFORM on SEQ EVENNESS
########################################################################
descr = "only SEQ PLATFORM on SEQ EVENNESS"

seqMethod_df <- genome_df %>%
  filter_small_method_classes() %>%
  handle_cov_outliers("SequencingMethod", "E_score") %>%
  select(E_score, SequencingMethod)

seqMethod_tree_fit <- get_tree_fit(seqMethod_df, "E_score")
seqMethod_tree_fit %>% collect_metrics()
seqMethod_tree <- extract_workflow(seqMethod_tree_fit)
seqMethod_tree

# Decision-Tree plot
plot_regression_tree(seqMethod_tree, descr)

# Variable-Importance plot: plotting importance scores for the model predictors
plot_variable_importance(seqMethod_tree, descr)

# RESULTS:
# When considering only sequencing methods with 5 or more observations and
# removing in-class outliers, there does appear a moderate importance
# on assembly method, with the moderate R^2 supporting this.


########################################################################
## Effects of only AMBIGUOUS NUCL on SEQ EVENNESS
########################################################################
descr = "only AMBIGUOUS NUCL on SEQ EVENNESS"

ambigNucl_df <- genome_df %>%
  filter_small_method_classes() %>%
  handle_outliers("E_score", 3, 0) %>%
  select(E_score, N_count)

ambigNucl_lm_fit <- linear_reg() %>%
  fit(E_score ~ ., data = ambigNucl_df)
tidy(ambigNucl_lm_fit)
summary(ambigNucl_lm_fit$fit)

ambigNucl_tree_fit <- get_tree_fit(ambigNucl_df, "E_score")
ambigNucl_tree_fit %>% collect_metrics()
ambigNucl_tree <- extract_workflow(ambigNucl_tree_fit)
ambigNucl_tree

# Decision-Tree plot
plot_regression_tree(ambigNucl_tree, descr)

# Variable-Importance plot: plotting importance scores for the model predictors
plot_variable_importance(ambigNucl_tree, descr)

# RESULTS:
# The linear regression result supports the previous Pearson metric that demonstrated
# statistical significance for ambiguous count's impact on evenness. The decision tree
# additionally supports this with importance given to the count.


########################################################################
## Effects of only IR MISMATCHES on SEQ EVENNESS
########################################################################
descr = "only IR MISMATCHES on SEQ EVENNESS"

IRmism_df <- genome_df %>%
  filter_small_method_classes() %>%
  handle_outliers("E_score", 3, 0) %>%
  select(E_score, IR_mismatches)

IRmism_lm_fit <- linear_reg() %>%
  fit(E_score ~ ., data = IRmism_df)
tidy(IRmism_lm_fit)
summary(IRmism_lm_fit$fit)

IRmism_tree_fit <- get_tree_fit(IRmism_df, "E_score")
IRmism_tree_fit %>% collect_metrics()
IRmism_tree <- extract_workflow(IRmism_tree_fit)
IRmism_tree

# Decision-Tree plot
plot_regression_tree(IRmism_tree, descr)

## DEBUGGING NECESSARY: ERROR OCCURS IN COMMAND HEREAFTER !

# Variable-Importance plot: plotting importance scores for the model predictors
plot_variable_importance(IRmism_tree, descr)

# RESULTS:
# The linear regression result supports the previous Pearson metric that demonstrated
# lack of statistical significance for mismatch's impact on evenness. The decision tree
# additionally supports this with no importance given to mismatches.

########################################################################
########################################################################

########################################################################
## Effects of READ LENGTH on SEQ EVENNESS
########################################################################
#descr = "READ LENGTH on SEQ EVENNESS"

#readLen_df <- genome_df %>%
#  filter_small_method_classes() %>%
#  handle_outliers("E_score", 3, 0) %>%
#  select(E_score, avgLength)

#readLen_tree_fit <- get_tree_fit(readLen_df, "E_score")
#readLen_tree_fit %>% collect_metrics()
#readLen_tree <- extract_workflow(readLen_tree_fit)
#readLen_tree

## Decision-Tree plot
#plot_regression_tree(readLen_tree, descr)

## Variable-Importance plot: plotting importance scores for the model predictors
#plot_variable_importance(readLen_tree, descr)

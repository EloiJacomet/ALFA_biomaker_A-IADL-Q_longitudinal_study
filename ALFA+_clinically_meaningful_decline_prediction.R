#
# ============================================================================= #
# ---- ALFA+ CLINICALLY MEANINGFUL FUNCTIONAL DECLINE PREDICTION ----
# ============================================================================= #
# Purpose:
#   - Create logistic regression models to:
#     - Test which fluid biomarkers best predict functional decline
#     - Test wheather adding demographic/genetic/cognitive variables improve accuracy
#     - Generate the best possible predictive model
#     - Evaluate the predictive capacity with rigurous metrics
# 
# Description:
#   This script aims to fit and analyze multiple univariate and multivariate
#   prediction models via Logistic Regression models to predict future clinically
#   meaningful functional decline using a validated -2.2 MIC in A-IADL-Q scores and
#   discover which fluid biomarkers exhibit the best prognostic.
#
# The Pipeline includes:
#   - Build univariate and multivariate prognostic models to predict who will decline based on the A-IADL-Q score.
#   - Test univariate prognostic accuracy of demographics and cognitive variables.
#   - Test univariate prognostic accuracy of individual plasma and CSF biomarkers.
#   - Compare clinical models with both SUCOG and PACC as cognitive assessments.
#   - Optimize models based on AUCs, DeLong's tests and likelihood ratio tests.
#   - Measure prediction capabilities of optimized models with ROC-AUC 10-fold CV.
#
# Inputs:
#   - Full ALFA+ dataset with a binary variable of decliners/non-decliners (-2.2 A-IADL-Q score V3-V1)
#
# Outputs:
#   - Model evaluations with AUCs, Bootstrapped 95% CI, NPV, PPV...
#   - ROC curves and forest plots for visualization of the results
#
# Reproducibility:
#   - Developed in R version 4.5.2
#   - Key Packages:
#     - dplyr (1.1.4)
#     - pROC (1.19.0.1)
#     - rsample (1.3.2)
#     - ggplot2 (4.0.1)
#     - tidytext (0.4.3)
#
# Code: Eloi Jacomet & Federica Anastasi
# Date: June 6th 2026 
# ============================================================================= #

rm(list = ls())

# ============================================================================= #
# 1. SETUP ----
# ============================================================================= #

## ---- 1.1 Load libraries ----

suppressPackageStartupMessages({
  library(dplyr)
  library(pROC)
  library(rsample)
  library(ggplot2)
  library(tidytext)
  library(here)
})

## ---- 1.2 Set working and output directories and working seed ----

wd <- here()

setwd(wd)

dir.create("results_functional_prediction")

set.seed(123)

## ---- 1.3 Load data ----

load("data/alfa_study_decliners_20260602.RData")

## ---- 1.4 Recode variables and clean the dataframe ----

# Keep only baseline variables for prediction
baseline <- alfa_study_decliners[alfa_study_decliners$visit == "V1", ]

cat("Baseline N:", nrow(baseline), "\n")

# Create binary yes/no categorical and 1/0 (for logistic regression) decliner variable
baseline$decliner_fct <- ifelse(baseline$iadl_decline == "Decliner", "Yes", "No")
baseline$decliner_fct[is.na(baseline$iadl_decline)] <- NA
baseline$decliner_fct <- factor(baseline$decliner_fct, levels = c("No", "Yes"))

baseline$decliner_bin <- ifelse(baseline$decliner_fct == "Yes", 1, 0)
baseline$decliner_bin[is.na(baseline$decliner_fct)] <- NA

analysis_df <- baseline[!is.na(baseline$decliner_bin), ]

cat("N after removing missing outcome:", nrow(analysis_df), "\n")
print(table(analysis_df$decliner_fct))

## ---- 1.5 Recode categorical predictors ----

analysis_df$APOE_binary <- ifelse(
  analysis_df$APOE_binary %in% c("0", 0), "Non-carrier",
  ifelse(analysis_df$APOE_binary %in% c("1", 1), "Carrier", as.character(analysis_df$APOE_binary))
)
analysis_df$APOE_binary <- factor(analysis_df$APOE_binary)

analysis_df$SCD_def <- factor(analysis_df$SCD_def)
analysis_df$A_bl <- factor(analysis_df$A_bl)
analysis_df$T_bl <- factor(analysis_df$T_bl)
analysis_df$N_bl <- factor(analysis_df$N_bl)

cat("\nA_bl distribution:\n")
print(table(analysis_df$A_bl, useNA = "ifany"))
cat("\nA_bl levels:\n")
print(levels(analysis_df$A_bl))

## ---- 1.6 Create biomarker ratios ----

# PLASMA Ab42/40
analysis_df$plasma_Ab42_Ab40_ratio <- ifelse(
  !is.na(analysis_df$BIO_Abeta1.42_PLASMA_Elecsys_UGOT_Roche) &
    !is.na(analysis_df$BIO_Abeta1.40_PLASMA_Elecsys_UGOT_Roche) &
    analysis_df$BIO_Abeta1.40_PLASMA_Elecsys_UGOT_Roche > 0,
  analysis_df$BIO_Abeta1.42_PLASMA_Elecsys_UGOT_Roche /
    analysis_df$BIO_Abeta1.40_PLASMA_Elecsys_UGOT_Roche,
  NA
)

# CSF Ab42/40
analysis_df$csf_Ab42_Ab40_ratio <- ifelse(
  !is.na(analysis_df$BIO_Abeta1.42_CSF_Elecsys_UGOT_Roche) &
    !is.na(analysis_df$BIO_Abeta1.40_CSF_Elecsys_UGOT_Roche) &
    analysis_df$BIO_Abeta1.40_CSF_Elecsys_UGOT_Roche > 0,
  analysis_df$BIO_Abeta1.42_CSF_Elecsys_UGOT_Roche /
    analysis_df$BIO_Abeta1.40_CSF_Elecsys_UGOT_Roche,
  NA
)

# PLASMA pTau181/Ab42
analysis_df$plasma_pTau181_Ab42_ratio <- ifelse(
  !is.na(analysis_df$BIO_pTau181_PLASMA_Elecsys_UGOT_Roche) &
    !is.na(analysis_df$BIO_Abeta1.42_PLASMA_Elecsys_UGOT_Roche) &
    analysis_df$BIO_Abeta1.42_PLASMA_Elecsys_UGOT_Roche > 0,
  analysis_df$BIO_pTau181_PLASMA_Elecsys_UGOT_Roche /
    analysis_df$BIO_Abeta1.42_PLASMA_Elecsys_UGOT_Roche,
  NA
)

# PLASMA pTau217/Ab42
analysis_df$plasma_pTau217_Ab42_ratio <- ifelse(
  !is.na(analysis_df$BIO_pTau217_PLASMA_Elecsys_UGOT_Roche) &
    !is.na(analysis_df$BIO_Abeta1.42_PLASMA_Elecsys_UGOT_Roche) &
    analysis_df$BIO_Abeta1.42_PLASMA_Elecsys_UGOT_Roche > 0,
  analysis_df$BIO_pTau217_PLASMA_Elecsys_UGOT_Roche /
    analysis_df$BIO_Abeta1.42_PLASMA_Elecsys_UGOT_Roche,
  NA
)

# CSF pTau181/Ab42
analysis_df$csf_pTau181_Ab42_ratio <- ifelse(
  !is.na(analysis_df$BIO_pTau181_CSF_Elecsys_UGOT_Roche) &
    !is.na(analysis_df$BIO_Abeta1.42_CSF_Elecsys_UGOT_Roche) &
    analysis_df$BIO_Abeta1.42_CSF_Elecsys_UGOT_Roche > 0,
  analysis_df$BIO_pTau181_CSF_Elecsys_UGOT_Roche /
    analysis_df$BIO_Abeta1.42_CSF_Elecsys_UGOT_Roche,
  NA
)

# CSF pTau217/Ab42
analysis_df$csf_pTau217_Ab42_ratio <- ifelse(
  !is.na(analysis_df$BIO_pTau217_CSF_MSD_LILLY_in.house) &
    !is.na(analysis_df$BIO_Abeta1.42_CSF_Elecsys_UGOT_Roche) &
    analysis_df$BIO_Abeta1.42_CSF_Elecsys_UGOT_Roche > 0,
  analysis_df$BIO_pTau217_CSF_MSD_LILLY_in.house /
    analysis_df$BIO_Abeta1.42_CSF_Elecsys_UGOT_Roche,
  NA
)

# CSF pTau231/Ab42
analysis_df$csf_pTau231_Ab42_ratio <- ifelse(
  !is.na(analysis_df$BIO_pTau231_CSF_ELISA_UGOT_ADx.NeuroSciences) &
    !is.na(analysis_df$BIO_Abeta1.42_CSF_Elecsys_UGOT_Roche) &
    analysis_df$BIO_Abeta1.42_CSF_Elecsys_UGOT_Roche > 0,
  analysis_df$BIO_pTau231_CSF_ELISA_UGOT_ADx.NeuroSciences /
    analysis_df$BIO_Abeta1.42_CSF_Elecsys_UGOT_Roche,
  NA
)

# CSF tTau/Ab42
analysis_df$csf_tTau_Ab42_ratio <- ifelse(
  !is.na(analysis_df$BIO_tTau_CSF_Elecsys_UGOT_Roche) &
    !is.na(analysis_df$BIO_Abeta1.42_CSF_Elecsys_UGOT_Roche) &
    analysis_df$BIO_Abeta1.42_CSF_Elecsys_UGOT_Roche > 0,
  analysis_df$BIO_tTau_CSF_Elecsys_UGOT_Roche /
    analysis_df$BIO_Abeta1.42_CSF_Elecsys_UGOT_Roche,
  NA
)

## ---- 1.7 Define variables for univariate analysis ----

uni_vars <- c(
  # demographics
  "Age_V1",
  "sex",
  "YearsEducation",
  "APOE_binary",
  
  # cognition / clinical
  "PACC",
  "FE_Total_SUCOG",
  "LE_Total_SUCOG",
  "ME_Total_SUCOG",
  "Total_SUCOG",
  "Frailty_Index",
  "SCD_def",
  "sucog_decr2yrs",
  
  # ATN
  "A_bl",
  "T_bl",
  "N_bl",
  
  # raw biomarkers
  "BIO_Abeta1.42_CSF_Elecsys_UGOT_Roche",
  "BIO_Abeta1.40_CSF_Elecsys_UGOT_Roche",
  "BIO_Abeta1.42_PLASMA_Elecsys_UGOT_Roche",
  "BIO_Abeta1.40_PLASMA_Elecsys_UGOT_Roche",
  "BIO_NFL_CSF_Elecsys_UGOT_Roche",
  "BIO_NFL_PLASMA_Elecsys_UGOT_Roche",
  "BIO_GFAP_CSF_Elecsys_UGOT_Roche",
  "BIO_GFAP_PLASMA_Elecsys_UGOT_Roche",
  "BIO_pTau181_CSF_Elecsys_UGOT_Roche",
  "BIO_pTau181_PLASMA_Elecsys_UGOT_Roche",
  "BIO_pTau217_CSF_MSD_LILLY_in.house",
  "BIO_pTau231_CSF_ELISA_UGOT_ADx.NeuroSciences",
  "BIO_tTau_CSF_Elecsys_UGOT_Roche",
  "BIO_pTau217_PLASMA_Elecsys_UGOT_Roche",
  
  # ratios
  "plasma_Ab42_Ab40_ratio",
  "csf_Ab42_Ab40_ratio",
  "plasma_pTau181_Ab42_ratio",
  "csf_pTau181_Ab42_ratio",
  "csf_pTau217_Ab42_ratio",
  "csf_pTau231_Ab42_ratio",
  "csf_tTau_Ab42_ratio",
  "plasma_pTau217_Ab42_ratio"
)

uni_vars <- uni_vars[uni_vars %in% names(analysis_df)]

cat("\nVariables to test univariately:\n")
print(uni_vars)

## ---- 1.8 Build one common case dataframe ----

# We remove all individuals with missing values of univariate variables
vars_needed <- c("decliner_bin", "decliner_fct", uni_vars)

complete_df <- analysis_df[complete.cases(analysis_df[, vars_needed]), ]

cat("\nComplete-case N:", nrow(complete_df), "\n")
print(table(complete_df$decliner_fct))
print(table(complete_df$decliner_bin, complete_df$A_bl))

## ---- 1.9 Create 3 stratified dataframes (global, A- and A+) ----

data_all <- complete_df

# Check coding of A_bl before filtering
cat("\nA_bl levels in complete-case dataset:\n")
print(levels(complete_df$A_bl))
print(table(complete_df$A_bl, useNA = "ifany"))

# Amyloid negative
data_Aneg <- complete_df[complete_df$A_bl == "A-", ]

# Amyloid positive
data_Apos <- complete_df[complete_df$A_bl == "A+", ]

cat("\nDataset sizes:\n")
cat("All:", nrow(data_all), "\n")  #215
cat("A- :", nrow(data_Aneg), "\n") #121
cat("A+ :", nrow(data_Apos), "\n") #94

cat("\nOutcome distribution by dataset:\n")
print(table(data_all$decliner_fct))
print(table(data_Aneg$decliner_fct))
print(table(data_Apos$decliner_fct))

cat("\nOutcome by amyloid group:\n")
print(table(complete_df$decliner_bin, complete_df$A_bl))

# ============================================================================= #
# 2. GENERAL UNIVARIATE PREDICTION MODELS ----
# ============================================================================= #

## ---- 2.1 Function to run univariate models ----

run_one_univariate <- function(df, var_name, dataset_name) {
  
  # Keep only the outcome variable and the outcome variable (1/0 for glm)
  sub_df <- df[, c("decliner_bin", var_name)]
  
  # Remove rows where the predictor is missing.
  sub_df <- sub_df[!is.na(sub_df[[var_name]]), ]
  
  # If no rows remain, stop and return nothing.
  if (nrow(sub_df) == 0) {
    return(NULL)
  }
  
  # Check that the predictor has at least 2 unique values (glm needs at least 2 values)
  n_unique <- dplyr::n_distinct(sub_df[[var_name]], na.rm = TRUE)
  if (n_unique < 2) {
    return(NULL)
  }
  
  # Create the model formula, for example
  formula_text <- paste("decliner_bin ~", var_name)
  uni_formula <- as.formula(formula_text)
  
  # Run a univariate logistic regression with outcome and variable of interest
  model <- try(glm(uni_formula, data = sub_df, family = binomial), silent = TRUE)
  
  # If the model fails, return nothing.
  if (inherits(model, "try-error")) {
    return(NULL)
  }
  
  # Extract coefficient table from the model summary.
  coef_mat <- summary(model)$coefficients
  
  # If there is no predictor coefficient, return nothing.
  if (nrow(coef_mat) < 2) {
    return(NULL)
  }
  
  # Predict probabilities from the logistic regression model.
  probs <- try(predict(model, type = "response"), silent = TRUE)
  if (inherits(probs, "try-error")) {
    return(NULL)
  }
  
  # Calculate ROC curve and AUC.
  roc_obj <- try(pROC::roc(sub_df$decliner_bin, probs, quiet = TRUE), silent = TRUE)
  if (inherits(roc_obj, "try-error")) {
    return(NULL)
  }
  
  auc_val <- as.numeric(pROC::auc(roc_obj))
  
  # Likelihood ratio test for the overall contribution of the predictor.
  model_lrt <- try(drop1(model, test = "Chisq"), silent = TRUE)
  
  overall_p <- NA
  if (!inherits(model_lrt, "try-error")) {
    if (nrow(model_lrt) >= 2) {
      overall_p <- model_lrt$`Pr(>Chi)`[2]
    }
  }
  
  # Intialize results dataframe
  results_var <- data.frame()
  
  # Loop over all predictor coefficients.
  for (i in 2:nrow(coef_mat)) {
    
    est <- coef_mat[i, "Estimate"]
    se  <- coef_mat[i, "Std. Error"]
    p   <- coef_mat[i, "Pr(>|z|)"]
    
    # Convert log-odds coefficient into odds ratio.
    OR <- exp(est)
    
    # Calculate 95% confidence interval for the odds ratio.
    CI_low <- exp(est - 1.96 * se)
    CI_high <- exp(est + 1.96 * se)
    
    # Store results for this predictor term.
    row_res <- data.frame(
      dataset = dataset_name,
      variable = var_name,
      term = rownames(coef_mat)[i],
      OR = OR,
      CI_low = CI_low,
      CI_high = CI_high,
      p_value = p,
      overall_model_p = overall_p,
      AUC = auc_val,
      N = nrow(sub_df)
    )
    
    results_var <- rbind(results_var, row_res)
  }
  
  return(results_var)
}

## ---- 2.2 Function to run prediction ----

run_univariate_dataset <- function(df, var_list, dataset_name) {
  
  # Create an empty object to store the results from all univariate models.
  all_results <- data.frame()
  
  # Print dataset information to monitor progress and sample composition.
  cat("\n========================================\n")
  cat("Running univariate models in:", dataset_name, "\n")
  
  # Print total sample size in this dataset.
  cat("N =", nrow(df), "\n")
  
  # Print the distribution of the outcome groups
  print(table(df$decliner_fct))
  
  cat("========================================\n")
  
  # Loop through all variables in the variable list.
  for (v in var_list) {
    
    # Print current variable name to track progress.
    cat("Running:", v, "\n")
    
    # Run the univariate logistic regression function for the current predictor.
    res <- run_one_univariate(df, v, dataset_name)
    
    # Append them to the final results table if it yielded results
    if (!is.null(res)) {
      all_results <- rbind(all_results, res)
    }
  }
  
  # Return a combined table containing the results
  return(all_results)
}

## ---- 2.3 Run univariate models for all previously selected variables ----

uni_all <- run_univariate_dataset(data_all, uni_vars, "All")
uni_Aneg <- run_univariate_dataset(data_Aneg, uni_vars, "A_negative")
uni_Apos <- run_univariate_dataset(data_Apos, uni_vars, "A_positive")

## ---- 2.4 Combine all results ----

uni_results <- rbind(uni_all, uni_Aneg, uni_Apos)

## ---- 2.5 Examine results and save complete csv data ----

cat("\nTop results by AUC:\n")
uni_results_auc <- uni_results[order(-uni_results$AUC), ]
print(head(uni_results_auc, 30))

cat("\nTop results by p-value:\n")
uni_results_p <- uni_results[order(uni_results$p_value), ]
print(head(uni_results_p, 30))

cat("\nTop AUC results - All:\n")
print(head(uni_results_auc[uni_results_auc$dataset == "All", ], 20))

cat("\nTop AUC results - A_negative:\n")
print(head(uni_results_auc[uni_results_auc$dataset == "A_negative", ], 20))

cat("\nTop AUC results - A_positive:\n")
print(head(uni_results_auc[uni_results_auc$dataset == "A_positive", ], 20))

univariate_results <- list(
  variables_tested = uni_vars,
  data_all_n = nrow(data_all),
  data_Aneg_n = nrow(data_Aneg),
  data_Apos_n = nrow(data_Apos),
  results_all = uni_all,
  results_Aneg = uni_Aneg,
  results_Apos = uni_Apos,
  results_combined = uni_results
)

write.csv(
  univariate_results$results_combined,
  file = "results_functional_prediction/univariate_results_all_datasets.csv",
  row.names = FALSE
)

# ============================================================================= #
# 3. BIOMARKER ONLY UNIVARIATE PREDICTION MODELS ----
# ============================================================================= #

## ---- 3.1 Extract only biomarkers and ration ----

biomarker_vars <- uni_vars[
  grepl("^BIO_", uni_vars) | grepl("_ratio$", uni_vars)
]

cat("\nBiomarker / ratio variables for ROC analysis:\n")
print(biomarker_vars)

## ---- 3.2 Function to run AUC and bootstrap 95% CI ----

run_one_biomarker_roc <- function(df, var_name, dataset_name, boot_n = 2000) {
  
  sub_df <- df[, c("decliner_bin", var_name)]
  
  # Remove missing outcome or missing predictor
  sub_df <- sub_df[
    !is.na(sub_df$decliner_bin) &
      !is.na(sub_df[[var_name]]),
  ]
  
  # Need at least some data
  if (nrow(sub_df) == 0) {
    return(NULL)
  }
  
  # Need both outcome classes
  if (length(unique(sub_df$decliner_bin)) < 2) {
    return(NULL)
  }
  
  # Predictor must vary
  if (dplyr::n_distinct(sub_df[[var_name]]) < 2) {
    return(NULL)
  }
  
  # ROC analysis
  roc_obj <- try(
    pROC::roc(
      response = sub_df$decliner_bin,
      predictor = sub_df[[var_name]],
      quiet = TRUE,
      direction = "auto"
    ),
    silent = TRUE
  )
  
  if (inherits(roc_obj, "try-error")) {
    return(NULL)
  }
  
  auc_val <- as.numeric(pROC::auc(roc_obj))
  
  # Bootstrap CI for AUC
  auc_ci <- try(
    pROC::ci.auc(
      roc_obj,
      method = "bootstrap",
      boot.n = boot_n,
      conf.level = 0.95
    ),
    silent = TRUE
  )
  
  if (inherits(auc_ci, "try-error")) {
    ci_low <- NA
    ci_high <- NA
  } else {
    ci_low <- as.numeric(auc_ci[1])
    ci_high <- as.numeric(auc_ci[3])
  }
  
  data.frame(
    dataset = dataset_name,
    variable = var_name,
    AUC = auc_val,
    CI_low = ci_low,
    CI_high = ci_high,
    N = nrow(sub_df),
    n_non_decliners = sum(sub_df$decliner_bin == 0),
    n_decliners = sum(sub_df$decliner_bin == 1)
  )
}

## ---- 3.3 Function to run biomarker-only prediction ----

run_biomarker_roc_dataset <- function(df, var_list, dataset_name, boot_n = 2000) {
  
  all_results <- data.frame()
  
  cat("\n========================================\n")
  cat("Running biomarker ROC analyses in:", dataset_name, "\n")
  cat("N =", nrow(df), "\n")
  print(table(df$decliner_fct))
  cat("========================================\n")
  
  for (v in var_list) {
    cat("Running ROC:", v, "\n")
    
    res <- run_one_biomarker_roc(
      df = df,
      var_name = v,
      dataset_name = dataset_name,
      boot_n = boot_n
    )
    
    if (!is.null(res)) {
      all_results <- rbind(all_results, res)
    }
  }
  
  return(all_results)
}

## ---- 3.4 Run univariate models for all 3 groups ----

roc_all <- run_biomarker_roc_dataset(data_all, biomarker_vars, "All")
roc_Aneg <- run_biomarker_roc_dataset(data_Aneg, biomarker_vars, "A_negative")
roc_Apos <- run_biomarker_roc_dataset(data_Apos, biomarker_vars, "A_positive")

## ---- 3.5 Combine results and save complete dataframe as csv ----

roc_biomarker_results <- rbind(roc_all, roc_Aneg, roc_Apos)

# Save table
write.csv(
  roc_biomarker_results,
  file = "results_functional_prediction/biomarker_ratio_roc_bootstrap_CI.csv",
  row.names = FALSE
)

## ---- 3.6 Forest plot of the biomarker-only univariate models ----

# Extract results into a new dataframe
roc_plot_df <- roc_biomarker_results

# Rename variables and stratas for prettier plots
variable_labels <- c(
  "BIO_Abeta1.42_PLASMA_Elecsys_UGOT_Roche" = "Plasma Aβ42",
  "BIO_pTau181_PLASMA_Elecsys_UGOT_Roche"   = "Plasma pTau181",
  "BIO_pTau181_CSF_Elecsys_UGOT_Roche"      = "CSF pTau181",
  "BIO_pTau217_PLASMA_Elecsys_UGOT_Roche"   = "Plasma pTau217",
  "BIO_pTau217_CSF_MSD_LILLY_in.house"      = "CSF pTau217",
  "BIO_pTau231_CSF_ELISA_UGOT_ADx.NeuroSciences" = "CSF pTau231",
  "BIO_tTau_CSF_Elecsys_UGOT_Roche"         = "CSF tTau",
  "BIO_GFAP_PLASMA_Elecsys_UGOT_Roche"      = "Plasma GFAP",
  "BIO_NFL_PLASMA_Elecsys_UGOT_Roche"       = "Plasma NfL",
  
  "plasma_pTau181_Ab42_ratio" = "Plasma pTau181/Aβ42",
  "plasma_pTau217_Ab42_ratio" = "Plasma pTau217/Aβ42",
  "csf_pTau231_Ab42_ratio"    = "CSF pTau231/Aβ42",
  "csf_tTau_Ab42_ratio"       = "CSF tTau/Aβ42",
  "csf_pTau181_Ab42_ratio"    = "CSF pTau181/Aβ42",
  "csf_pTau217_Ab42_ratio"    = "CSF pTau217/Aβ42",
  
  "BIO_Abeta1.42_CSF_Elecsys_UGOT_Roche"    = "CSF Aβ42",
  "BIO_Abeta1.40_PLASMA_Elecsys_UGOT_Roche" = "Plasma Aβ40",
  "csf_Ab42_Ab40_ratio"                     = "CSF Aβ42/Aβ40",
  "plasma_Ab42_Ab40_ratio"                  = "Plasma Aβ42/Aβ40",
  "BIO_Abeta1.40_CSF_Elecsys_UGOT_Roche"    = "CSF Aβ40",
  
  "BIO_GFAP_CSF_Elecsys_UGOT_Roche"         = "CSF GFAP",
  "BIO_NFL_CSF_Elecsys_UGOT_Roche"          = "CSF NfL"
)

roc_plot_df$variable <- dplyr::recode(
  roc_plot_df$variable,
  !!!variable_labels
)

roc_plot_df$dataset <- factor(
  roc_plot_df$dataset,
  levels = c("All", "A_negative", "A_positive"),
  labels = c("All", "A-", "A+")
)

# Create plot
p_auc <- ggplot(
  roc_plot_df,
  aes(
    x = AUC,
    y = tidytext::reorder_within(variable, AUC, dataset)
  )
) +
  geom_vline(xintercept = 0.80, linetype = "dashed") +
  geom_errorbarh(
    aes(xmin = CI_low, xmax = CI_high),
    height = 0.25
  ) +
  geom_point(size = 2.5) +
  tidytext::scale_y_reordered() +
  facet_wrap(~ dataset, scales = "free_y") +
  labs(
    title = "Biomarker and ratio ROC performance",
    subtitle = "Univariate AUC with 95% bootstrap confidence intervals",
    x = "AUC",
    y = NULL
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 8)
  )

print(p_auc)

# Save plot
ggsave(
  filename = "results_functional_prediction/biomarker_ratio_roc_bootstrap_CI.png",
  plot = p_auc,
  width = 10,
  height = 5,
  dpi = 300
)

# ============================================================================= #
# 4. BASE CLINICAL MODELS + COVARIATES: PACC VS. SUCOG ----
# ============================================================================= #

##  ---- 4.1 Define base clinical model formulas ----

formula_pacc <- decliner_bin ~ Age_V1 + sex + APOE_binary + PACC

formula_sucog <- decliner_bin ~ Age_V1 + sex + APOE_binary + Total_SUCOG


## ---- 4.2 Function to run multivariate models and extract AUC + 95% CI  ----

run_base_model_roc <- function(df,
                               formula_obj,
                               model_name,
                               dataset_name,
                               boot_n = 2000) {
  
  # Identify variables used in the model
  vars_needed <- all.vars(formula_obj)
  
  # Keep only needed variables
  sub_df <- df[, vars_needed]
  
  # Remove rows with missing values in any model variable
  sub_df <- sub_df[complete.cases(sub_df), ]
  
  # Stop if no data remain
  if (nrow(sub_df) == 0) {
    return(NULL)
  }
  
  # Need both decliners and non-decliners
  if (length(unique(sub_df$decliner_bin)) < 2) {
    return(NULL)
  }
  
  # Fit logistic regression model
  model <- try(
    glm(
      formula_obj,
      data = sub_df,
      family = binomial
    ),
    silent = TRUE
  )
  
  # Stop if model fails
  if (inherits(model, "try-error")) {
    return(NULL)
  }
  
  # Predicted probability of being a decliner
  probs <- try(
    predict(model, type = "response"),
    silent = TRUE
  )
  
  if (inherits(probs, "try-error")) {
    return(NULL)
  }
  
  # ROC curve based on predicted probabilities
  roc_obj <- try(
    pROC::roc(
      response = sub_df$decliner_bin,
      predictor = probs,
      quiet = TRUE
    ),
    silent = TRUE
  )
  
  if (inherits(roc_obj, "try-error")) {
    return(NULL)
  }
  
  # AUC
  auc_val <- as.numeric(pROC::auc(roc_obj))
  
  # Bootstrap 95% CI for AUC
  auc_ci <- try(
    pROC::ci.auc(
      roc_obj,
      method = "bootstrap",
      boot.n = boot_n,
      conf.level = 0.95
    ),
    silent = TRUE
  )
  
  if (inherits(auc_ci, "try-error")) {
    ci_low <- NA
    ci_high <- NA
  } else {
    ci_low <- as.numeric(auc_ci[1])
    ci_high <- as.numeric(auc_ci[3])
  }
  
  # Return model performance summary
  data.frame(
    dataset = dataset_name,
    model = model_name,
    AUC = auc_val,
    CI_low = ci_low,
    CI_high = ci_high,
    N = nrow(sub_df),
    n_non_decliners = sum(sub_df$decliner_bin == 0),
    n_decliners = sum(sub_df$decliner_bin == 1)
  )
}

## ---- 4.3 Run prediction models for PACC and SUCOG in all strata  ----

base_model_results <- rbind(
  
  # PACC-based model
  run_base_model_roc(
    df = data_all,
    formula_obj = formula_pacc,
    model_name = "Age + Sex + APOE + PACC",
    dataset_name = "All"
  ),
  
  run_base_model_roc(
    df = data_Aneg,
    formula_obj = formula_pacc,
    model_name = "Age + Sex + APOE + PACC",
    dataset_name = "A_negative"
  ),
  
  run_base_model_roc(
    df = data_Apos,
    formula_obj = formula_pacc,
    model_name = "Age + Sex + APOE + PACC",
    dataset_name = "A_positive"
  ),
  
  
  # SUCOG-based model
  run_base_model_roc(
    df = data_all,
    formula_obj = formula_sucog,
    model_name = "Age + Sex + APOE + SUCOG",
    dataset_name = "All"
  ),
  
  run_base_model_roc(
    df = data_Aneg,
    formula_obj = formula_sucog,
    model_name = "Age + Sex + APOE + SUCOG",
    dataset_name = "A_negative"
  ),
  
  run_base_model_roc(
    df = data_Apos,
    formula_obj = formula_sucog,
    model_name = "Age + Sex + APOE + SUCOG",
    dataset_name = "A_positive"
  )
)

## ---- 4.4 Print and save results  ----

cat("\nBase clinical model results:\n")
print(base_model_results)

write.csv(
  base_model_results,
  file = "results_functional_prediction/base_clinical_model_PACC_vs_SUCOG_auc_bootstrap_CI.csv",
  row.names = FALSE
)

## ---- 4.5 Plot cognitive multivariate models  ----

base_plot_df <- base_model_results

base_plot_df$dataset <- factor(
  base_plot_df$dataset,
  levels = c("All", "A_negative", "A_positive"),
  labels = c("All", "A-", "A+")
)

base_plot_df$model <- factor(
  base_plot_df$model,
  levels = c(
    "Age + Sex + APOE + PACC",
    "Age + Sex + APOE + SUCOG"
  )
)

p_base <- ggplot(
  base_plot_df,
  aes(
    x = AUC,
    y = model
  )
) +
  geom_vline(xintercept = 0.80, linetype = "dashed") +
  geom_errorbarh(
    aes(xmin = CI_low, xmax = CI_high),
    height = 0.20
  ) +
  geom_point(size = 3) +
  facet_wrap(~ dataset) +
  labs(
    title = "Base clinical model performance",
    subtitle = "AUC with 95% bootstrap confidence intervals",
    x = "AUC",
    y = NULL
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 9)
  )

print(p_base)

ggsave(
  filename = "results_functional_prediction/base_clinical_model_PACC_vs_SUCOG_auc_plot.png",
  plot = p_base,
  width = 10,
  height = 5,
  dpi = 300
)

## ---- 4.6 Statistical comparison between groups with DeLong's Test ----

# Keep only variables needed for both models
vars_compare <- c(
  "decliner_bin",
  "Age_V1",
  "sex",
  "APOE_binary",
  "PACC",
  "Total_SUCOG"
)

compare_df <- data_Apos[, vars_compare]

# Remove missing values
compare_df <- compare_df[complete.cases(compare_df), ]

cat("\nA+ comparison sample size:\n")
print(nrow(compare_df))

cat("\nOutcome distribution:\n")
print(table(compare_df$decliner_fct))

# Fit PACC model
model_pacc <- glm(
  decliner_bin ~ Age_V1 + sex + APOE_binary + PACC,
  data = compare_df,
  family = binomial
)

# Predicted probabilities
probs_pacc <- predict(model_pacc, type = "response")

# ROC curve
roc_pacc <- pROC::roc(
  response = compare_df$decliner_bin,
  predictor = probs_pacc,
  quiet = TRUE
)

# Fit SUCOG model
model_sucog <- glm(
  decliner_bin ~ Age_V1 + sex + APOE_binary + Total_SUCOG,
  data = compare_df,
  family = binomial
)

# Predicted probabilities
probs_sucog <- predict(model_sucog, type = "response")

# ROC curve
roc_sucog <- pROC::roc(
  response = compare_df$decliner_bin,
  predictor = probs_sucog,
  quiet = TRUE
)

# Compare ROC curves using DeLong test
roc_comparison <- pROC::roc.test(
  roc_pacc,
  roc_sucog,
  method = "delong",
  paired = TRUE
)

cat("\n========================================\n")
cat("DeLong test: PACC vs SUCOG in A+\n")
cat("========================================\n")

print(roc_comparison)

# Print AUCs clearly
cat("\nPACC model AUC:\n")
print(as.numeric(pROC::auc(roc_pacc)))

cat("\nSUCOG model AUC:\n")
print(as.numeric(pROC::auc(roc_sucog)))

# ============================================================================= #
# 5. ADDED VALUE OF INDIVIDUAL BIOMARKERS BEYOND CLINICAL MODELS ----
# ============================================================================= #

## ---- 5.1 Function to run base vs. base + biomarker and compare ----

run_added_value_one_biomarker <- function(df,
                                          biomarker,
                                          dataset_name,
                                          base_formula,
                                          base_model_name,
                                          boot_n = 2000) {
  
  vars_needed <- unique(c(
    all.vars(base_formula),
    biomarker
  ))
  
  sub_df <- df[, vars_needed]
  sub_df <- sub_df[complete.cases(sub_df), ]
  
  if (nrow(sub_df) == 0) {
    return(NULL)
  }
  
  if (length(unique(sub_df$decliner_bin)) < 2) {
    return(NULL)
  }
  
  if (dplyr::n_distinct(sub_df[[biomarker]]) < 2) {
    return(NULL)
  }
  
  extended_formula <- as.formula(
    paste(deparse(base_formula), "+", biomarker)
  )
  
  # Fit base model
  model_base <- try(
    glm(
      base_formula,
      data = sub_df,
      family = binomial
    ),
    silent = TRUE
  )
  
  if (inherits(model_base, "try-error")) {
    return(NULL)
  }
  
  # Fit extended model
  model_extended <- try(
    glm(
      extended_formula,
      data = sub_df,
      family = binomial
    ),
    silent = TRUE
  )
  
  if (inherits(model_extended, "try-error")) {
    return(NULL)
  }
  
  # Predicted probabilities
  probs_base <- predict(model_base, type = "response")
  probs_extended <- predict(model_extended, type = "response")
  
  # ROC curves
  roc_base <- try(
    pROC::roc(
      response = sub_df$decliner_bin,
      predictor = probs_base,
      quiet = TRUE
    ),
    silent = TRUE
  )
  
  roc_extended <- try(
    pROC::roc(
      response = sub_df$decliner_bin,
      predictor = probs_extended,
      quiet = TRUE
    ),
    silent = TRUE
  )
  
  if (inherits(roc_base, "try-error") | inherits(roc_extended, "try-error")) {
    return(NULL)
  }
  
  auc_base <- as.numeric(pROC::auc(roc_base))
  auc_extended <- as.numeric(pROC::auc(roc_extended))
  delta_auc <- auc_extended - auc_base
  
  # Bootstrap CI for extended model AUC
  auc_ci <- try(
    pROC::ci.auc(
      roc_extended,
      method = "bootstrap",
      boot.n = boot_n,
      conf.level = 0.95
    ),
    silent = TRUE
  )
  
  if (inherits(auc_ci, "try-error")) {
    ci_low <- NA
    ci_high <- NA
  } else {
    ci_low <- as.numeric(auc_ci[1])
    ci_high <- as.numeric(auc_ci[3])
  }
  
  # DeLong test: base AUC vs extended AUC
  delong_test <- try(
    pROC::roc.test(
      roc_base,
      roc_extended,
      method = "delong",
      paired = TRUE
    ),
    silent = TRUE
  )
  
  if (inherits(delong_test, "try-error")) {
    delong_p <- NA
  } else {
    delong_p <- delong_test$p.value
  }
  
  # Likelihood ratio test: nested model comparison
  lrt_test <- try(
    anova(
      model_base,
      model_extended,
      test = "Chisq"
    ),
    silent = TRUE
  )
  
  if (inherits(lrt_test, "try-error")) {
    lrt_p <- NA
  } else {
    lrt_p <- lrt_test$`Pr(>Chi)`[2]
  }
  
  data.frame(
    dataset = dataset_name,
    base_model = base_model_name,
    biomarker = biomarker,
    model = paste0(base_model_name, " + ", biomarker),
    AUC_base = auc_base,
    AUC_extended = auc_extended,
    CI_low = ci_low,
    CI_high = ci_high,
    delta_AUC = delta_auc,
    delong_p = delong_p,
    lrt_p = lrt_p,
    N = nrow(sub_df),
    n_non_decliners = sum(sub_df$decliner_bin == 0),
    n_decliners = sum(sub_df$decliner_bin == 1)
  )
}

## ---- 5.2 Function to run the added value for each biomarker ----

run_added_value_dataset <- function(df,
                                    biomarker_list,
                                    dataset_name,
                                    base_formula,
                                    base_model_name,
                                    boot_n = 2000) {
  
  all_results <- data.frame()
  
  cat("\n========================================\n")
  cat("Running added-value models in:", dataset_name, "\n")
  cat("Base model:", base_model_name, "\n")
  cat("N =", nrow(df), "\n")
  print(table(df$decliner_fct))
  cat("========================================\n")
  
  for (b in biomarker_list) {
    
    cat("Running:", base_model_name, "+", b, "\n")
    
    res <- run_added_value_one_biomarker(
      df = df,
      biomarker = b,
      dataset_name = dataset_name,
      base_formula = base_formula,
      base_model_name = base_model_name,
      boot_n = boot_n
    )
    
    if (!is.null(res)) {
      all_results <- rbind(all_results, res)
    }
  }
  
  return(all_results)
}

## ---- 5.3 Run added biomarker value for both PACC and SUCOG in all strata ----

formula_pacc <- decliner_bin ~ Age_V1 + sex + APOE_binary + PACC
formula_sucog <- decliner_bin ~ Age_V1 + sex + APOE_binary + Total_SUCOG

# PACC base model + biomarkers
added_pacc_all <- run_added_value_dataset(
  data_all,
  biomarker_vars,
  "All",
  formula_pacc,
  "Age + Sex + APOE + PACC"
)

added_pacc_Aneg <- run_added_value_dataset(
  data_Aneg,
  biomarker_vars,
  "A_negative",
  formula_pacc,
  "Age + Sex + APOE + PACC"
)

added_pacc_Apos <- run_added_value_dataset(
  data_Apos,
  biomarker_vars,
  "A_positive",
  formula_pacc,
  "Age + Sex + APOE + PACC"
)

added_pacc_results <- rbind(
  added_pacc_all,
  added_pacc_Aneg,
  added_pacc_Apos
)


# SUCOG base model + biomarkers
added_sucog_all <- run_added_value_dataset(
  data_all,
  biomarker_vars,
  "All",
  formula_sucog,
  "Age + Sex + APOE + SUCOG"
)

added_sucog_Aneg <- run_added_value_dataset(
  data_Aneg,
  biomarker_vars,
  "A_negative",
  formula_sucog,
  "Age + Sex + APOE + SUCOG"
)

added_sucog_Apos <- run_added_value_dataset(
  data_Apos,
  biomarker_vars,
  "A_positive",
  formula_sucog,
  "Age + Sex + APOE + SUCOG"
)

added_sucog_results <- rbind(
  added_sucog_all,
  added_sucog_Aneg,
  added_sucog_Apos
)

## ---- 5.4 Save all dataframes as csv ----

write.csv(
  added_pacc_results,
  file = "results_functional_prediction/added_value_PACC_base_plus_biomarkers.csv",
  row.names = FALSE
)

write.csv(
  added_sucog_results,
  file = "results_functional_prediction/added_value_SUCOG_base_plus_biomarkers.csv",
  row.names = FALSE
)

## ---- 5.5 Plot added biomarker for both PACC and SUCOG ----

plot_added_value_results <- function(results_df,
                                     plot_title,
                                     output_file) {
  
  plot_df <- results_df
  
  plot_df$dataset <- factor(
    plot_df$dataset,
    levels = c("All", "A_negative", "A_positive"),
    labels = c("All", "A-", "A+")
  )
  
  # Significance annotation:
  # *  = biomarker significantly improves model by LRT
  # ** = biomarker significantly improves both LRT and DeLong
  #
  # LRT tests whether the biomarker adds model information.
  # DeLong tests whether AUC significantly improves.
  
  plot_df$sig_label <- ifelse(
    !is.na(plot_df$lrt_p) & plot_df$lrt_p < 0.05 &
      !is.na(plot_df$delong_p) & plot_df$delong_p < 0.05,
    "**",
    ifelse(
      !is.na(plot_df$lrt_p) & plot_df$lrt_p < 0.05,
      "*",
      ""
    )
  )
  
  p <- ggplot(
    plot_df,
    aes(
      x = AUC_extended,
      y = tidytext::reorder_within(biomarker, AUC_extended, dataset)
    )
  ) +
    geom_vline(xintercept = 0.80, linetype = "dashed") +
    geom_errorbarh(
      aes(xmin = CI_low, xmax = CI_high),
      height = 0.25
    ) +
    geom_point(size = 2.5) +
    geom_text(
      aes(label = sig_label),
      nudge_x = 0.015,
      size = 5
    ) +
    tidytext::scale_y_reordered() +
    facet_wrap(~ dataset, scales = "free_y") +
    labs(
      title = plot_title,
      subtitle = "* LRT p < 0.05; ** LRT p < 0.05 and DeLong p < 0.05",
      x = "AUC of clinical model + biomarker",
      y = NULL
    ) +
    theme_bw() +
    theme(
      strip.text = element_text(face = "bold"),
      axis.text.y = element_text(size = 8)
    )
  
  print(p)
  
  ggsave(
    filename = output_file,
    plot = p,
    width = 12,
    height = 8,
    dpi = 300
  )
  
  return(p)
}

# Plot: PACC base + biomarkers
p_added_pacc <- plot_added_value_results(
  results_df = added_pacc_results,
  plot_title = "Added value of biomarkers beyond PACC clinical model",
  output_file = "results_functional_prediction/added_value_PACC_plus_biomarkers_auc_plot.png"
)


# Plot: SUCOG base + biomarkers
 p_added_sucog <- plot_added_value_results(
  results_df = added_sucog_results,
  plot_title = "Added value of biomarkers beyond SUCOG clinical model",
  output_file = "results_functional_prediction/added_value_SUCOG_plus_biomarkers_auc_plot.png"
)

# ============================================================================= #
# 6. A+ FOCUSED ONE AND TWO-BIOMARKERS MODEL OPTIMIZATION ----
# ============================================================================= #

## ---- 6.1 Biomarker and participant selection ----
 
# Focus only on amyloid-positive participants
ad_df <- data_Apos
 
# Candidate biomarkers selected from previous results
candidate_biomarkers <- c(
  "BIO_NFL_PLASMA_Elecsys_UGOT_Roche",
  "BIO_pTau217_CSF_MSD_LILLY_in.house",
  "BIO_NFL_CSF_Elecsys_UGOT_Roche",
  "BIO_pTau231_CSF_ELISA_UGOT_ADx.NeuroSciences",
  "BIO_pTau181_PLASMA_Elecsys_UGOT_Roche",
  "plasma_pTau181_Ab42_ratio",
  "csf_pTau217_Ab42_ratio",
  "BIO_pTau217_PLASMA_Elecsys_UGOT_Roche",
  "plasma_pTau217_Ab42_ratio"
)
 
candidate_biomarkers <- candidate_biomarkers[
  candidate_biomarkers %in% names(ad_df)
]
 
## ---- 6.2 Create base formulas ----
 
base_formulas <- list(
  PACC_base = decliner_bin ~ Age_V1 + sex + APOE_binary + PACC,
  SUCOG_base = decliner_bin ~ Age_V1 + sex + APOE_binary + Total_SUCOG
)
 
## ---- 6.3 Function to calculate evaluation metrics ----
 
get_classification_metrics <- function(outcome, probs, threshold) {
   
   pred_class <- ifelse(probs >= threshold, 1, 0)
   
   TP <- sum(pred_class == 1 & outcome == 1)
   TN <- sum(pred_class == 0 & outcome == 0)
   FP <- sum(pred_class == 1 & outcome == 0)
   FN <- sum(pred_class == 0 & outcome == 1)
   
   sensitivity <- TP / (TP + FN)
   specificity <- TN / (TN + FP)
   PPV <- TP / (TP + FP)
   NPV <- TN / (TN + FN)
   accuracy <- (TP + TN) / length(outcome)
   
   data.frame(
     threshold = threshold,
     sensitivity = sensitivity,
     specificity = specificity,
     PPV = PPV,
     NPV = NPV,
     accuracy = accuracy,
     TP = TP,
     TN = TN,
     FP = FP,
     FN = FN
   )
}
 
## ---- 6.4 Function to run prediction for candidate models ----
 
evaluate_candidate_model <- function(df, formula_obj, model_name, boot_n = 2000) {
   
   vars_needed <- all.vars(formula_obj)
   sub_df <- df[, vars_needed]
   sub_df <- sub_df[complete.cases(sub_df), ]
   
   if (nrow(sub_df) == 0) return(NULL)
   if (length(unique(sub_df$decliner_bin)) < 2) return(NULL)
   
   model <- try(
     glm(formula_obj, data = sub_df, family = binomial),
     silent = TRUE
   )
   
   if (inherits(model, "try-error")) return(NULL)
   
   probs <- predict(model, type = "response")
   
   roc_obj <- try(
     pROC::roc(
       response = sub_df$decliner_bin,
       predictor = probs,
       quiet = TRUE
     ),
     silent = TRUE
   )
   
   if (inherits(roc_obj, "try-error")) return(NULL)
   
   auc_val <- as.numeric(pROC::auc(roc_obj))
   
   auc_ci <- try(
     pROC::ci.auc(
       roc_obj,
       method = "bootstrap",
       boot.n = boot_n,
       conf.level = 0.95
     ),
     silent = TRUE
   )
   
   if (inherits(auc_ci, "try-error")) {
     ci_low <- NA
     ci_high <- NA
   } else {
     ci_low <- as.numeric(auc_ci[1])
     ci_high <- as.numeric(auc_ci[3])
   }
   
   # Optimal threshold by Youden index
   best_coords <- pROC::coords(
     roc_obj,
     x = "best",
     best.method = "youden",
     ret = c("threshold", "sensitivity", "specificity"),
     transpose = FALSE
   )
   
   threshold <- as.numeric(best_coords["threshold"])
   
   metrics <- get_classification_metrics(
     outcome = sub_df$decliner_bin,
     probs = probs,
     threshold = threshold
   )
   
   data.frame(
     model = model_name,
     formula = paste(deparse(formula_obj), collapse = " "),
     AUC = auc_val,
     CI_low = ci_low,
     CI_high = ci_high,
     N = nrow(sub_df),
     n_non_decliners = sum(sub_df$decliner_bin == 0),
     n_decliners = sum(sub_df$decliner_bin == 1),
     metrics
   )
}

## ---- 6.5 Generate models with all single and double biomarkers combinations ----

candidate_results <- data.frame()

for (base_name in names(base_formulas)) {
  
  base_formula <- base_formulas[[base_name]]
  
  # Base model alone
  res_base <- evaluate_candidate_model(
    df = ad_df,
    formula_obj = base_formula,
    model_name = base_name
  )
  
  candidate_results <- rbind(candidate_results, res_base)
  
  # Base + one biomarker
  for (b in candidate_biomarkers) {
    
    formula_1bio <- update(
      base_formula,
      paste(". ~ . +", b)
    )
    
    res_1bio <- evaluate_candidate_model(
      df = ad_df,
      formula_obj = formula_1bio,
      model_name = paste0(base_name, " + ", b)
    )
    
    candidate_results <- rbind(candidate_results, res_1bio)
  }
  
  # Base + two biomarkers
  biomarker_pairs <- combn(candidate_biomarkers, 2, simplify = FALSE)
  
  for (pair in biomarker_pairs) {
    
    formula_2bio <- update(
      base_formula,
      paste(". ~ . +", paste(pair, collapse = " + "))
    )
    
    res_2bio <- evaluate_candidate_model(
      df = ad_df,
      formula_obj = formula_2bio,
      model_name = paste0(base_name, " + ", paste(pair, collapse = " + "))
    )
    
    candidate_results <- rbind(candidate_results, res_2bio)
  }
}

## ---- 6.6 Order results by AUC and save results ----
 
candidate_results <- candidate_results[order(-candidate_results$AUC), ]

cat("\nTop A+ candidate models:\n")
print(head(candidate_results, 20))

write.csv(
  candidate_results,
  file = "results_functional_prediction/Apos_candidate_models_AUC_PPV_NPV.csv",
  row.names = FALSE
)

## ---- 6.7 Plot top models ----

top_models_plot <- candidate_results[1:min(15, nrow(candidate_results)), ]

top_models_plot$model <- factor(
  top_models_plot$model,
  levels = rev(top_models_plot$model)
)

p_top_ad <- ggplot(
  top_models_plot,
  aes(
    x = AUC,
    y = model
  )
) +
  geom_vline(xintercept = 0.80, linetype = "dashed") +
  geom_errorbarh(
    aes(xmin = CI_low, xmax = CI_high),
    height = 0.25
  ) +
  geom_point(size = 3) +
  labs(
    title = "Top A+ candidate prediction models",
    subtitle = "AUC with 95% bootstrap CI; PPV/NPV calculated at Youden-optimal threshold",
    x = "AUC",
    y = NULL
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8)
  )

print(p_top_ad)

ggsave(
  filename = "results_functional_prediction/Apos_top_candidate_models_AUC_plot.png",
  plot = p_top_ad,
  width = 12,
  height = 7,
  dpi = 300
) 
 
# ============================================================================= #
# 7. LEAVE-ONE-OUT SIMPLIFICATION ----
# ============================================================================= #

## ---- 7.1 Select candidate models to simplify ----

simplification_models <- list(
  
  # Model 1: clinically scalable plasma model
  SUCOG_plasma_NfL = list(
    full = decliner_bin ~ Age_V1 + sex + APOE_binary + Total_SUCOG +
      BIO_NFL_PLASMA_Elecsys_UGOT_Roche
  ),
  
  # Model 2: strong plasma model
  SUCOG_plasma_NfL_plasma_pTau217_ratio = list(
    full = decliner_bin ~ Age_V1 + sex + APOE_binary + Total_SUCOG +
      BIO_NFL_PLASMA_Elecsys_UGOT_Roche +
      plasma_pTau217_Ab42_ratio
  ),
  
  # Model 3: best-performing tau/neurodegeneration model
  SUCOG_plasma_NfL_CSF_pTau217 = list(
    full = decliner_bin ~ Age_V1 + sex + APOE_binary + Total_SUCOG +
      BIO_NFL_PLASMA_Elecsys_UGOT_Roche +
      BIO_pTau217_CSF_MSD_LILLY_in.house
  ),
  
  # Model 4: best clinically sensible model
  SUCOG_plasma_NfL_plasma_pTau217 = list(
    full = decliner_bin ~ Age_V1 + sex + APOE_binary + Total_SUCOG +
      BIO_NFL_PLASMA_Elecsys_UGOT_Roche +
      BIO_pTau217_PLASMA_Elecsys_UGOT_Roche
  )
)

## ---- 7.2 Leave-one-out formula ----

remove_predictor_from_formula <- function(formula_obj, predictor_to_remove) {
  
  outcome <- all.vars(formula_obj)[1]
  predictors <- attr(terms(formula_obj), "term.labels")
  
  predictors_reduced <- predictors[predictors != predictor_to_remove]
  
  new_formula <- as.formula(
    paste(outcome, "~", paste(predictors_reduced, collapse = " + "))
  )
  
  return(new_formula)
}

## ---- 7.3 Fit all combinations for a specific model ----

evaluate_model_for_simplification <- function(df,
                                              formula_obj,
                                              model_name,
                                              model_version,
                                              boot_n = 2000) {
  
  vars_needed <- all.vars(formula_obj)
  sub_df <- df[, vars_needed]
  sub_df <- sub_df[complete.cases(sub_df), ]
  
  if (nrow(sub_df) == 0) return(NULL)
  if (length(unique(sub_df$decliner_bin)) < 2) return(NULL)
  
  model <- try(
    glm(
      formula_obj,
      data = sub_df,
      family = binomial
    ),
    silent = TRUE
  )
  
  if (inherits(model, "try-error")) return(NULL)
  
  probs <- try(
    predict(model, type = "response"),
    silent = TRUE
  )
  
  if (inherits(probs, "try-error")) return(NULL)
  
  roc_obj <- try(
    pROC::roc(
      response = sub_df$decliner_bin,
      predictor = probs,
      quiet = TRUE
    ),
    silent = TRUE
  )
  
  if (inherits(roc_obj, "try-error")) return(NULL)
  
  auc_val <- as.numeric(pROC::auc(roc_obj))
  
  auc_ci <- try(
    pROC::ci.auc(
      roc_obj,
      method = "bootstrap",
      boot.n = boot_n,
      conf.level = 0.95
    ),
    silent = TRUE
  )
  
  if (inherits(auc_ci, "try-error")) {
    ci_low <- NA
    ci_high <- NA
  } else {
    ci_low <- as.numeric(auc_ci[1])
    ci_high <- as.numeric(auc_ci[3])
  }
  
  data.frame(
    model_name = model_name,
    model_version = model_version,
    formula = deparse(formula_obj),
    AUC = auc_val,
    CI_low = ci_low,
    CI_high = ci_high,
    N = nrow(sub_df),
    n_non_decliners = sum(sub_df$decliner_bin == 0),
    n_decliners = sum(sub_df$decliner_bin == 1),
    AIC = AIC(model),
    BIC = BIC(model)
  )
}

## ---- 7.4 Function to compare full model with reduced ones ----

compare_simplified_models <- function(df,
                                      full_formula,
                                      model_name,
                                      predictors_to_test = c(
                                        "Age_V1",
                                        "sex",
                                        "APOE_binary"
                                      ),
                                      boot_n = 2000) {
  
  results <- data.frame()
  
  # Evaluate full model
  full_res <- evaluate_model_for_simplification(
    df = df,
    formula_obj = full_formula,
    model_name = model_name,
    model_version = "Full model",
    boot_n = boot_n
  )
  
  if (is.null(full_res)) return(NULL)
  
  # Add comparison columns immediately so all rows have same columns
  full_res$delta_AUC_vs_full <- 0
  full_res$delong_p_vs_full <- NA
  full_res$lrt_p_vs_full <- NA
  
  results <- rbind(results, full_res)
  
  # Fit full model once for paired comparisons
  vars_full <- all.vars(full_formula)
  full_df <- df[, vars_full]
  full_df <- full_df[complete.cases(full_df), ]
  
  full_model <- glm(
    full_formula,
    data = full_df,
    family = binomial
  )
  
  full_probs <- predict(full_model, type = "response")
  
  full_roc <- pROC::roc(
    response = full_df$decliner_bin,
    predictor = full_probs,
    quiet = TRUE
  )
  
  full_auc <- as.numeric(pROC::auc(full_roc))
  
  # Loop over variables to remove
  for (pred in predictors_to_test) {
    
    reduced_formula <- remove_predictor_from_formula(
      formula_obj = full_formula,
      predictor_to_remove = pred
    )
    
    # Use the full model complete-case dataset for paired comparison
    paired_df <- full_df
    
    reduced_model <- try(
      glm(
        reduced_formula,
        data = paired_df,
        family = binomial
      ),
      silent = TRUE
    )
    
    if (inherits(reduced_model, "try-error")) next
    
    reduced_probs <- predict(reduced_model, type = "response")
    
    reduced_roc <- pROC::roc(
      response = paired_df$decliner_bin,
      predictor = reduced_probs,
      quiet = TRUE
    )
    
    reduced_auc <- as.numeric(pROC::auc(reduced_roc))
    
    auc_ci <- try(
      pROC::ci.auc(
        reduced_roc,
        method = "bootstrap",
        boot.n = boot_n,
        conf.level = 0.95
      ),
      silent = TRUE
    )
    
    if (inherits(auc_ci, "try-error")) {
      ci_low <- NA
      ci_high <- NA
    } else {
      ci_low <- as.numeric(auc_ci[1])
      ci_high <- as.numeric(auc_ci[3])
    }
    
    # DeLong test: full vs reduced AUC
    delong_test <- try(
      pROC::roc.test(
        full_roc,
        reduced_roc,
        method = "delong",
        paired = TRUE
      ),
      silent = TRUE
    )
    
    if (inherits(delong_test, "try-error")) {
      delong_p <- NA
    } else {
      delong_p <- delong_test$p.value
    }
    
    # LRT: reduced vs full model
    lrt_test <- try(
      anova(
        reduced_model,
        full_model,
        test = "Chisq"
      ),
      silent = TRUE
    )
    
    if (inherits(lrt_test, "try-error")) {
      lrt_p <- NA
    } else {
      lrt_p <- lrt_test$`Pr(>Chi)`[2]
    }
    
    reduced_res <- data.frame(
      model_name = model_name,
      model_version = paste0("Without ", pred),
      formula = deparse(reduced_formula),
      AUC = reduced_auc,
      CI_low = ci_low,
      CI_high = ci_high,
      N = nrow(paired_df),
      n_non_decliners = sum(paired_df$decliner_bin == 0),
      n_decliners = sum(paired_df$decliner_bin == 1),
      AIC = AIC(reduced_model),
      BIC = BIC(reduced_model),
      delta_AUC_vs_full = reduced_auc - full_auc,
      delong_p_vs_full = delong_p,
      lrt_p_vs_full = lrt_p
    )
    
    results <- rbind(results, reduced_res)
  }
  
  return(results)
}

## ---- 7.5 Run simplified comparison for all candidate models ----

simplification_results <- data.frame()

for (m in names(simplification_models)) {
  
  cat("\n========================================\n")
  cat("Testing simplification for:", m, "\n")
  cat("========================================\n")
  
  res <- compare_simplified_models(
    df = data_Apos,
    full_formula = simplification_models[[m]]$full,
    model_name = m,
    predictors_to_test = c("Age_V1", "sex", "APOE_binary"),
    boot_n = 2000
  )
  
  if (!is.null(res)) {
    simplification_results <- rbind(simplification_results, res)
  }
}

cat("\nSimplification results:\n")
print(simplification_results)

# Save results
write.csv(
  simplification_results,
  file = "results_functional_prediction/Apos_model_simplification_results.csv",
  row.names = FALSE
)

## ---- 7.6 Plot simplification results ----

simplification_plot_df <- simplification_results

simplification_plot_df$model_version <- factor(
  simplification_plot_df$model_version,
  levels = c(
    "Full model",
    "Without Age_V1",
    "Without sex",
    "Without APOE_binary"
  )
)

p_simplify <- ggplot(
  simplification_plot_df,
  aes(
    x = AUC,
    y = model_version
  )
) +
  geom_vline(xintercept = 0.80, linetype = "dashed") +
  geom_errorbarh(
    aes(xmin = CI_low, xmax = CI_high),
    height = 0.20
  ) +
  geom_point(size = 3) +
  facet_wrap(~ model_name, scales = "free_y") +
  labs(
    title = "Can the A+ prediction models be simplified?",
    subtitle = "Effect of removing Age, Sex, or APOE from selected top models",
    x = "AUC",
    y = NULL
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 9)
  )

print(p_simplify)

ggsave(
  filename = "results_functional_prediction/Apos_model_simplification_auc_plot.png",
  plot = p_simplify,
  width = 13,
  height = 6,
  dpi = 300
) 
 
# ============================================================================= #
# 8. 10-FOLD CROSS-VALIDATION OF SELECTED CANDIDATE MODLES ----
# ============================================================================= #

## ---- 8.1 Define final SUCOG candidate models ----

cv_df <- data_Apos

cv_models <- list(
  
  SUCOG_only =
    decliner_bin ~ Total_SUCOG,
  
  SUCOG_plasma_NfL =
    decliner_bin ~ Total_SUCOG +
    BIO_NFL_PLASMA_Elecsys_UGOT_Roche,
  
  SUCOG_plasma_NfL_CSF_pTau217 =
    decliner_bin ~ Total_SUCOG +
    BIO_NFL_PLASMA_Elecsys_UGOT_Roche +
    BIO_pTau217_CSF_MSD_LILLY_in.house,
  
  SUCOG_plasma_NfL_plasma_pTau217 =
    decliner_bin ~ Total_SUCOG +
    BIO_NFL_PLASMA_Elecsys_UGOT_Roche +
    BIO_pTau217_PLASMA_Elecsys_UGOT_Roche,
  
  SUCOG_plasma_NfL_adjusted =
    decliner_bin ~ Age_V1 + sex + APOE_binary +
    Total_SUCOG +
    BIO_NFL_PLASMA_Elecsys_UGOT_Roche,
  
  SUCOG_plasma_NfL_CSF_pTau217_adjusted =
    decliner_bin ~ Age_V1 + sex + APOE_binary +
    Total_SUCOG +
    BIO_NFL_PLASMA_Elecsys_UGOT_Roche +
    BIO_pTau217_CSF_MSD_LILLY_in.house,
  
  SUCOG_plasma_NfL_plasma_pTau217_adjusted =
    decliner_bin ~ Age_V1 + sex + APOE_binary +
    BIO_NFL_PLASMA_Elecsys_UGOT_Roche +
    BIO_pTau217_PLASMA_Elecsys_UGOT_Roche
)

## ---- 8.2 Function for SUCOG cross-validation of the models ----

run_repeated_cv_auc <- function(df,
                                formula_obj,
                                model_name,
                                k = 10,
                                repeats = 50) {
  
  vars_needed <- all.vars(formula_obj)
  sub_df <- df[, vars_needed]
  sub_df <- sub_df[complete.cases(sub_df), ]
  
  if (nrow(sub_df) == 0) return(NULL)
  if (length(unique(sub_df$decliner_bin)) < 2) return(NULL)
  
  all_preds <- data.frame()
  fold_aucs <- c()
  
  for (r in 1:repeats) {
    
    # Stratified fold assignment
    fold_id <- rep(NA, nrow(sub_df))
    
    for (class_value in unique(sub_df$decliner_bin)) {
      
      class_idx <- which(sub_df$decliner_bin == class_value)
      class_idx <- sample(class_idx)
      
      fold_id[class_idx] <- rep(
        1:k,
        length.out = length(class_idx)
      )
    }
    
    for (fold in 1:k) {
      
      train_df <- sub_df[fold_id != fold, ]
      test_df  <- sub_df[fold_id == fold, ]
      
      # Need both outcome classes in training and test
      if (length(unique(train_df$decliner_bin)) < 2) next
      if (length(unique(test_df$decliner_bin)) < 2) next
      
      model <- try(
        glm(
          formula_obj,
          data = train_df,
          family = binomial
        ),
        silent = TRUE
      )
      
      if (inherits(model, "try-error")) next
      
      probs <- try(
        predict(
          model,
          newdata = test_df,
          type = "response"
        ),
        silent = TRUE
      )
      
      if (inherits(probs, "try-error")) next
      
      roc_fold <- try(
        pROC::roc(
          response = test_df$decliner_bin,
          predictor = probs,
          quiet = TRUE
        ),
        silent = TRUE
      )
      
      if (!inherits(roc_fold, "try-error")) {
        fold_aucs <- c(fold_aucs, as.numeric(pROC::auc(roc_fold)))
      }
      
      all_preds <- rbind(
        all_preds,
        data.frame(
          model = model_name,
          repeat_id = r,
          fold = fold,
          observed = test_df$decliner_bin,
          predicted_prob = probs
        )
      )
    }
  }
  
  if (nrow(all_preds) == 0) return(NULL)
  
  # Overall pooled CV AUC across all held-out predictions
  pooled_roc <- pROC::roc(
    response = all_preds$observed,
    predictor = all_preds$predicted_prob,
    quiet = TRUE
  )
  
  pooled_auc <- as.numeric(pROC::auc(pooled_roc))
  
  data.frame(
    model = model_name,
    pooled_cv_AUC = pooled_auc,
    mean_fold_AUC = mean(fold_aucs, na.rm = TRUE),
    sd_fold_AUC = sd(fold_aucs, na.rm = TRUE),
    median_fold_AUC = median(fold_aucs, na.rm = TRUE),
    min_fold_AUC = min(fold_aucs, na.rm = TRUE),
    max_fold_AUC = max(fold_aucs, na.rm = TRUE),
    N = nrow(sub_df),
    n_repeats = repeats,
    n_folds = k,
    n_valid_folds = length(fold_aucs)
  )
}

## ---- 8.3 Run 10-fold cross-validation for each SUCOG  model ----

cv_results <- data.frame()

for (m in names(cv_models)) {
  
  cat("\n========================================\n")
  cat("Running repeated 10-fold CV for:", m, "\n")
  cat("========================================\n")
  
  res <- run_repeated_cv_auc(
    df = cv_df,
    formula_obj = cv_models[[m]],
    model_name = m,
    k = 10,
    repeats = 50
  )
  
  if (!is.null(res)) {
    cv_results <- rbind(cv_results, res)
  }
}

cv_results <- cv_results[order(-cv_results$pooled_cv_AUC), ]

cat("\nRepeated 10-fold CV results:\n")
print(cv_results)

# Save results
write.csv(
  cv_results,
  file = "results_functional_decline/repeated_10foldCV_results_no_caret.csv",
  row.names = FALSE
)

## ---- 8.4 Forest plot of SUCOG candidate mdoels after cross-validation ----

cv_results$model <- factor(
  cv_results$model,
  levels = rev(cv_results$model)
)

p_cv <- ggplot(
  cv_results,
  aes(
    x = pooled_cv_AUC,
    y = model
  )
) +
  geom_vline(xintercept = 0.80, linetype = "dashed") +
  geom_errorbarh(
    aes(
      xmin = mean_fold_AUC - sd_fold_AUC,
      xmax = mean_fold_AUC + sd_fold_AUC
    ),
    height = 0.20
  ) +
  geom_point(size = 3) +
  labs(
    title = "Repeated 10-fold cross-validated AUC",
    subtitle = "Point = pooled held-out AUC; bars = mean fold AUC ± SD",
    x = "Cross-validated AUC",
    y = NULL
  ) +
  theme_bw()

print(p_cv)

ggsave(
  filename = "results_functional_prediction/repeated_10foldCV_AUC_plot_no_caret.png",
  plot = p_cv,
  width = 10,
  height = 5,
  dpi = 300
) 

## ---- 8.5 Define final PACC candidate models ----

pacc_cv_df <- data_Apos

pacc_cv_models <- list(
  
  PACC_only =
    decliner_bin ~ PACC,
  
  PACC_plasma_NfL =
    decliner_bin ~ PACC +
    BIO_NFL_PLASMA_Elecsys_UGOT_Roche,
  
  PACC_plasma_NfL_CSF_pTau217 =
    decliner_bin ~ PACC +
    BIO_NFL_PLASMA_Elecsys_UGOT_Roche +
    BIO_pTau217_CSF_MSD_LILLY_in.house,
  
  PACC_plasma_NfL_adjusted =
    decliner_bin ~ Age_V1 + sex + APOE_binary +
    PACC +
    BIO_NFL_PLASMA_Elecsys_UGOT_Roche,
  
  PACC_plasma_NfL_CSF_pTau217_adjusted =
    decliner_bin ~ Age_V1 + sex + APOE_binary +
    PACC +
    BIO_NFL_PLASMA_Elecsys_UGOT_Roche +
    BIO_pTau217_CSF_MSD_LILLY_in.house,
  
  PACC_plasma_NfL_plasma_pTau217 =
    decliner_bin ~ PACC +
    BIO_NFL_PLASMA_Elecsys_UGOT_Roche +
    BIO_pTau217_PLASMA_Elecsys_UGOT_Roche,
  
  PACC_plasma_NfL_plasma_pTau217_adjusted =
    decliner_bin ~ Age_V1 + sex + APOE_binary +
    PACC +
    BIO_NFL_PLASMA_Elecsys_UGOT_Roche +
    BIO_pTau217_PLASMA_Elecsys_UGOT_Roche
)

## ---- 8.6 Run 10-fold cross-validation for each PACC model ----

pacc_cv_results <- data.frame()

for (m in names(pacc_cv_models)) {
  
  cat("\n========================================\n")
  cat("Running repeated 10-fold CV for:", m, "\n")
  cat("========================================\n")
  
  res <- run_repeated_cv_auc(
    df = pacc_cv_df,
    formula_obj = pacc_cv_models[[m]],
    model_name = m,
    k = 10,
    repeats = 50
  )
  
  if (!is.null(res)) {
    pacc_cv_results <- rbind(pacc_cv_results, res)
  }
}

# Rank by pooled CV AUC
pacc_cv_results <- pacc_cv_results[
  order(-pacc_cv_results$pooled_cv_AUC),
]

cat("\nRepeated 10-fold CV results - PACC models:\n")
print(pacc_cv_results)

write.csv(
  pacc_cv_results,
  file = "results_functional_prediction/repeated_10foldCV_PACC_models.csv",
  row.names = FALSE
) 
 
## ---- 8.7 Forest plot of PACC candidate mdoels after cross-validation ----

pacc_cv_results$model <- factor(
  pacc_cv_results$model,
  levels = rev(pacc_cv_results$model)
)

p_pacc_cv <- ggplot(
  pacc_cv_results,
  aes(
    x = pooled_cv_AUC,
    y = model
  )
) +
  geom_vline(xintercept = 0.80, linetype = "dashed") +
  geom_errorbarh(
    aes(
      xmin = mean_fold_AUC - sd_fold_AUC,
      xmax = mean_fold_AUC + sd_fold_AUC
    ),
    height = 0.20
  ) +
  geom_point(size = 3) +
  labs(
    title = "Repeated 10-fold cross-validated AUC - PACC models",
    subtitle = "Point = pooled held-out AUC; bars = mean fold AUC ± SD",
    x = "Cross-validated AUC",
    y = NULL
  ) +
  theme_bw()

print(p_pacc_cv)

ggsave(
  filename = "results_functional_prediction/repeated_10foldCV_PACC_AUC_plot.png",
  plot = p_pacc_cv,
  width = 10,
  height = 5,
  dpi = 300
) 

# ============================================================================= #
# 9. FINAL MODEL SELECTION, CROSS-VALIDATION, PPV/NPV AND ROC PANELS ----
# ============================================================================= #

## ---- 9.1 Function for repeated 10-fold CV with stored predictions ----

run_repeated_cv_with_predictions <- function(df,
                                             formula_obj,
                                             model_name,
                                             k = 10,
                                             repeats = 50) {
  
  vars_needed <- all.vars(formula_obj)
  sub_df <- df[, vars_needed]
  sub_df <- sub_df[complete.cases(sub_df), ]
  
  if (nrow(sub_df) == 0) return(NULL)
  if (length(unique(sub_df$decliner_bin)) < 2) return(NULL)
  
  all_preds <- data.frame()
  
  for (r in 1:repeats) {
    
    # Stratified fold assignment
    fold_id <- rep(NA, nrow(sub_df))
    
    for (class_value in unique(sub_df$decliner_bin)) {
      class_idx <- which(sub_df$decliner_bin == class_value)
      class_idx <- sample(class_idx)
      
      fold_id[class_idx] <- rep(
        1:k,
        length.out = length(class_idx)
      )
    }
    
    for (fold in 1:k) {
      
      train_df <- sub_df[fold_id != fold, ]
      test_df  <- sub_df[fold_id == fold, ]
      
      if (length(unique(train_df$decliner_bin)) < 2) next
      if (length(unique(test_df$decliner_bin)) < 2) next
      
      model <- try(
        glm(
          formula_obj,
          data = train_df,
          family = binomial
        ),
        silent = TRUE
      )
      
      if (inherits(model, "try-error")) next
      
      probs <- try(
        predict(
          model,
          newdata = test_df,
          type = "response"
        ),
        silent = TRUE
      )
      
      if (inherits(probs, "try-error")) next
      
      all_preds <- rbind(
        all_preds,
        data.frame(
          model = model_name,
          repeat_id = r,
          fold = fold,
          observed = test_df$decliner_bin,
          predicted_prob = probs
        )
      )
    }
  }
  
  if (nrow(all_preds) == 0) return(NULL)
  
  return(all_preds)
}

## ---- 9.2 Define final selected models ----

final_sucog_models <- list(
  
  SUCOG_base_adjusted =
    decliner_bin ~ Total_SUCOG + Age_V1 + sex + APOE_binary,
  
  SUCOG_plasma_NfL_adjusted =
    decliner_bin ~ Total_SUCOG + Age_V1 + sex + APOE_binary +
    BIO_NFL_PLASMA_Elecsys_UGOT_Roche,
  
  SUCOG_plasma_NfL_CSF_pTau217_adjusted =
    decliner_bin ~ Total_SUCOG + Age_V1 + sex + APOE_binary +
    BIO_NFL_PLASMA_Elecsys_UGOT_Roche +
    BIO_pTau217_CSF_MSD_LILLY_in.house,
  
  SUCOG_plasma_NfL_plasma_pTau217_adjusted =
    decliner_bin ~ Total_SUCOG + Age_V1 + sex + APOE_binary +
    BIO_NFL_PLASMA_Elecsys_UGOT_Roche +
    BIO_pTau217_PLASMA_Elecsys_UGOT_Roche
)

## ---- 9.3 Run CV for each model ----

sucog_cv_predictions <- data.frame()

for (m in names(final_sucog_models)) {
  
  cat("\n========================================\n")
  cat("Running CV predictions for:", m, "\n")
  cat("========================================\n")
  
  pred_m <- run_repeated_cv_with_predictions(
    df = data_Apos,
    formula_obj = final_sucog_models[[m]],
    model_name = m,
    k = 10,
    repeats = 50
  )
  
  if (!is.null(pred_m)) {
    sucog_cv_predictions <- rbind(
      sucog_cv_predictions,
      pred_m
    )
  }
}

## ---- 9.4 Calculate AUC, NPV, PPV, optimal threshold and other metrics ----

sucog_cv_metrics <- data.frame()
sucog_cv_rocs <- list()

for (m in unique(sucog_cv_predictions$model)) {
  
  pred_df <- sucog_cv_predictions[
    sucog_cv_predictions$model == m,
  ]
  
  roc_obj <- pROC::roc(
    response = pred_df$observed,
    predictor = pred_df$predicted_prob,
    quiet = TRUE
  )
  
  sucog_cv_rocs[[m]] <- roc_obj
  
  auc_val <- as.numeric(pROC::auc(roc_obj))
  
  # Optimal threshold by Youden index
  best_coords <- pROC::coords(
    roc_obj,
    x = "best",
    best.method = "youden",
    ret = c("threshold", "sensitivity", "specificity"),
    transpose = FALSE
  )
  
  threshold <- best_coords[, "threshold"]
  
  pred_class <- ifelse(
    pred_df$predicted_prob >= threshold,
    1,
    0
  )
  
  TP <- sum(pred_class == 1 & pred_df$observed == 1)
  TN <- sum(pred_class == 0 & pred_df$observed == 0)
  FP <- sum(pred_class == 1 & pred_df$observed == 0)
  FN <- sum(pred_class == 0 & pred_df$observed == 1)
  
  sensitivity <- TP / (TP + FN)
  specificity <- TN / (TN + FP)
  PPV <- TP / (TP + FP)
  NPV <- TN / (TN + FN)
  accuracy <- (TP + TN) / nrow(pred_df)
  
  sucog_cv_metrics <- rbind(
    sucog_cv_metrics,
    data.frame(
      model = m,
      AUC = auc_val,
      threshold = threshold,
      sensitivity = sensitivity,
      specificity = specificity,
      PPV = PPV,
      NPV = NPV,
      accuracy = accuracy,
      TP = TP,
      TN = TN,
      FP = FP,
      FN = FN,
      total_predictions = nrow(pred_df)
    )
  )
}


cat("\nFinal SUCOG CV-based classification metrics:\n")
print(sucog_cv_metrics)


write.csv(
  sucog_cv_metrics,
  file = "results_functional_prediction/final_SUCOG_models_adjusted_CV_PPV_NPV.csv",
  row.names = FALSE
)

write.csv(
  sucog_cv_predictions,
  file = "results_functional_prediction/final_SUCOG_models_adjusted_CV_predictions.csv",
  row.names = FALSE
)

## ---- 9.5 Forest plot of CV AUCs with 95% CI  ----

roc_plot_data <- data.frame()

for (m in names(sucog_cv_rocs)) {
  
  roc_obj <- sucog_cv_rocs[[m]]
  
  roc_df <- data.frame(
    model = m,
    specificity = roc_obj$specificities,
    sensitivity = roc_obj$sensitivities
  )
  
  roc_plot_data <- rbind(
    roc_plot_data,
    roc_df
  )
}


# Add AUC values to model labels
auc_labels <- sucog_cv_metrics[, c("model", "AUC")]

roc_plot_data <- merge(
  roc_plot_data,
  auc_labels,
  by = "model"
)

roc_plot_data$model_label <- paste0(
  roc_plot_data$model,
  "\nAUC = ",
  round(roc_plot_data$AUC, 3)
)

## ---- 9.6 Plot ROC curves with CV AUCs  ----

p_sucog_roc_panel <- ggplot(
  roc_plot_data,
  aes(
    x = 1 - specificity,
    y = sensitivity
  )
) +
  geom_line(linewidth = 1) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed"
  ) +
  facet_wrap(~ model_label) +
  labs(
    title = "Cross-validated ROC curves for final SUCOG models",
    subtitle = "Repeated 10-fold CV held-out predictions",
    x = "1 - Specificity",
    y = "Sensitivity"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold")
  )

print(p_sucog_roc_panel)

ggsave(
  filename = "results_funcional_prediction/final_SUCOG_models_CV_adjusted_ROC_panel.png",
  plot = p_sucog_roc_panel,
  width = 12,
  height = 8,
  dpi = 300
)



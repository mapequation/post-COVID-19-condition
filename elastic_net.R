# Load required libraries

library("org.Hs.eg.db")
library(missMethyl)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

library(glmnet)
library(tidyverse)
library(pscl)
library(effectsize)
library(dplyr)

symptoms <- "other"

############# Function to run elastic net with alpha tuning #############

run_elastic_net_with_alpha_tuning <- function(data, symptoms, alpha_grid = seq(0.1, 0.9, by = 0.1)) {
  # Prepare response
  successes <- data[,paste0(symptoms,".symptoms.positive")]
  failures <- data[,paste0(symptoms,".symptoms.total")] - data[,paste0(symptoms,".symptoms.positive")]
  y_binomial <- cbind(successes, failures)
  
  # Extract CpG feature columns
  cpg_cols <- grep("^cg", colnames(data), value = TRUE)
  X <- as.matrix(data[, cpg_cols])
  
  # Initialize result storage
  best_model <- NULL
  best_alpha <- NULL
  best_lambda <- NULL
  lowest_cvm <- Inf
  
  # Search over alpha grid
  for (alpha_val in alpha_grid) {
    set.seed(123)
    cv_fit <- cv.glmnet(X, y_binomial, alpha = alpha_val, family = "binomial", standardize = TRUE)
    
    if (min(cv_fit$cvm) < lowest_cvm) {
      lowest_cvm <- min(cv_fit$cvm)
      best_model <- cv_fit
      best_alpha <- alpha_val
      best_lambda <- cv_fit$lambda.min
    }
  }
  
  # Extract coefficients at best lambda
  coefs <- coef(best_model, s = best_lambda)
  coefs_df <- as.data.frame(as.matrix(coefs))
  coefs_df <- tibble(
    feature = rownames(coefs_df),
    coefficient = coefs_df[, 1]
  ) %>%
    filter(feature != "(Intercept)", coefficient != 0) %>%
    mutate(abs_coef = abs(coefficient)) %>%
    arrange(desc(abs_coef))  # Rank by absolute coefficient
  
  return(list(
    ranked_features = coefs_df,
    features = coefs_df$feature,
    best_alpha = best_alpha,
    best_lambda = best_lambda
  ))
}

############# Stability Selection Function #############

stability_selection <- function(data, alpha_val = 0.5, n_iter = 100, subsample_frac = 0.8) {
  successes <- data[,paste0(symptoms,".symptoms.positive")]
  failures <- data[,paste0(symptoms,".symptoms.total")] - data[,paste0(symptoms,".symptoms.positive")]
  y_binomial <- cbind(successes, failures)
  
  cpg_cols <- grep("^cg", colnames(data), value = TRUE)
  X <- as.matrix(data[, cpg_cols])
  
  feature_counts <- setNames(rep(0, length(cpg_cols)), cpg_cols)
  
  set.seed(42)
  for (i in 1:n_iter) {
    idx <- sample(1:nrow(X), size = floor(nrow(X) * subsample_frac), replace = TRUE)
    X_sub <- X[idx, ]
    y_sub <- y_binomial[idx, ]
    
    fit <- cv.glmnet(X_sub, y_sub, alpha = alpha_val, family = "binomial", standardize = TRUE)
    coefs <- coef(fit, s = fit$lambda.min)
    nonzero <- rownames(coefs)[which(coefs != 0)]
    nonzero <- setdiff(nonzero, "(Intercept)")
    
    feature_counts[nonzero] <- feature_counts[nonzero] + 1
  }
  
  stability_df <- tibble(
    feature = names(feature_counts),
    selection_frequency = feature_counts / n_iter
  ) %>%
    arrange(desc(selection_frequency))
  
  return(stability_df)
}

############# Cohen's d and 95% CI #############

calculate_cohens_d_all_cpgs <- function(data, symptom_positive_col = "neurological.symptoms.positive") {
  # Create binary symptom group: presence vs absence
  data <- data %>%
    mutate(symptom_group = ifelse(.data[[symptom_positive_col]] > 0, "positive", "negative")) %>%
    mutate(symptom_group = factor(symptom_group, levels = c("positive", "negative")))
  
  cpg_cols <- grep("^cg", colnames(data), value = TRUE)
  
  results <- list()
  
  for (cg in cpg_cols) {
    # Prepare data frame for effectsize
    df <- data %>%
      select(all_of(c(cg, "symptom_group"))) %>%
      rename(methylation = all_of(cg))
    
    # Only calculate if both groups have >1 sample
    group_counts <- table(df$symptom_group)
    if (all(group_counts >= 2)) {
      d_res <- tryCatch(
        cohens_d(methylation ~ symptom_group, data = df, ci = 0.95),
        error = function(e) NULL
      )
      
      if (!is.null(d_res)) {
        results[[cg]] <- tibble(
          feature = cg,
          cohens_d = d_res$Cohens_d,
          ci_lower = d_res$CI_low,
          ci_upper = d_res$CI_high
        )
      }
    }
  }
  
  # Combine all results into one dataframe
  results_df <- bind_rows(results) %>%
    arrange(desc(abs(cohens_d)))
  
  return(results_df)
}


# Load data
main_dir <- "COVUM/elastic_net/"
data.3months <- read.csv(paste0(main_dir, "glmnet_3months.csv"))
data.6months <- read.csv(paste0(main_dir, "glmnet_6months.csv"))
data.12months <- read.csv(paste0(main_dir, "glmnet_12months.csv"))

# Run models with alpha tuning
result_3m <- run_elastic_net_with_alpha_tuning(data.3months, symptoms)
result_6m <- run_elastic_net_with_alpha_tuning(data.6months, symptoms)
result_12m <- run_elastic_net_with_alpha_tuning(data.12months, symptoms)

# Find intersection of selected CpGs
common_cpgs <- Reduce(intersect, list(result_3m$features, result_6m$features, result_12m$features))

# Print results
cat("Common CpGs across all time points:\n")
print(common_cpgs)

cat("\nBest alpha values:\n")
cat("3 months:", result_3m$best_alpha, "\n")
cat("6 months:", result_6m$best_alpha, "\n")
cat("12 months:", result_12m$best_alpha, "\n")

# Run stability selection on each time point
stab_3m <- stability_selection(data.3months, alpha_val = result_3m$best_alpha)
stab_6m <- stability_selection(data.6months, alpha_val = result_6m$best_alpha)
stab_12m <- stability_selection(data.12months, alpha_val = result_12m$best_alpha)

# Cohen's d and 95% CI
cohen_d_results.3m <- calculate_cohens_d_all_cpgs(data.3months, paste0(symptoms,".symptoms.positive"))
cohen_d_results.6m <- calculate_cohens_d_all_cpgs(data.6months, paste0(symptoms,".symptoms.positive"))
cohen_d_results.12m <- calculate_cohens_d_all_cpgs(data.12months, paste0(symptoms,".symptoms.positive"))

stability.3m <- c()
stability.6m <- c()
stability.12m <- c()
cohens_d.3m <- c()
cohens_d.6m <- c()
cohens_d.12m <- c()
ci_lower.3m <- c()
ci_lower.6m <- c()
ci_lower.12m <- c()
ci_upper.3m <- c()
ci_upper.6m <- c()
ci_upper.12m <- c()
N <- length(common_cpgs)
for (ix in 1:N){
  cpgs <- common_cpgs[ix]
  stability.3m[ix] <- stab_3m$selection_frequency[stab_3m$feature==cpgs]
  stability.6m[ix] <- stab_6m$selection_frequency[stab_6m$feature==cpgs]
  stability.12m[ix] <- stab_12m$selection_frequency[stab_12m$feature==cpgs]
  cohens_d.3m[ix] <- cohen_d_results.3m$cohens_d[cohen_d_results.3m$feature==cpgs]
  cohens_d.6m[ix] <- cohen_d_results.6m$cohens_d[cohen_d_results.6m$feature==cpgs]
  cohens_d.12m[ix] <- cohen_d_results.12m$cohens_d[cohen_d_results.12m$feature==cpgs]
  ci_lower.3m[ix] <- cohen_d_results.3m$ci_lower[cohen_d_results.3m$feature==cpgs]
  ci_lower.6m[ix] <- cohen_d_results.6m$ci_lower[cohen_d_results.6m$feature==cpgs]
  ci_lower.12m[ix] <- cohen_d_results.12m$ci_lower[cohen_d_results.12m$feature==cpgs]
  ci_upper.3m[ix] <- cohen_d_results.3m$ci_upper[cohen_d_results.3m$feature==cpgs]
  ci_upper.6m[ix] <- cohen_d_results.6m$ci_upper[cohen_d_results.6m$feature==cpgs]
  ci_upper.12m[ix] <- cohen_d_results.12m$ci_upper[cohen_d_results.12m$feature==cpgs]
}
results <- data.frame(
  feature=common_cpgs,
  stability.3m = stability.3m,
  stability.6m = stability.6m,
  stability.12m = stability.12m,
  cohens_d.3m = cohens_d.3m,
  cohens_d.6m = cohens_d.6m,  
  cohens_d.12m = cohens_d.12m,  
  ci_lower.3m = ci_lower.3m,  
  ci_upper.3m = ci_upper.3m,  
  ci_lower.6m = ci_lower.6m,  
  ci_upper.6m = ci_upper.6m,  
  ci_lower.12m = ci_lower.12m,
  ci_upper.12m = ci_upper.12m
  )

print(results)
